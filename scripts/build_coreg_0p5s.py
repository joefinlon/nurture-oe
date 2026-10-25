"""0.5 s x 15 m co-registration of the HALO cloud-file lidar profiles and CloudCube (any segment).

Target grid: the 0.5 s profiles of the HALO cloud file (``<date>_cloud_data.mat``) x HALO 15 m
``z`` (MSL). Lidar fields are native 0.5 s (no averaging): beta_p and delta_p at 532 nm, the
1064 nm fields and the phase mask ``cp`` (= C&S ``CloudMask/Type``). Units are inferred, not
declared (``io/halo_cloud.py``) [GAP].

Radar: each CloudCube gate goes to the 0.5 s window ``[t - 0.25 s, t + 0.25 s)`` containing the
time the lidar footprint reaches it (``coreg.cloudcube_to_halo_grid(assign="footprint")``);
linear-Z mean over the 4-6 rays; gates on the HALO MSL datum [DECISION D23]. The lidar is taken
as body-fixed like the radar (HALO ``surface/offset_angle`` = cos(pitch) cos(roll)), so the
relative geometric offset is zero and this equals time matching. An empirical ~1 s radar-lidar
lag is unexplained and not applied (``TIME_OFFSET_S`` = 0) [GAP G27], [GAP G23]. Each CloudCube
granule is aggregated on its own and every profile belongs to one granule, so no pixel mixes Ka
calibrations [GAP G24]. Gas-corrected Z and DWR are display-only.

Radar pixels are censored means [OPEN O21]; about one gate in three is duplicated on the 15 m
grid [OPEN O29]. Whether 0.5 s rows become y rows is S2's choice [OPEN O36]; beta_p rows are
valid only where the HSRL inversion holds [DECISION D21].

Usage: ``python scripts/build_coreg_0p5s.py [--start ISO --end ISO]`` (default 15:30-16:00 UTC).
Writes ``<cache_dir>/coreg_cloud_<YYYYMMDD>_<HHMM>-<HHMM>_0p5s.nc``.
"""

import argparse
import datetime as dt
import sys
import warnings

import numpy as np
import xarray as xr

import nurture_oe
from nurture_oe import coreg
from nurture_oe.config import cache_dir, cloudcube_granules, example_file
from nurture_oe.io.cloudcube import BANDS, read_cloudcube
from nurture_oe.io.halo_cloud import read_halo_cloud

SEGMENT = ("2026-02-17T15:30:00", "2026-02-17T16:00:00")
PAD = np.timedelta64(6, "s")
LIDAR_POINTING = "body_fixed"
TIME_OFFSET_S = 0.0
# output name -> cloud-file name (one place to change when the CF-compliant file arrives)
LIDAR_FIELDS = {"bsc532": "bsc", "dep532": "dep", "bsc1064": "bsc1064", "dep1064": "dep1064", "phase": "cp"}


def out_path(segment):
    t0, t1 = (np.datetime64(s, "s").astype(dt.datetime) for s in segment)
    return cache_dir() / f"coreg_cloud_{t0:%Y%m%d}_{t0:%H%M}-{t1:%H%M}_0p5s.nc"


def build(segment=SEGMENT):
    t0, t1 = (np.datetime64(s, "ns") for s in segment)
    path = example_file("halo_cloud")
    lid_pad = read_halo_cloud(path, variables=list(LIDAR_FIELDS.values()), time_range=(t0 - PAD, t1 + PAD))
    sfc = read_halo_cloud(path, group="surface", variables=["altSurface"], time_range=(t0 - PAD, t1 + PAD))
    win = coreg.subprofile_windows(lid_pad.time.values)
    keep = (win.t_lo.values >= t0) & (win.t_hi.values <= t1)
    win, lid, sfc = win.isel(time=keep), lid_pad.isel(profile=keep), sfc.isel(profile=keep)
    n = win.sizes["time"]

    grans = cloudcube_granules(t0, t1)
    gran_idx = np.full(n, -1, dtype="int8")
    zmin = {b: np.full(n, np.nan) for b in BANDS}
    parts = []
    for g, (gpath, gs, ge) in enumerate(grans):
        mine = (win.t_lo.values >= gs) & (win.t_hi.values <= ge)
        if not mine.any():
            continue
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            cc = read_cloudcube(gpath, time_range=(max(t0, gs) - PAD, min(t1, ge) + PAD))
        extra = {f"zgc_{b}": coreg.gas_corrected_display(cc, b) for b in ("w", "g")}
        # padded nav: footprint rays precede their profile by up to ~2.5 s
        grid = coreg.cloudcube_to_halo_grid(cc, lid_pad, win.isel(time=mine), extra=extra, assign="footprint",
                                            lidar_pointing=LIDAR_POINTING, time_offset_s=TIME_OFFSET_S)
        gran_idx[mine] = g
        for b in BANDS:
            zmin[b][mine] = cc.attrs.get(f"zmin_1km_dbz_{b}", np.nan)
        parts.append(grid)
    if not parts:
        raise RuntimeError(f"no CloudCube granule covers {segment}")
    radar = xr.concat(parts, dim="time", combine_attrs="drop_conflicts").reindex(time=win.time.values)
    assignment = parts[0].attrs["gate_assignment"]
    with np.errstate(invalid="ignore"):
        dwr = radar["zgc_w"] - radar["zgc_g"]
    dwr.attrs = {"long_name": "DWR W-G from gas-corrected-only (display) Z", "units": "dB"}

    out = radar.assign(dwr_wg_gc=dwr)
    for name, src in LIDAR_FIELDS.items():
        out[name] = (("time", "z"), lid[src].values, {**lid[src].attrs, "source_name": src})
    typ = lid[LIDAR_FIELDS["phase"]].values
    out["liquid"] = (("time", "z"), (typ == 0).astype("int8"), {"long_name": "cloud-file cp == 0 (WATER), 0.5 s"})
    out["ice"] = (("time", "z"), np.isin(typ, (1, 2)).astype("int8"), {"long_name": "cloud-file cp in {1, 2}, 0.5 s"})
    out["surface_alt_lidar"] = ("time", sfc["altSurface"].values, {**sfc["altSurface"].attrs, "long_name": "cloud-file surface/altSurface"})
    out["gps_alt"] = ("time", lid["gps_alt"].values, lid["gps_alt"].attrs)
    out["cc_granule"] = ("time", gran_idx, {"long_name": "index into global attribute source_cloudcube; -1 none"})
    for b in BANDS:
        out[f"zmin_1km_dbz_{b}"] = ("time", zmin[b], {"long_name": f"granule radar_sensitivity_{b}_dbz", "units": "dBZ"})
    out.attrs = {
        "title": f"HALO cloud-file lidar + CloudCube on the 0.5 s x HALO 15 m grid, {segment[0]} to {segment[1]} UTC",
        "source_halo_cloud": lid.attrs["source_file"],
        "source_cloudcube": [p.name for p, _, _ in grans],
        "nurture_oe_version": nurture_oe.__version__,
        "created": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/build_coreg_0p5s.py",
        "window_rule": win.attrs["rule"], "window_half_width_s": float(win.attrs["half_width_s"]),
        "aircraft_alt_reference": "halo (cloud-file nav/gps_alt, MSL) [DECISION D23]",
        "gate_assignment": assignment, "min_range_m": parts[0].attrs["min_range_m"],
        "lidar_units": "inferred, not declared by the file; see *_inferred attributes [GAP]",
        "gas_correction": "zgc_* = zm + 2*int katten dr, display only (appendix A §A.4.1)",
        "note": "Radar values are censored means over 4-6 rays per pixel [OPEN O21]; about one gate in three is "
                "duplicated on the 15 m grid [OPEN O29]. Lidar values are native 0.5 s and unscreened. Row "
                "construction from this file is S2's choice [OPEN O36]; beta_p rows only where HSRL-valid [DECISION D21].",
    }
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--start", default=SEGMENT[0])
    ap.add_argument("--end", default=SEGMENT[1])
    a = ap.parse_args(argv)
    seg = (a.start, a.end)
    ds = build(seg)
    path = out_path(seg)
    enc = {v: {"zlib": True, "complevel": 4} for v in ds.data_vars if ds[v].ndim == 2}
    ds.to_netcdf(path, encoding=enc)
    print(f"wrote {path} {dict(ds.sizes)}")
    print("rays per pixel (5/50/95 %):", np.nanpercentile(np.where(ds.n_rays > 0, ds.n_rays, np.nan), [5, 50, 95]))


if __name__ == "__main__":
    sys.exit(main())
