"""Task (d): co-register HALO + CloudCube for 2026-02-17 15:30-16:00 UTC and run the Ka coverage check.

Writes ``<cache_dir>/coreg_20260217_1530-1600.nc`` and ``docs/reports/ka_coverage_*.csv``.
"""

import datetime as dt
import sys

import numpy as np
import pandas as pd
import xarray as xr

import nurture_oe
from nurture_oe import coreg
from nurture_oe.checks.ka_coverage import ka_coverage
from nurture_oe.config import REPO_ROOT, cache_dir, example_file
from nurture_oe.io.cloudcube import read_cloudcube
from nurture_oe.io.halo_aerowv import read_halo_aerowv
from nurture_oe.io.halo_cloudsurface import read_halo_cloudsurface

SEGMENT = ("2026-02-17T15:30:00", "2026-02-17T16:00:00")
PAD = np.timedelta64(6, "s")
CLEARANCES_M = (0, 60, 120, 250, 500, 1000)
HALO_VARS = ["532_total_attn_bsc", "532_dep", "gps_alt", "Temperature", "Aerosol_ID"]


def build(segment=SEGMENT):
    t0, t1 = (np.datetime64(s, "ns") for s in segment)
    halo = read_halo_aerowv(example_file("halo_aerowv"), variables=HALO_VARS, time_range=segment)
    # keep only HALO profiles whose whole window lies inside the segment
    halo = halo.sel(time=(halo.time - coreg.HALF_WINDOW >= t0) & (halo.time + coreg.HALF_WINDOW <= t1))
    win = coreg.halo_windows(halo.time.values)
    cs = read_halo_cloudsurface(example_file("halo_cloudsurface"), variables=["Type", "altSurface"],
                                time_range=(t0 - PAD, t1 + PAD))
    cc_path = example_file("cloudcube_1521")
    cc = read_cloudcube(cc_path, time_range=(t0 - PAD, t1 + PAD))
    extra = {f"zgc_{b}": coreg.gas_corrected_display(cc, b) for b in ("w", "g")}
    grid = coreg.cloudcube_to_halo_grid(cc, halo, win, extra=extra)
    fr = coreg.cloudmask_fractions(cs, win)
    pm = coreg.phase_masks(fr)
    with np.errstate(invalid="ignore"):
        dwr = grid["zgc_w"] - grid["zgc_g"]
    dwr.attrs = {"long_name": "DWR W-G from gas-corrected-only (display) Z", "units": "dB"}
    out = xr.merge([halo, grid.drop_vars("time").assign_coords(time=halo.time),
                    fr.drop_vars("time").assign_coords(time=halo.time),
                    pm.drop_vars("time").assign_coords(time=halo.time).astype("int8")],
                   combine_attrs="drop")
    out["dwr_wg_gc"] = dwr.drop_vars("time").assign_coords(time=halo.time)
    sfc = np.nan_to_num(cs["altSurface"].values)[:, None]
    s, _ = coreg.window_sums(sfc, cs["time"].values, win)
    k, _ = coreg.window_sums(np.isfinite(cs["altSurface"].values)[:, None].astype(float), cs["time"].values, win)
    with np.errstate(invalid="ignore", divide="ignore"):
        out["surface_alt_lidar"] = ("time", (s / k)[:, 0], {"long_name": "window mean of C&S altSurface", "units": "m"})
    out.attrs = {
        "title": "HALO + CloudCube co-registered on the HALO grid, 2026-02-17 15:30-16:00 UTC (STM segment)",
        "source_halo_aerowv": halo.attrs["source_file"], "source_halo_cloudsurface": cs.attrs["source_file"],
        "source_cloudcube": cc.attrs["source_file"], "nurture_oe_version": nurture_oe.__version__,
        "created": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "window_rule": win.attrs["rule"], "aircraft_alt_reference": grid.attrs["aircraft_alt_reference"],
        "gate_assignment": grid.attrs["gate_assignment"], "min_range_m": grid.attrs["min_range_m"],
        **{k: v for k, v in cc.attrs.items() if k.startswith(("zmin_1km_dbz_", "radar_sensitivity_"))},
        "phase_rule": "liquid: class 0 in >=50% of sub-profiles; ice: classes 1-2 in >=50% (appendix A §A.10.5)",
        "gas_correction": "zgc_* = zm + 2*int katten dr, display only (appendix A §A.4.1)",
    }
    return out


def main():
    ds = build()
    path = cache_dir() / "coreg_20260217_1530-1600.nc"
    ds.to_netcdf(path)
    print(f"wrote {path} {dict(ds.sizes)}")
    ice = ds["ice"].values.astype(bool)
    summary, by_range = ka_coverage(ds, ice, ds["Temperature"].values)
    rep = REPO_ROOT / "docs" / "reports"
    hag = ds["z"].values[None, :] - ds["surface_alt_lidar"].values[:, None]
    sweep = []
    for clr in CLEARANCES_M:
        sm, _ = ka_coverage(ds, ice & (hag > clr))
        sweep.append(sm.assign(clearance_above_surface_m=clr))
    sweep = pd.concat(sweep, ignore_index=True)
    sweep.to_csv(rep / "ka_coverage_clearance_20260217_1530-1600.csv", index=False)
    print(sweep.to_string(index=False))
    summary.to_csv(rep / "ka_coverage_summary_20260217_1530-1600.csv", index=False)
    by_range.to_csv(rep / "ka_coverage_by_range_20260217_1530-1600.csv", index=False)
    print(summary.to_string(index=False))
    print(by_range.to_string(index=False))


if __name__ == "__main__":
    sys.exit(main())
