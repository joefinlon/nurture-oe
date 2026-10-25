"""Checks of the HALO 0.5 s cloud file (``Cloud_Plots/20260217_cloud_data.mat``, delivered 2026-10-23).

Read-only. Establishes what the interim file (no metadata) contains by comparison with the other
02-17 products, and sizes the radar-lidar footprint offset at 0.5 s [GAP G27], [OPEN O36]:

1. header dump -> ``docs/reports/headers/20260217_cloud_data.txt``;
2. grids, time and nav vs C&S and Aerosol/WV;
3. ``cloudphase/cp`` vs C&S ``CloudMask/Type`` over the whole flight; class labels in the PI's figure;
4. units / meaning of ``bsc``, ``bsc1064``, ``dep``, ``dep1064`` from 10 s window means ([t-5, t+5),
   fixed; no offset sweep) vs Aerosol/WV on 15:30-16:00;
5. statistics of undocumented variables;
6. radar gate offsets from the aircraft nadir in the 0.5 s file (``scripts/build_coreg_0p5s.py``
   must have run); ``surface/offset_angle`` vs cos(pitch) cos(roll) (lidar pointing model);
7. lag diagnostic: Spearman(ln beta_p, Z_W) in ice pixels vs radar-lidar time lag by altitude band,
   against the lags a nadir or a body-fixed lidar would imply.

Writes ``docs/reports/check_halo_cloud.json`` and ``docs/reports/figures/check_halo_cloud_lag.png``.
"""

import json
import re
import struct
import sys
import warnings
import zlib

import h5py
import numpy as np
import scipy.io as sio
import xarray as xr
from scipy.stats import spearmanr

from nurture_oe import coreg
from nurture_oe.config import REPO_ROOT, cache_dir, example_file
from nurture_oe.io.cloudcube import read_cloudcube
from nurture_oe.io.halo_aerowv import read_halo_aerowv
from nurture_oe.io.halo_cloud import read_halo_cloud

SEGMENT = (np.datetime64("2026-02-17T15:30:00", "ns"), np.datetime64("2026-02-17T16:00:00", "ns"))
PAD = np.timedelta64(10, "s")
LAGS_S = np.round(np.arange(-3.0, 1.0 + 1e-9, 0.25), 2)
BANDS_KM = ((1, 2), (2, 3), (3, 4), (4, 5))
REPORTS = REPO_ROOT / "docs" / "reports"


def pct(a, q=(1, 5, 25, 50, 75, 95, 99)):
    a = np.asarray(a, float)
    a = a[np.isfinite(a)]
    return {str(k): float(v) for k, v in zip(q, np.percentile(a, q))} if a.size else {}


def dump_header(path, out):
    lines = [f"# {path.name}  (h5py walk; MATLAB attrs shown)", f"userblock {h5py.File(path, 'r').userblock_size}"]
    with h5py.File(path, "r") as f:
        def visit(name, o):
            attrs = {k: (o.attrs[k].decode() if isinstance(o.attrs[k], bytes) else o.attrs[k])
                     for k in o.attrs if k != "MATLAB_fields"}
            if isinstance(o, h5py.Dataset):
                lines.append(f"D /{name}  shape={o.shape} dtype={o.dtype} chunks={o.chunks} "
                             f"compression={o.compression}  attrs={attrs}")
            else:
                lines.append(f"G /{name}  attrs={attrs}")
        f.visititems(visit)
    out.write_text("\n".join(lines) + "\n")


def fig_strings(path, words=("WATER", "MIX", "ICE", "SPEC", "AER.", "NLM")):
    """Order of first appearance of the class labels in a MATLAB .fig (v5) graphics workspace."""
    ws = sio.loadmat(path)["__function_workspace__"].tobytes()
    chunks, off = [], 8
    while off + 8 <= len(ws):
        typ, n = struct.unpack("<II", ws[off:off + 8])
        if typ >> 16:
            off += 8
            continue
        data = ws[off + 8:off + 8 + n]
        chunks.append(zlib.decompress(data) if typ == 15 else data)
        off += 8 + n + (0 if typ == 15 else (8 - n % 8) % 8)
    blob = b"".join(chunks)
    pos = {w: blob.find(w.encode()) for w in words}
    found = {w: p for w, p in pos.items() if p >= 0}
    return {"found": sorted(found, key=found.get), "missing": [w for w in words if w not in found],
            "title": [s.decode() for s in set(re.findall(rb"HALO [ -~]{3,40}", blob))]}


def win_mean(v, t, w, weights=None):
    x = v if weights is None else v * weights
    ok = np.isfinite(x)
    s, _ = coreg.window_sums(np.where(ok, x, 0.0), t, w)
    n, _ = coreg.window_sums(ok.astype(float), t, w)
    with np.errstate(invalid="ignore", divide="ignore"):
        return s / n


def unit_checks(cl, halo):
    w = coreg.halo_windows(halo.time.values)
    t = cl.time.values
    out = {}
    for a, b in (("bsc", "532_bsc"), ("bsc1064", "1064_bsc")):
        m, hb = win_mean(cl[a].values, t, w), halo[b].values
        for thr in (1e-4, 1e-2):
            ok = np.isfinite(m) & np.isfinite(hb) & (hb > thr)
            out[f"{a}_over_{b}_x1e3_hbsc>{thr}"] = {"n": int(ok.sum()), **pct(m[ok] / hb[ok] * 1e3, (5, 25, 50, 75, 95))}
    bsr = halo["532_bsr"].values
    for a, bp, lam in (("dep", "bsc", "532"), ("dep1064", "bsc1064", "1064")):
        dw = win_mean(cl[a].values, t, w, cl[bp].values) / win_mean(cl[bp].values, t, w)   # beta_p-weighted
        for lo, hi in ((0.5, 2), (2, 10), (10, np.inf)):
            ok = (np.isfinite(dw) & (bsr >= lo) & (bsr < hi) & np.isfinite(halo[f"{lam}_dep"].values)
                  & np.isfinite(halo[f"{lam}_aer_dep"].values))
            out[f"{a}_bsr[{lo},{hi})"] = {
                "n": int(ok.sum()),
                "median(dep - {}_dep)".format(lam): float(np.median(dw[ok] - halo[f"{lam}_dep"].values[ok])),
                "median(dep - {}_aer_dep)".format(lam): float(np.median(dw[ok] - halo[f"{lam}_aer_dep"].values[ok]))}
    ice = np.isin(cl.cp.values, (1, 2)) & (cl.bsc.values > 1e-5)
    out["ice_color_ratio_bsc1064_over_bsc_0p5s"] = pct(cl.bsc1064.values[ice] / cl.bsc.values[ice], (5, 25, 50, 75, 95))
    with np.errstate(invalid="ignore", divide="ignore"):
        cr = halo["1064_bsc"].values / halo["532_bsc"].values
    hice = (halo["Aerosol_ID"].values == 1) & (halo["532_bsc"].values > 1e-2)
    out["ice_color_ratio_aerowv_1064_bsc_over_532_bsc_10s"] = pct(cr[hice], (5, 25, 50, 75, 95))
    return out


def lag_diagnostic(cl_pad, win, cc):
    """Spearman(ln beta_p, Z_W) over ice pixels vs radar-lidar lag, by 1 km altitude band.

    The radar is re-aggregated by ray time with its rays shifted by ``lag`` (radar time - lidar
    time), so lags are not limited to the 0.5 s profile step. The vertex of a parabola through
    the three best lags gives the peak. Also returns the median offset a nadir-pointing lidar
    would imply (``footprint_time(lidar_pointing="nadir")``) for comparison.
    """
    cl = cl_pad.sel(profile=np.isin(cl_pad.time.values, win.time.values))
    alt = coreg.aircraft_alt_from_halo(cc, cl_pad)
    with np.errstate(invalid="ignore", divide="ignore"):
        lid = np.log(cl.bsc.values)
    ice = np.isin(cl.cp.values, (1, 2))
    zkm = cl.z.values / 1e3
    curves = {b: [] for b in BANDS_KM}
    for lag in LAGS_S:
        shifted = cc.assign_coords(time=cc.time.values - np.timedelta64(int(round(lag * 1e9)), "ns"))
        g = coreg.cloudcube_to_halo_grid(shifted, cl, win, fields=("zm_w",),
                                         aircraft_alt=alt.assign_coords(time=shifted.time.values))
        for b in BANDS_KM:
            m = (zkm >= b[0]) & (zkm < b[1])
            x, y = lid[:, m], g.zm_w.values[:, m]
            ok = ice[:, m] & np.isfinite(x) & np.isfinite(y)
            curves[b].append(float(spearmanr(x[ok], y[ok])[0]))
    nadir = coreg.cloudcube_to_halo_grid(cc, cl_pad, win, fields=("zm_w",), assign="footprint", lidar_pointing="nadir")
    out = {}
    for b in BANDS_KM:
        r = np.array(curves[b])
        i = int(np.nanargmax(r))
        vtx = LAGS_S[i]
        if 0 < i < r.size - 1:
            y0, y1, y2 = r[i - 1:i + 2]
            vtx = LAGS_S[i] + 0.5 * (LAGS_S[1] - LAGS_S[0]) * (y0 - y2) / (y0 - 2 * y1 + y2)
        m = (zkm >= b[0]) & (zkm < b[1])
        out[f"{b[0]}-{b[1]} km"] = {
            "curve": dict(zip(map(float, LAGS_S), map(float, r))), "peak_lag_s": float(vtx), "peak_r": float(r[i]),
            "n_ice_pixels": int(ice[:, m].sum()),
            "depth_below_aircraft_m": float(np.nanmedian(cl.gps_alt.values) - zkm[m].mean() * 1e3),
            "nadir_lidar_expected_lag_s": -float(np.nanmedian(nadir.footprint_dt_mean.values[:, m]))}
    return out


def lag_figure(lag, out):
    import ultraplot as uplt
    blue, orange, gray = "#1c5cab", "#eb6834", "#898781"
    fig, axs = uplt.subplots(ncols=2, refwidth=3.2, refaspect=1.2, share=False)
    ramp = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]
    for color, (band, res) in zip(ramp, lag.items()):
        c = res["curve"]
        axs[0].plot([float(k) for k in c], np.array(list(c.values())) - res["peak_r"], lw=1.4, color=color, label=band)
    axs[0].axvline(0, color=gray, lw=0.8)
    axs[0].format(xlabel="radar time − lidar time (s)", ylabel=r"Spearman($\ln\beta_p$, $Z_W$) − peak",
                  title="ice pixels, by altitude band")
    axs[0].set_xticks(np.arange(-3, 1.01, 1.0))
    axs[0].legend(loc="lr", ncols=1, frame=False)
    d = np.array([r["depth_below_aircraft_m"] for r in lag.values()]) / 1e3
    axs[1].plot(d, [r["peak_lag_s"] for r in lag.values()], color=blue, marker="o", label="observed peak")
    axs[1].plot(d, [r["nadir_lidar_expected_lag_s"] for r in lag.values()], color=orange, marker="s", ls="--",
                label="expected, nadir lidar")
    axs[1].plot(d, np.zeros_like(d), color=gray, lw=0.8, label="expected, body-fixed lidar")
    axs[1].format(xlabel="depth below aircraft (km)", ylabel="lag (s)", title="peak lag vs depth",
                  xlim=(d.min() - 0.3, d.max() + 0.3), ylim=(-2.4, 0.2))
    axs[1].legend(loc="b", ncols=2, frame=False)
    fig.format(suptitle="HALO 0.5 s lidar vs CloudCube W, 2026-02-17 15:30–16:00")
    fig.savefig(out, dpi=170)


def offset_angle_check(path):
    s = read_halo_cloud(path, group="surface", variables=["offset_angle", "pitch", "roll"])
    oa = s["offset_angle"].values
    model = np.cos(np.deg2rad(s["pitch"].values)) * np.cos(np.deg2rad(s["roll"].values))
    ok = np.isfinite(oa) & np.isfinite(model)
    turn = ok & (np.abs(s["roll"].values) > 10)
    return {"n": int(ok.sum()), "max_abs(offset_angle - cos(pitch)cos(roll))": float(np.abs(oa[ok] - model[ok]).max()),
            "n_turn_profiles_roll>10deg": int(turn.sum()),
            "acos_deg": pct(np.rad2deg(np.arccos(np.clip(oa[ok], -1, 1))), (1, 50, 99))}


def main():
    path = example_file("halo_cloud")
    res = {"file": path.name, "size_bytes": path.stat().st_size}
    dump_header(path, REPORTS / "headers" / "20260217_cloud_data.txt")

    with h5py.File(path, "r") as f, h5py.File(example_file("halo_cloudsurface"), "r") as c, \
            h5py.File(example_file("halo_aerowv"), "r") as h:
        dt_s = (f["data/nav/gps_time"][:, 0] - c["Nav/gps_time"][:, 0]) * 3600
        res["time_vs_cs"] = {"n_profiles": int(dt_s.size), "n_differ": int((dt_s != 0).sum()),
                             "max_abs_s": float(np.abs(dt_s).max())}
        res["nav_identical_to_cs"] = {k: bool(np.array_equal(f[f"data/nav/{k}"][...], c[f"Nav/{k}"][...], equal_nan=True))
                                      for k in ("gps_alt", "gps_lat", "gps_lon")}
        res["cloudphase_alt_equals_aerowv_z"] = bool(np.array_equal(f["data/cloudphase/alt"][0], h["z"][:]))
        res["cprofile_alt_equals_cs_wcp_altitude"] = bool(np.allclose(
            f["data/cprofile/alt"][:, 0], c["DataProducts/WaterCloudProfile/Altitude"][:, 0]))
        res["cloud_IAB_equals_cs_IAB532"] = bool(np.array_equal(f["data/cloud/IAB"][...], c["DataProducts/Cloud/IAB532"][...],
                                                                equal_nan=True))
        same, nan_same, n = 0, 0, 0
        for i in range(0, dt_s.size, 4000):
            a, b = f["data/cloudphase/cp"][i:i + 4000], c["DataProducts/CloudMask/Type"][i:i + 4000]
            both = np.isfinite(a) & np.isfinite(b)
            same += int((a[both] == b[both]).sum())
            n += int(both.sum())
            nan_same += int((np.isfinite(a) == np.isfinite(b)).sum())
        res["cp_vs_cs_Type_whole_flight"] = {"pixels_both_finite": n, "equal": same,
                                             "nan_pattern_equal_frac": nan_same / f["data/cloudphase/cp"].size}
    res["cloudphase_fig_labels"] = fig_strings(path.parent / "20260217_F1_cloudphase.fig")

    t0, t1 = SEGMENT
    cl = read_halo_cloud(path, time_range=(t0 - PAD, t1 + PAD))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        halo = read_halo_aerowv(example_file("halo_aerowv"), groups=("DataProducts",),
                                variables=["532_bsc", "1064_bsc", "532_dep", "532_aer_dep", "1064_dep",
                                           "1064_aer_dep", "532_bsr", "Aerosol_ID"], time_range=SEGMENT)
    halo = halo.sel(time=(halo.time - coreg.HALF_WINDOW >= t0) & (halo.time + coreg.HALF_WINDOW <= t1))
    res["units_10s"] = unit_checks(cl, halo)

    seg = cl.sel(profile=(cl.time >= t0) & (cl.time < t1))
    res["undocumented_cloudphase"] = {v: pct(seg[v].values) for v in ("S0", "cT2p", "X1064m_syn")}
    sf = read_halo_cloud(path, group="surface", variables=["offset_angle", "pitch", "roll"], time_range=SEGMENT)
    res["undocumented_surface"] = {v: pct(sf[v].values) for v in ("offset_angle", "pitch", "roll")}
    res["offset_angle_check"] = offset_angle_check(path)
    res["dep_532_unscreened"] = {"frac_abs_gt_1": float(np.mean(np.abs(seg.dep.values[np.isfinite(seg.dep.values)]) > 1))}

    fp = xr.open_dataset(cache_dir() / "coreg_cloud_20260217_1530-1600_0p5s.nc")
    zkm = fp.z.values / 1e3
    res["radar_gate_offsets_from_nadir"] = {}
    for zk in (1, 2, 4, 6, 8):
        j = int(np.argmin(np.abs(zkm - zk)))
        res["radar_gate_offsets_from_nadir"][f"{zk} km"] = {
            "along_m": pct(fp.along_track_offset_mean.values[:, j], (5, 50, 95)),
            "cross_m": pct(fp.cross_track_offset_mean.values[:, j], (5, 50, 95)),
            "n_rays": pct(fp.n_rays.values[:, j], (5, 50, 95))}

    win = coreg.subprofile_windows(cl.time.values)
    win = win.isel(time=(win.t_lo.values >= t0) & (win.t_hi.values <= t1))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        cc = read_cloudcube(example_file("cloudcube_1521"), time_range=(t0 - PAD, t1 + PAD), bands=("w",), profiles=("zm",))
    lag = lag_diagnostic(cl, win, cc)
    res["lag_diagnostic"] = lag
    lag_figure(lag, REPORTS / "figures" / "check_halo_cloud_lag.png")

    (REPORTS / "check_halo_cloud.json").write_text(json.dumps(res, indent=1))
    print(json.dumps({k: res[k] for k in ("time_vs_cs", "cp_vs_cs_Type_whole_flight", "cloudphase_fig_labels")}, indent=1))
    print(res["offset_angle_check"])
    for band, r in lag.items():
        print(band, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items() if k != "curve"})


if __name__ == "__main__":
    sys.exit(main())
