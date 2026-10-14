"""Task (a) data-level checks on the 2026-02-17 example flight.

Prints results and writes ``docs/reports/task_a_checks.json`` (numbers quoted in
``docs/reports/task_a_readers.md``). Read-only on the data store.
"""

import json
import sys
import warnings

import numpy as np

from nurture_oe import coreg
from nurture_oe.checks.ka_coverage import zmin_dbz
from nurture_oe.config import REPO_ROOT, example_file
from nurture_oe.io.cloudcube import read_cloudcube
from nurture_oe.io.halo_aerowv import read_halo_aerowv
from nurture_oe.io.halo_cloudsurface import read_halo_cloudsurface
from nurture_oe.io.metnav import read_metnav

# Variables named in appendix A §A.10.2 (01-24 header dump, summarized). Not a full list.
PK_0124_NAMED = {
    "532_total_attn_bsc", "1064_total_attn_bsc", "532_bsc", "532_bsr", "532_bsc_cloud_screened",
    "532_bsr_cloud_screened", "1064_bsc", "1064_bsr", "1064_bsc_cloud_screened", "1064_bsr_cloud_screened",
    "532_ext", "532_Sa", "1064_ext", "532_extinction_from_backscatter", "532_AOT_hi", "532_AOT_hi_above_cloud",
    "532_AOT_hi_col", "532_dep", "1064_dep", "532_dep_cloud_screened", "1064_dep_cloud_screened",
    "532_aer_dep", "1064_aer_dep", "532_aer_dep_cloud_screened", "1064_aer_dep_cloud_screened",
    "WVD_1064_532", "Angstrom_Spherical", "Aerosol_ID", "cloud_top_height", "cloud_top_height_935",
    "h2o_mmr_v", "RHi_v", "RHw_v", "PWV_v", "PWV_bot_altitude_v", "PWV_top_altitude_v", "range_cell_v",
    "Temperature", "Pressure", "Number_Density", "Relative_Humidity", "U", "V", "imu_pitch", "imu_roll",
    "TrueHeading", "DEM_altitude", "532_ext_unc", "532_Sa_unc", "h2o_mmr_v_unc", "WVD_1064_532_unc",
    "532_bsc_cloud_screened_unc", "532_dep_cloud_screened_unc", "1064_dep_cloud_screened_unc",
}


def fnum(x, nd=4):
    return None if x is None or not np.isfinite(x) else round(float(x), nd)


def aerowv_checks(a, cs):
    out = {}
    names = set(a.data_vars)
    out["n_profiles"], out["n_z"] = a.sizes["time"], a.sizes["z"]
    out["z_ascending"] = bool(np.all(np.diff(a.z.values) > 0))
    out["z_first_last_step"] = [fnum(a.z.values[0], 2), fnum(a.z.values[-1], 2), fnum(np.median(np.diff(a.z.values)), 4)]
    out["vars_not_named_in_PK_A10_2"] = sorted(names - PK_0124_NAMED)
    out["PK_A10_2_named_not_in_file"] = sorted(PK_0124_NAMED - names)
    out["State_Type"] = float(a["State_Type"])
    out["State_units"] = {v: a[v].attrs.get("units") for v in
                          ("Temperature", "Pressure", "Number_Density", "Relative_Humidity", "U", "V")}
    out["State_Pressure_range"] = [fnum(np.nanmin(a.Pressure), 4), fnum(np.nanmax(a.Pressure), 4)]
    out["time_spacing_s_pct_0_5_50_95_100"] = np.percentile(np.diff(a.time_native.values), [0, 5, 50, 95, 100]).round(2).tolist()

    # cloud masks from C&S (appendix A §A.10.5 rule) and PK regressions
    win = coreg.halo_windows(a.time.values)
    fr = coreg.cloudmask_fractions(cs, win)
    pm = coreg.phase_masks(fr)
    below = (a.gps_alt.values[:, None] - a.z.values[None, :]) > 2500
    liq, ice = pm.liquid.values & below, pm.ice.values & below
    out["n_sub_pct_0_1_5_50_100"] = np.percentile(fr.n_sub.values, [0, 1, 5, 50, 100]).tolist()
    out["n_liquid_gt2p5km_below"], out["n_ice_gt2p5km_below"] = int(liq.sum()), int(ice.sum())
    fill = {}
    for v in ("532_total_attn_bsc", "532_dep", "532_bsc", "532_bsr", "1064_total_attn_bsc", "1064_dep",
              "532_aer_dep", "1064_aer_dep", "532_bsc_cloud_screened", "532_ext", "532_Sa",
              "532_extinction_from_backscatter", "Aerosol_ID", "h2o_mmr_v"):
        x = np.isfinite(a[v].values)
        fill[v] = [fnum(x[liq].mean(), 3), fnum(x[ice].mean(), 3)]
    out["fill_fraction_liquid_ice"] = fill
    pairs = [("532_bsc_cloud_screened", "532_bsc_cloud_screened_unc"), ("532_dep_cloud_screened", "532_dep_cloud_screened_unc"),
             ("1064_bsc_cloud_screened", "1064_bsc_cloud_screened_unc"), ("1064_dep_cloud_screened", "1064_dep_cloud_screened_unc"),
             ("532_aer_dep_cloud_screened", "532_aer_dep_cloud_screened_unc"), ("1064_aer_dep_cloud_screened", "1064_aer_dep_cloud_screened_unc")]
    out["unc_mask_equals_screened_mask"] = {u: bool(np.array_equal(np.isfinite(a[p].values), np.isfinite(a[u].values)))
                                            for p, u in pairs}
    cloud = pm.cloud.values
    out["screened_unc_finite_in_cloud"] = {u: fnum(np.isfinite(a[u].values)[cloud].mean(), 3) for _, u in pairs}
    out["aerosol_id_eq1_vs_cs_ice"] = {
        "P(C&S ice | Aerosol_ID==1)": fnum(pm.ice.values[a.Aerosol_ID.values == 1].mean(), 3),
        "P(Aerosol_ID==1 | C&S ice, >2.5 km below)": fnum((a.Aerosol_ID.values == 1)[ice].mean(), 3),
        "Aerosol_ID classes present": sorted(int(c) for c in np.unique(a.Aerosol_ID.values[np.isfinite(a.Aerosol_ID.values)])),
    }

    # NaN census by region
    dist_below = a.gps_alt.values[:, None] - a.z.values[None, :]
    sfc = a.DEM_altitude.values[:, None]
    regions = {
        "above_aircraft": dist_below < 0,
        "0-0.5km_below": (dist_below >= 0) & (dist_below < 500),
        "0.5-2.5km_below": (dist_below >= 500) & (dist_below < 2500),
        "below_DEM": (a.z.values[None, :] < sfc) & (dist_below >= 2500),
        "rest": (dist_below >= 2500) & (a.z.values[None, :] >= sfc),
    }
    census = {}
    for v in a.data_vars:
        if a[v].dims != ("time", "z"):
            continue
        x = np.isfinite(a[v].values)
        census[v] = {k: fnum(x[m].mean(), 3) for k, m in regions.items()}
    out["finite_fraction_by_region"] = census

    # h2o_mmr_v truncation relative to cloud_top_height_935
    h = np.isfinite(a.h2o_mmr_v.values)
    cth = a.cloud_top_height_935.values
    gaps = []
    for i in np.where(np.isfinite(cth))[0]:
        col = h[i] & (a.z.values > cth[i])
        if col.any():
            gaps.append(a.z.values[col].min() - cth[i])
    out["h2o_lowest_valid_above_cth935_m_pct_5_50_95"] = np.percentile(gaps, [5, 50, 95]).round(1).tolist() if gaps else None
    out["h2o_n_profiles_with_cth935"] = int(np.isfinite(cth).sum())
    return out


def attitude_checks(a, m):
    """Are HALO imu_pitch / imu_roll samples at t or 10 s means? Compare with 1 Hz MetNav. [G21], [O28]"""
    out = {}
    tm = m.time.values
    win = coreg.halo_windows(a.time.values)

    def at_t(x):
        ok = np.isfinite(x)
        return np.interp(a.time.values.astype(np.int64), tm[ok].astype(np.int64), x[ok])

    def mean10(x):
        ok = np.isfinite(x)
        s, _ = coreg.window_sums(np.nan_to_num(x)[:, None], tm, win)
        k, _ = coreg.window_sums(ok[:, None].astype(float), tm, win)
        with np.errstate(invalid="ignore", divide="ignore"):
            return (s / k)[:, 0]

    turning = np.abs(at_t(m["Roll_Angle"].values)) > 5
    out["n_profiles_turning(|roll|>5deg)"] = int(turning.sum())
    for hv, mv in (("imu_pitch", "Pitch_Angle"), ("imu_roll", "Roll_Angle")):
        x, h = m[mv].values, a[hv].values
        res = {}
        for name, ref in (("sampled_at_t", at_t(x)), ("mean_[t-5,t+5)", mean10(x))):
            d = h - ref
            res[name] = {"rms_all": fnum(np.sqrt(np.nanmean(d ** 2)), 3),
                         "rms_turning": fnum(np.sqrt(np.nanmean(d[turning] ** 2)), 3),
                         "median_diff": fnum(np.nanmedian(d), 3)}
        out[hv] = res
        out[f"{hv}_pct_1_50_99"] = np.nanpercentile(h, [1, 50, 99]).round(2).tolist()
    return out


def cloudcube_checks(path, a, cs):
    out = {}
    cc = read_cloudcube(path, profiles=())
    r = cc.range.values
    out["file"] = path.name
    out["n_time"], out["n_range"] = cc.sizes["time"], cc.sizes["range"]
    out["dt_s"] = np.percentile(np.diff(cc.time_native.values), [0, 50, 100]).round(5).tolist()
    out["gate_spacing_m"] = fnum(np.median(np.diff(r)), 4)
    out["range_first_last"] = [fnum(r[0], 2), fnum(r[-1], 2)]
    out["zmin_1km_dbz"] = {b: cc.attrs[f"zmin_1km_dbz_{b}"] for b in ("ka", "w", "g")}
    tr = cc.time_native.values
    ta = a.time_native.values
    sel = (ta > tr[0]) & (ta < tr[-1])
    out["halo_gps_alt_minus_cc_alt_m_pct_1_50_99"] = np.nanpercentile(
        a.gps_alt.values[sel] - np.interp(ta[sel], tr, cc["alt"].values), [1, 50, 99]).round(1).tolist()
    cth = coreg.cloudcube_cos_offnadir(cc).values
    th = np.degrees(np.arccos(np.clip(cth, -1, 1)))
    out["offnadir_deg_pct_1_50_99"] = np.nanpercentile(th, [1, 50, 99]).round(2).tolist()
    out["abs(offnadir - incidence_angle)_deg_pct_50_99"] = np.nanpercentile(
        np.abs(th - cc["incidence_angle"].values), [50, 99]).round(3).tolist()
    out["landmask_counts"] = {int(k): int(v) for k, v in zip(*np.unique(cc["landmask"].values, return_counts=True))}
    return out


def cloudcube_profile_checks(path, a, cs, segment):
    """Threshold floor vs Z_min(r), katten shape, surface altitude: on a time slice."""
    out = {"segment": list(segment)}
    cc = read_cloudcube(path, time_range=segment)
    r = cc.range.values
    for b in ("ka", "w", "g"):
        z = cc[f"zm_{b}"].values
        z1 = cc.attrs[f"zmin_1km_dbz_{b}"]
        res = {"finite_fraction": fnum(np.isfinite(z).mean(), 4)}
        for lo, hi in ((500, 3000), (3000, 6000), (6000, 9000), (9000, 11500)):
            mm = (r >= lo) & (r < hi)
            zz = z[:, mm]
            ok = np.isfinite(zz)
            if ok.sum() < 100:
                continue
            d = zz[ok] - zmin_dbz(np.broadcast_to(r[mm], zz.shape)[ok], z1)
            res[f"zm-Zmin_attr p0.1/p1/p50 @ {lo/1e3:.1f}-{hi/1e3:.1f} km"] = np.percentile(d, [0.1, 1, 50]).round(1).tolist()
        k = cc[f"katten_{b}"].values
        res["katten_finite_fraction"] = fnum(np.isfinite(k).mean(), 3)
        res["katten_dB_per_km_pct_1_50_99"] = np.nanpercentile(k, [1, 50, 99]).round(4).tolist()
        res["katten_finite_at_r<0"] = fnum(np.isfinite(k[:, r < 0]).mean(), 3)
        out[b] = res
    # surface peak altitude vs lidar surface
    cth = coreg.cloudcube_cos_offnadir(cc).values
    alt_h = coreg.aircraft_alt_from_halo(cc, a).values
    zw = np.where(r[None, :] < 5000, -np.inf, np.nan_to_num(cc.zm_w.values, nan=-np.inf))
    rp = r[np.argmax(zw, axis=1)]
    alts = np.interp(cc.time.values.astype(np.int64), cs.time.values.astype(np.int64), cs.altSurface.values)
    out["W_surface_peak_minus_CS_altSurface_m_pct_5_25_50_75_95"] = {
        "HALO MSL reference": np.nanpercentile(alt_h - rp * cth - alts, [5, 25, 50, 75, 95]).round(1).tolist(),
        "CloudCube ellipsoid alt": np.nanpercentile(cc["alt"].values - rp * cth - alts, [5, 25, 50, 75, 95]).round(1).tolist(),
    }
    return out


def metnav_checks(m):
    dt = np.diff(m.time_native.values)
    return {"n": m.sizes["time"], "t_first_last": [str(m.time.values[0]), str(m.time.values[-1])],
            "dt_s_unique_counts": {str(k): int(v) for k, v in zip(*np.unique(dt, return_counts=True))},
            "variables": {v: m[v].attrs["units"] for v in m.data_vars if not v.endswith(("_llod", "_ulod"))},
            "missing_fraction": {v: fnum(np.isnan(m[v].values).mean(), 3) for v in m.data_vars
                                 if not v.endswith(("_llod", "_ulod"))},
            "llod_ulod_flag_vars": [v for v in m.data_vars if v.endswith(("_llod", "_ulod"))]}


def main():
    results = {}
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        a = read_halo_aerowv(example_file("halo_aerowv"))
        results["reader_warnings"] = sorted({str(x.message) for x in w})
    cs = read_halo_cloudsurface(example_file("halo_cloudsurface"), variables=["Type", "altSurface"])
    m = read_metnav(example_file("metnav"))
    results["halo_aerowv"] = aerowv_checks(a, cs)
    results["halo_cloudsurface"] = {"n_profiles": cs.sizes["profile"], "n_z": cs.sizes["z"],
                                    "z_identical_to_aerowv": bool(np.array_equal(cs.z.values, a.z.values)),
                                    "Type_classes": {str(int(k)): int(v) for k, v in zip(
                                        *np.unique(cs.Type.values[np.isfinite(cs.Type.values)], return_counts=True))},
                                    "dt_s_pct_0_1_50_99_100": np.percentile(np.diff(cs.time_native.values * 3600),
                                                                            [0, 1, 50, 99, 100]).round(3).tolist()}
    results["attitude_vs_metnav"] = attitude_checks(a, m)
    results["metnav"] = metnav_checks(m)
    results["cloudcube"] = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for name in ("cloudcube_1513", "cloudcube_1521", "cloudcube_1900"):
            results["cloudcube"][name] = cloudcube_checks(example_file(name), a, cs)
        results["cloudcube_profiles_1521"] = cloudcube_profile_checks(
            example_file("cloudcube_1521"), a, cs, ("2026-02-17T15:30:00", "2026-02-17T16:00:00"))
        results["cloudcube_profiles_1900"] = cloudcube_profile_checks(
            example_file("cloudcube_1900"), a, cs, ("2026-02-17T19:00:00", "2026-02-17T19:25:00"))
    out = REPO_ROOT / "docs" / "reports" / "task_a_checks.json"
    out.write_text(json.dumps(results, indent=1, default=str))
    print(json.dumps(results, indent=1, default=str))


if __name__ == "__main__":
    sys.exit(main())
