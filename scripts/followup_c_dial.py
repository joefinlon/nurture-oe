"""Follow-up §A.8(5c): in-ice DIAL q_v and the h2o_mmr_v ice fill fraction [OPEN O32].

1. Reconciles the ``h2o_mmr_v`` fill fraction in ice (0.932 on 2026-10-14 vs 0.97 recorded
   2026-10-06) by trying windows and pixel definitions.
2. In C&S ice pixels (classes 1-2 in >= 50 % of sub-profiles, > 2.5 km below the aircraft),
   computes RH_i from ``h2o_mmr_v`` with ``/State`` T and p (atm) and compares it with
   saturation and with MERRA-2 ``/State/Relative_Humidity``, by depth below the cloud-layer top
   and by relative ``h2o_mmr_v_unc``.

Writes ``docs/reports/followup_c_dial.json`` and ``figures/followup_c_dial.png``.
"""

import json
import sys
import warnings

import numpy as np

from nurture_oe import coreg
from nurture_oe.config import REPO_ROOT, example_file
from nurture_oe.io.halo_aerowv import read_halo_aerowv
from nurture_oe.io.halo_cloudsurface import read_halo_cloudsurface
from nurture_oe.thermo import esat_ice, esat_liquid, rh_ice

VARS = ["h2o_mmr_v", "h2o_mmr_v_unc", "range_cell_v", "RHi_v", "Aerosol_ID", "532_total_attn_bsc",
        "Temperature", "Pressure", "Relative_Humidity", "gps_alt"]
DEPTH_BINS = [0, 150, 300, 600, 1000, 2000, 4000, 12000]
UNC_BINS = [0, 0.05, 0.10, 0.20, 0.50, np.inf]
PCT = [5, 25, 50, 75, 95]


def dist(x):
    x = np.asarray(x)[np.isfinite(x)]
    if x.size == 0:
        return {"n": 0}
    return {"n": int(x.size), "pct_5_25_50_75_95": np.percentile(x, PCT).round(1).tolist(),
            "frac_lt_90": round(float((x < 90).mean()), 3), "frac_90_110": round(float(((x >= 90) & (x <= 110)).mean()), 3),
            "frac_gt_110": round(float((x > 110).mean()), 3)}


def layer_top_depth(cloud, z):
    """Depth (m) of each cloudy bin below the top of its contiguous cloud layer; NaN elsewhere."""
    depth = np.full(cloud.shape, np.nan)
    dz = np.median(np.diff(z))
    for i in range(cloud.shape[0]):
        top = None
        for j in range(cloud.shape[1] - 1, -1, -1):   # z ascending: scan downward
            if cloud[i, j]:
                if top is None:
                    top = z[j] + dz / 2
                depth[i, j] = top - z[j]
            else:
                top = None
    return depth


def main():
    warnings.simplefilter("ignore")
    a = read_halo_aerowv(example_file("halo_aerowv"), variables=VARS)
    cs = read_halo_cloudsurface(example_file("halo_cloudsurface"), variables=["Type", "altSurface"])
    z = a.z.values
    below = (a.gps_alt.values[:, None] - z[None, :]) > 2500
    h = np.isfinite(a.h2o_mmr_v.values)
    out = {"definition": "ice = C&S classes 1-2 in >= 50 % of sub-profiles of [t-5 s, t+5 s), > 2.5 km below aircraft"}

    # 1. fill-fraction reconciliation
    rec = {}
    for name, shift in (("[t-5,t+5)", 0), ("[t,t+10)", 5), ("[t-10,t)", -5)):
        w = coreg.halo_windows(a.time.values + np.timedelta64(int(shift * 1e9), "ns"))
        f = coreg.cloudmask_fractions(cs, w)
        ice_w = coreg.phase_masks(f).ice.values & below
        rec[f"C&S ice, window {name}"] = (int(ice_w.sum()), round(float(h[ice_w].mean()), 3))
        rec[f"C&S ice class 2 only, window {name}"] = (int((ice_w & (f.frac_2.values >= .5)).sum()),
                                                       round(float(h[(f.frac_2.values >= .5) & below].mean()), 3))
        rec[f"C&S ice in all 20 sub-profiles, window {name}"] = (
            int(((f.frac_1 + f.frac_2).values >= .999).sum()), round(float(h[((f.frac_1 + f.frac_2).values >= .999) & below].mean()), 3))
    win = coreg.halo_windows(a.time.values)
    fr = coreg.cloudmask_fractions(cs, win)
    pm = coreg.phase_masks(fr)
    ice = pm.ice.values & below
    k, n = coreg.window_sums(np.isfinite(cs["altSurface"].values)[:, None].astype(float), cs["time"].values, win)
    sfc_seen = (k[:, 0] > 0)[:, None]
    for nm, m in (("C&S ice, lidar surface seen in >= 1 sub-profile", ice & sfc_seen),
                  ("C&S ice, lidar surface never seen", ice & ~sfc_seen),
                  ("C&S ice & 532_total_attn_bsc finite", ice & np.isfinite(a["532_total_attn_bsc"].values)),
                  ("HALO Aerosol_ID == 1 (ice), > 2.5 km below", (a.Aerosol_ID.values == 1) & below),
                  ("C&S ice & Aerosol_ID == 1", ice & (a.Aerosol_ID.values == 1))):
        rec[nm] = (int(m.sum()), round(float(h[m].mean()), 3))
    hour = (a.time.values.astype("datetime64[h]") - a.time.values.astype("datetime64[D]")).astype(int)
    out["fill_by_hour_UTC"] = {int(hh): (int((ice & (hour == hh)[:, None]).sum()), round(float(h[ice & (hour == hh)[:, None]].mean()), 3))
                               for hh in np.unique(hour) if (ice & (hour == hh)[:, None]).any()}
    out["fill_reconciliation (n_pixels, fill)"] = rec

    # 2. RH_i in ice
    t, p = a.Temperature.values, a.Pressure.values
    rhi = rh_ice(a.h2o_mmr_v.values, p, t)
    rhm = a.Relative_Humidity.values
    rhm_as_ice = rhm                                   # MERRA-2 RH read as over ice (see report)
    rhm_w_to_ice = rhm * esat_liquid(t) / esat_ice(t)  # alternative: MERRA-2 RH over water
    rel_unc = a.h2o_mmr_v_unc.values / a.h2o_mmr_v.values
    depth = layer_top_depth(pm.cloud.values, z)
    m = ice & np.isfinite(rhi)
    pure = m & (fr.frac_0.values == 0)
    out["check_vs_product_RHi_v_ratio_pct_1_50_99"] = np.nanpercentile((rhi / a.RHi_v.values)[m], [1, 50, 99]).round(4).tolist()
    out["T_K_in_ice_pct_5_50_95"] = np.percentile(t[m], [5, 50, 95]).round(1).tolist()
    out["dlnEsi_dT_per_K_at_median_T"] = round(float(np.log(esat_ice(np.median(t[m]) + 0.5) / esat_ice(np.median(t[m]) - 0.5))), 3)
    out["MERRA2_RH_max_by_T"] = {f"{lo}-{lo + 10} K": [round(float(np.nanmax(rhm[(t >= lo) & (t < lo + 10)])), 1),
                                                      round(float(100 * esat_ice(lo + 5) / esat_liquid(lo + 5)), 1)]
                                 for lo in range(220, 270, 10)}
    out["RHi_DIAL_all_ice"] = dist(rhi[m])
    rhw = rhi * esat_ice(t) / esat_liquid(t)       # same vapour pressure, over supercooled water
    out["RHw_DIAL_all_ice"] = dist(rhw[m])
    out["RHi_at_water_saturation_pct_5_50_95"] = np.percentile((100 * esat_liquid(t) / esat_ice(t))[m], [5, 50, 95]).round(1).tolist()
    out["frac_RHw_gt_105"] = round(float((rhw[m] > 105).mean()), 3)
    out["RHi_DIAL_ice_no_liquid_subprofiles"] = dist(rhi[pure])
    out["RHi_MERRA2_as_ice_same_pixels"] = dist(rhm_as_ice[m])
    out["RHi_MERRA2_if_over_water_same_pixels"] = dist(rhm_w_to_ice[m])
    d = rhi - rhm_as_ice
    out["DIAL_minus_MERRA2_pct_points"] = {"n": int(np.isfinite(d[m]).sum()), "pct_5_25_50_75_95": np.nanpercentile(d[m], PCT).round(1).tolist(),
                                           "corr": round(float(np.corrcoef(rhi[m & np.isfinite(rhm)], rhm[m & np.isfinite(rhm)])[0, 1]), 3)}
    out["by_depth_below_layer_top_m"] = {}
    for lo, hi in zip(DEPTH_BINS[:-1], DEPTH_BINS[1:]):
        mm = m & (depth >= lo) & (depth < hi)
        out["by_depth_below_layer_top_m"][f"{lo}-{hi}"] = {
            "DIAL": dist(rhi[mm]), "MERRA2": dist(rhm_as_ice[mm]), "DIAL_RHw": dist(rhw[mm]),
            "fill_fraction": round(float(h[ice & (depth >= lo) & (depth < hi)].mean()), 3),
            "median_range_cell_v_m": round(float(np.nanmedian(a.range_cell_v.values[mm])), 0) if mm.any() else None}
    out["by_relative_unc"] = {}
    for lo, hi in zip(UNC_BINS[:-1], UNC_BINS[1:]):
        mm = m & (rel_unc >= lo) & (rel_unc < hi)
        out["by_relative_unc"][f"{lo:.2f}-{hi:.2f}"] = {"DIAL": dist(rhi[mm]), "MERRA2": dist(rhm_as_ice[mm])}
    # control: the same comparison in clear air (C&S class 4 in >= 90 % of sub-profiles), T < 0 °C
    clear = (fr.frac_4.values >= 0.9) & below & (t < 273.15) & np.isfinite(rhi) & np.isfinite(rhm)
    above_ice = np.zeros_like(clear)
    nbin = int(round(500 / np.median(np.diff(z))))
    for s_ in range(1, nbin + 1):
        above_ice[:, s_:] |= ice[:, :-s_]
    out["control_clear_air"] = {
        "all_clear_T<0C": {"DIAL": dist(rhi[clear]), "MERRA2": dist(rhm[clear]),
                           "DIAL_minus_MERRA2_pct_5_50_95": np.percentile((rhi - rhm)[clear], [5, 50, 95]).round(1).tolist()},
        "clear_within_500m_above_ice": {"DIAL": dist(rhi[clear & above_ice]), "MERRA2": dist(rhm[clear & above_ice]),
                                        "DIAL_minus_MERRA2_pct_5_50_95": np.percentile((rhi - rhm)[clear & above_ice], [5, 50, 95]).round(1).tolist()}}
    out["rel_unc_in_ice_pct_5_50_95"] = np.nanpercentile(rel_unc[m], [5, 50, 95]).round(3).tolist()
    (REPO_ROOT / "docs" / "reports" / "followup_c_dial.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))
    figure(rhi, rhm_as_ice, rhm_w_to_ice, depth, rel_unc, m)


def figure(rhi, rhm, rhw, depth, rel_unc, m):
    import ultraplot as uplt
    fig, axs = uplt.subplots(ncols=3, refwidth=3.0, refaspect=1.1, share=False)
    bins = np.arange(0, 201, 5)
    ax = axs[0]
    ax.hist(rhi[m], bins=bins, color="#1c5cab", alpha=0.6, label="DIAL h2o_mmr_v")
    ax.hist(rhm[m], bins=bins, color="#eb6834", alpha=0.5, label="MERRA-2 RH (as ice)")
    ax.hist(rhw[m], bins=bins, histtype="step", color="#7a2c10", label="MERRA-2 RH if over water")
    ax.axvline(100, color="#0b0b0b", lw=0.8)
    ax.format(xlabel="RH$_i$ (%)", ylabel="ice pixels", title="RH$_i$ in C&S ice")
    ax.legend(loc="ur", ncols=1, fontsize=7)
    for ax, x, edges, lab in ((axs[1], depth, DEPTH_BINS, "depth below layer top (m)"),
                              (axs[2], rel_unc, [0, 0.05, 0.1, 0.2, 0.5, 1.0], "relative h2o_mmr_v_unc")):
        ctr, q = [], []
        for lo, hi in zip(edges[:-1], edges[1:]):
            mm = m & (x >= lo) & (x < hi)
            if mm.sum() < 20:
                continue
            ctr.append(0.5 * (lo + hi))
            q.append([np.percentile(rhi[mm], PCT), np.nanpercentile(rhm[mm], PCT)])
        q = np.array(q)
        for j, c, nm in ((0, "#1c5cab", "DIAL"), (1, "#eb6834", "MERRA-2")):
            ax.fill_between(ctr, q[:, j, 1], q[:, j, 3], color=c, alpha=0.25)
            ax.plot(ctr, q[:, j, 2], color=c, marker="o", ms=3, label=f"{nm} median, IQR")
        ax.axhline(100, color="#0b0b0b", lw=0.8)
        ax.format(xlabel=lab, ylabel="RH$_i$ (%)", xscale="log" if "depth" in lab else "linear")
        ax.legend(loc="ur", ncols=1, fontsize=7)
    axs.format(suptitle="HALO DIAL RH$_i$ in lidar-classified ice, 2026-02-17 (T, p: MERRA-2 /State)")
    p = REPO_ROOT / "docs" / "reports" / "figures" / "followup_c_dial.png"
    fig.savefig(p, dpi=150)
    print(f"wrote {p}")


if __name__ == "__main__":
    sys.exit(main())
