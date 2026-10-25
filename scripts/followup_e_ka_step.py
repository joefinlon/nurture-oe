"""Follow-up §A.8(5e): Ka calibration vs sensitivity across the 18:59:37 / 19:00:00 boundary [GAP G24], [OPEN O15].

Native rays and gates within +-1, +-2, +-3 and +-6 min of the CloudCube granule boundary (no
co-registration; the three bands share the range grid). Step = median(after) - median(before);
90 % intervals from a block bootstrap over 10 s blocks of rays.

1. DWR_Ka-W (display gas correction with ``katten``) in Rayleigh-like ice: Ka, W, G detected,
   |DWR_W-G| <= 1 dB, Z_W below a cap, Ka above a COMMON threshold (the less sensitive
   granule's Z_min(r) + 3 dB), T < -5 °C (HALO /State), a clearance above the surface echo.
   Reported as "not testable" where either side has < 50 gates.
2. Surface echo per band: peak ``zm`` and the contiguous extent of detected gates below and
   above the peak (below the surface there is no atmosphere, so the extent below the peak
   measures the range response).
3. sigma0_Ka - sigma0_W and sigma0_G - sigma0_W (control) over land (``landmask`` 3),
   ``incidence_angle`` 1.8-2.0°, raw and with two-way ``katten`` gas correction to the surface.
   No liquid-free columns exist within +-6 min (C&S class 0 in every column), so all rays are used.

Writes ``docs/reports/followup_e_ka_step.json`` and ``figures/followup_e_ka_step.png``.
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

BOUNDARY = np.datetime64("2026-02-17T18:59:50", "ns")
SPANS_MIN = (1, 2, 3, 6)
BLOCK_S = 10.0
N_BOOT = 2000
DWR_WG_MAX = 1.0
T_MAX_K = 268.15
INC_RANGE = (1.80, 2.00)
CLEARANCES_M = (500.0, 1500.0)
ZW_MAX_DBZ = (-5.0, 100.0)


def load_side(name, t0, t1, halo):
    cc = read_cloudcube(example_file(name), time_range=(t0, t1))
    r = cc.range.values
    zgc = {b: coreg.gas_corrected_display(cc, b).values for b in ("ka", "w", "g")}
    # surface range per ray: W peak beyond 5 km
    zw = np.where(r[None, :] < 5000, -np.inf, np.nan_to_num(cc.zm_w.values, nan=-np.inf))
    isfc = np.argmax(zw, axis=1)
    r_sfc = r[isfc]
    alt = coreg.aircraft_alt_from_halo(cc, halo).values
    cth = coreg.cloudcube_cos_offnadir(cc).values
    zgate = alt[:, None] - r[None, :] * cth[:, None]
    # temperature at gate altitude from the nearest HALO profile
    k = np.clip(np.searchsorted(halo.time.values, cc.time.values), 0, halo.sizes["time"] - 1)
    tk = np.stack([np.interp(zgate[i], halo.z.values, halo.Temperature.values[k[i]], left=np.nan, right=np.nan)
                   for i in range(zgate.shape[0])])
    # two-way gas to the surface gate (dB), per band
    dr = np.median(np.diff(r)) / 1000.0
    gas_sfc = {}
    for b in ("ka", "w", "g"):
        kk = np.nan_to_num(cc[f"katten_{b}"].values) * (r[None, :] >= 0)
        cum = np.cumsum(kk, axis=1) * dr
        gas_sfc[b] = 2 * cum[np.arange(cum.shape[0]), isfc]
    return {"cc": cc, "r": r, "zgc": zgc, "r_sfc": r_sfc, "tk": tk, "gas_sfc": gas_sfc,
            "t": cc.time.values, "zmin_ka": cc.attrs["zmin_1km_dbz_ka"]}


def blocks(t):
    return ((t - t[0]) / np.timedelta64(int(BLOCK_S * 1e9), "ns")).astype(int)


def boot_step(xb, bb, xa, ba, rng):
    """Median(after) - median(before) with block-bootstrap 5-95 % interval."""
    ub, ua = np.unique(bb), np.unique(ba)
    gb = [xb[bb == u] for u in ub]
    ga = [xa[ba == u] for u in ua]
    steps = []
    for _ in range(N_BOOT):
        sb = np.concatenate([gb[i] for i in rng.integers(0, len(gb), len(gb))])
        sa = np.concatenate([ga[i] for i in rng.integers(0, len(ga), len(ga))])
        steps.append(np.median(sa) - np.median(sb))
    return {"median_before": round(float(np.median(xb)), 2), "median_after": round(float(np.median(xa)), 2),
            "step": round(float(np.median(xa) - np.median(xb)), 2),
            "step_boot_p5_p95": np.percentile(steps, [5, 95]).round(2).tolist(),
            "n_before": int(xb.size), "n_after": int(xa.size), "blocks_before": len(gb), "blocks_after": len(ga)}


def rayleigh_dwr(side, zthr_1km, margin, clearance_m, zw_max):
    """DWR_Ka-W at gates of Rayleigh-like ice: weak Z_W, |DWR_W-G| <= 1 dB, Ka above a common threshold."""
    zk, zw, zg = side["zgc"]["ka"], side["zgc"]["w"], side["zgc"]["g"]
    r = side["r"][None, :]
    with np.errstate(invalid="ignore", divide="ignore"):
        ok = (np.isfinite(zk) & np.isfinite(zw) & np.isfinite(zg) & (np.abs(zw - zg) <= DWR_WG_MAX) & (zw < zw_max)
              & (side["cc"].zm_ka.values > zmin_dbz(r, zthr_1km) + margin) & (side["tk"] < T_MAX_K)
              & (r > 500) & (r < side["r_sfc"][:, None] - clearance_m))
    rows = np.nonzero(ok)[0]
    return (zk - zw)[ok], blocks(side["t"])[rows]


def surface_echo(side, band):
    """Per ray: surface-peak zm, and contiguous echo extent above / below the peak (m)."""
    z = side["cc"][f"zm_{band}"].values
    r = side["r"]
    zz = np.where(r[None, :] < 5000, -np.inf, np.nan_to_num(z, nan=-np.inf))
    j = np.argmax(zz, axis=1)
    f = np.isfinite(z)
    up = np.full(z.shape[0], np.nan)
    dn = np.full(z.shape[0], np.nan)
    for i in range(z.shape[0]):
        if not f[i, j[i]]:
            continue
        k = j[i]
        while k > 0 and f[i, k - 1]:
            k -= 1
        up[i] = r[j[i]] - r[k]
        k = j[i]
        while k < r.size - 1 and f[i, k + 1]:
            k += 1
        dn[i] = r[k] - r[j[i]]
    return z[np.arange(z.shape[0]), j], up, dn


def sigma0_rays(side):
    cc = side["cc"]
    inc = cc.incidence_angle.values
    m = (cc.landmask.values == 3) & (inc >= INC_RANGE[0]) & (inc <= INC_RANGE[1])
    out = {}
    for a, b in (("ka", "w"), ("g", "w")):
        raw = cc[f"s0_{a}"].values - cc[f"s0_{b}"].values
        gc = raw + side["gas_sfc"][a] - side["gas_sfc"][b]
        ok = m & np.isfinite(raw)
        out[f"{a}-{b}"] = (raw[ok], gc[ok], blocks(side["t"])[ok])
    out["incidence_median"] = round(float(np.median(inc[m])), 2)
    return out


def main():
    warnings.simplefilter("ignore")
    rng = np.random.default_rng(20261020)
    halo = read_halo_aerowv(example_file("halo_aerowv"), variables=["gps_alt", "Temperature"],
                            time_range=("2026-02-17T18:20", "2026-02-17T19:27"))
    res = {"boundary": "18:59:37.9 (last ray, 15:21 granule) / 19:00:00.1 (first ray, 19:00 granule)",
           "selection": {"DWR_WG_abs_max_dB": DWR_WG_MAX, "T_max_K": T_MAX_K, "landmask": 3, "incidence_deg": INC_RANGE,
                         "block_s": BLOCK_S, "n_boot": N_BOOT}}
    series = {}
    for span in SPANS_MIN:
        d = np.timedelta64(span * 60, "s")
        before = load_side("cloudcube_1521", BOUNDARY - d, BOUNDARY, halo)
        after = load_side("cloudcube_1900", BOUNDARY, BOUNDARY + d, halo)
        key = f"+-{span} min"
        out = res[key] = {"rays_before_after": [int(before["t"].size), int(after["t"].size)]}
        zthr = max(before["zmin_ka"], after["zmin_ka"])
        for clr in CLEARANCES_M:
            for zw_max in ZW_MAX_DBZ:
                xb, bb = rayleigh_dwr(before, zthr, 3.0, clr, zw_max)
                xa, ba = rayleigh_dwr(after, zthr, 3.0, clr, zw_max)
                k2 = f"DWR_Ka-W Rayleigh-like ice (Z_W < {zw_max} dBZ, >= {clr} m above sfc, Ka > Zmin({zthr}) + 3 dB)"
                out[k2] = boot_step(xb, bb, xa, ba, rng) if min(xb.size, xa.size) >= 50 else {
                    "n_before": int(xb.size), "n_after": int(xa.size), "note": "not testable (< 50 gates on a side)"}
        for b in ("ka", "w", "g"):
            pb, ub, db = surface_echo(before, b)
            pa, ua, da = surface_echo(after, b)
            bb_, ba_ = blocks(before["t"]), blocks(after["t"])
            okb, oka = np.isfinite(pb) & np.isfinite(db), np.isfinite(pa) & np.isfinite(da)
            out[f"surface peak zm_{b} (dBZ)"] = boot_step(pb[okb], bb_[okb], pa[oka], ba_[oka], rng)
            out[f"surface echo extent below peak, {b} (m)"] = boot_step(db[okb], bb_[okb], da[oka], ba_[oka], rng)
            out[f"surface echo extent above peak, {b} (m)"] = boot_step(ub[okb], bb_[okb], ua[oka], ba_[oka], rng)
            if span == max(SPANS_MIN) and b == "ka":
                series["below"] = [(before["t"][okb], bb_[okb], db[okb]), (after["t"][oka], ba_[oka], da[oka])]
        pk = {}
        for sd, nm in ((before, "before"), (after, "after")):
            p_ka, _, _ = surface_echo(sd, "ka")
            p_w, _, _ = surface_echo(sd, "w")
            pk[nm] = p_ka - p_w
        okb, oka = np.isfinite(pk["before"]), np.isfinite(pk["after"])
        out["surface peak zm_ka - zm_w (dB)"] = boot_step(pk["before"][okb], blocks(before["t"])[okb], pk["after"][oka], blocks(after["t"])[oka], rng)
        sb, sa = sigma0_rays(before), sigma0_rays(after)
        out["incidence_median_before_after"] = [sb["incidence_median"], sa["incidence_median"]]
        for comb in ("ka-w", "g-w"):
            rb, gb_, bb = sb[comb]
            ra, ga_, ba = sa[comb]
            out[f"sigma0 {comb}, raw"] = boot_step(rb, bb, ra, ba, rng)
            out[f"sigma0 {comb}, gas-corrected"] = boot_step(gb_, bb, ga_, ba, rng)
            if span == max(SPANS_MIN):
                series[f"s0 {comb}"] = [(before["t"], bb, gb_), (after["t"], ba, ga_)]
    late = load_side("cloudcube_1900", np.datetime64("2026-02-17T19:20:00", "ns"), np.datetime64("2026-02-17T19:26:00", "ns"), halo)
    r = late["r"][None, :]
    high = (r > 500) & (r < late["r_sfc"][:, None] - 1500)
    res["19:20-19:26 check"] = {
        "ka_extent_below_surface_peak_m_median": float(np.nanmedian(surface_echo(late, "ka")[2])),
        "w_extent_below_surface_peak_m_median": float(np.nanmedian(surface_echo(late, "w")[2])),
        "detections_gt_1.5km_above_surface": {b: int((np.isfinite(late["cc"][f"zm_{b}"].values) & high).sum()) for b in ("ka", "w", "g")},
        "surface_peak_zm_median (attenuation check)": {b: round(float(np.nanmedian(surface_echo(late, b)[0])), 1) for b in ("ka", "w", "g")}}
    (REPO_ROOT / "docs" / "reports" / "followup_e_ka_step.json").write_text(json.dumps(res, indent=1, default=str))
    print(json.dumps(res, indent=1, default=str))
    figure(series)


def figure(series):
    import matplotlib.dates as mdates
    import ultraplot as uplt
    fig, axs = uplt.subplots(nrows=len(series), refwidth=6.5, refaspect=4.0, sharex=True, sharey=False)
    labels = {"below": "Ka surface echo: contiguous extent below the surface peak (m) — range response",
              "s0 ka-w": r"$\sigma_0^{Ka}-\sigma_0^{W}$ (dB), land, gas-corrected",
              "s0 g-w": r"$\sigma_0^{G}-\sigma_0^{W}$ (dB), land, gas-corrected (control)"}
    colors = ("#3987e5", "#eb6834")
    for ax, (k, parts) in zip(axs, series.items()):
        for (t, blk, x), c in zip(parts, colors):
            u = np.unique(blk)
            tc = np.array([t[0] + np.timedelta64(int((i + 0.5) * BLOCK_S * 1e9), "ns") for i in u])
            med = np.array([np.median(x[blk == i]) for i in u])
            ax.plot(tc, med, ".", color=c, ms=4)
            ax.hlines(np.median(x), tc[0], tc[-1], color=c, lw=1.5)
        ax.axvline(BOUNDARY, color="#0b0b0b", lw=0.8)
        ax.format(ylabel="m" if k == "below" else "dB", title=labels[k], titleloc="l")
    axs[-1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs.format(xlabel="UTC, 2026-02-17", suptitle="Ka across the 18:59:37 / 19:00:00 granule boundary (10 s block medians; lines: side medians)")
    p = REPO_ROOT / "docs" / "reports" / "figures" / "followup_e_ka_step.png"
    fig.savefig(p, dpi=150)
    print(f"wrote {p}")


if __name__ == "__main__":
    sys.exit(main())
