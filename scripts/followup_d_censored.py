"""Follow-up §A.8(5d): censored-mean effect in co-registered radar rows [OPEN O21].

On 2026-02-17 15:30-16:00 (10 s x 15 m grid of task (d)):

1. Descriptive: detected fraction f of valid pixels per band. ``<band>_all`` in the cache is the
   detecting-ray mean *without* the f >= 0.5 threshold, so it equals ``<band>`` wherever that is
   valid; the non-detecting rays carry no value.
2. Bounds on real partial pixels: the all-ray linear mean lies between f*m_det (non-detecting
   rays at Z = 0) and f*m_det + (1 - f)*Z_floor(r) (non-detecting rays at the floor). Bias of the
   censored mean = 10 log10(m_det / m_true). Floor: Z_min(r) from the granule attribute plus the
   2026-10-14 offset of the delivered floor (0.1th percentile, 6-9 km, 15:30-16:00) [GAP G24].
3. Synthetic censoring: in pixels where every ray detects (true mean known), censor rays below
   Z_min(r) + delta for delta = -12 ... +12 dB, recompute the detecting-ray mean, and bin the bias
   by the resulting detected fraction. Same for DWR_W-G with both bands censored at the same delta.

Writes ``docs/reports/followup_d_censored.json``, ``followup_d_censored_*.csv`` and
``figures/followup_d_censored.png``.
"""

import json
import sys
import warnings

import numpy as np
import pandas as pd

from nurture_oe import coreg
from nurture_oe.checks.ka_coverage import zmin_dbz
from nurture_oe.config import REPO_ROOT, cache_dir, example_file
from nurture_oe.io.cloudcube import read_cloudcube
from nurture_oe.io.halo_aerowv import read_halo_aerowv

SEGMENT = ("2026-02-17T15:30:00", "2026-02-17T16:00:00")
PAD = np.timedelta64(6, "s")
BANDS = ("ka", "w", "g")
FLOOR_OFFSET_DB = {"ka": -9.0, "w": -10.9, "g": -7.0}   # task (a), 6-9 km, 15:30-16:00 (appendix A §A.4)
DELTAS = np.arange(-12.0, 12.1, 1.5)
F_BINS = np.round(np.arange(0.0, 1.0001, 0.1), 2)
PCT = [5, 25, 50, 75, 95]


def rays_on_grid(cc, halo, win):
    """Per-ray dBZ on the HALO z grid (nearest gate, HALO MSL datum), plus slant range."""
    alt = coreg.aircraft_alt_from_halo(cc, halo).values
    cth = coreg.cloudcube_cos_offnadir(cc).values
    idx, r_need = coreg.gate_index_on_grid(cc, halo.z.values, alt, cth)
    valid = idx >= 0
    safe = np.where(valid, idx, 0)
    out = {}
    for b in BANDS:
        v = np.take_along_axis(cc[f"zm_{b}"].values, safe, axis=1)
        out[b] = np.where(valid, v, np.nan)
    return out, np.where(valid, r_need, np.nan), valid


def wsum(x, t, win):
    return coreg.window_sums(x, t, win)[0]


def binned(f, x, label):
    rows = []
    for lo, hi in zip(F_BINS[:-1], F_BINS[1:]):
        m = (f > lo) & (f <= hi) & np.isfinite(x)
        if m.sum() < 30:
            continue
        rows.append({"f_bin": f"({lo:.1f}, {hi:.1f}]", "n": int(m.sum()),
                     **{f"{label}_p{p}": round(float(np.percentile(x[m], p)), 2) for p in PCT}})
    return pd.DataFrame(rows)


def main():
    warnings.simplefilter("ignore")
    t0, t1 = (np.datetime64(s, "ns") for s in SEGMENT)
    halo = read_halo_aerowv(example_file("halo_aerowv"), variables=["gps_alt"], time_range=SEGMENT)
    halo = halo.sel(time=(halo.time - coreg.HALF_WINDOW >= t0) & (halo.time + coreg.HALF_WINDOW <= t1))
    win = coreg.halo_windows(halo.time.values)
    cc = read_cloudcube(example_file("cloudcube_1521"), time_range=(t0 - PAD, t1 + PAD))
    t = cc.time.values
    ray, r, valid = rays_on_grid(cc, halo, win)
    nvalid = wsum(valid.astype(float), t, win)
    rmean = wsum(np.nan_to_num(r), t, win) / np.where(nvalid > 0, nvalid, np.nan)
    ds = __import__("xarray").open_dataset(cache_dir() / "coreg_20260217_1530-1600.nc")
    out = {"segment": SEGMENT, "n_profiles": int(halo.sizes["time"]),
           "note_all": "<band>_all = detecting-ray mean without the f >= 0.5 threshold (not an all-ray mean)"}

    # 1. descriptive + 2. bounds
    desc, bounds = {}, []
    lin = {}
    for b in BANDS:
        det = np.isfinite(ray[b])
        nd = wsum(det.astype(float), t, win)
        s = wsum(np.where(det, 10 ** (ray[b] / 10), 0.0), t, win)
        with np.errstate(invalid="ignore", divide="ignore"):
            f = nd / nvalid
            mdet = s / nd
        lin[b] = (f, mdet)
        same = float(np.nanmax(np.abs(10 * np.log10(mdet[f >= 0.5]) - ds[f"zm_{b}"].values[f >= 0.5])))
        v = f >= 0.5
        desc[b] = {"max_abs_diff_vs_2026-10-14_cache_dB": round(same, 3),
                   "n_valid(f>=0.5)": int(v.sum()), "frac_valid_with_f<1": round(float((v & (f < 1)).sum() / v.sum()), 3),
                   "frac_valid_with_f<0.9": round(float((v & (f < 0.9)).sum() / v.sum()), 3),
                   "n_partial(0<f<0.5)": int(((f > 0) & (f < 0.5)).sum()),
                   "zm_minus_zm_all_where_valid_max_abs_dB": float(np.nanmax(np.abs(ds[f"zm_{b}"] - ds[f"zm_{b}_all"]).values))}
        floor = 10 ** ((zmin_dbz(rmean, ds.attrs[f"zmin_1km_dbz_{b}"]) + FLOOR_OFFSET_DB[b]) / 10)
        with np.errstate(invalid="ignore", divide="ignore"):
            bias_hi = -10 * np.log10(f)                                    # non-detects at Z = 0
            bias_lo = 10 * np.log10(mdet / (f * mdet + (1 - f) * floor))  # non-detects at the floor
        df = binned(f, bias_hi, "bias_upper_dB").merge(binned(f, bias_lo, "bias_lower_dB"), on=["f_bin", "n"], how="outer")
        df.insert(0, "band", b)
        bounds.append(df)
        lin[b] = (f, mdet, bias_lo, bias_hi)
    out["descriptive"] = desc
    bounds = pd.concat(bounds, ignore_index=True)
    # DWR_W-G bounds on real pixels (both bands valid, at least one partial)
    fw, _, low_w, up_w = lin["w"]
    fg, _, low_g, up_g = lin["g"]
    both = (fw >= 0.5) & (fg >= 0.5) & ((fw < 1) | (fg < 1))
    fmin_real = np.minimum(fw, fg)
    dwr_bounds = binned(fmin_real[both], (low_w - up_g)[both], "dwr_bias_min_dB").merge(
        binned(fmin_real[both], (up_w - low_g)[both], "dwr_bias_max_dB"), on=["f_bin", "n"])
    out["n_WG_valid"] = int(((fw >= 0.5) & (fg >= 0.5)).sum())
    out["n_WG_valid_with_a_partial_band"] = int(both.sum())
    out["dwr_bounds_real_pixels_by_min_f"] = dwr_bounds.to_dict(orient="records")
    print(dwr_bounds.to_string(index=False))

    # 3. synthetic censoring on fully detected pixels
    syn, syn_dwr = [], []
    true = {}
    for b in BANDS:
        f, mdet = lin[b][:2]
        true[b] = (f >= 0.999) & (nvalid > 0)
    zr = {b: zmin_dbz(r, cc.attrs[f"zmin_1km_dbz_{b}"]) for b in BANDS}
    cens = {}
    for d in DELTAS:
        for b in BANDS:
            keep = np.isfinite(ray[b]) & (ray[b] > zr[b] + d)
            nd = wsum(keep.astype(float), t, win)
            s = wsum(np.where(keep, 10 ** (ray[b] / 10), 0.0), t, win)
            with np.errstate(invalid="ignore", divide="ignore"):
                fp = nd / nvalid
                bias = 10 * np.log10((s / nd) / lin[b][1])
            m = true[b] & (fp > 0) & (fp < 1)
            syn.append(pd.DataFrame({"band": b, "delta_dB": d, "f": fp[m], "bias_dB": bias[m]}))
            cens[b] = (fp, bias)
        m = true["w"] & true["g"] & (cens["w"][0] > 0) & (cens["g"][0] > 0) & ((cens["w"][0] < 1) | (cens["g"][0] < 1))
        syn_dwr.append(pd.DataFrame({"delta_dB": d, "f_w": cens["w"][0][m], "f_g": cens["g"][0][m],
                                     "dwr_bias_dB": (cens["w"][1] - cens["g"][1])[m]}))
    syn = pd.concat(syn, ignore_index=True)
    syn_dwr = pd.concat(syn_dwr, ignore_index=True)
    out["n_fully_detected_pixels"] = {b: int(true[b].sum()) for b in BANDS}
    tab = pd.concat([binned(g.f.values, g.bias_dB.values, "bias_dB").assign(band=b) for b, g in syn.groupby("band")])
    dw = syn_dwr
    fmin = np.minimum(dw.f_w.values, dw.f_g.values)
    tab_dwr = {"by_min_f": binned(fmin, dw.dwr_bias_dB.values, "dwr_bias_dB"),
               "by_f_w (f_g >= 0.9)": binned(dw.f_w.values[dw.f_g >= 0.9], dw.dwr_bias_dB.values[dw.f_g >= 0.9], "dwr_bias_dB"),
               "by_f_g (f_w >= 0.9)": binned(dw.f_g.values[dw.f_w >= 0.9], dw.dwr_bias_dB.values[dw.f_w >= 0.9], "dwr_bias_dB")}
    rep = REPO_ROOT / "docs" / "reports"
    bounds.to_csv(rep / "followup_d_censored_bounds.csv", index=False)
    tab.to_csv(rep / "followup_d_censored_synthetic.csv", index=False)
    pd.concat([v.assign(selection=k) for k, v in tab_dwr.items()]).to_csv(rep / "followup_d_censored_synthetic_dwr.csv", index=False)
    out["bounds_real_partial_pixels"] = bounds.to_dict(orient="records")
    out["synthetic"] = tab.to_dict(orient="records")
    out["synthetic_dwr"] = {k: v.to_dict(orient="records") for k, v in tab_dwr.items()}
    (rep / "followup_d_censored.json").write_text(json.dumps(out, indent=1, default=str))
    print(json.dumps(desc, indent=1))
    print(bounds.to_string(index=False))
    print(tab.to_string(index=False))
    for k, v in tab_dwr.items():
        print(k)
        print(v.to_string(index=False))
    figure(bounds, tab, tab_dwr)


def figure(bounds, tab, tab_dwr):
    import ultraplot as uplt
    colors = {"ka": "#1baf7a", "w": "#3987e5", "g": "#eb6834"}
    fig, axs = uplt.subplots(ncols=3, refwidth=3.0, refaspect=1.15, share=False)

    def ctr(s):
        lo, hi = s.strip("(]").split(",")
        return 0.5 * (float(lo) + float(hi))

    ax = axs[0]
    for b, g in bounds.groupby("band"):
        x = g.f_bin.map(ctr)
        ax.fill_between(x, g.bias_lower_dB_p50, g.bias_upper_dB_p50, color=colors[b], alpha=0.3, label="_nolegend_")
        ax.plot(x, g.bias_upper_dB_p50, color=colors[b], lw=1, label=f"{b.upper()}: median upper / lower bound")
        ax.plot(x, g.bias_lower_dB_p50, color=colors[b], lw=1, ls="--", label="_nolegend_")
    ax.format(xlabel="detected fraction f", ylabel="censored-mean bias (dB)", title="Real partial pixels: bounds")
    ax.legend(loc="ur", ncols=1, fontsize=7)
    ax = axs[1]
    for b, g in tab.groupby("band"):
        x = g.f_bin.map(ctr)
        ax.fill_between(x, g.bias_dB_p25, g.bias_dB_p75, color=colors[b], alpha=0.25, label="_nolegend_")
        ax.plot(x, g.bias_dB_p50, color=colors[b], marker="o", ms=3, label=b.upper())
    ax.format(xlabel="detected fraction f", ylabel="bias (dB)", title="Synthetic censoring (median, IQR)")
    ax.legend(loc="ur", ncols=1, fontsize=7)
    ax = axs[2]
    for (k, g), c in zip(tab_dwr.items(), ("#0b0b0b", "#3987e5", "#eb6834")):
        x = g.f_bin.map(ctr)
        ax.plot(x, g.dwr_bias_dB_p50, color=c, marker="o", ms=3, label=k)
        if k == "by_min_f":
            ax.fill_between(x, g.dwr_bias_dB_p25, g.dwr_bias_dB_p75, color=c, alpha=0.15, label="_nolegend_")
    ax.axhline(0, color="#898781", lw=0.8)
    ax.format(xlabel="detected fraction", ylabel="DWR$_{W-G}$ bias (dB)", title="Synthetic: DWR$_{W-G}$")
    ax.legend(loc="ur", ncols=1, fontsize=7)
    axs.format(suptitle="Censored (detecting-ray) window means, 2026-02-17 15:30–16:00, 10 s × 15 m")
    p = REPO_ROOT / "docs" / "reports" / "figures" / "followup_d_censored.png"
    fig.savefig(p, dpi=150)
    print(f"wrote {p}")


if __name__ == "__main__":
    sys.exit(main())
