"""Follow-up §A.8(5b): HALO 10 s window offset sweep and window bookkeeping [GAP G23].

Correlation of window-mean C&S ``IAB532`` with column-integrated HALO ``532_total_attn_bsc``
as a function of the window-centre offset (10 s window ``[t + off - 5, t + off + 5)``,
off = -10 ... +10 s in 0.5 s steps), on the 273 "heterogeneous" profiles of task (a).

Heterogeneous (task (a) definition, 2026-10-14): with C&S ``IAB532`` NaN -> 0, the std/mean of
``IAB532`` in [t - 10, t + 10) exceeds 0.5 and the mean exceeds 0.005.

Writes ``docs/reports/followup_b_window.json``, ``followup_b_window_sweep.csv`` and
``figures/followup_b_window.png``.
"""

import json
import sys

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

from nurture_oe.config import REPO_ROOT, example_file
from nurture_oe.io.halo_aerowv import read_halo_aerowv
from nurture_oe.io.halo_cloudsurface import read_halo_cloudsurface

OFFSETS = np.round(np.arange(-10.0, 10.0 + 1e-9, 0.5), 1)
HALF = 5.0
N_BOOT = 2000
BLOCK = 10
SEGMENT = (15.5 * 3600, 16.0 * 3600)


def load():
    a = read_halo_aerowv(example_file("halo_aerowv"), variables=["532_total_attn_bsc"])
    c = read_halo_cloudsurface(example_file("halo_cloudsurface"), variables=["IAB532"])
    ta = a.time_native.values.astype(float)
    tc = c.time_native.values.astype(float) * 3600.0
    dz = np.gradient(a.z.values) / 1000.0
    col = np.nansum(a["532_total_attn_bsc"].values * dz, axis=1)
    iab = np.nan_to_num(c["IAB532"].values)
    return ta, tc, col, iab


def heterogeneous(ta, tc, iab):
    cs = np.concatenate([[0], np.cumsum(iab)])
    cs2 = np.concatenate([[0], np.cumsum(iab ** 2)])
    i0, i1 = np.searchsorted(tc, ta - 10), np.searchsorted(tc, ta + 10)
    n = i1 - i0
    mu = (cs[i1] - cs[i0]) / n
    sd = np.sqrt(np.maximum((cs2[i1] - cs2[i0]) / n - mu ** 2, 0))
    return (sd / np.maximum(mu, 1e-4) > 0.5) & (mu > 0.005)


def window_mean(ta, tc, iab, off):
    cs = np.concatenate([[0], np.cumsum(iab)])
    j0 = np.searchsorted(tc, ta + off - HALF)
    j1 = np.searchsorted(tc, ta + off + HALF)
    return (cs[j1] - cs[j0]) / np.maximum(j1 - j0, 1)


def robust_std(x):
    q = np.percentile(x, [25, 75])
    return (q[1] - q[0]) / 1.349


def sweep(ta, tc, col, iab, het):
    m = np.stack([window_mean(ta, tc, iab, o) for o in OFFSETS])          # (n_off, n_prof)
    fixed = het & (col > 0.003) & np.all(m > 0.003, axis=0)
    rows = []
    for k, o in enumerate(OFFSETS):
        var = het & (m[k] > 0.003) & (col > 0.003)
        rows.append({"offset_s": o,
                     "corr_variable_set": np.corrcoef(m[k][var], col[var])[0, 1], "n_variable_set": int(var.sum()),
                     "corr_fixed_set": np.corrcoef(m[k][fixed], col[fixed])[0, 1],
                     "spearman_fixed_set": spearmanr(m[k][fixed], col[fixed])[0],
                     "robust_std_logratio_fixed_set": robust_std(np.log(m[k][fixed] / col[fixed]))})
    return pd.DataFrame(rows), m, fixed


def corr_rows(m, col, idx):
    x = m[:, idx]
    y = col[idx]
    xm = x - x.mean(axis=1, keepdims=True)
    ym = y - y.mean()
    return (xm @ ym) / np.sqrt((xm ** 2).sum(axis=1) * (ym ** 2).sum())


def bootstrap(m, col, fixed, rng):
    """Peak-offset distribution: i.i.d. and moving-block (``BLOCK`` consecutive profiles)."""
    idx = np.where(fixed)[0]
    n = idx.size
    peaks = {"iid": [], f"block{BLOCK}": [], "iid_spearman": []}
    rm = np.apply_along_axis(rankdata, 1, m[:, idx])
    curves = []
    for _ in range(N_BOOT):
        s = idx[rng.integers(0, n, n)]
        c = corr_rows(m, col, s)
        curves.append(c)
        peaks["iid"].append(OFFSETS[np.argmax(c)])
        j = rng.integers(0, n, n)
        rs = rm[:, j]
        rc = rankdata(col[idx][j])
        peaks["iid_spearman"].append(OFFSETS[np.argmax(corr_rows(rs, rc, np.arange(n)))])
        starts = rng.integers(0, n - BLOCK + 1, int(np.ceil(n / BLOCK)))
        sb = idx[(starts[:, None] + np.arange(BLOCK)[None, :]).ravel()[:n]]
        peaks[f"block{BLOCK}"].append(OFFSETS[np.argmax(corr_rows(m, col, sb))])
    return {k: np.array(v) for k, v in peaks.items()}, np.array(curves)


def widths(curve):
    """Full width at half drop (peak to the curve minimum over the sweep), in s."""
    k = int(np.argmax(curve))
    half = 0.5 * (curve[k] + curve.min())
    lo = k
    while lo > 0 and curve[lo - 1] >= half:
        lo -= 1
    hi = k
    while hi < curve.size - 1 and curve[hi + 1] >= half:
        hi += 1
    return OFFSETS[lo], OFFSETS[hi]


def bookkeeping(ta, tc, sel=None):
    """Sub-profile counts per [t-5, t+5) window and overlaps between consecutive windows."""
    if sel is not None:
        keep = (ta - HALF >= sel[0]) & (ta + HALF <= sel[1])
        ta = ta[keep]
    j0 = np.searchsorted(tc, ta - HALF)
    j1 = np.searchsorted(tc, ta + HALF)
    n = j1 - j0
    distinct = np.array([np.unique(tc[a:b]).size for a, b in zip(j0, j1)])
    sp = np.diff(ta)
    overlap = np.maximum(0, j1[:-1] - j0[1:])           # sub-profiles shared by consecutive windows
    lo, hi = (ta[0] - HALF, ta[-1] + HALF)
    in_span = (tc >= lo) & (tc < hi)
    cover = np.zeros(tc.size, int)
    for a, b in zip(j0, j1):
        cover[a:b] += 1
    vc = pd.Series(n).value_counts().sort_index()
    return {"n_profiles": int(ta.size),
            "n_sub_counts": {int(k): int(v) for k, v in vc.items()},
            "n_windows_not_20": int((n != 20).sum()),
            "n_windows_not_20_distinct": int((distinct != 20).sum()),
            "n_windows_with_duplicate_times": int((distinct != n).sum()),
            "spacing_s_pct_0_1_50_99_100": np.percentile(sp, [0, 1, 50, 99, 100]).round(3).tolist(),
            "n_consecutive_pairs_spacing_lt_10s": int((sp < 10 - 1e-6).sum()),
            "n_consecutive_pairs_sharing_subprofiles": int((overlap > 0).sum()),
            "n_subprofiles_shared(in two windows)": int((cover[in_span] == 2).sum()),
            "max_windows_per_subprofile": int(cover.max()),
            "n_subprofiles_in_no_window(within span)": int((cover[in_span] == 0).sum()),
            "n_subprofiles_in_span": int(in_span.sum())}


def main():
    ta, tc, col, iab = load()
    het = heterogeneous(ta, tc, iab)
    df, m, fixed = sweep(ta, tc, col, iab, het)
    rng = np.random.default_rng(20261020)
    peaks, curves = bootstrap(m, col, fixed, rng)
    df["corr_fixed_boot_p05"] = np.percentile(curves, 5, axis=0)
    df["corr_fixed_boot_p95"] = np.percentile(curves, 95, axis=0)
    rep = REPO_ROOT / "docs" / "reports"
    df.round(4).to_csv(rep / "followup_b_window_sweep.csv", index=False)
    cf, cv = df.corr_fixed_set.values, df.corr_variable_set.values
    out = {
        "n_heterogeneous": int(het.sum()), "n_fixed_set": int(fixed.sum()),
        "peak_offset_s_fixed_set": float(OFFSETS[np.argmax(cf)]), "peak_corr_fixed_set": round(float(cf.max()), 3),
        "peak_offset_s_variable_set": float(OFFSETS[np.argmax(cv)]), "peak_corr_variable_set": round(float(cv.max()), 3),
        "corr_at_0_variable_set": round(float(cv[OFFSETS == 0][0]), 3),
        "full_width_half_drop_s_fixed_set": list(map(float, widths(cf))),
        "curve_min_over_sweep_fixed_set": round(float(cf.min()), 3),
        "offsets_within_0.01_of_peak_fixed_set": OFFSETS[cf >= cf.max() - 0.01].tolist(),
"peak_offset_s_spearman_fixed_set": float(OFFSETS[np.argmax(df.spearman_fixed_set.values)]),
        "peak_spearman_fixed_set": round(float(df.spearman_fixed_set.max()), 3),
        "spearman_at_0_fixed_set": round(float(df.spearman_fixed_set.values[OFFSETS == 0][0]), 3),
        "full_width_half_drop_s_spearman": list(map(float, widths(df.spearman_fixed_set.values))),
        "min_robust_std_offset_s": float(OFFSETS[np.argmin(df.robust_std_logratio_fixed_set.values)]),
        "bootstrap": {k: {"n": N_BOOT, "peak_offset_pct_2.5_50_97.5": np.percentile(v, [2.5, 50, 97.5]).tolist(),
                          "fraction_peak_in_[-1,+1]": round(float(np.mean(np.abs(v) <= 1.0)), 3)}
                      for k, v in peaks.items()},
        "bookkeeping_flight": bookkeeping(ta, tc),
        "bookkeeping_1530_1600": bookkeeping(ta, tc, SEGMENT),
    }
    (rep / "followup_b_window.json").write_text(json.dumps(out, indent=1))
    print(df.round(3).to_string(index=False))
    print(json.dumps(out, indent=1))
    figure(df, peaks)


def figure(df, peaks):
    import ultraplot as uplt
    fig, axs = uplt.subplots(ncols=2, refwidth=3.6, refaspect=1.3, share=False)
    ax = axs[0]
    ax.fill_between(df.offset_s, df.corr_fixed_boot_p05, df.corr_fixed_boot_p95, color="#cde2fb", label="bootstrap 5–95 %")
    ax.plot(df.offset_s, df.corr_fixed_set, color="#1c5cab", lw=1.5, label=f"fixed set")
    ax.plot(df.offset_s, df.corr_variable_set, color="#eb6834", lw=1, ls="--", label="variable set (task a)")
    ax.plot(df.offset_s, df.spearman_fixed_set, color="#13805a", lw=1.2, label="Spearman, fixed set")
    ax.axvline(0, color="#898781", lw=0.8)
    ax.format(xlabel="window-centre offset (s)", ylabel="corr(windowed IAB532, column β'532)",
              title="Correlation vs offset, window [t+off−5, t+off+5)")
    ax.legend(loc="lc", ncols=1, fontsize=7)
    ax = axs[1]
    bins = np.arange(-10.25, 10.3, 0.5)
    for (k, v), c in zip(peaks.items(), ("#1c5cab", "#eb6834", "#13805a")):
        ax.hist(v, bins=bins, color=c, alpha=0.5, label=k)
    ax.format(xlabel="peak offset (s)", ylabel="bootstrap replicates", title="Peak location (bootstrap)")
    ax.legend(loc="ur", ncols=1, fontsize=7)
    p = REPO_ROOT / "docs" / "reports" / "figures" / "followup_b_window.png"
    fig.savefig(p, dpi=150)
    print(f"wrote {p}")


if __name__ == "__main__":
    sys.exit(main())
