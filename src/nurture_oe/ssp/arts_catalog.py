"""Catalog of the ARTS microwave single-scattering database v1.1 (``SSD/``), task (c).

Groundwork for ``[OPEN O27]`` (index §5.2, §10 S1). Read-only, plain netCDF4: each size file
holds one group per (frequency, temperature) with ``SingleScatteringData``, ``ShapeData`` and
``CalculationData``. In ``SSD/`` a missing (frequency, temperature) group means the property
was not computed for that size (DDA size-parameter limit); nothing is extrapolated there
(the extensions are in ``StandardHabits/``, which is not used here).

Reference values ``README_TABLE`` are transcribed from the database ReadMe (Ekelund et al.
2021, Tables 1-2): Id -> (name, Dmax range um, Dveq range um, n sizes, alpha, beta, software).
The ReadMe notes that alpha and beta "are calculated for the upper end of d_max".
"""

import re
from pathlib import Path

import netCDF4
import numpy as np
import pandas as pd

GROUP_RE = re.compile(r"Freq(\d+\.\d+)GHz_T(\d+\.\d+)K")
SIZE_RE = re.compile(r"Dmax(\d+)um_Dveq(\d+)um_Mass([\d.eE+-]+)kg\.nc$")

# Ka-band grid frequencies (CloudCube's exact Ka frequency is unknown [GAP G1]), W, and the two
# grid frequencies bracketing CloudCube G band (~239 GHz, exact value unknown [GAP G1]).
TARGET_FREQS_GHZ = (31.3, 31.5, 35.6, 94.1, 228.0, 247.2)

README_TABLE = {
    9: ("Plate type 1", (13, 10000), (10, 2596), 45, 0.76, 2.48, "RSP"),
    7: ("Column type 1", (14, 10000), (10, 1815), 45, 0.037, 2.05, "RSP"),
    16: ("Thin Plate", (25, 5059), (10, 2000), 35, 30, 3.00, "RSP"),
    15: ("Thick Plate", (16, 3246), (10, 2000), 35, 110, 3.00, "RSP"),
    12: ("Block Column", (13, 2632), (10, 2000), 35, 210, 3.00, "RSP"),
    13: ("Short Column", (17, 3303), (10, 2000), 34, 110, 3.00, "RSP"),
    14: ("Long Column", (24, 4835), (10, 2000), 35, 34, 3.00, "RSP"),
    3: ("Sector Snowflake", (20, 12000), (20, 1415), 34, 0.00081, 1.44, "RSP"),
    24: ("Ice Sphere (Mie)", (1, 50000), (1, 50000), 200, 480, 3.00, "Mie"),
    27: ("ICON Cloud Ice", (13, 10000), (10, 2929), 45, 1.6, 2.56, "RC"),
    31: ("GEM Cloud Ice", (10, 3088), (10, 3000), 45, 440, 3.00, "RC"),
    6: ("6-Bullet Rosette", (16, 10000), (10, 2371), 45, 0.48, 2.42, "RSP"),
    2: ("5-Bullet Rosette", (17, 10000), (10, 2231), 45, 0.4, 2.43, "SFTK"),
    10: ("Perpendicular 4-Bullet Rosette", (18, 10000), (10, 2071), 45, 0.32, 2.43, "SFTK"),
    11: ("Flat 4-Bullet Rosette", (18, 10000), (10, 2071), 45, 0.32, 2.43, "SFTK"),
    4: ("Perpendicular 3-Bullet Rosette", (19, 10000), (10, 2137), 45, 0.44, 2.47, "SFTK"),
    5: ("Flat 3-Bullet Rosette", (20, 10000), (10, 1882), 45, 0.2, 2.43, "SFTK"),
    36: ("Brown-Francis 6-Column Rosette", (172, 10000), (139, 1829), 32, 0.019, 1.9, "SFTK"),
    1: ("Evans Snow Aggregate", (32, 11755), (50, 2506), 35, 0.20, 2.39, "ESP"),
    26: ("Tyynela Dendrite Aggregate", (595, 20826), (228, 3328), 35, 0.10, 2.25, "ESP"),
    8: ("8-Column Aggregate", (19, 9714), (10, 5000), 39, 65, 3.00, "RSP"),
    17: ("Small Column Aggregate", (105, 3855), (37, 738), 35, 0.14, 2.45, "SFTK"),
    18: ("Large Column Aggregate", (368, 19981), (128, 3021), 35, 0.25, 2.43, "SFTK"),
    21: ("Small Block Aggregate", (100, 7328), (72, 1665), 35, 0.21, 2.33, "SFTK"),
    22: ("Large Block Aggregate", (349, 21875), (253, 4607), 35, 0.35, 2.27, "SFTK"),
    19: ("Small Plate Aggregate", (99, 7054), (53, 1376), 35, 0.077, 2.25, "SFTK"),
    20: ("Large Plate Aggregate", (349, 22860), (197, 4563), 34, 0.21, 2.26, "SFTK"),
    30: ("ICON Hail", (120, 5349), (94, 5000), 35, 380, 2.99, "RC"),
    28: ("ICON Snow", (120, 20000), (94, 3219), 35, 0.031, 1.95, "RC"),
    34: ("GEM Hail", (120, 5031), (94, 5000), 35, 540, 3.02, "RC"),
    32: ("GEM Snow", (170, 10459), (94, 5000), 35, 24, 2.86, "RC"),
    23: ("Spherical Graupel", (622, 9744), (454, 5293), 30, 13, 2.69, "SFTK"),
    29: ("ICON Graupel", (170, 6658), (94, 5000), 35, 390, 3.13, "RC"),
    33: ("GEM Graupel", (120, 6597), (94, 5000), 35, 170, 2.96, "RC"),
    25: ("Liquid Sphere (Mie)", (1, 50000), (1, 50000), 200, 480, 3.00, "Mie"),
}

# TAMUice2016 habit folders whose *names* match ARTS habits. Flag only: geometric equivalence is
# unverified and is an [OPEN O27] question for a later session.
TAMU_NAME_MATCH = {
    9: "plate", 16: "plate", 15: "plate",
    7: "single_column", 12: "single_column", 13: "single_column", 14: "single_column",
    2: "SBR (bullet rosette)", 4: "SBR (bullet rosette)", 5: "SBR (bullet rosette)", 6: "SBR (bullet rosette)",
    10: "SBR (bullet rosette)", 11: "SBR (bullet rosette)",
    8: "8_columns",
    19: "5_plates / 10_plates (plate aggregates)", 20: "5_plates / 10_plates (plate aggregates)",
}


def habit_dirs(root):
    """Habit folders under ``<root>/SSD/TotallyRandom`` as (phase, class, riming, path), sorted by Id."""
    base = Path(root) / "SSD" / "TotallyRandom"
    out = []
    for p in sorted(base.glob("*/*/*/*_Id*")):
        phase, klass, riming = p.relative_to(base).parts[:3]
        out.append((phase, klass, riming, p))
    return sorted(out, key=lambda x: habit_id(x[3]))


def habit_id(path):
    return int(Path(path).name.rsplit("_Id", 1)[1])


def read_summary(habit_dir):
    """``DataSummary.txt`` as a dict (numbers converted where possible)."""
    out = {}
    for line in (Path(habit_dir) / "DataSummary.txt").read_text().splitlines():
        if "=" not in line:
            continue
        k, v = line.split("=", 1)
        try:
            out[k.strip()] = float(v)
        except ValueError:
            out[k.strip()] = v.strip()
    return out


def size_table(habit_dir):
    """One row per size file: shape data and the (frequency, temperature) groups present.

    Returns
    -------
    pandas.DataFrame
        ``file``, ``dmax_m``, ``dveq_m``, ``mass_kg``, ``d_aeq_m`` (``diameter_area_eq_aerodynamical``),
        ``orient_type``, ``freqs_ghz`` (sorted tuple), ``temps_k`` (sorted tuple), ``groups`` (set of
        (f, T)), ``n_freq``, ``fmax_ghz``.
    """
    rows = []
    for f in sorted((Path(habit_dir) / "TotallyRandom").glob("*.nc")):
        with netCDF4.Dataset(f) as nc:
            groups = {(float(m.group(1)), float(m.group(2))) for m in map(GROUP_RE.fullmatch, nc.groups) if m}
            g0 = nc.groups[next(iter(nc.groups))]
            sd = g0["ShapeData"]
            val = {k: float(sd[k][...]) if k in sd.variables else np.nan
                   for k in ("diameter_max", "diameter_vol_eq", "mass", "diameter_area_eq_aerodynamical")}
            orient = g0["SingleScatteringData"].getncattr("orient_type")
        freqs = tuple(sorted({g[0] for g in groups}))
        rows.append({"file": f.name, "dmax_m": val["diameter_max"], "dveq_m": val["diameter_vol_eq"],
                     "mass_kg": val["mass"], "d_aeq_m": val["diameter_area_eq_aerodynamical"],
                     "orient_type": orient, "freqs_ghz": freqs, "temps_k": tuple(sorted({g[1] for g in groups})),
                     "groups": groups, "n_freq": len(freqs), "fmax_ghz": max(freqs) if freqs else np.nan})
    return pd.DataFrame(rows).sort_values("dveq_m", ignore_index=True)


def fit_power(d, y, select=None):
    """Least-squares fit of ``y = a d**b`` in log space; returns (a, b, rms of log10 residual)."""
    d, y = np.asarray(d, float), np.asarray(y, float)
    ok = (d > 0) & (y > 0) & np.isfinite(d) & np.isfinite(y)
    if select is not None:
        ok &= select
    if ok.sum() < 2:
        return np.nan, np.nan, np.nan
    b, la = np.polyfit(np.log10(d[ok]), np.log10(y[ok]), 1)
    res = np.log10(y[ok]) - (la + b * np.log10(d[ok]))
    return 10 ** la, b, float(np.sqrt(np.mean(res ** 2)))


def largest_with(sizes, f_ghz, temps=None):
    """Largest Dveq (and its Dmax) with computed data at ``f_ghz`` (all of ``temps`` if given)."""
    def has(groups):
        if temps is None:
            return any(abs(g[0] - f_ghz) < 1e-3 for g in groups)
        return all(any(abs(g[0] - f_ghz) < 1e-3 and abs(g[1] - t) < 1e-3 for g in groups) for t in temps)
    m = sizes["groups"].map(has).to_numpy()
    if not m.any():
        return np.nan, np.nan, 0
    i = np.flatnonzero(m)[-1]
    return sizes["dveq_m"].iloc[i], sizes["dmax_m"].iloc[i], int(m.sum())


def catalog_row(phase, klass, riming, habit_dir, sizes, summary, target_freqs=TARGET_FREQS_GHZ):
    hid = habit_id(habit_dir)
    ref = README_TABLE.get(hid)
    all_t = sorted({t for g in sizes["groups"] for (_, t) in g})
    all_f = sorted({f for g in sizes["groups"] for (f, _) in g})
    # does every size have every temperature at each of its frequencies?
    t_complete = all(len({t for (f2, t) in g if f2 == f}) == len(all_t) for g in sizes["groups"] for f in {f for f, _ in g})
    a_full, b_full, r_full = fit_power(sizes.dmax_m, sizes.mass_kg)
    upper = sizes.dmax_m.to_numpy() >= np.sqrt(sizes.dmax_m.min() * sizes.dmax_m.max())   # upper half in log Dmax
    a_up, b_up, r_up = fit_power(sizes.dmax_m, sizes.mass_kg, upper)
    area = np.pi / 4 * sizes.d_aeq_m.to_numpy() ** 2
    a_area, b_area, r_area = fit_power(sizes.dmax_m, area)
    row = {
        "habit_id": hid, "habit_folder": Path(habit_dir).name, "readme_name": ref[0] if ref else "",
        "phase": phase, "class": f"{klass}/{riming}", "orientation": ";".join(sorted(set(sizes.orient_type))),
        "status": summary.get("STATUS", ""), "n_sizes": len(sizes), "readme_n_sizes": ref[3] if ref else np.nan,
        "dmax_min_um": sizes.dmax_m.min() * 1e6, "dmax_max_um": sizes.dmax_m.max() * 1e6,
        "dveq_min_um": sizes.dveq_m.min() * 1e6, "dveq_max_um": sizes.dveq_m.max() * 1e6,
        "mass_alpha_fit_all": a_full, "mass_beta_fit_all": b_full, "mass_fit_rms_log10_all": r_full,
        "mass_alpha_fit_upper": a_up, "mass_beta_fit_upper": b_up,
        "mass_alpha_summary": summary.get("ALPHA", np.nan), "mass_beta_summary": summary.get("BETA", np.nan),
        "mass_alpha_readme": ref[4] if ref else np.nan, "mass_beta_readme": ref[5] if ref else np.nan,
        "area_from_d_aeq_present": bool(np.isfinite(sizes.d_aeq_m).any()),
        "area_gamma_fit": a_area, "area_sigma_fit": b_area, "area_fit_rms_log10": r_area,
        "n_freq_grid": len(all_f), "freq_min_ghz": min(all_f), "freq_max_ghz": max(all_f),
        "n_freq_per_size_min": int(sizes.n_freq.min()), "n_freq_per_size_max": int(sizes.n_freq.max()),
        "n_sizes_full_freq_grid": int((sizes.n_freq == len(all_f)).sum()),
        "temps_k": ";".join(f"{t:g}" for t in all_t), "all_temps_at_every_freq_size": t_complete,
        "tamu2016_name_match (flag only, O27)": TAMU_NAME_MATCH.get(hid, ""),
        "software": ref[6] if ref else "",
    }
    for f in target_freqs:
        dv, dm, n = largest_with(sizes, f)
        dv_all, _, _ = largest_with(sizes, f, all_t)
        row[f"dveq_max_um_at_{f:g}GHz"] = dv * 1e6
        row[f"dmax_um_at_{f:g}GHz"] = dm * 1e6
        row[f"n_sizes_at_{f:g}GHz"] = n
        row[f"dveq_max_um_at_{f:g}GHz_all_T"] = dv_all * 1e6
    both = sizes["groups"].map(lambda g: any(abs(x[0] - 228.0) < 1e-3 for x in g) and any(abs(x[0] - 247.2) < 1e-3 for x in g))
    row["dveq_max_um_with_228_and_247.2GHz"] = sizes.dveq_m[both].max() * 1e6 if both.any() else np.nan
    return row


def build_catalog(root, target_freqs=TARGET_FREQS_GHZ):
    """Habit catalog (one row per habit) and per-size coverage table for ``<root>/SSD/TotallyRandom``."""
    rows, size_rows = [], []
    for phase, klass, riming, d in habit_dirs(root):
        sizes = size_table(d)
        summary = read_summary(d)
        rows.append(catalog_row(phase, klass, riming, d, sizes, summary, target_freqs))
        s = sizes.drop(columns="groups").assign(habit_id=habit_id(d), habit_folder=d.name)
        for f in target_freqs:
            s[f"has_{f:g}GHz"] = sizes["groups"].map(lambda g, f=f: any(abs(x[0] - f) < 1e-3 for x in g))
        s["freqs_ghz"] = s["freqs_ghz"].map(lambda t: ";".join(f"{x:g}" for x in t))
        s["temps_k"] = s["temps_k"].map(lambda t: ";".join(f"{x:g}" for x in t))
        size_rows.append(s)
    return pd.DataFrame(rows), pd.concat(size_rows, ignore_index=True)
