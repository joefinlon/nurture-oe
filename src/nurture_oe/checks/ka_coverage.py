"""Ka coverage check for ``[OPEN O15]`` condition 3 (appendix A §A.4(i), §A.8(4)).

Fraction of ice pixels with W and G detected where Ka also detects, on the co-registered
HALO grid. Ka "detects" where its window-mean Z exceeds

    Z_min^Ka(r) = Z_1km + 20 log10(r / 1 km)   [dBZ],

with r the slant range from the aircraft (from ``alt``). ``Z_1km`` comes from the granule's
``radar_sensitivity_ka_dbz`` attribute (file governs); PK's -19.26 is reported alongside.
"""

import numpy as np
import pandas as pd

from ..io.cloudcube import PK_ZMIN_1KM_DBZ


def zmin_dbz(r_m, zmin_1km_dbz):
    """Range-dependent minimum detectable reflectivity, ``Z_1km + 20 log10(r/1 km)``."""
    with np.errstate(divide="ignore", invalid="ignore"):
        return zmin_1km_dbz + 20.0 * np.log10(np.asarray(r_m) / 1000.0)


def ka_detects(grid, zmin_1km_dbz=None, min_detect_frac=0.5):
    """Boolean Ka detection on the co-registered grid.

    Parameters
    ----------
    grid : xarray.Dataset
        Output of ``coreg.cloudcube_to_halo_grid`` (needs ``zm_ka``, ``zm_ka_all``,
        ``zm_ka_detect_frac``, ``range_mean``).
    zmin_1km_dbz : float or None
        ``None`` -> "finite" definition only (detect fraction >= ``min_detect_frac``).
        Otherwise also require window-mean Z_Ka > Z_min^Ka(range_mean).
    """
    finite = grid["zm_ka_detect_frac"].values >= min_detect_frac
    if zmin_1km_dbz is None:
        return finite
    return finite & (grid["zm_ka_all"].values > zmin_dbz(grid["range_mean"].values, zmin_1km_dbz))


def ka_coverage(grid, ice, temperature_k=None, min_detect_frac=0.5, range_bins_km=(0, 2, 4, 6, 8, 10, 12, 14)):
    """Ka coverage tables.

    Parameters
    ----------
    grid : xarray.Dataset
        Co-registered CloudCube fields on the HALO grid.
    ice : array-like of bool (time, z)
        Ice-pixel mask (C&S classes 1-2 in >= 50 % of sub-profiles).
    temperature_k : array-like (time, z), optional
        HALO ``/State/Temperature`` (K) for the T < 0 °C subset (O15 condition 2).

    Returns
    -------
    summary : pandas.DataFrame
        One row per (definition, subset): n_denominator, n_ka, fraction.
    by_range : pandas.DataFrame
        Fraction per slant-range bin for each definition (all ice pixels).
    """
    wg = (grid["zm_w_detect_frac"].values >= min_detect_frac) & (grid["zm_g_detect_frac"].values >= min_detect_frac)
    denom = np.asarray(ice, bool) & wg & np.isfinite(grid["range_mean"].values)
    file_zmin = grid.attrs.get("zmin_1km_dbz_ka")
    defs = {f"Z>Zmin, file attr ({file_zmin:+.2f} dBZ)": file_zmin,
            f"Z>Zmin, PK ({PK_ZMIN_1KM_DBZ['ka']:+.2f} dBZ)": PK_ZMIN_1KM_DBZ["ka"],
            "finite zm_ka only": None}
    subsets = {"all ice, W&G detected": denom}
    if temperature_k is not None:
        subsets["... and T < 0 °C"] = denom & (np.asarray(temperature_k) < 273.15)
    rows, rng_rows = [], []
    r_km = grid["range_mean"].values / 1000.0
    for dname, z1 in defs.items():
        det = ka_detects(grid, z1, min_detect_frac)
        for sname, m in subsets.items():
            n, k = int(m.sum()), int((m & det).sum())
            rows.append({"definition": dname, "subset": sname, "n_ice_WG": n, "n_ka": k,
                         "fraction": k / n if n else np.nan})
        for lo, hi in zip(range_bins_km[:-1], range_bins_km[1:]):
            m = denom & (r_km >= lo) & (r_km < hi)
            n, k = int(m.sum()), int((m & det).sum())
            rng_rows.append({"definition": dname, "range_km": f"{lo}-{hi}", "n_ice_WG": n, "n_ka": k,
                             "fraction": k / n if n else np.nan})
    return pd.DataFrame(rows), pd.DataFrame(rng_rows)
