"""Moist thermodynamics for diagnostics (not yet the F(x) gas model).

Saturation vapour pressures from Murphy & Koop (2005, QJRMS 131, 1539-1565), Eqs. 7 (ice)
and 10 (liquid, valid 123-332 K, incl. supercooled).
"""

import numpy as np

EPS = 0.621981  # M_w / M_d
ATM_PA = 101325.0


def esat_ice(t_k):
    """Saturation vapour pressure over ice (Pa), Murphy & Koop (2005) Eq. 7."""
    t = np.asarray(t_k, float)
    return np.exp(9.550426 - 5723.265 / t + 3.53068 * np.log(t) - 0.00728332 * t)


def esat_liquid(t_k):
    """Saturation vapour pressure over (supercooled) liquid water (Pa), Murphy & Koop (2005) Eq. 10."""
    t = np.asarray(t_k, float)
    return np.exp(54.842763 - 6763.22 / t - 4.210 * np.log(t) + 0.000367 * t
                  + np.tanh(0.0415 * (t - 218.8)) * (53.878 - 1331.22 / t - 9.44523 * np.log(t) + 0.014025 * t))


def vapour_pressure(mmr_kgkg, p_pa):
    """Vapour pressure (Pa) from water-vapour mass mixing ratio (kg/kg, per dry air) and pressure."""
    w = np.asarray(mmr_kgkg, float)
    return w * np.asarray(p_pa, float) / (EPS + w)


def rh_ice(mmr_gkg, p_atm, t_k):
    """Relative humidity over ice (%) from ``h2o_mmr_v`` (g/kg) and HALO ``/State`` p (atm), T (K)."""
    e = vapour_pressure(np.asarray(mmr_gkg, float) * 1e-3, np.asarray(p_atm, float) * ATM_PA)
    return 100.0 * e / esat_ice(t_k)
