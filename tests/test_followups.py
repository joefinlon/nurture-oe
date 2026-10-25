"""Tests for code added with the appendix A §A.8(5) follow-ups (2026-10-20)."""

import warnings

import numpy as np
import pytest

from nurture_oe import coreg
from nurture_oe.config import example_file
from nurture_oe.thermo import esat_ice, esat_liquid, rh_ice


def _times(seconds):
    return np.datetime64("2026-02-17T15:30:00", "ns") + (np.asarray(seconds) * 1e9).astype("timedelta64[ns]")


def test_esat_triple_point():
    assert esat_ice(273.16) == pytest.approx(611.657, abs=0.01)
    assert esat_liquid(273.16) == pytest.approx(611.657, abs=0.05)
    assert esat_ice(253.15) < esat_liquid(253.15)


def test_rh_ice_saturated_air():
    t, p_atm = 250.0, 0.5
    e = esat_ice(t)
    w = 0.621981 * e / (p_atm * 101325.0 - e)          # mixing ratio at ice saturation
    assert rh_ice(w * 1e3, p_atm, t) == pytest.approx(100.0, rel=1e-6)


def test_subprofile_windows_tile():
    t = _times(np.arange(10) * 0.5)
    w = coreg.subprofile_windows(t)
    assert np.all(w.t_hi.values[:-1] == w.t_lo.values[1:])   # no gaps, no overlap
    assert w.attrs["half_width_s"] == pytest.approx(0.25)


def test_geoid_matches_cloudcube_alt_conversion():
    """CloudCube alt = MetNav MSL + N(EGM96) to centimetres (follow-up (a), 19:00 granule)."""
    pytest.importorskip("pyproj")
    from nurture_oe.geoid import geoid_undulation
    from nurture_oe.io.cloudcube import read_cloudcube
    from nurture_oe.io.metnav import read_metnav
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cc = read_cloudcube(example_file("cloudcube_1900"), profiles=(), time_range=("2026-02-17T19:00", "2026-02-17T19:05"))
    m = read_metnav(example_file("metnav"))
    tm = m.time.values.astype("datetime64[ns]").astype(np.int64)
    tc = cc.time.values.astype("datetime64[ns]").astype(np.int64)
    msl = np.interp(tc, tm, m["GPS_Altitude_MSL"].values)
    n96 = geoid_undulation(cc["lat"].values, cc["lon"].values, "EGM96")
    res = msl - cc["alt"].values + n96
    assert np.nanpercentile(np.abs(res), 99) < 0.15
    assert 10 < np.nanmedian(-n96) < 30


def test_arts_fit_power_and_coverage():
    import pandas as pd
    from nurture_oe.ssp.arts_catalog import fit_power, largest_with
    d = np.array([1e-4, 1e-3, 1e-2])
    a, b, rms = fit_power(d, 0.2 * d ** 2.4)
    assert a == pytest.approx(0.2) and b == pytest.approx(2.4) and rms < 1e-12
    sizes = pd.DataFrame({"dveq_m": [1e-4, 1e-3, 2e-3], "dmax_m": [2e-4, 2e-3, 4e-3],
                          "groups": [{(247.2, 230.0), (94.1, 230.0)}, {(247.2, 230.0), (94.1, 230.0)}, {(94.1, 230.0)}]})
    assert largest_with(sizes, 94.1)[0] == 2e-3
    assert largest_with(sizes, 247.2)[0] == 1e-3


def test_arts_catalog_one_habit():
    from nurture_oe.config import load_paths
    from nurture_oe.ssp.arts_catalog import habit_dirs, read_summary, size_table
    root = load_paths()["arts_ssdb"]
    d = [h for h in habit_dirs(root) if h[3].name == "5bulletRosette_Id2"][0][3]
    s = size_table(d)
    summ = read_summary(d)
    assert len(s) == summ["N_SIZE"] == 45
    assert s.dmax_m.max() == pytest.approx(summ["DMAX_MAX"], rel=1e-4)
    assert set(s.orient_type) == {"totally_random"}
    assert s.n_freq.max() == 34
