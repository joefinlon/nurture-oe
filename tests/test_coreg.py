"""Co-registration: synthetic unit tests plus regressions against PK numbers."""

import numpy as np
import pytest
import xarray as xr

from nurture_oe import coreg
from nurture_oe.checks.ka_coverage import zmin_dbz


def _times(seconds):
    return np.datetime64("2026-02-17T15:30:00", "ns") + (np.asarray(seconds) * 1e9).astype("timedelta64[ns]")


def test_window_sums_centered():
    w = coreg.halo_windows(_times([10.0, 20.0]))
    t = _times(np.arange(0, 30, 0.5))
    s, n = coreg.window_sums(np.ones((t.size, 1)), t, w)
    assert n.tolist() == [20, 20]          # [5, 15) and [15, 25) at 0.5 s
    assert s[:, 0].tolist() == [20, 20]


def test_zmin_matches_pk():
    assert zmin_dbz(1000.0, -19.26) == pytest.approx(-19.26)
    assert zmin_dbz(6000.0, -19.26) == pytest.approx(-3.7, abs=0.01)   # appendix A §A.4(i)


def _synthetic_cc(n=4, alt=10000.0):
    rng = np.arange(-100.0, 12000.0, 20.0)
    t = _times(np.arange(n) * 0.1)
    up = coreg.local_up_ecef(np.full(n, 54.0), np.full(n, -70.0))
    zm = np.full((n, rng.size), np.nan)
    zm[0, :], zm[1, :] = 0.0, 10.0      # two detecting rays, two non-detecting
    ds = xr.Dataset(
        {"zm_w": (("time", "range"), zm), "katten_w": (("time", "range"), np.full((n, rng.size), 0.5)),
         "alt": ("time", np.full(n, alt)), "lat": ("time", np.full(n, 54.0)), "lon": ("time", np.full(n, -70.0)),
         "look_vector_ecef": (("time", "nxyz"), -up)},
        coords={"time": t, "range": rng}, attrs={"zmin_1km_dbz_ka": -16.32})
    return ds


def test_gate_index_nadir():
    cc = _synthetic_cc()
    z = np.array([10000.0, 9000.0, 5000.0, 11000.0])
    idx, r = coreg.gate_index_on_grid(cc, z, cc["alt"].values, np.ones(4))
    assert r[0].tolist() == [0.0, 1000.0, 5000.0, -1000.0]
    assert cc.range.values[idx[0, 1]] == 1000.0
    assert idx[0, 3] == -1      # above the aircraft: excluded (min_range=0)


def test_gas_correction_constant_k():
    cc = _synthetic_cc()
    zgc = coreg.gas_corrected_display(cc, "w")
    i = int(np.argmin(np.abs(cc.range.values - 5000.0)))
    assert float(zgc[0, i]) == pytest.approx(0.0 + 2 * 0.5 * 5.0, abs=1e-6)   # 2 * k * r_km


def test_linear_mean_and_detect_fraction():
    cc = _synthetic_cc()
    halo = xr.Dataset(coords={"time": _times([0.15]), "z": np.array([5000.0])})
    w = coreg.halo_windows(halo.time.values)
    g = coreg.cloudcube_to_halo_grid(cc, halo, w, fields=("zm_w",), aircraft_alt="cloudcube")
    assert int(g.n_rays[0]) == 4
    assert float(g.zm_w_detect_frac[0, 0]) == pytest.approx(0.5)
    assert float(g.zm_w[0, 0]) == pytest.approx(10 * np.log10((1 + 10) / 2))


@pytest.mark.slow
def test_pk_phase_pixel_counts(aerowv, cs_type):
    """Appendix A §A.10.5: 10,395 liquid / 62,567 ice pixels > 2.5 km below the aircraft (1 % tolerance)."""
    w = coreg.halo_windows(aerowv.time.values)
    fr = coreg.cloudmask_fractions(cs_type, w)
    pm = coreg.phase_masks(fr)
    below = (aerowv.gps_alt.values[:, None] - aerowv.z.values[None, :]) > 2500
    assert int((pm.liquid.values & below).sum()) == pytest.approx(10395, rel=0.01)
    assert int((pm.ice.values & below).sum()) == pytest.approx(62567, rel=0.01)
    assert np.percentile(fr.n_sub.values, 5) == 20


@pytest.mark.slow
def test_pk_fill_fractions(aerowv, cs_type):
    """Appendix A §A.10.5(b): unscreened 532 0.91 liquid / 0.80 ice; screened 0.055 / 0.50."""
    w = coreg.halo_windows(aerowv.time.values)
    pm = coreg.phase_masks(coreg.cloudmask_fractions(cs_type, w))
    below = (aerowv.gps_alt.values[:, None] - aerowv.z.values[None, :]) > 2500
    liq, ice = pm.liquid.values & below, pm.ice.values & below
    f = np.isfinite(aerowv["532_total_attn_bsc"].values)
    assert f[liq].mean() == pytest.approx(0.91, abs=0.01) and f[ice].mean() == pytest.approx(0.80, abs=0.01)
    s = np.isfinite(aerowv["532_bsc_cloud_screened"].values)
    assert s[liq].mean() == pytest.approx(0.055, abs=0.005) and s[ice].mean() == pytest.approx(0.50, abs=0.01)
    assert np.array_equal(s, np.isfinite(aerowv["532_bsc_cloud_screened_unc"].values))


@pytest.mark.slow
def test_segment_coregistration_geometry(aerowv, cs_type):
    """15:30-16:00: ~100 rays per window; W surface peak within ~1 gate of the lidar surface."""
    from nurture_oe.config import example_file
    from nurture_oe.io.cloudcube import read_cloudcube
    seg = ("2026-02-17T15:30:00", "2026-02-17T16:00:00")
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cc = read_cloudcube(example_file("cloudcube_1521"), time_range=seg, profiles=("zm",), bands=("w",))
    halo = aerowv.sel(time=slice(np.datetime64(seg[0]) + coreg.HALF_WINDOW, np.datetime64(seg[1]) - coreg.HALF_WINDOW))
    w = coreg.halo_windows(halo.time.values)
    g = coreg.cloudcube_to_halo_grid(cc, halo, w, fields=("zm_w",))
    assert np.percentile(g.n_rays.values, 5) >= 99
    r = cc.range.values
    zw = np.where(r[None, :] < 5000, -np.inf, np.nan_to_num(cc.zm_w.values, nan=-np.inf))
    peak_alt = coreg.aircraft_alt_from_halo(cc, aerowv).values - r[np.argmax(zw, 1)] * coreg.cloudcube_cos_offnadir(cc).values
    sfc = np.interp(cc.time.values.astype(np.int64), cs_type.time.values.astype(np.int64), cs_type.altSurface.values)
    assert abs(np.nanmedian(peak_alt - sfc)) < 20.0
