"""Solver `pvc`: analytical A2-A6, numerical form N3, regression R2-R5, frozen fixture (appendix C §C.4.3)."""
import numpy as np
import pytest
from scipy.integrate import quad

from lidar_ms import Instrument, SolverConfig, simulate
from lidar_ms import grid as lgrid
from lidar_ms import mixing, pvc
from lidar_ms.anisotropy import half_level_factors
from lidar_ms.compat import read_ms_input
from lidar_ms.single import single
from msref import FROZEN_DIR, REF_DIR, REF_PAIRS, load_ref

NATIVE = SolverConfig(solver="pvc")
COMPAT = SolverConfig(solver="pvc", compat="multiscatter-1.2.11")
SINGLE = SolverConfig(solver="single")


def case(name, **kw):
    return read_ms_input(REF_DIR / f"{name}.in", **kw)


def run(c, config=NATIVE, instrument=None, optics=None):
    return simulate(c.grid, instrument or c.instrument, optics or c.optics, config)


def one_way_tau(c):
    a = c.optics.alpha[0, :, :3].sum(-1)
    dr = c.grid.dr[0]
    return np.cumsum(a * dr) - 0.5 * a * dr, a


# --- analytical ------------------------------------------------------------------------------------
def test_a2_vanishing_footprint():
    c = case("halo_liquid_30")
    tiny = Instrument(532e-9, 1e-9, 1e-9)
    np.testing.assert_allclose(run(c, NATIVE, tiny).total[0], run(c, SINGLE, tiny).total[0], rtol=1e-6)


def test_a3_wide_field_of_view():
    c = case("halo_liquid_30_iso")
    wide = Instrument(532e-9, c.instrument.rho_tr, 0.5)
    out = run(c, NATIVE, wide)
    a, dr = c.optics.alpha[0], c.grid.dr[0]
    ac, an = a[:, :3].sum(-1), a[:, 3]           # all particulate classes are lobe-forming (a = 8 um)
    tp = (ac + 2 * an) * dr
    T = np.exp(-np.concatenate([[0.0], np.cumsum(tp)[:-1]]))
    expect = c.optics.beta[0].sum(-1) * T * (-np.expm1(-tp)) / tp
    np.testing.assert_allclose(out.total[0], expect, rtol=1e-12)


def test_a4_no_lobe_equals_single():
    c = case("halo_liquid_30")
    opt = c.optics.replace(a_eq=np.full(c.optics.shape, 1e-7))   # Theta = 1.7 rad > theta_max
    a = run(c, NATIVE, optics=opt)
    b = run(c, SINGLE, optics=opt)
    np.testing.assert_array_equal(a.total_by_class, b.total_by_class)


def _a5_integral(c, d):
    inst = c.instrument
    a = c.optics.alpha[0, :, :3].sum(-1)
    i0 = int(np.argmax(a > 0))
    redge = c.grid.r_half[0, i0]
    th2 = (inst.wavelength / (np.pi * c.optics.a_eq[0, i0, 2])) ** 2
    rt2, rf2 = inst.rho_tr ** 2, inst.rho_fov ** 2
    R = redge + d
    f = lambda x: 1 - np.exp(-rf2 * R * R / (rt2 * R * R + th2 * (d - x) ** 2))
    return inst.kappa * a[i0] * quad(f, 0, d, epsabs=0, epsrel=1e-12)[0]


def test_a5_thin_slab():
    """M at the half-levels inside `halo_thin_slab` against the double-scattering integral (6 %)."""
    c = case("halo_thin_slab")
    out = run(c)
    a = c.optics.alpha[0, :, :3].sum(-1)
    idx = np.flatnonzero(a > 0)
    redge = c.grid.r_half[0, idx[0]]
    rel = [abs(out.ms_factor[0, i + 1] / _a5_integral(c, c.grid.r_half[0, i + 1] - redge) - 1) for i in idx]
    assert max(rel) <= 0.06, max(rel)


def _state_case(q_scale=1.0, rho_fov=0.5e-3):
    N, B = 80, 1
    r = (9000.0 + 15.0 * np.arange(N))[:, None]
    dr = lgrid.bin_widths(r.T).T
    rh = lgrid.half_levels(r.T).T
    ac = np.where((np.arange(N) > 20) & (np.arange(N) < 50), 0.02, 0.0)[:, None]
    an = np.full((N, B), 1e-5)
    q = q_scale * ac * (532e-9 / (np.pi * 8e-6)) ** 2
    tu = np.exp(-np.concatenate([np.zeros((1, B)), np.cumsum(2 * (ac + an) * dr, axis=0)]))
    inst = Instrument.from_full_angles(532e-9, 0.8e-3, 2 * rho_fov)
    return pvc.propagate(dr, rh, ac, an, q, tu, inst.rho_tr, inst.rho_fov, inst.kappa, return_state=True), rh, inst


def test_a6_moment_bookkeeping():
    (Mt, g2, valid, st), rh, inst = _state_case(q_scale=0.0)
    np.testing.assert_allclose(st["s2"][1:, 0], inst.rho_tr ** 2 * rh[1:, 0] ** 2, rtol=1e-12)
    (Mt, g2, valid, st), rh, inst = _state_case(q_scale=1.0, rho_fov=0.5)        # no clip
    v = valid[:, 0] & ~st["clipped"][:, 0]
    assert v.sum() > 20
    tu = np.exp(-np.concatenate([[0.0], np.cumsum(2 * (np.where((np.arange(80) > 20) & (np.arange(80) < 50), 0.02, 0.0) + 1e-5) * 15.0)]))
    T, Ta, Tb = st["T"][v, 0], st["Ta"][v, 0], st["Tb"][v, 0]
    np.testing.assert_allclose(T, tu[v] + Ta + Tb, rtol=1e-12)
    lhs = T * st["s2"][v, 0]
    rhs = tu[v] * st["s2u"][v, 0] + Ta * st["s2a"][v, 0] + Tb * st["s2b"][v, 0]
    np.testing.assert_allclose(lhs, rhs, rtol=1e-12)


# --- numerical form ---------------------------------------------------------------------------------
def _internals(c, scale=1.0):
    o = c.optics
    alpha = np.moveaxis(o.alpha * np.array([scale, scale, scale, 1.0]), 0, -1)
    beta = np.moveaxis(o.beta, 0, -1)
    r = c.grid.r.T
    dr = lgrid.bin_widths(r.T).T
    beta_ss, tu, tau = single(alpha.sum(1), beta, dr)
    th2 = mixing.lobe_variance(np.moveaxis(o.a_eq, 0, -1), c.instrument.wavelength)
    lobe = mixing.lobe_mask(th2, np.array([True, True, True, False]), 0.1)
    ac, an, q = mixing.pvc_bin(alpha, th2, lobe)
    Mt, g2, valid = pvc.propagate(dr, lgrid.half_levels(r.T).T, ac, an, q, tu, c.instrument.rho_tr,
                                  c.instrument.rho_fov, c.instrument.kappa)
    Fm, Fp, _ = half_level_factors(g2, th2, o.descriptors, o.kinds)
    return beta, tu, Mt, Fm, Fp, tau


@pytest.mark.parametrize("name,scale", [("i3rc_lidar_exp5", 1.0), ("halo_liquid_30", 1.0), ("halo_liquid_30", 100.0)])
def test_n3_log_mean_form(name, scale):
    c = case(name)
    beta, tu, Mt, Fm, Fp, tau = _internals(c, scale)
    lm = pvc.bin_means(beta, tu, Mt, Fm, Fp)
    assert np.all(np.isfinite(lm))
    direct = pvc.bin_means_direct(beta, tu, Mt, Fm, Fp, tau)
    ok = np.isfinite(direct) & (beta > 0)
    assert ok.any()
    np.testing.assert_allclose(lm[ok], direct[ok], rtol=1e-10)
    if scale > 1:
        assert np.real(tu[-1]).max() == 0.0            # T_u underflows


def test_n3_package_finite():
    for name in ("i3rc_lidar_exp5", "i3rc_lidar_exp6"):
        for cfg in (NATIVE, COMPAT):
            assert np.all(np.isfinite(run(case(name), cfg).total))


# --- regression -------------------------------------------------------------------------------------
PVC_REFS = [(nm, t) for nm, t in REF_PAIRS if t in ("pvc", "pvc_hsrl")]


@pytest.mark.parametrize("name,tag", PVC_REFS)
def test_r2_compat_vs_multiscatter(name, tag):
    c = case(name)
    out = run(c, COMPAT)
    tb = out.total_by_class[0]
    assert np.all(np.isfinite(tb))
    ref = load_ref(name, tag).data
    if tag == "pvc_hsrl":
        pairs = [(tb[:, :3].sum(-1), ref[:, 4]), (tb[:, 3], ref[:, 5])]
    else:
        pairs = [(tb.sum(-1), ref[:, 4])]
    for got, want in pairs:
        ok = np.isfinite(want)
        np.testing.assert_allclose(got[ok], want[ok], rtol=1e-7, atol=0)


def test_r3_native_vs_compat():
    c = case("halo_liquid_30")
    nat, cmp_ = run(c, NATIVE).total_by_class[0], run(c, COMPAT).total_by_class[0]
    cloud = c.optics.alpha[0, :, 0] > 0
    np.testing.assert_allclose(nat[cloud, 0], cmp_[cloud, 0], rtol=0.01)
    above = c.z > 1500
    # "identical": the two modes use different (algebraically equal) single-scatter arithmetic
    np.testing.assert_allclose(nat[above, 3], cmp_[above, 3], rtol=1e-12)
    assert np.max(np.abs(nat[cloud, 3] / cmp_[cloud, 3] - 1)) > 0.5          # Dep2 is present


@pytest.mark.parametrize("name,tol1,tol3", [("halo_liquid_30_iso", 0.05, 0.10), ("halo_liquid_10_iso", 0.05, 0.10)])
def test_r4_native_vs_h06(name, tol1, tol3):
    c = case(name)
    got = run(c).total[0]
    ref = load_ref(name, "h06").data[:, 4]
    tau, a = one_way_tau(c)
    rel = np.abs(got / ref - 1)
    assert rel[(a > 0) & (tau <= 1)].max() <= tol1
    assert rel[(tau > 1) & (tau <= 3)].max() <= tol3


def test_r5_explicit2_vs_thin_slab_integral():
    """M implied by the C code's exact double scattering vs the bin mean of the A5 integral (0.2 %)."""
    c = case("halo_thin_slab")
    e2 = load_ref("halo_thin_slab", "explicit2").data[:, 4]
    s1 = load_ref("halo_thin_slab", "single").data[:, 4]
    a = c.optics.alpha[0, :, :3].sum(-1)
    idx = np.flatnonzero(a > 0)
    rh = c.grid.r_half[0]
    redge = rh[idx[0]]
    for i in idx:
        lo, hi = rh[i] - redge, rh[i + 1] - redge
        mbar = quad(lambda d: _a5_integral(c, d), lo, hi, epsrel=1e-10)[0] / (hi - lo)
        assert abs((e2[i] / s1[i] - 1) / mbar - 1) <= 2e-3


# --- frozen fixture (§C.2.4, regenerated at M4) --------------------------------------------------
def test_fixture_halo_liquid_30():
    fz = np.loadtxt(FROZEN_DIR / "halo_liquid_30_pvc_native.csv", delimiter=",")
    c = case("halo_liquid_30")
    out = run(c)
    gate = fz[:, 0].astype(int)
    ratio = out.total_by_class[0, gate, 0] / out.beta_ss[0, gate, 0]
    np.testing.assert_allclose(ratio, fz[:, 3], rtol=1e-10)
    np.testing.assert_allclose(out.ms_factor[0, gate + 1], fz[:, 4], rtol=1e-10)
    provisional = {7.5: 1.14, 37.5: 2.21, 67.5: 4.53, 97.5: 9.32}
    for depth, val in provisional.items():
        assert fz[fz[:, 1] == depth, 3][0] == pytest.approx(val, abs=0.005)
