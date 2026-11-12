"""Depolarization closures: K1-K4, K6 and J4 (appendix C §C.4.3). K5 waits for the Mie fixture (M9)."""
import numpy as np
import pytest

from lidar_ms import (Instrument, NearBackscatter, OpticalProfile, OpticalTangent, RangeGrid, SolverConfig,
                      check_jvp, simulate)
from lidar_ms.anisotropy import class_factor, p_droplet, p_yang, table_factor
from lidar_ms.compat import read_ms_input
from lidar_ms.constants import FLAG_BITS
from msref import FROZEN_DIR, REF_DIR

HALO = Instrument.from_full_angles(532e-9, 0.8e-3, 1.0e-3)
PVC = SolverConfig(solver="pvc")
X_DROP = 10.0 * (np.arange(301) / 300.0) ** 2              # gamma / Theta
SPHERE = NearBackscatter("gamma_over_theta", X_DROP, p_droplet(X_DROP), np.zeros(X_DROP.size))
X_ICE = np.linspace(0.0, 0.5, 501)                          # gamma, rad


def ice_table(delta):
    c0 = delta / (1 + delta)
    p = p_yang(X_ICE)
    return NearBackscatter("gamma", X_ICE, (1 - c0) * p, c0 * p)


def case(name):
    return read_ms_input(REF_DIR / f"{name}.in")


def ice_over_liquid(delta_ice=0.4):
    """`halo_ice_over_liquid`: the ice layer is the file's isotropic ('other') particulate class."""
    c = case("halo_ice_over_liquid")
    d = np.array(c.optics.delta)
    d[..., 2] = delta_ice
    return c, c.optics.replace(delta=d)


def delta_acc(out, j=0):
    dr = np.diff(out.r[0]).mean()
    ip = np.cumsum(out.beta_perp[0, :, j]) * dr
    ia = np.cumsum(out.beta_par[0, :, j]) * dr
    with np.errstate(invalid="ignore", divide="ignore"):
        return ip / ia


# --- K1, K2 ---------------------------------------------------------------------------------------
def test_k1_c1_identity():
    c = case("halo_liquid_30")
    out = simulate(c.grid, c.instrument, c.optics, SolverConfig(solver="pvc", closure="hu"))
    cl = c.optics.alpha[0, :, 0] > 0
    a_s = out.a_s[0, cl]
    np.testing.assert_allclose(delta_acc(out)[cl], (1 - np.sqrt(a_s)) / (1 + np.sqrt(a_s)), rtol=1e-12)


def test_k2_c1_halo_liquid_30():
    c = case("halo_liquid_30")
    out = simulate(c.grid, c.instrument, c.optics, PVC)
    assert list(out.closure_used) == ["hu"]
    cl = np.flatnonzero(c.optics.alpha[0, :, 0] > 0)
    with np.errstate(invalid="ignore"):
        d = out.beta_perp[0, :, 0] / out.beta_par[0, :, 0]
    assert out.a_s[0, cl[-1]] == pytest.approx(0.635, abs=0.005)
    depth = 1500.0 - c.z + 7.5
    for dep, val in {7.5: 0.034, 37.5: 0.180, 67.5: 0.293, 97.5: 0.360}.items():
        assert d[np.isclose(depth, dep)][0] == pytest.approx(val, abs=0.005)
    fz = np.loadtxt(FROZEN_DIR / "halo_liquid_30_c1_native.csv", delimiter=",")
    gate = fz[:, 0].astype(int)
    np.testing.assert_allclose(d[gate], fz[:, 2], rtol=1e-10)
    np.testing.assert_allclose(out.a_s[0, gate], fz[:, 3], rtol=1e-10)


# --- K3 ----------------------------------------------------------------------------------------------
def test_k3_c1_under_ice():
    c, opt = ice_over_liquid()
    out = simulate(c.grid, c.instrument, opt, SolverConfig(solver="pvc", closure="hu"))
    liq = np.flatnonzero(opt.alpha[0, :, 0] > 0)
    bit3 = 1 << FLAG_BITS["C1_UNDER_OTHER_LOBE"]
    assert np.all(out.flags[0, liq] & bit3)
    assert not np.any(out.flags[0, opt.alpha[0, :, 0] == 0] & bit3)
    assert delta_acc(out)[liq[0]] == pytest.approx(0.33, abs=0.01)


def test_k3_auto_selection():
    c, opt = ice_over_liquid()
    tables = {"droplet": SPHERE, "isotropic": ice_table(0.4)}
    with_t = simulate(c.grid, c.instrument, opt, PVC, near=tables)
    without = simulate(c.grid, c.instrument, opt, PVC)
    assert list(with_t.closure_used) == ["near_backscatter"]
    assert list(without.closure_used) == ["hu"]
    liq = opt.alpha[0, :, 0] > 0
    assert np.all(without.flags[0, liq] & (1 << FLAG_BITS["C1_UNDER_OTHER_LOBE"]))
    assert not np.any(with_t.flags[0] & (1 << FLAG_BITS["C1_UNDER_OTHER_LOBE"]))
    # liquid only -> C1; no liquid and no tables -> single-scatter split
    assert list(simulate(*case_parts("halo_liquid_30"), PVC).closure_used) == ["hu"]
    noliq = opt.replace(alpha=opt.alpha * np.array([0, 1, 1, 1]), beta=opt.beta * np.array([0, 1, 1, 1]))
    assert list(simulate(c.grid, c.instrument, noliq, PVC).closure_used) == ["none"]
    assert list(simulate(c.grid, c.instrument, noliq, PVC, near={"isotropic": ice_table(0.4)}).closure_used) \
        == ["near_backscatter"]


def case_parts(name):
    c = case(name)
    return c.grid, c.instrument, c.optics


# --- K4 ----------------------------------------------------------------------------------------------
def test_k4_components_sum_to_total():
    tab = ice_table(0.4)
    g2 = np.logspace(-8, -1, 30)
    tot, _ = class_factor("table", "ice", g2, 1.0, tab, None, "total")
    par, _ = class_factor("table", "ice", g2, 1.0, tab, None, "par")
    perp, _ = class_factor("table", "ice", g2, 1.0, tab, None, "perp")
    np.testing.assert_allclose(par + perp, tot, rtol=1e-10)


def test_k4_vanishing_ms_returns_single_scatter_delta():
    c, opt = ice_over_liquid(0.4)
    tiny = Instrument(532e-9, 1e-9, 1e-9)
    out = simulate(c.grid, tiny, opt, SolverConfig(solver="pvc", closure="near_backscatter"),
                   near={"droplet": SPHERE, "isotropic": ice_table(0.4)})
    ice = opt.alpha[0, :, 2] > 0
    tb = out.total_by_class[0, ice, 2]
    np.testing.assert_allclose(out.beta_perp[0, ice, 2] / tb, 0.4 / 1.4, rtol=1e-10)
    liq = opt.alpha[0, :, 0] > 0
    np.testing.assert_allclose(out.beta_perp[0, liq, 0], 0.0, atol=1e-10 * tb.max())


def test_k4_sphere_zero_cross_pol_in_single_scattering():
    c, opt = ice_over_liquid(0.4)
    out = simulate(c.grid, c.instrument, opt, SolverConfig(solver="single", closure="near_backscatter"),
                   near={"droplet": SPHERE, "isotropic": ice_table(0.4)})
    assert np.all(out.beta_perp[0, :, 0] == 0)


def test_k4_rayleigh_table():
    g = np.linspace(0.0, 0.3, 3001)
    cg = np.cos(g)
    c_par = (2.25 * (1 + cg ** 2) + 1.5 * cg) / 6.0         # (C.7) with Rayleigh P11, P22, P33 at pi - gamma
    c_perp = (1 - cg) ** 2 / 8.0
    NearBackscatter("gamma", g, c_par, c_perp)
    g2 = np.logspace(-8, -3, 40)
    F, _ = table_factor(g2, g, np.broadcast_to(c_perp, (g2.size, g.size)))
    assert np.all(F <= 1e-6)


# --- K5 ----------------------------------------------------------------------------------------------
@pytest.mark.skip(reason="needs the Mie fixture (tests/lidar_ms/mie_fixture.py, milestone M9); "
                         "ordering conflict recorded in docs/findings_for_pk.md")
def test_k5_c2_under_ice():
    pass


# --- K6 ----------------------------------------------------------------------------------------------
def _liquid(segs, n=200, r0=9000.0, desc=("droplet",), a_eq=8e-6, extra=None):
    J = len(desc) + 1
    sh = (1, n, J)
    a = np.zeros(sh)
    for j, (i0, i1, al) in segs:
        a[0, i0:i1, j] = al
    a[0, :, -1] = 1e-5
    b = a / 18.8
    b[..., -1] = 1e-5 * 3 / (8 * np.pi)
    kinds = tuple(extra or ["liquid"] * len(desc)) + ("molecular",)
    ae = np.full(sh, a_eq) if np.isscalar(a_eq) else np.broadcast_to(np.asarray(a_eq), sh).copy()
    d = np.zeros(sh)
    return RangeGrid(r0 + 15.0 * np.arange(n)), OpticalProfile(
        kinds=kinds, descriptors=tuple(desc) + ("isotropic",), names=tuple(f"c{j}" for j in range(J)),
        alpha=a, omega=np.ones(sh), g=np.full(sh, 0.85), a_eq=ae, beta=b, delta=d)


def _flags(grid, opt, config=PVC, near=None, inst=HALO):
    return int(np.bitwise_or.reduce(simulate(grid, inst, opt, config, near=near).flags[0]))


def _k6_cases():
    rise = NearBackscatter("gamma", np.array([0, 1e-3, 1.0]), np.array([1.0, 8, 8]), np.zeros(3))
    spike = NearBackscatter("gamma", np.array([0, 1e-4, 1.0]), np.array([1.0, 0, 0]), np.zeros(3))
    short = NearBackscatter("gamma", np.linspace(0, 0.005, 11), 0.7 * p_yang(np.linspace(0, 0.005, 11)),
                            0.3 * p_yang(np.linspace(0, 0.005, 11)))
    none = SolverConfig(solver="pvc", closure="none")
    return {
        # p > 1 away from pi drives A_s below 0.144 (not reachable with closed-form descriptors)
        "C1_DELTA_RANGE": lambda: _flags(*_liquid([(0, (20, 40, 0.03))], desc=("table",)), near={"c0": rise}),
        "C1_SMALL_FOOTPRINT": lambda: _flags(*_liquid([(0, (20, 40, 0.03))], r0=1700.0)),
        # second liquid layer whose near-180 spike makes F ~ 0: local single-scatter fraction ~ 1
        "C1_NEGATIVE": lambda: _flags(*_liquid([(0, (10, 20, 0.03)), (1, (140, 150, 0.03))],
                                               desc=("droplet", "table"), n=260), near={"c1": spike}),
        "C1_UNDER_OTHER_LOBE": lambda: _flags(*_liquid([(0, (40, 60, 0.03)), (1, (10, 30, 0.001))],
                                                       desc=("droplet", "isotropic"), extra=["liquid", "ice"],
                                                       a_eq=np.array([8e-6, 50e-6, 1e-5]))),
        "LOBE_MIXTURE": lambda: _flags(*_liquid([(0, (20, 40, 0.005)), (1, (20, 40, 0.001))],
                                                desc=("droplet", "isotropic"), extra=["liquid", "ice"],
                                                a_eq=np.array([8e-6, 100e-6, 1e-5])), config=none),
        "TABLE_EXTRAPOLATED": lambda: _flags(*_liquid([(0, (20, 60, 0.002))], desc=("table",), extra=["ice"],
                                                      a_eq=50e-6), config=none, near={"c0": short}),
        "EXTINGUISHED": lambda: _flags(*_liquid([(0, (20, 120, 0.02))], desc=("isotropic",), extra=["ice"],
                                                a_eq=50e-6), config=none),
    }


@pytest.mark.parametrize("name", list(_k6_cases()))
def test_k6_each_flag_alone(name):
    assert _k6_cases()[name]() == 1 << FLAG_BITS[name]


def test_k6_clean_case_raises_nothing():
    c = case("halo_liquid_30")
    assert not simulate(c.grid, c.instrument, c.optics, PVC).flags.any()


# --- J4 ----------------------------------------------------------------------------------------------
def _tangents(opt, K=3, seed=0):
    rng = np.random.default_rng(seed)
    P, N, J = opt.shape
    return OpticalTangent(alpha=np.real(opt.alpha)[:, None] * rng.standard_normal((P, K, N, J)),
                          a_eq=np.real(opt.a_eq)[:, None] * rng.standard_normal((P, K, N, J)),
                          beta=np.real(opt.beta)[:, None] * rng.standard_normal((P, K, N, J)),
                          delta=0.05 * rng.standard_normal((P, K, N, J)) * (np.real(opt.delta)[:, None] > 0))


def test_j4_through_c1():
    c = case("halo_liquid_30")
    err = check_jvp(c.grid, c.instrument, c.optics, _tangents(c.optics), SolverConfig(solver="pvc", closure="hu"))
    assert err <= 1e-6, err


def test_j4_through_c2():
    c, opt = ice_over_liquid(0.4)
    err = check_jvp(c.grid, c.instrument, opt, _tangents(opt), SolverConfig(solver="pvc", closure="near_backscatter"),
                    near={"droplet": SPHERE, "isotropic": ice_table(0.4)})
    assert err <= 1e-6, err


def test_diagnose_closures_returns_the_other_split():
    c, opt = ice_over_liquid(0.4)
    tables = {"droplet": SPHERE, "isotropic": ice_table(0.4)}
    out = simulate(c.grid, c.instrument, opt, SolverConfig(solver="pvc", diagnose_closures=True), near=tables)
    hu = simulate(c.grid, c.instrument, opt, SolverConfig(solver="pvc", closure="hu"), near=tables)
    assert list(out.closure_used) == ["near_backscatter"] and list(out.diagnostic["closure"]) == ["hu"]
    np.testing.assert_allclose(out.diagnostic["beta_perp"], hu.beta_perp, rtol=1e-14)
