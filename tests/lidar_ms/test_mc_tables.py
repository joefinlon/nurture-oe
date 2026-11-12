"""Solver `mc` building blocks: X1-X4 (appendix C §C.4.3, milestone M9)."""
import numpy as np
import pytest

import mie_fixture as mf
from lidar_ms.mc import fixtures as fx
from lidar_ms.mc import geometry as geo
from lidar_ms.mc import tables as tb

TH = fx.angle_grid(1e-4)


def _sample_cos(p11, n=1_000_000, seed=1):
    table = tb.inverse_cdf(TH, p11)
    rng = np.random.default_rng(seed)
    c = np.cos(tb.sample_theta(table, None, rng.random(n)))
    return c.mean(), c.std(ddof=1) / np.sqrt(n)


# --- X1 ------------------------------------------------------------------------------------------
@pytest.mark.parametrize("g", [0.0, -0.3, 0.5, 0.85, 0.95])
def test_x1_hg_inverse_cdf_closed_form(g):
    table = tb.inverse_cdf(TH, fx.henyey_greenstein_p11(TH, g))
    u = np.linspace(0, 1, table.size)
    assert np.max(np.abs(np.cos(table) - fx.hg_inverse_cdf(u, g))) <= 1e-4


def _x1_fixtures():
    lobe = fx.gaussian_lobe(TH, 0.0212)["P11"]
    return {
        "hg_0.85": (fx.henyey_greenstein(TH, 0.85)["P11"], 0.85),
        "hg_-0.3": (fx.henyey_greenstein(TH, -0.3)["P11"], -0.3),
        "dhg": (fx.double_henyey_greenstein(TH, 0.9, 0.9, -0.5)["P11"], 0.9 * 0.9 - 0.1 * 0.5),
        "isotropic": (fx.isotropic(TH)["P11"], 0.0),
        "rayleigh": (fx.rayleigh(TH)["P11"], 0.0),
        "gaussian_lobe": (lobe, fx.asymmetry(TH, lobe)),
    }


@pytest.mark.parametrize("name", list(_x1_fixtures()))
def test_x1_sampled_asymmetry(name):
    p11, g = _x1_fixtures()[name]
    assert fx.sphere_integral(TH, p11) == pytest.approx(2.0, rel=1e-3)
    m, se = _sample_cos(p11)
    assert abs(m - g) <= 3 * se, (m, g, se)


def test_x1_mie_sampled_asymmetry():
    d = mf.droplet_optics()
    table = tb.inverse_cdf(d["theta"], d["P11"])
    rng = np.random.default_rng(2)
    c = np.cos(tb.sample_theta(table, None, rng.random(1_000_000)))
    se = c.std(ddof=1) / 1e3
    # g from the series; the tabulated matrix differs from it by its angular quadrature only
    assert abs(c.mean() - d["g"]) <= 3 * se + 2e-4, (c.mean(), d["g"], se)


# --- X2 ------------------------------------------------------------------------------------------
@pytest.mark.parametrize("m,x", [(1.334 + 0j, 0.5), (1.334 + 0j, 12.0), (1.334 + 1e-3j, 60.0),
                                 (1.55 + 0.01j, 5.213), (1.334 + 0j, 150.0)])
def test_x2_optical_theorem_and_normalization(m, x):
    o = mf.single(m, x, np.array([0.0]))
    assert o["qext"] == pytest.approx(4.0 / x ** 2 * np.real(o["S1"][0]), rel=1e-6)
    assert o["S1"][0] == pytest.approx(o["S2"][0], rel=1e-12)
    nq = 2 * mf.n_stop(x) + 20
    mu, w = np.polynomial.legendre.leggauss(nq)
    s = mf.single(m, x, np.arccos(mu))
    i = np.abs(s["S1"]) ** 2 + np.abs(s["S2"]) ** 2
    assert np.sum(w * i) / x ** 2 == pytest.approx(o["qsca"], rel=1e-6)
    assert np.sum(w * i * mu) / np.sum(w * i) == pytest.approx(o["g"], rel=1e-6)
    if m.imag == 0:
        assert o["qsca"] == pytest.approx(o["qext"], rel=1e-10)


def test_x2_rayleigh_limit():
    m, x = 1.334 + 0j, 1e-3
    th = np.linspace(0, np.pi, 91)
    o = mf.single(m, x, th)
    K = (m * m - 1) / (m * m + 2)
    assert o["qsca"] == pytest.approx(8.0 / 3.0 * x ** 4 * abs(K) ** 2, rel=1e-4)
    i1, i2 = np.abs(o["S1"]) ** 2, np.abs(o["S2"]) ** 2
    p11 = 0.5 * (i1 + i2)
    np.testing.assert_allclose(p11 / p11[0], 0.5 * (1 + np.cos(th) ** 2), rtol=1e-5)
    np.testing.assert_allclose(0.5 * (i2 - i1) / p11, -np.sin(th) ** 2 / (1 + np.cos(th) ** 2), atol=1e-5)


def test_x2_extinction_efficiency_tends_to_two():
    m = 1.334 + 0j
    dev = []
    for x0 in (30.0, 300.0, 3000.0):
        x = x0 * (1 + np.linspace(-0.1, 0.1, 401))           # average over ripple
        a, b = mf.coefficients(m, x)
        qe, _, _ = mf.efficiencies(a, b, x)
        dev.append(abs(qe.mean() - 2.0))
    assert dev[0] > dev[1] > dev[2]
    assert dev[2] < 0.02, dev


@pytest.mark.parametrize("a_eq", [5e-6, 8e-6, 12e-6])
def test_x2_water_cloud_lidar_ratio(a_eq):
    d = mf.droplet_optics(a_eq=a_eq)
    assert 18.0 <= d["lidar_ratio"] <= 20.0, d["lidar_ratio"]
    assert d["omega"] == pytest.approx(1.0, abs=1e-12)
    assert fx.sphere_integral(d["theta"], d["P11"]) == pytest.approx(2.0, rel=1e-12)
    assert d["a_eq"] == pytest.approx(a_eq, rel=1e-4)


# --- X3 ------------------------------------------------------------------------------------------
def test_x3_rotation_inverse():
    psi = np.random.default_rng(3).uniform(-np.pi, np.pi, 1000)
    prod = geo.rotation_matrix(psi) @ geo.rotation_matrix(-psi)
    np.testing.assert_allclose(prod, np.broadcast_to(np.eye(4), prod.shape), atol=1e-10)
    S = np.random.default_rng(4).standard_normal((4, 1000))
    c2, s2 = np.cos(2 * psi), np.sin(2 * psi)
    np.testing.assert_allclose(geo.rotate_stokes(S, c2, s2), np.einsum("nij,jn->in", geo.rotation_matrix(psi), S),
                               atol=1e-12)


def test_x3_orthonormal_over_many_scatterings():
    rng = np.random.default_rng(5)
    n = 2000
    u = np.zeros((3, n))
    u[2] = 1.0
    e = np.zeros((3, n))
    e[0] = 1.0
    worst = 0.0
    for k in range(1000):
        if k % 2:
            v = geo.new_direction(u, e, rng.uniform(-1, 1, n), rng.uniform(0, 2 * np.pi, n))
        else:                                               # arbitrary target direction (DDIS branch, LE)
            v = geo.normalize(rng.standard_normal((3, n)))
        _, _, c2, s2, l_out = geo.scattering_plane(u, e, v)
        u, e = v, l_out
        worst = max(worst, np.max(np.abs(geo.norm(u) - 1)), np.max(np.abs(geo.norm(e) - 1)),
                    np.max(np.abs(geo.dot(u, e))), np.max(np.abs(c2 * c2 + s2 * s2 - 1)))
    assert worst <= 1e-10, worst


def test_x3_azimuth_equals_rotation_angle():
    """Sampling about u at azimuth phi from e rotates the reference axis by psi = phi (§C.2.7)."""
    rng = np.random.default_rng(6)
    n = 1000
    u = geo.normalize(rng.standard_normal((3, n)))
    e = geo.perpendicular_axis(rng.standard_normal((3, n)), u)
    phi = rng.uniform(-np.pi, np.pi, n)
    v = geo.new_direction(u, e, rng.uniform(-0.99, 0.99, n), phi)
    _, _, c2, s2, _ = geo.scattering_plane(u, e, v)
    np.testing.assert_allclose(c2, np.cos(2 * phi), atol=1e-10)
    np.testing.assert_allclose(s2, np.sin(2 * phi), atol=1e-10)


# --- X4 ------------------------------------------------------------------------------------------
def _layered():
    r_half = 1000.0 + 15.0 * np.arange(41)
    a = np.zeros((40, 2))
    a[5:15, 0] = 0.02
    a[22:30, 0] = 0.005
    a[:, 1] = 1e-5
    return tb.build_medium(r_half, a, np.ones_like(a), np.zeros((40, 2), int))


def test_x4_tau_inversion_round_trip():
    med = _layered()
    tau = np.random.default_rng(7).uniform(0, med.tau_nodes[-1], 100_000)
    np.testing.assert_allclose(tb.tau_of_z(med, tb.z_of_tau(med, tau)), tau, rtol=1e-12, atol=1e-12 * tau.max())
    assert np.all(np.diff(med.tau_nodes) > 0)
    z = np.random.default_rng(8).uniform(med.r_half[0], med.z_max, 100_000)
    np.testing.assert_allclose(tb.z_of_tau(med, tb.tau_of_z(med, z)), z, rtol=1e-12)


def _optical_path(med, z0, z1, mu):
    lo, hi = np.minimum(z0, z1), np.maximum(z0, z1)
    top, bot = med.z_nodes[:-1], med.z_nodes[1:]
    ov = np.clip(np.minimum(hi[:, None], bot[None]) - np.maximum(lo[:, None], top[None]), 0, None)
    return (ov * med.alpha[None]).sum(axis=1) / np.abs(mu)


def test_x4_free_paths_are_exponential():
    """Sampled optical paths through the layered medium follow P(t > c) = exp(-c)."""
    med = _layered()
    rng = np.random.default_rng(9)
    n = 200_000
    z0 = rng.uniform(med.r_half[0], med.z_max, n)
    mu = rng.uniform(-1, 1, n)
    mu = np.where(np.abs(mu) < 1e-3, 1e-3, mu)
    tau = -np.log(rng.random(n))
    target = tb.tau_of_z(med, z0) + mu * tau
    out = (target < 0) | (target > med.tau_nodes[-1])
    z1 = tb.z_of_tau(med, np.clip(target, 0, med.tau_nodes[-1]))
    t = np.where(out, np.inf, _optical_path(med, z0, z1, mu))
    t_exit = np.where(mu > 0, med.tau_nodes[-1] - tb.tau_of_z(med, z0), tb.tau_of_z(med, z0)) / np.abs(mu)
    for c in (0.05, 0.2, 0.5, 1.0, 2.0):
        reach = t_exit > c                                   # paths that can be observed beyond c
        p = np.mean(t[reach] > c)
        se = np.sqrt(p * (1 - p) / reach.sum())
        assert abs(p - np.exp(-c)) <= 3 * se, (c, p, np.exp(-c), se)
    np.testing.assert_allclose(t[~out], tau[~out], rtol=1e-9, atol=1e-12)
