"""Solver `mc`, scalar mode: X4 (transport), X5-X9 (appendix C §C.4.3, milestone M10).

Monte Carlo tests carry @pytest.mark.slow. Seeds are fixed (McConfig.seed = 0 unless stated).
"""
import os

import numpy as np
import pytest
from scipy.integrate import quad

from lidar_ms import Instrument, McConfig, OpticalProfile, RangeGrid, SolverConfig, simulate
from lidar_ms import grid as lgrid
from lidar_ms.mc import driver as dv
from lidar_ms.mc import fixtures as fx
from lidar_ms.mc import tables as tb
from lidar_ms.mc import transport as tr
from mc_cases import HALO, consistent_beta, lobe_case, mc_config, phase, tau_mid
from test_mc_tables import _layered, _optical_path
from test_pvc import _a5_integral, case

pytestmark = pytest.mark.slow
WORKERS = min(8, os.cpu_count() or 1)


def hg_case(n=40, g=0.85, alpha=0.01, cloud=(5, 30)):
    """HG cloud (not lobe-forming: a_eq 1 um gives Theta = 0.17 rad) over Rayleigh air, HALO geometry."""
    r = 9000.0 + 15.0 * np.arange(n)
    th = fx.angle_grid(1e-3)
    ph = phase(th, [fx.henyey_greenstein(th, g), fx.rayleigh(th)], np.tile([0, 1], (n, 1)))
    a = np.zeros((1, n, 2))
    a[0, cloud[0]:cloud[1], 0] = alpha
    a[0, :, 1] = 1e-5
    gg = np.zeros_like(a)
    gg[..., 0] = g
    ae = np.full_like(a, 1e-6)
    opt = OpticalProfile(kinds=("liquid", "molecular"), descriptors=("isotropic", "isotropic"), alpha=a,
                         omega=np.ones_like(a), g=gg, a_eq=ae, beta=np.zeros_like(a), delta=np.zeros_like(a))
    return RangeGrid(r), consistent_beta(opt, ph), ph


# --- X4 (transport) ------------------------------------------------------------------------------
def test_x4_transport_free_path_exponential():
    med = _layered()
    rng = np.random.default_rng(10)
    n = 200_000
    z0 = rng.uniform(med.r_half[0], med.z_max, n)
    uz = rng.uniform(-1, 1, n)
    uz = np.where(np.abs(uz) < 1e-3, 1e-3, uz)
    s, exits = tr.free_path(med, z0, uz, rng)
    t = np.where(exits, np.inf, _optical_path(med, z0, z0 + s * uz, uz))
    t_exit = np.where(uz > 0, med.tau_nodes[-1] - tb.tau_of_z(med, z0), tb.tau_of_z(med, z0)) / np.abs(uz)
    for c in (0.05, 0.2, 0.5, 1.0, 2.0):
        reach = t_exit > c
        p = np.mean(t[reach] > c)
        se = np.sqrt(p * (1 - p) / reach.sum())
        assert abs(p - np.exp(-c)) <= 3 * se, (c, p, se)


# --- X5 --------------------------------------------------------------------------------------------
def test_x5_order1_equals_single_scattering():
    grid, opt, ph = hg_case()
    out = simulate(grid, HALO, opt, mc_config(400_000), phase=ph)
    ref = simulate(grid, HALO, opt, SolverConfig(solver="single"))
    tau = tau_mid(grid, opt)
    o1 = out.beta_order[0, :, :, 0].sum(-1)
    run = dv.run_profile(grid.r[0], opt.alpha[0], opt.omega[0], opt.a_eq[0], opt.kinds, ph, HALO,
                         McConfig(n_photons=400_000, chunk=50_000, polarized=False), 0.1)
    se = run["batches"][:, :, :, 0].sum(axis=(2, 3)).std(axis=0, ddof=1) / np.sqrt(20)
    np.testing.assert_allclose(run["beta"][:, :, 0].sum(axis=(1, 2)), o1, rtol=1e-12)
    m = tau <= 3
    dev = np.abs(o1 - ref.beta_ss[0].sum(-1))[m]
    assert np.all(dev <= 3 * se[m]), np.max(dev / se[m])
    np.testing.assert_allclose(out.beta_ss[0], out.beta_order[0, :, :, 0])


# --- X6 --------------------------------------------------------------------------------------------
def test_x6_time_of_flight():
    grid, opt, ph = hg_case()
    r_half = lgrid.half_levels(grid.r)[0]
    med = tb.build_medium(r_half, opt.alpha[0], opt.omega[0], ph.type_index)
    tabs = dv.build_tables(ph, [0, 1])
    rec = []
    N, J = med.sigma.shape[0] - 1, med.sigma.shape[1]
    dv.run_chunk(50_000, med, tabs, HALO, McConfig(polarized=False), np.random.default_rng(11),
                 np.zeros((N, J, 3, 2)), np.zeros(N), record=rec)
    order = np.concatenate([g["order"] for g in rec])
    r_app = np.concatenate([g["r_app"] for g in rec])
    dist = np.concatenate([g["dist"] for g in rec])
    bins = np.concatenate([g["bin"] for g in rec])
    assert np.all(r_app >= dist)
    one = order == 1
    assert one.sum() > 10_000 and (order > 1).sum() > 10_000
    np.testing.assert_array_equal(bins[one], np.searchsorted(r_half, dist[one], side="right") - 1)


# --- X7 --------------------------------------------------------------------------------------------
def test_x7_ddis_unbiased():
    """epsilon = 0 and 0.1 agree per bin and order. Without DDIS the local estimate is heavy-tailed
    (forward peak seen from below), so the batch errors need 4e6 packets to be trustworthy."""
    grid, opt, ph = hg_case()
    runs = [dv.run_profile(grid.r[0], opt.alpha[0], opt.omega[0], opt.a_eq[0], opt.kinds, ph, HALO,
                           McConfig(n_photons=4_000_000, chunk=100_000, polarized=False, epsilon_ddis=eps, seed=s,
                                    workers=WORKERS), 0.1) for eps, s in ((0.0, 1), (0.1, 2))]
    m = tau_mid(grid, opt) <= 3
    tot = [r["beta"].sum(axis=(1, 2, 3)) for r in runs]
    se = [r["batches"].sum(axis=(2, 3, 4)).std(axis=0, ddof=1) / np.sqrt(20) for r in runs]
    dev = np.abs(tot[0] - tot[1])[m]
    lim = 3 * np.hypot(se[0], se[1])[m]
    assert np.all(dev <= lim), np.max(dev / lim)
    for o in range(3):                                      # each order, summed over the tau <= 3 bins
        tot = [r["beta"][m][:, :, o].sum() for r in runs]
        se = [r["batches"][:, m][..., o, :].sum(axis=(1, 2, 3)).std(ddof=1) / np.sqrt(20) for r in runs]
        assert abs(tot[0] - tot[1]) <= 3 * np.hypot(*se), (o, tot, se)


@pytest.mark.parametrize("eps", [0.1, 0.5])
def test_x7_mean_weight_factor_is_one(eps):
    grid, inst, opt, ph, _ = lobe_case("halo_thin_slab", iso_fraction=0.02)
    med = tb.build_medium(lgrid.half_levels(grid.r)[0], opt.alpha[0], opt.omega[0], ph.type_index)
    tabs = dv.build_tables(ph, [0, 1])
    rng = np.random.default_rng(12)
    n = 400_000
    z = rng.uniform(med.r_half[20], med.r_half[40], n)
    x = np.stack([rng.normal(0, 2, n), rng.normal(0, 2, n), z])
    u = tr.launch(n, 0.3, rng)["u"]                         # broad set of directions
    u = np.where(rng.random(n) < 0.5, u, -u)
    p = dict(x=x, u=u, e=tr.geo.perpendicular_axis(np.array([1.0, 0, 0]), u), S=None, W=np.ones(n),
             L=z, order=np.ones(n, np.int64))
    layer = tb.layer_of_z(med, z)
    d = -x / tr.geo.norm(x)
    _, w = tr.redirect(p, med, tabs, layer, d, eps, rng)
    assert abs(w.mean() - 1) <= 3 * w.std(ddof=1) / np.sqrt(n), (w.mean(), w.std())


# --- X8 --------------------------------------------------------------------------------------------
def test_x8_thin_slab_double_scattering():
    """Order-2 / order-1 per slab bin vs the bin mean of the A5 integral.

    epsilon_ddis = 0.5 and 4e8 packets: the return-lobe half of the double scattering arrives through
    rare DDIS events (a packet redirected towards the detector must collide again inside the slab);
    see docs/findings_for_pk.md (2026-11-12)."""
    grid, inst, opt, ph, _ = lobe_case("halo_thin_slab", iso_fraction=0.02)
    c = case("halo_thin_slab")
    a = c.optics.alpha[0, :, :3].sum(-1)
    idx = np.flatnonzero(a > 0)
    redge = c.grid.r_half[0, idx[0]]
    mbar = []
    for i in idx:
        lo, hi = c.grid.r_half[0, i] - redge, c.grid.r_half[0, i + 1] - redge
        mbar.append(quad(lambda d: _a5_integral(c, d), lo, hi, epsrel=1e-10)[0] / (hi - lo))
    mbar = np.array(mbar)
    mc = McConfig(n_photons=400_000_000, chunk=200_000, polarized=False, epsilon_ddis=0.5, workers=WORKERS)
    run = dv.run_profile(grid.r[0], opt.alpha[0], opt.omega[0], opt.a_eq[0], opt.kinds, ph, inst, mc, 0.1)
    b = run["batches"][:, idx].sum(axis=(2, 4))                 # (batch, bin, order)
    ratio = run["beta"][idx].sum(axis=(1, 3))
    ratio = ratio[:, 1] / ratio[:, 0]
    se = (b[:, :, 1] / b[:, :, 0]).std(axis=0, ddof=1) / np.sqrt(b.shape[0])
    ok = (np.abs(ratio - mbar) <= 3 * se) | (np.abs(ratio / mbar - 1) <= 0.03)
    assert ok.all(), (ratio / mbar - 1)[~ok]
    tot = run["beta"][idx].sum(axis=(1, 3))
    assert tot[:, 1].sum() / (tot[:, 0] * mbar).sum() == pytest.approx(1.0, abs=0.01)


# --- X9 --------------------------------------------------------------------------------------------
def test_x9_wide_field_of_view():
    """Total return vs A3 near cloud top (tau <= 1), with the lobe fixture whose non-lobe light is 2 %
    isotropic and the rest absorbed (A3's exact limit). With a fully isotropic non-lobe half the Monte
    Carlo exceeds A3 by wide-angle light (docs/findings_for_pk.md, 2026-11-12)."""
    wide = Instrument(532e-9, HALO.rho_tr, 0.5)
    grid, _, opt, ph, _ = lobe_case("halo_liquid_30_iso", instrument=wide, iso_fraction=0.02)
    run = dv.run_profile(grid.r[0], opt.alpha[0], opt.omega[0], opt.a_eq[0], opt.kinds, ph, wide,
                         McConfig(n_photons=1_000_000, chunk=50_000, polarized=False), 0.1)
    al, dr = opt.alpha[0], grid.dr[0]
    ac, an = al[:, :3].sum(-1), al[:, 3]
    tp = (ac + 2 * an) * dr
    T = np.exp(-np.concatenate([[0.0], np.cumsum(tp)[:-1]]))
    expect = opt.beta[0].sum(-1) * T * (-np.expm1(-tp)) / tp
    tot = run["beta"].sum(axis=(1, 2, 3))
    se = run["batches"].sum(axis=(2, 3, 4)).std(axis=0, ddof=1) / np.sqrt(20)
    m = (ac > 0) & (tau_mid(grid, opt) <= 1)
    assert m.sum() >= 2
    ok = (np.abs(tot - expect) <= 3 * se) | (np.abs(tot / expect - 1) <= 0.05)
    assert ok[m].all(), (tot / expect - 1)[m]


# --- plumbing --------------------------------------------------------------------------------------
def test_mc_returns_and_reproducibility():
    grid, opt, ph = hg_case(n=20, cloud=(5, 12))
    cfg = mc_config(40_000, chunk=10_000)
    a = simulate(grid, HALO, opt, cfg, phase=ph)
    b = simulate(grid, HALO, opt, SolverConfig(solver="mc", mc=McConfig(n_photons=40_000, chunk=10_000,
                                                                          polarized=False, workers=2)), phase=ph)
    np.testing.assert_array_equal(a.beta_par, b.beta_par)
    assert list(a.closure_used) == ["mc_scalar"] and np.all(a.beta_perp == 0)
    assert a.beta_order.shape == (1, 20, 2, 3) and a.stderr_par.shape == (1, 20, 2)
    np.testing.assert_allclose(a.beta_order.sum(-1), a.beta_par, rtol=1e-12)
    assert np.all(np.isnan(a.ms_factor)) and np.all(a.beta_wide == 0)
    with pytest.raises(ValueError, match="PhaseTables"):
        simulate(grid, HALO, opt, cfg)
    th = fx.angle_grid(5e-3)
    coarse = phase(th, [fx.isotropic(th), fx.rayleigh(th)], ph.type_index)
    lobed = opt.replace(a_eq=np.full(opt.shape, 8e-6))
    with pytest.raises(ValueError, match="forward lobe"):
        simulate(grid, HALO, lobed, cfg, phase=coarse)
    with pytest.warns(UserWarning, match="differs"):
        simulate(grid, HALO, opt.replace(beta=opt.beta * 1.1), cfg, phase=ph)
