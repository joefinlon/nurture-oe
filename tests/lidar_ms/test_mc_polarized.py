"""Solver `mc`, polarized mode: P1-P4 (appendix C §C.2.7, §C.4.3, milestone M11).

P2 compares with an independent Jones-vector quadrature (`jones_double.py`); "P2-Mie" repeats it for
one Mie sphere size (an extra check of the polarized Mie path). P4 needs the cases of Hu et al. (2006,
Table 1), which are not in Project Knowledge: it is skipped, and a labelled proxy at the shipped
CALIPSO-like fixture geometry records how far the Monte Carlo lies from Hu's Eq. 3 there.
"""
import os
from functools import lru_cache

import numpy as np
import pytest

import mie_fixture as mf
from jones_double import double_scatter, tabulated_amplitude
from lidar_ms import Instrument, McConfig, OpticalProfile, PhaseTables, RangeGrid, simulate
from lidar_ms.mc import fixtures as fx
from lidar_ms.mc.driver import run_profile
from mc_cases import HALO, consistent_beta, mc_config, phase
from msref import load_input

pytestmark = pytest.mark.slow
WORKERS = min(8, os.cpu_count() or 1)
DELTA_ICE = 0.4


def droplet_ice_case(n=30):
    """HALO geometry: Mie droplets (a_eq 8 um), a non-lobe depolarizing HG class (delta = 0.4), Rayleigh air."""
    d = mf.droplet_optics()
    th = d["theta"]
    mats = [{k: d[k] for k in fx.ELEMENTS}, fx.henyey_greenstein(th, 0.75, delta_pi=DELTA_ICE), fx.rayleigh(th)]
    r = 9000.0 + 15.0 * np.arange(n)
    a = np.zeros((1, n, 3))
    a[0, 8:20, 0] = 0.02
    a[0, 4:12, 1] = 0.002
    a[0, :, 2] = 1e-5
    ae = np.full_like(a, 8e-6)
    ae[..., 1] = 1e-6                                       # Theta = 0.17 rad: not lobe-forming
    opt = OpticalProfile(kinds=("liquid", "ice", "molecular"), descriptors=("droplet", "isotropic", "isotropic"),
                         alpha=a, omega=np.ones_like(a), g=np.zeros_like(a), a_eq=ae, beta=np.zeros_like(a),
                         delta=np.zeros_like(a))
    ph = PhaseTables(theta=th, type_index=np.tile([0, 1, 2], (n, 1)), **fx.stack(*mats))
    return RangeGrid(r), consistent_beta(opt, ph), ph


def rayleigh_slab():
    r = 505.0 + 10.0 * np.arange(150)
    a = np.zeros((r.size, 2))
    a[(r > 1000) & (r < 1300), 1] = 3e-4
    th = fx.angle_grid(1e-3)
    return r, a, phase(th, [fx.isotropic(th), fx.rayleigh(th)], np.tile([0, 1], (r.size, 1)))


PENCIL = Instrument(532e-9, 1e-6, 0.1)


def _order2_totals(run, dr):
    tot = (run["batches"][:, :, :, 1, :] * dr).sum(axis=(1, 2))      # (batch, pol)
    return tot.mean(axis=0), tot.std(axis=0, ddof=1) / np.sqrt(tot.shape[0]), tot


def _check_order2(run, dr, co, cx):
    m, se, tot = _order2_totals(run, dr)
    assert abs(m[0] - co) <= 3 * se[0], (m[0] / co - 1, se[0] / co)
    assert abs(m[1] - cx) <= 3 * se[1], (m[1] / cx - 1, se[1] / cx)
    ratio = tot[:, 1] / tot[:, 0]
    assert abs(m[1] / m[0] - cx / co) <= 3 * ratio.std(ddof=1) / np.sqrt(ratio.size)


# --- P1 ----------------------------------------------------------------------------------------------
def test_p1_order1_single_scatter_depolarization():
    grid, opt, ph = droplet_ice_case()
    out = simulate(grid, HALO, opt, mc_config(200_000, polarized=True), phase=ph)
    run = run_profile(grid.r[0], opt.alpha[0], opt.omega[0], opt.a_eq[0], opt.kinds, ph, HALO,
                      McConfig(n_photons=200_000, chunk=50_000), 0.1)
    o1 = run["beta"][:, :, 0, :]                                      # (N, J, pol)
    ice = opt.alpha[0, :, 1] > 0
    np.testing.assert_allclose(o1[ice, 1, 1] / o1[ice, 1, 0], DELTA_ICE, rtol=1e-10)
    liq = opt.alpha[0, :, 0] > 0
    assert np.all(o1[liq, 0, 0] > 0)
    assert np.all(np.abs(o1[liq, 0, 1]) <= 1e-12 * o1[liq, 0, 0])   # spheres: zero
    assert np.all(np.abs(o1[:, 2, 1]) <= 1e-12 * o1[:, 2, 0])       # Rayleigh without anisotropy: zero
    assert list(out.closure_used) == ["mc"]
    assert np.any(out.beta_perp[0, liq, 0] > 0)                      # multiple scattering depolarizes droplets
    np.testing.assert_allclose(out.beta_order.sum(-1), out.total_by_class, rtol=1e-12)


# --- P2 ----------------------------------------------------------------------------------------------
def test_p2_rayleigh_double_scattering():
    r, a, ph = rayleigh_slab()
    co, cx = double_scatter(1000.0, 1300.0, 3e-4, 0.1)
    run = run_profile(r, a, np.ones_like(a), np.zeros_like(a), ("other", "molecular"), ph, PENCIL,
                      McConfig(n_photons=8_000_000, chunk=200_000, workers=WORKERS), 0.1)
    _check_order2(run, 10.0, co, cx)


def test_p2_mie_sphere_double_scattering():
    """Extra check: one Mie sphere size (a = 2.5 um), amplitudes S1, S2 in the Jones quadrature."""
    a_eq = 2.5e-6
    th = np.linspace(0, np.pi, 20001)
    o = mf.single(mf.M_WATER_532, 2 * np.pi / 532e-9 * a_eq, th)
    sc = np.sqrt(2.0 / fx.sphere_integral(th, 0.5 * (np.abs(o["S1"]) ** 2 + np.abs(o["S2"]) ** 2)))
    s1, s2 = o["S1"] * sc, o["S2"] * sc
    m = dict(P11=0.5 * (abs(s1) ** 2 + abs(s2) ** 2), P12=0.5 * (abs(s2) ** 2 - abs(s1) ** 2),
             P33=np.real(s2 * np.conj(s1)), P34=np.imag(s2 * np.conj(s1)))
    m["P22"], m["P44"] = m["P11"], m["P33"]
    co, cx = double_scatter(1000.0, 1300.0, 3e-4, 0.1, amplitude=tabulated_amplitude(th, s1, s2))
    r, a, _ = rayleigh_slab()
    a = a[:, ::-1].copy()                                             # the slab is the liquid class
    ph = phase(th, [m, fx.rayleigh(th)], np.tile([0, 1], (r.size, 1)))
    run = run_profile(r, a, np.ones_like(a), np.full_like(a, a_eq), ("liquid", "molecular"), ph, PENCIL,
                      McConfig(n_photons=40_000_000, chunk=200_000, workers=WORKERS, epsilon_ddis=0.5), 0.1)
    _check_order2(run, 10.0, co, cx)


# --- P3 ----------------------------------------------------------------------------------------------
def _paired(args, mc):
    runs = [run_profile(*args, mc, 0.1, pol_angle=ang) for ang in (0.0, np.pi / 4)]
    diff = runs[0]["batches"] - runs[1]["batches"]
    return runs, diff


@pytest.mark.parametrize("which", ["rayleigh_slab", "droplets_halo"])
def test_p3_rotated_polarization(which):
    if which == "rayleigh_slab":
        r, a, ph = rayleigh_slab()
        args = (r, a, np.ones_like(a), np.zeros_like(a), ("other", "molecular"), ph, PENCIL)
        mc = McConfig(n_photons=4_000_000, chunk=200_000, workers=WORKERS)
    else:
        grid, opt, ph = droplet_ice_case()
        args = (grid.r[0], opt.alpha[0], opt.omega[0], opt.a_eq[0], opt.kinds, ph, HALO)
        mc = McConfig(n_photons=1_000_000, chunk=100_000, workers=WORKERS)
    runs, diff = _paired(args, mc)
    nb = diff.shape[0]
    scale = np.abs(runs[0]["batches"]).sum(axis=0) / nb
    for axes in ((1, 2), (0, 1)):                            # of (bin, class, order, pol)
        ax = tuple(x + 1 for x in axes)
        d = diff.sum(axis=ax)
        m = d.mean(axis=0)
        se = d.std(axis=0, ddof=1) / np.sqrt(nb)
        floor = 1e-12 * scale.sum(axis=axes).sum(axis=-1, keepdims=True)   # order 1: identical up to rounding
        lim = family_limit(m.size, nb) * se + floor
        assert np.all(np.abs(m) <= lim), (axes, np.max(np.abs(m) / lim))
    assert runs[0]["beta"][..., 1].sum() > 0


def family_limit(n_tests, n_batches):
    """Multiplier on the batch SE with the same family-wise false-alarm rate (0.27 %) as one 3-SE test,
    for n_tests comparisons (Bonferroni, Student t with n_batches - 1 dof)."""
    from scipy.stats import norm, t
    p = 2 * norm.sf(3.0) / n_tests
    return t.isf(p / 2, n_batches - 1)


# --- P4 ----------------------------------------------------------------------------------------------
def hu_eq3(d):
    """Hu et al. (2006) Eq. 3 as given in appendix C §C.2.6."""
    return 0.999 - 3.906 * d + 6.263 * d ** 2 - 3.554 * d ** 3


@pytest.mark.skip(reason="P4 needs the cases of Hu et al. (2006, Table 1), not in Project Knowledge; "
                         "recorded in docs/findings_for_pk.md (2026-11-12)")
def test_p4_hu_table1():
    pass


@lru_cache(maxsize=1)
def _p4_proxy():
    """A_s / Eq. 3 - 1 and delta_acc in the cloud bins, Mie droplets (a_eq 8 um, 20 km-1) at the shipped
    CALIPSO-like fixture geometry (700 km, rho_tr 50 urad, rho_fov 65 urad); not F-Lidar [GAP G6]."""
    hdr = load_input("calipso_cirrus_aerosol_wa_example")
    inst = Instrument(532e-9, hdr.rho_tr, hdr.rho_fov[0])
    d = mf.droplet_optics()
    th = d["theta"]
    N, dr = 60, 15.0
    ph = PhaseTables(theta=th, type_index=np.tile([0, 1], (N, 1)),
                     **fx.stack({k: d[k] for k in fx.ELEMENTS}, fx.rayleigh(th)))
    r = hdr.altitude - 2000.0 + dr * np.arange(N)
    a = np.zeros((N, 2))
    a[5:, 0] = 0.02
    a[:, 1] = 1e-5
    run = run_profile(r, a, np.ones_like(a), np.full_like(a, 8e-6), ("liquid", "molecular"), ph, inst,
                      McConfig(n_photons=4_000_000, chunk=50_000, epsilon_ddis=0.5, workers=WORKERS), 0.1)
    b = run["beta"][5:, 0]                                            # cloud bins: (N, order, pol)
    i_t = np.cumsum(b.sum(axis=(1, 2)) * dr)
    i_s = np.cumsum(b[:, 0].sum(-1) * dr)
    dacc = np.cumsum(b[..., 1].sum(-1) * dr) / np.cumsum(b[..., 0].sum(-1) * dr)
    return i_s / i_t / hu_eq3(dacc) - 1, dacc


def test_p4_proxy_offset_recorded():
    """The proxy's offset from Eq. 3 as recorded in docs/findings_for_pk.md, so that a change is noticed."""
    rel, dacc = _p4_proxy()
    assert np.all(np.isfinite(rel)) and np.all(dacc < 0.45)
    deep = dacc >= 0.1
    assert deep.sum() >= 40
    assert np.all((rel[deep] > -0.10) & (rel[deep] < -0.04)), rel[deep]


@pytest.mark.xfail(strict=True, reason="P4 proxy (not P4): at the CALIPSO-like fixture geometry A_s lies 4-8 % below "
                                       "Hu's Eq. 3 (8 um droplets, 20 km-1); docs/findings_for_pk.md (2026-11-12)")
def test_p4_proxy_within_two_percent():
    rel, _ = _p4_proxy()
    assert np.max(np.abs(rel)) <= 0.02, rel
