"""Analytical test A1 and regression R1 (appendix C §C.4.3)."""
import numpy as np
import pytest

from lidar_ms import Instrument, OpticalProfile, RangeGrid, SolverConfig, simulate
from lidar_ms.compat import read_ms_input
from msref import REF_DIR, REF_PAIRS, load_ref

HALO = Instrument.from_full_angles(532e-9, 0.8e-3, 1.0e-3)


def slab(alpha=0.02, beta=1e-3, n=60, dr=15.0, r0=9000.0):
    r = r0 + dr * np.arange(n)
    shape = (1, n, 2)
    a = np.zeros(shape)
    a[0, :, 0] = alpha
    b = np.zeros(shape)
    b[0, :, 0] = beta
    opt = OpticalProfile(kinds=("liquid", "molecular"), descriptors=("droplet", "isotropic"), alpha=a,
                         omega=np.ones(shape), g=np.zeros(shape), a_eq=np.full(shape, 8e-6), beta=b,
                         delta=np.zeros(shape))
    return RangeGrid(r), opt


def test_a1_beer_lambert():
    alpha, beta, dr = 0.02, 1e-3, 15.0
    g, opt = slab(alpha, beta, dr=dr)
    out = simulate(g, HALO, opt, SolverConfig(solver="single"))
    i = np.arange(60)
    expect = beta * np.exp(-2 * alpha * dr * i) * (1 - np.exp(-2 * alpha * dr)) / (2 * alpha * dr)
    np.testing.assert_allclose(out.total_by_class[0, :, 0], expect, rtol=1e-12)
    np.testing.assert_allclose(out.beta_ss[0, :, 0], expect, rtol=1e-12)
    assert np.all(out.ms_factor == 0)


SINGLE = [nm for nm, t in REF_PAIRS if t == "single"]


@pytest.mark.parametrize("name", SINGLE)
def test_r1_single_vs_multiscatter(name):
    case = read_ms_input(REF_DIR / f"{name}.in")
    ref = load_ref(name, "single").data[:, 4]
    out = simulate(case.grid, case.instrument, case.optics,
                   SolverConfig(solver="single", compat="multiscatter-1.2.11"))
    got = out.total[0]
    assert np.all(np.isfinite(got))
    ok = np.isfinite(ref)
    np.testing.assert_allclose(got[ok], ref[ok], rtol=1e-7, atol=0)
