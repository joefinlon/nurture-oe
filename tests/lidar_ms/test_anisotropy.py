"""Numerical-form tests N1, N2 (appendix C §C.4.3)."""
import numpy as np
import pytest
from scipy.integrate import quad

from lidar_ms.anisotropy import factor_droplet, factor_yang, p_droplet, p_yang, table_factor

GAMMA2 = np.logspace(-6, -1, 10)        # rad^2: ten tabulated coangle variances
THETA2 = (532e-9 / (np.pi * 8e-6)) ** 2


def eq34(p, g2):
    """Direct quadrature of H08 Eq. 34: int_0^inf p(sqrt(g2 t)) exp(-t) dt."""
    val, _ = quad(lambda t: p(np.sqrt(g2 * t)) * np.exp(-t), 0, np.inf, epsabs=1e-13, epsrel=1e-12, limit=400)
    return val


def test_n1_zero_coangle():
    assert factor_droplet(0.0, THETA2) == pytest.approx(1.0, abs=1e-15)
    assert factor_yang(np.array(0.0)) == pytest.approx(1.0, abs=1e-15)


@pytest.mark.parametrize("g2", GAMMA2)
def test_n1_droplet_vs_quadrature(g2):
    ref = eq34(lambda gam: p_droplet(gam / np.sqrt(THETA2)), g2)
    assert factor_droplet(g2, THETA2) == pytest.approx(ref, rel=1e-8)


@pytest.mark.parametrize("g2", GAMMA2)
def test_n1_yang_vs_quadrature(g2):
    ref = eq34(p_yang, g2)
    assert factor_yang(np.array(g2)) == pytest.approx(ref, rel=1e-8)


def test_n1_complex_safe():
    h = 1e-30
    g2 = 1e-3
    for f, df in ((lambda z: factor_yang(z), None), (lambda z: factor_droplet(z, THETA2), None)):
        cs = np.imag(f(np.array(g2 + 1j * h))) / h
        fd = (f(np.array(g2 * (1 + 1e-6))) - f(np.array(g2 * (1 - 1e-6)))) / (2e-6 * g2)
        assert cs == pytest.approx(np.real(fd), rel=1e-6)


def test_n2_table_reproduces_eq36():
    # Node placement is not fixed by appendix C: quadratic spacing in gamma / Theta (dense near the
    # peak) gives 4.4e-5; uniform spacing over 0-10 gives 4.4e-4.
    x = 10.0 * (np.arange(301) / 300.0) ** 2    # gamma / Theta
    p = p_droplet(x)
    ratio = np.logspace(-4, 3, 200)             # gamma2 / Theta^2
    F, tail = table_factor(ratio * THETA2, x, np.broadcast_to(p, (ratio.size, x.size)), THETA2)
    err = np.max(np.abs(F - factor_droplet(ratio * THETA2, THETA2)))
    assert err <= 2e-4, err


def test_table_limits():
    x = np.linspace(0.0, 0.2, 21)
    p = np.linspace(1.0, 0.5, 21)
    F0, _ = table_factor(np.array([0.0]), x, p[None])
    assert F0[0] == 1.0
    Finf, tail = table_factor(np.array([1e6]), x, p[None])
    assert Finf[0] == pytest.approx(0.5, abs=1e-6) and tail[0] > 0.99
