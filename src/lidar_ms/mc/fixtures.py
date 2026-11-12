"""Phase-matrix fixtures for tests of solver `mc` (appendix C §C.2.2: test fixtures only).

Henyey-Greenstein, double Henyey-Greenstein, isotropic, Rayleigh and the Gaussian forward-lobe
fixture of Asm4. Every function returns a dict of the six elements P11, P12, P22, P33, P34, P44 on
the angle grid `theta` [rad], normalized so that int P11 dOmega = 4 pi. None of these is a
microphysical model: Mie (liquid) and TAMUice2016 (ice) enter only through Stage 1 [DECISION D18].
"""
import numpy as np

ELEMENTS = ("P11", "P12", "P22", "P33", "P34", "P44")


def angle_grid(theta_min_step, n_tail=2000, n_lobe=400, lobe=0.3):
    """Angle grid from 0 to pi: uniform steps `theta_min_step` up to the lobe region, then geometric to
    `lobe` rad, then uniform in theta to pi. First spacing equals `theta_min_step`."""
    h = float(theta_min_step)
    head = h * np.arange(11)
    mid = np.geomspace(head[-1], lobe, n_lobe)[1:] if lobe > head[-1] else np.array([])
    start = mid[-1] if mid.size else head[-1]
    tail = np.linspace(start, np.pi, n_tail)[1:]
    return np.concatenate([head, mid, tail])


def normalize_p11(theta, p11):
    """Scale p11 so that int_0^pi p11 sin(theta) dtheta = 2 (trapezoid on a refined grid)."""
    return p11 * 2.0 / sphere_integral(theta, p11)


def sphere_integral(theta, f):
    """int_0^pi f sin(theta) dtheta, with f linear in theta between nodes (exact sin weights)."""
    t0, t1 = theta[:-1], theta[1:]
    f0, f1 = f[..., :-1], f[..., 1:]
    h = t1 - t0
    # int_{t0}^{t1} (f0 + (f1 - f0)(t - t0)/h) sin t dt, closed form
    c0, c1, s0, s1 = np.cos(t0), np.cos(t1), np.sin(t0), np.sin(t1)
    i0 = c0 - c1                                     # int sin
    i1 = (s1 - s0 - h * c1) / np.where(h > 0, h, 1)  # int (t - t0)/h sin
    return np.sum(f0 * i0 + (f1 - f0) * i1, axis=-1)


def scalar_matrix(theta, p11, delta_pi=0.0):
    """Complete a P11 with polarization elements of a simple valid family.

    With s = sin^2(theta/2) and d = (1 - delta_pi)/(1 + delta_pi):
    P22 = P11 (1 - (1 - d) s), P33 = P11 (1 - (1 + d) s), P44 = P11 (1 - 2 d s), P12 = P34 = 0.
    Forward (s = 0) the matrix is diagonal; at pi, P33 = -P22 and P44 = P11 - 2 P22 (randomly
    oriented particles), and (P11 - P22)/(P11 + P22) = delta_pi. delta_pi = 0 is sphere-like.
    """
    s = np.sin(0.5 * theta) ** 2
    d = (1.0 - delta_pi) / (1.0 + delta_pi)
    z = np.zeros_like(p11)
    return dict(P11=p11, P12=z, P22=p11 * (1 - (1 - d) * s), P33=p11 * (1 - (1 + d) * s), P34=z.copy(),
                P44=p11 * (1 - 2 * d * s))


def henyey_greenstein_p11(theta, g):
    mu = np.cos(theta)
    return (1 - g * g) / (1 + g * g - 2 * g * mu) ** 1.5


def henyey_greenstein(theta, g, delta_pi=0.0):
    return scalar_matrix(theta, henyey_greenstein_p11(theta, g), delta_pi)


def hg_inverse_cdf(u, g):
    """Closed-form inverse CDF of Henyey-Greenstein in mu = cos(theta); u in [0, 1] from theta = 0."""
    u = np.asarray(u, dtype=float)
    if abs(g) < 1e-8:
        return 1.0 - 2.0 * u
    # CDF from mu = 1 downwards: u = (1 - g^2)/(2g) [1/(1 - g) - 1/t], t = sqrt(1 + g^2 - 2g mu)
    t = 1.0 / (1.0 / (1.0 - g) - 2.0 * g * u / (1.0 - g * g))
    return (1 + g * g - t * t) / (2 * g)


def double_henyey_greenstein(theta, f, g1, g2, delta_pi=0.0):
    p = f * henyey_greenstein_p11(theta, g1) + (1 - f) * henyey_greenstein_p11(theta, g2)
    return scalar_matrix(theta, p, delta_pi)


def isotropic(theta, delta_pi=0.0):
    return scalar_matrix(theta, np.ones_like(theta), delta_pi)


def rayleigh(theta):
    """Rayleigh scattering matrix (no molecular anisotropy)."""
    mu = np.cos(theta)
    p11 = 0.75 * (1 + mu * mu)
    return dict(P11=p11, P12=-0.75 * (1 - mu * mu), P22=p11.copy(), P33=1.5 * mu, P34=np.zeros_like(mu),
                P44=1.5 * mu)


def gaussian_lobe(theta, lobe_width, lobe_fraction=0.5, delta_pi=0.0):
    """Gaussian forward lobe of 1/e half-width `lobe_width` holding `lobe_fraction` of the scattered
    energy (Asm4: 1/2), the rest isotropic. Near 180 deg the phase function is flat (descriptor
    "isotropic")."""
    lobe = np.exp(-(theta / lobe_width) ** 2)
    lobe = lobe * 2.0 / sphere_integral(theta, lobe)
    return scalar_matrix(theta, lobe_fraction * lobe + (1 - lobe_fraction), delta_pi)


def asymmetry(theta, p11):
    """g = <cos theta> of a P11 normalized to int P11 sin dtheta = 2."""
    return sphere_integral(theta, p11 * np.cos(theta)) / sphere_integral(theta, p11)


def stack(*mats):
    """Stack fixtures into (T, A) arrays per element."""
    return {k: np.stack([m[k] for m in mats]) for k in ELEMENTS}
