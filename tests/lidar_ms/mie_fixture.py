"""Test-only Lorenz-Mie routine for `lidar_ms` [DECISION D32].

Bohren & Huffman (1983) series (BHMIE): logarithmic derivative by downward recurrence, Riccati-Bessel
functions by upward recurrence, angular functions pi_n, tau_n. Vectorized over radii (coefficients)
and angles (matrix products). Never imported by package code; Mie optics enter `lidar_ms` only
through Stage 1 in production [DECISION D18].

Phase-matrix convention (Bohren & Huffman, Mishchenko): P11 = (|S1|^2 + |S2|^2)/2,
P12 = (|S2|^2 - |S1|^2)/2, P33 = P44 = Re(S2 S1*), P34 = Im(S2 S1*), P22 = P11; all normalized by
4 pi / (k^2 C_sca) so that int P11 dOmega = 4 pi. Rayleigh limit: P12 = -3/4 sin^2, P33 = 3/2 cos.

The refractive index of water used by the tests (1.334 + 0i at 532 nm; Hale & Querry 1973 give
n = 1.334 there and k ~ 1e-9) is a fixture value passed explicitly, not a database specification.
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np

from lidar_ms import NearBackscatter, PhaseTables
from lidar_ms.mc.fixtures import sphere_integral

M_WATER_532 = 1.334 + 0.0j


def n_stop(x):
    return int(np.max(np.asarray(x) + 4.0 * np.asarray(x) ** (1.0 / 3.0) + 2.0))


def coefficients(m, x):
    """Mie coefficients a_n, b_n, shape (R, nmax), zero beyond each radius's n_stop.

    Parameters
    ----------
    m : complex refractive index (relative)
    x : (R,) size parameters
    """
    x = np.atleast_1d(np.asarray(x, dtype=float))
    nmax = n_stop(x)
    stops = (x + 4.0 * x ** (1.0 / 3.0) + 2.0).astype(int)
    mx = m * x
    nmx = int(max(nmax, np.max(np.abs(mx)))) + 16
    D = np.zeros((x.size, nmx + 1), dtype=complex)
    for n in range(nmx, 0, -1):
        D[:, n - 1] = n / mx - 1.0 / (D[:, n] + n / mx)
    psi0, psi1 = np.cos(x), np.sin(x)
    chi0, chi1 = -np.sin(x), np.cos(x)
    xi1 = psi1 - 1j * chi1
    a = np.zeros((x.size, nmax), dtype=complex)
    b = np.zeros((x.size, nmax), dtype=complex)
    with np.errstate(over="ignore", invalid="ignore"):
        for n in range(1, nmax + 1):
            psi = (2 * n - 1) / x * psi1 - psi0
            chi = (2 * n - 1) / x * chi1 - chi0
            xi = psi - 1j * chi
            Dn = D[:, n]
            ta = (Dn / m + n / x)
            tb = (m * Dn + n / x)
            an = (ta * psi - psi1) / (ta * xi - xi1)
            bn = (tb * psi - psi1) / (tb * xi - xi1)
            keep = n <= stops
            a[:, n - 1] = np.where(keep, an, 0.0)
            b[:, n - 1] = np.where(keep, bn, 0.0)
            psi0, psi1, chi0, chi1, xi1 = psi1, psi, chi1, chi, xi
    return a, b


def angular_functions(mu, nmax):
    """pi_n(mu), tau_n(mu) for n = 1..nmax, shape (nmax, A)."""
    mu = np.asarray(mu, dtype=float)
    pi = np.zeros((nmax, mu.size))
    tau = np.zeros((nmax, mu.size))
    p0, p1 = np.zeros_like(mu), np.ones_like(mu)
    for n in range(1, nmax + 1):
        if n > 1:
            p0, p1 = p1, ((2 * n - 1) / (n - 1)) * mu * p1 - (n / (n - 1)) * p0
        pi[n - 1] = p1
        tau[n - 1] = n * mu * p1 - (n + 1) * p0
    return pi, tau


def amplitudes(a, b, mu):
    """S1, S2 (R, A)."""
    nmax = a.shape[1]
    n = np.arange(1, nmax + 1)
    fn = (2 * n + 1) / (n * (n + 1))
    pi, tau = angular_functions(mu, nmax)
    return (a * fn) @ pi + (b * fn) @ tau, (a * fn) @ tau + (b * fn) @ pi


def efficiencies(a, b, x):
    """Q_ext, Q_sca, g (R,) from the series."""
    nmax = a.shape[1]
    n = np.arange(1, nmax + 1)
    qe = 2.0 / x ** 2 * np.sum((2 * n + 1) * np.real(a + b), axis=1)
    qs = 2.0 / x ** 2 * np.sum((2 * n + 1) * (np.abs(a) ** 2 + np.abs(b) ** 2), axis=1)
    t1 = np.sum((n[:-1] * (n[:-1] + 2) / (n[:-1] + 1)) *
                np.real(a[:, :-1] * np.conj(a[:, 1:]) + b[:, :-1] * np.conj(b[:, 1:])), axis=1)
    t2 = np.sum((2 * n + 1) / (n * (n + 1)) * np.real(a * np.conj(b)), axis=1)
    g = 4.0 / (x ** 2 * qs) * (t1 + t2)
    return qe, qs, g


def single(m, x, theta):
    """One sphere: dict with S1, S2 (A,), qext, qsca, g."""
    a, b = coefficients(m, [x])
    s1, s2 = amplitudes(a, b, np.cos(theta))
    qe, qs, g = efficiencies(a, b, np.array([x]))
    return dict(S1=s1[0], S2=s2[0], qext=qe[0], qsca=qs[0], g=g[0], a=a[0], b=b[0])


def gamma_weights(radii, r_eff, v_eff):
    """Hansen & Travis (1974) gamma distribution n(r) ~ r^((1-3v)/v) exp(-r/(r_eff v)), times dr."""
    lw = (1 - 3 * v_eff) / v_eff * np.log(radii) - radii / (r_eff * v_eff)
    w = np.exp(lw - lw.max()) * np.gradient(radii)
    return w / w.sum()


def distribution(m, wavelength, radii, weights, theta, block=1500):
    """Size-averaged optics (radii processed in blocks to bound memory).

    Returns
    -------
    dict: P11, P12, P22, P33, P34, P44 (A,) normalized by the series C_sca (int P11 dOmega = 4 pi up to
    the angular quadrature); c_ext, c_sca [m2, weight-averaged]; omega, g, lidar_ratio [sr] (exact,
    from the series at pi), a_eq [m] (pi a_eq^2 = <G>), r_eff [m].
    """
    k = 2 * np.pi / wavelength
    G = np.pi * radii ** 2
    acc = np.zeros((4, theta.size))
    c_ext = c_sca = g_sum = 0.0
    s_back = 0.0
    for i0 in range(0, radii.size, block):
        sl = slice(i0, i0 + block)
        x = k * radii[sl]
        a, b = coefficients(m, x)
        qe, qs, g = efficiencies(a, b, x)
        s1, s2 = amplitudes(a, b, np.cos(theta))
        w = weights[sl]
        c_ext += np.sum(w * qe * G[sl])
        c_sca += np.sum(w * qs * G[sl])
        g_sum += np.sum(w * qs * G[sl] * g)
        i1, i2 = np.abs(s1) ** 2, np.abs(s2) ** 2
        s21 = s2 * np.conj(s1)
        acc += np.stack([w @ (0.5 * (i1 + i2)), w @ (0.5 * (i2 - i1)), w @ np.real(s21), w @ np.imag(s21)])
        n = np.arange(1, a.shape[1] + 1)
        sb = np.abs(np.sum((2 * n + 1) * (-1.0) ** n * (a - b), axis=1)) ** 2 / 4.0   # S11(pi)
        s_back += np.sum(w * sb)
    norm = 4 * np.pi / (k ** 2 * c_sca)
    P11, P12, P33, P34 = norm * acc
    omega = c_sca / c_ext
    return dict(P11=P11, P12=P12, P22=P11.copy(), P33=P33, P34=P34, P44=P33.copy(), c_ext=c_ext, c_sca=c_sca,
                omega=omega, g=g_sum / c_sca, lidar_ratio=k ** 2 * c_ext / s_back,
                a_eq=np.sqrt(np.sum(weights * radii ** 2)),
                r_eff=np.sum(weights * radii ** 3) / np.sum(weights * radii ** 2))


def droplet_radii(r_eff, v_eff, wavelength, dx=0.1, span=(0.05, 3.5)):
    """Radius nodes with size-parameter step dx over span x r_eff."""
    k = 2 * np.pi / wavelength
    x = np.arange(span[0] * r_eff * k, span[1] * r_eff * k, dx)
    return x / k


def droplet_angle_grid(lobe_width, n_fine=400, n_tail=1600):
    """Angle grid: steps of lobe_width / 20 through the forward lobe, uniform in theta beyond, and a
    fine block near 180 deg for the near-backscatter tables."""
    h = lobe_width / 20.0
    head = np.arange(0.0, 8 * lobe_width, h)
    back = np.pi - np.arange(0.0, 0.35, 5e-4)[::-1]
    mid = np.linspace(head[-1], back[0], n_tail)[1:-1]
    return np.concatenate([head, mid, back])


@lru_cache(maxsize=8)
def droplet_optics(a_eq=8e-6, v_eff=0.1, wavelength=532e-9, m=M_WATER_532, dx=0.02):
    """Gamma-distributed water droplets whose equivalent-area radius equals `a_eq`.

    The Hansen-Travis distribution with effective variance v has <r^2> = r_eff^2 (1 - v)(1 - 2v), so
    r_eff = a_eq / sqrt((1 - v)(1 - 2v)).

    The phase-matrix elements are rescaled together so that int P11 sin dtheta = 2 exactly on the
    returned grid (P11 linear in theta); the change is a few 1e-4 and keeps Monte Carlo tables
    consistent. `lidar_ratio` is the exact series value.
    """
    r_eff = a_eq / np.sqrt((1 - v_eff) * (1 - 2 * v_eff))
    radii = droplet_radii(r_eff, v_eff, wavelength, dx=dx)
    w = gamma_weights(radii, r_eff, v_eff)
    theta = droplet_angle_grid(wavelength / (np.pi * a_eq))
    out = distribution(m, wavelength, radii, w, theta)
    scale = 2.0 / sphere_integral(theta, out["P11"])
    for key in ("P11", "P12", "P22", "P33", "P34", "P44"):
        out[key] = out[key] * scale
    out["theta"] = theta
    out["lobe_width"] = wavelength / (np.pi * a_eq)
    return out


def near_backscatter_table(optics, x_max=10.0, n=301):
    """NearBackscatter in gamma / Theta (droplets, H08) from a Mie phase matrix, Eq. (C.7)."""
    theta, lobe = optics["theta"], optics["lobe_width"]
    x = x_max * (np.arange(n) / (n - 1)) ** 2
    th = np.pi - x * lobe
    rev = theta[::-1]
    f = {k: np.interp(-th, -rev, optics[k][::-1]) for k in ("P11", "P22", "P33")}
    p0 = optics["P11"][-1]
    c_par = (2 * f["P11"] + f["P22"] - f["P33"]) / (4 * p0)
    c_perp = (2 * f["P11"] - f["P22"] + f["P33"]) / (4 * p0)
    return NearBackscatter("gamma_over_theta", x, c_par, c_perp)


def phase_tables(mats, type_index):
    """PhaseTables from a list of element dicts on a common theta grid."""
    keys = ("P11", "P12", "P22", "P33", "P34", "P44")
    return PhaseTables(theta=mats[0]["theta"], type_index=np.asarray(type_index),
                       **{k: np.stack([mm[k] for mm in mats]) for k in keys})
