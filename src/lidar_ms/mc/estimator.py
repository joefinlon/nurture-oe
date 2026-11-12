"""Local estimate, field-of-view test, time of flight and tallies for solver `mc` (appendix C §C.2.7).

At every collision, before redirection, each class j contributes

    Delta_j = W (alpha_j omega_j / alpha) (1 / 4 pi) L(psi2) P_j(Theta_d) L(psi) S
              exp(-tau_z(z) d / z) r_app^2 / d^2  1[sqrt(x^2 + y^2) <= rho_fov z],   r_app = (L + d) / 2,

tallied in the range bin containing r_app, by class, scattering order (1, 2, >= 3) and polarization
(co: (I + Q)/2, cross: (I - Q)/2 with Q referred to the analyzer). Scalar mode tallies I as co-polarized.
"""
import numpy as np

from . import geometry as geo
from . import tables as tb
from .transport import apply_matrix

N_ORDERS = 3


def detector_direction(x):
    """d = -x / |x| and d = |x|."""
    dist = geo.norm(x)
    return -x / dist, dist


def fov_weight(x, rho_fov, receiver):
    """Receiver acceptance of light arriving from x: top hat (rho <= rho_fov z) or Gaussian."""
    rho2 = x[0] ** 2 + x[1] ** 2
    lim2 = (rho_fov * x[2]) ** 2
    if receiver == "gaussian":
        return np.exp(-rho2 / lim2)
    return (rho2 <= lim2).astype(float)


def analyzer(d, pol_angle):
    """Analyzer axis p (x-hat projected perpendicular to d) and the double angle (cos 2eta, sin 2eta) of
    the transmitted polarization direction relative to p, about d (eta = 0 for pol_angle = 0)."""
    p = geo.perpendicular_axis(np.array([1.0, 0.0, 0.0]), d)
    if pol_angle == 0.0:
        return p, None
    a = geo.perpendicular_axis(np.array([np.cos(pol_angle), np.sin(pol_angle), 0.0]), d)
    c = geo.dot(p, a)
    s = geo.dot(geo.cross(p, a), d)
    return p, geo.double_angle(c, s)


def local_estimate(p, med, tabs, layer, rho_fov, receiver, polarized=False, pol_angle=0.0):
    """Contributions of the packets at their current collision.

    Returns
    -------
    d : (3, n) detector direction (also used by DDIS)
    bin_ : (n,) int  range bin of r_app (-1 or N where outside the grid)
    co, cross : (J, n) contributions (cross is None in scalar mode)
    r_app, dist : (n,)
    """
    x, u, e = p["x"], p["u"], p["e"]
    d, dist = detector_direction(x)
    r_app = 0.5 * (np.maximum(p["L"], dist) + dist)          # L >= |x| (triangle inequality); clip rounding
    bin_ = np.searchsorted(med.r_half, r_app, side="right") - 1
    cos_t, sin_t, c2, s2, l_out = geo.scattering_plane(u, e, d)
    th = np.arctan2(sin_t, cos_t)
    i, f = tb.interp_index(tabs.theta, th)
    z = x[2]
    base = (p["W"] / (4 * np.pi * med.alpha[layer]) * np.exp(-tb.tau_of_z(med, z) * dist / z)
            * (r_app / dist) ** 2 * fov_weight(x, rho_fov, receiver))
    sig = med.sigma[layer]
    J = sig.shape[1]
    co = np.zeros((J, layer.size))
    cross = np.zeros((J, layer.size)) if polarized else None
    if polarized:
        S = geo.rotate_stokes(p["S"], c2, s2)
        pa, eta = analyzer(d, pol_angle)
        a2c, a2s = geo.analyzer_rotation(l_out, d, pa)
    for j in range(J):
        sj = sig[:, j]
        if not np.any(sj > 0):
            continue
        if polarized:
            out = geo.rotate_stokes(apply_matrix(tb.lookup(tabs.P, med.types[layer, j], i, f), S), a2c, a2s)
            q = out[1] if eta is None else eta[0] * out[1] + eta[1] * out[2]
            co[j] = base * sj * 0.5 * (out[0] + q)
            cross[j] = base * sj * 0.5 * (out[0] - q)
        else:
            co[j] = base * sj * tb.lookup(tabs.P[0], med.types[layer, j], i, f)
    return d, bin_, co, cross, r_app, dist


def tally(acc, peak, bin_, order, co, cross):
    """Add contributions to acc (N, J, 3, 2) and track the largest single contribution per bin, peak (N,)."""
    N, J = acc.shape[:2]
    ok = (bin_ >= 0) & (bin_ < N)
    if not ok.any():
        return
    k = bin_[ok]
    o = np.minimum(order[ok], N_ORDERS) - 1
    base = (k * J)[None, :] + np.arange(J)[:, None]          # (J, m)
    idx = (base * N_ORDERS + o[None, :]) * 2
    flat = acc.reshape(-1)
    flat += np.bincount(idx.ravel(), weights=co[:, ok].ravel(), minlength=flat.size)
    total = co[:, ok].sum(axis=0)
    if cross is not None:
        flat += np.bincount((idx + 1).ravel(), weights=cross[:, ok].ravel(), minlength=flat.size)
        total = total + cross[:, ok].sum(axis=0)
    np.maximum.at(peak, k, np.abs(total))
