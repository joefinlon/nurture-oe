"""Photon transport for solver `mc`: launch, free path, redirection with detector directional
importance sampling (DDIS), Russian roulette (appendix C §C.2.7; Buras & Mayer 2011, Eqs. 5, 6, 10).

Packets are a structure of arrays (a dict): x, u, e (3, n); W, L (n,); order (n,) int; S (4, n) in
polarized mode, None in scalar mode. Pure functions of arrays; the random generator is passed in.
"""
import numpy as np

from . import geometry as geo
from . import tables as tb

HORIZONTAL = 1e-9            # |u_z| below this: path from tau / alpha (§C.2.7)
ROULETTE_SURVIVAL = 0.1


def launch(n, rho_tr, rng, polarized=False, pol_angle=0.0):
    """n packets at the origin. u ~ (theta_x, theta_y, 1), theta_x, theta_y ~ N(0, rho_tr^2 / 2);
    e = x-hat projected perpendicular to u; S = (1, cos 2chi, sin 2chi, 0) with chi = pol_angle."""
    th = rng.normal(0.0, rho_tr / np.sqrt(2.0), (2, n))
    u = geo.normalize(np.stack([th[0], th[1], np.ones(n)]))
    e = geo.perpendicular_axis(np.array([1.0, 0.0, 0.0]), u)
    S = None
    if polarized:
        S = np.zeros((4, n))
        S[0] = 1.0
        S[1] = np.cos(2 * pol_angle)
        S[2] = np.sin(2 * pol_angle)
    return dict(x=np.zeros((3, n)), u=u, e=e, S=S, W=np.ones(n), L=np.zeros(n), order=np.zeros(n, np.int64))


def compact(p, keep):
    """Keep the packets where `keep` is True."""
    return {k: (None if v is None else v[..., keep]) for k, v in p.items()}


def free_path(med, z, uz, rng):
    """Geometric path s to the next collision and the packets that leave [0, z_max] instead.

    tau = -ln xi; tau_z(z_new) = tau_z(z) + u_z tau; s = (z_new - z) / u_z, or tau / alpha if |u_z| < 1e-9.
    """
    tau = -np.log1p(-rng.random(z.size))
    target = tb.tau_of_z(med, z) + uz * tau
    exits = (target < 0.0) | (target > med.tau_nodes[-1])
    z_new = tb.z_of_tau(med, np.clip(target, 0.0, med.tau_nodes[-1]))
    flat = np.abs(uz) < HORIZONTAL
    s = np.where(flat, tau / med.alpha[tb.layer_of_z(med, z)], (z_new - z) / np.where(flat, 1.0, uz))
    return s, exits | (flat & ((z < 0) | (z > med.z_max)))


def move(p, med, rng):
    """Advance every packet to its next collision; returns the surviving packets.

    Dropped: packets leaving [0, z_max], and packets with (L + |x|)/2 beyond the last range bin.
    """
    s, exits = free_path(med, p["x"][2], p["u"][2], rng)
    x = p["x"] + s * p["u"]
    L = p["L"] + s
    keep = ~exits & ((L + geo.norm(x)) * 0.5 <= med.r_half[-1])
    p = dict(p, x=x, L=L, order=p["order"] + 1)
    return compact(p, keep)


def mixture(med, tabs, layer, i, f, polarized):
    """Scattering-weighted phase matrix at the interpolation point (i, f).

    Returns P_mix (6, n) or (1, n) (sum_j sigma_j P_j / sigma) and sigma = sum_j sigma_j (n,)."""
    sig = med.sigma[layer]                                  # (n, J)
    tot = sig.sum(axis=1)
    rows = tabs.P if polarized else tabs.P[:1]
    acc = np.zeros((rows.shape[0], layer.size))
    for j in range(sig.shape[1]):
        if not np.any(sig[:, j] > 0):
            continue
        acc += sig[:, j] * tb.lookup(rows, med.types[layer, j], i, f)
    return acc / np.where(tot > 0, tot, 1.0), tot


def apply_matrix(P, S):
    """Mueller matrix of randomly oriented particles (P11, P12, P22, P33, P34, P44) times S (4, n)."""
    P11, P12, P22, P33, P34, P44 = P
    return np.stack([P11 * S[0] + P12 * S[1], P12 * S[0] + P22 * S[1],
                     P33 * S[2] + P34 * S[3], -P34 * S[2] + P44 * S[3]])


def redirect(p, med, tabs, layer, d, epsilon, rng, polarized=False):
    """Implicit capture and new direction with DDIS (§C.2.7). Returns (packets, w).

    With probability epsilon the direction is drawn from p_ddis about the detector direction d;
    otherwise class j ~ alpha_j omega_j and theta ~ P11_j about u. In both branches
    w = p_mix(mu) / ((1 - eps) p_mix(mu) + eps p_ddis(mu')); in polarized mode
    S <- P_mix(mu) L(psi) S / p_mix(mu), W <- W w S_I, S <- S / S_I.
    """
    n = layer.size
    u, e = p["u"], p["e"]
    sig = med.sigma[layer]
    tot = sig.sum(axis=1)
    W = p["W"] * tot / med.alpha[layer]
    # class and angle from P11 about u
    cum = np.cumsum(sig, axis=1)
    j = np.minimum((cum < (rng.random(n) * tot)[:, None]).sum(axis=1), sig.shape[1] - 1)
    types = med.types[layer, j]
    theta = tb.sample_theta(tabs.theta_of_u, types, rng.random(n))
    phi = rng.uniform(0.0, 2 * np.pi, n)
    v = geo.new_direction(u, e, np.cos(theta), phi)
    if epsilon > 0:
        dd = rng.random(n) < epsilon
        if dd.any():
            th2 = tb.sample_theta(tabs.ddis_theta_of_u, None, rng.random(dd.sum()))
            ph2 = rng.uniform(0.0, 2 * np.pi, dd.sum())
            dsub = d[:, dd]
            a, b = geo.orthonormal_pair(dsub)
            v[:, dd] = np.cos(th2) * dsub + np.sin(th2) * (np.cos(ph2) * a + np.sin(ph2) * b)
    cos_t, sin_t, c2, s2, l_out = geo.scattering_plane(u, e, v)
    th = np.arctan2(sin_t, cos_t)
    i, f = tb.interp_index(tabs.theta, th)
    Pm, _ = mixture(med, tabs, layer, i, f, polarized)
    p_mix = Pm[0]
    if epsilon > 0:
        th_d = np.arctan2(geo.norm(geo.cross(v, d)), geo.dot(v, d))
        p_dd = np.interp(th_d, tabs.theta, tabs.ddis_p)
        w = p_mix / ((1.0 - epsilon) * p_mix + epsilon * p_dd)
    else:
        w = np.ones(n)
    W = W * w
    S = p["S"]
    if polarized:
        S = apply_matrix(Pm, geo.rotate_stokes(S, c2, s2)) / p_mix
        s_i = S[0]
        W = W * np.maximum(s_i, 0.0)                         # S_I = 0: fully blocked (e.g. Rayleigh at 90 deg)
        S = np.where(s_i > 0, S / np.where(s_i > 0, s_i, 1.0), np.array([1.0, 0.0, 0.0, 0.0])[:, None])
    return dict(p, u=v, e=l_out, W=W, S=S), w


def roulette(p, w_min, rng):
    """Russian roulette: packets with W < w_min survive with probability 0.1 and W <- 10 W."""
    low = p["W"] < w_min
    if not low.any():
        return p
    live = rng.random(low.size) < ROULETTE_SURVIVAL
    W = np.where(low & live, p["W"] / ROULETTE_SURVIVAL, p["W"])
    return compact(dict(p, W=W), ~low | live)
