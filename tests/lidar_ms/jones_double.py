"""Direct quadrature of second-order polarized backscatter by a slab of identical scatterers
(tests P2 and P2-Mie, appendix C §C.2.7).

Independent of the Monte Carlo's Stokes rotations: the field is propagated as a Jones vector with the
amplitude matrix of the scatterer in its scattering plane (Bohren & Huffman 1983, §3.2),

    E' = S2(theta) (E . l_in) l_out + S1(theta) (E . n) n,   n = u x u' / |u x u'|, l = n x u,

which is exact for fully polarized light scattered by one sphere size (Rayleigh: S1 = sqrt(3/2),
S2 = sqrt(3/2) cos theta). Amplitudes are normalized so that (|S1|^2 + |S2|^2)/2 = P11 with
int P11 dOmega = 4 pi. Pencil beam along z with E0 = x-hat; homogeneous slab [zt, zb] with extinction
alpha (omega = 1); top-hat receiver (accepts sqrt(x^2 + y^2) <= rho_fov z); no time gate (the caller's
grid must extend beyond every r_app). Returns the expected order-2 co- and cross-polarized local-estimate
tallies per launched packet, i.e. sum_k beta_k dr_k / kappa.
"""
import numpy as np


def rayleigh_amplitude(mu):
    c = np.sqrt(1.5)
    return np.full(np.shape(mu), c, dtype=complex), c * np.asarray(mu, dtype=complex)


def tabulated_amplitude(theta, s1, s2):
    """Amplitude function from tables on a fine theta grid (complex linear interpolation)."""
    def amp(mu):
        th = np.arccos(np.clip(mu, -1, 1))
        f = lambda s: np.interp(th, theta, s.real) + 1j * np.interp(th, theta, s.imag)
        return f(s1), f(s2)
    return amp


def _gl(a, b, n):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (b - a) * x + 0.5 * (b + a), 0.5 * (b - a) * w


def _dot(a, b):
    return (a * b).sum(0)


def _cross(a, b):
    return np.stack([a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]])


def _scatter(E, u, v, amplitude):
    """Jones field after scattering from direction u into v."""
    nv = _cross(u, v)
    nn = np.sqrt(_dot(nv, nv))
    n = nv / np.where(nn > 0, nn, 1.0)
    s1, s2 = amplitude(_dot(u, v))
    l_in, l_out = _cross(n, u), _cross(n, v)
    return s2 * _dot(E, l_in) * l_out + s1 * _dot(E, n) * n


def double_scatter(zt, zb, alpha, rho_fov, amplitude=rayleigh_amplitude, n_z=24, n_theta=1200, n_phi=24,
                   n_s=24, theta_nodes=None):
    co = cross = 0.0
    zs, wz = _gl(zt, zb, n_z)
    if theta_nodes is None:
        th, wth = _gl(0.0, np.pi, n_theta)
    else:
        th, wth = theta_nodes
    phi = 2 * np.pi * np.arange(n_phi) / n_phi
    wphi = 2 * np.pi / n_phi
    xs, ws = np.polynomial.legendre.leggauss(n_s)
    ct, st = np.cos(th), np.sin(th)
    shape = (th.size, phi.size, n_s)
    u0 = np.array([0.0, 0.0, 1.0])[:, None, None, None] * np.ones(shape)
    E0 = np.array([1.0, 0.0, 0.0])[:, None, None, None] * np.ones(shape)
    u1 = np.stack([st[:, None] * np.cos(phi)[None, :], st[:, None] * np.sin(phi)[None, :],
                   np.broadcast_to(ct[:, None], (th.size, phi.size))])[..., None] * np.ones(shape)
    E1 = _scatter(E0, u0, u1, amplitude)
    for z1, w1 in zip(zs, wz):
        smax = np.where(ct > 0, (zb - z1) / np.where(ct > 0, ct, 1),
                        np.where(ct < 0, (z1 - zt) / np.where(ct < 0, -ct, 1), np.inf))
        den = st - rho_fov * ct
        smax = np.where(den > 0, np.minimum(smax, rho_fov * z1 / np.where(den > 0, den, 1)), smax)
        s = (0.5 * smax[:, None] * (xs[None, :] + 1))[:, None, :]                     # (T, 1, S)
        wsx = (0.5 * smax[:, None] * ws[None, :])[:, None, :]
        x2 = u1 * s
        x2[2] += z1
        d2 = np.sqrt(_dot(x2, x2))
        d = -x2 / d2
        E2 = _scatter(E1, u1, d, amplitude)
        p = np.array([1.0, 0, 0])[:, None, None, None] - d[0] * d
        p = p / np.sqrt(_dot(p, p))
        q = _cross(d, p)
        ic = np.abs(_dot(E2, p)) ** 2
        ix = np.abs(_dot(E2, q)) ** 2
        z2 = x2[2]
        r_app = 0.5 * (z1 + s + d2)
        att = np.exp(-alpha * (z1 - zt)) * np.exp(-alpha * s) * np.exp(-alpha * (z2 - zt) * d2 / z2)
        k = alpha * alpha / (4 * np.pi) ** 2 * att * (r_app / d2) ** 2
        wt = w1 * (wth * st)[:, None, None] * wphi * wsx
        co += np.sum(wt * k * ic)
        cross += np.sum(wt * k * ix)
    return co, cross
