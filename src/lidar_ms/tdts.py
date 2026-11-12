"""Solver `pvc+tdts`, wide-angle part: time-dependent two-stream (Hogan & Battaglia 2008),
appendix C §C.2.5 (a)-(g).

Diffuse energy travelling away from (I+) and toward (I-) the instrument is stepped in time on a
regular grid, Delta t = Delta r / c, for 2 N_out steps. The stencil is applied to the whole active
range of bins at once and vectorized over the batch axis B. Kernels accept float64 or complex128.

Internal arrays are padded with one ghost row at each end (index 0 and N+1) where I = V = 0.
"""
import numpy as np

from . import grid as _grid
from .constants import MU1

ONE_OVER_FOUR_PI = 1.0 / (4.0 * np.pi)
D_FACTOR = 4.0 / 3.0
COMPAT_D_FACTOR = 1.33333333          # MS_D_FACTOR in ms.h


def _safe(x, cond):
    return np.where(cond, x, 1.0)


# --- (a) sources ---------------------------------------------------------------------------------
def source_terms(alpha, omega, g, lobe, dr, kappa, first_bin=None):
    """Native source from the quasi-direct beam, §C.2.5 (a).

    Parameters
    ----------
    alpha, omega, g : (N, J, B) per class
    lobe : (N, J, B) bool, lobe-forming classes (L)
    dr : (B,) regular spacing
    kappa : float
    first_bin : (B,) int, optional
        Sources above this bin are zero (Dep8: start at the first cloudy bin).

    Returns
    -------
    S_plus, S_minus : (N, B)
    """
    zero = np.zeros((), dtype=alpha.dtype)
    aL = np.where(lobe, alpha, zero).sum(axis=1)
    aN = np.where(lobe, zero, alpha).sum(axis=1)
    sig = alpha * omega
    sigN = np.where(lobe, zero, sig).sum(axis=1)
    sigNg = np.where(lobe, zero, sig * g).sum(axis=1)
    d = (0.5 * aL + aN) * dr
    cum = np.concatenate([np.zeros_like(d[:1]), np.cumsum(d, axis=0)[:-1]], axis=0)
    pos = np.real(d) > 0
    Td = np.exp(-cum) * np.where(pos, -np.expm1(-_safe(d, pos)) / _safe(d, pos), 1.0)
    sw = np.where(lobe, alpha * (omega - 0.5), zero).sum(axis=1) + sigN
    swg = np.where(lobe, alpha * (omega * g - 0.5), zero).sum(axis=1) + sigNg
    ok = np.real(sw) > 0
    gw = np.where(ok, swg / _safe(sw, ok), 0.0)
    E = np.where(ok, kappa * Td * sw * dr, 0.0)
    if first_bin is not None:
        E = np.where(np.arange(E.shape[0])[:, None] >= first_bin[None, :], E, 0.0)
    fwd = np.real(gw) >= 1.0 / (3.0 * MU1)
    sp = np.where(fwd, E, 0.5 * E * (1.0 + 3.0 * gw * MU1))
    sm = np.where(fwd, 0.0, 0.5 * E * (1.0 - 3.0 * gw * MU1))
    return sp, sm


def source_width(r, dr, q, rho_tr):
    """Lateral variance of the quasi-direct beam at bin centres, §C.2.5 (a), native (exact; Dep5).

    Step 1 of §C.2.4 run one-way (q -> q/2, no equivalent medium), initialized at r_0.
    r : (N, B) ranges; dr : (B,); q : (N, B) = sum over lobe-forming classes of alpha_j Theta_j^2.
    """
    N, B = r.shape
    rt2 = rho_tr ** 2
    out = np.empty(r.shape, dtype=np.result_type(q, float))
    s2 = rt2 * r[0] ** 2 + 0.0 * q[0]
    c = rt2 * r[0] + 0.0 * q[0]
    z2 = rt2 + 0.0 * q[0]
    out[0] = s2
    for i in range(N - 1):
        qh = 0.5 * q[i]
        s2 = s2 + z2 * dr * dr + 2.0 * c * dr + qh * dr ** 3 / 3.0
        c = c + z2 * dr + 0.5 * qh * dr * dr
        z2 = z2 + qh * dr
        out[i + 1] = s2
    return out


# --- (b) scaled medium -------------------------------------------------------------------------
def scaled_medium(alpha, omega, g, lobe):
    """Delta-Eddington scaling of lobe-forming classes only, §C.2.5 (b).

    Returns alpha', omega', g' : (N, B) and the H1 class shares : (N, J, B).
    """
    zero = np.zeros((), dtype=alpha.dtype)
    sig = alpha * omega
    a_p = np.where(lobe, alpha * (1.0 - omega * g * g), alpha).sum(axis=1)
    sj = np.where(lobe, sig * (1.0 - g * g), sig)                       # scaled scattering per class
    s_p = sj.sum(axis=1)
    sg_p = np.where(lobe, sig * g * (1.0 - g), sig * g).sum(axis=1)
    ok_a = np.real(a_p) > 0
    ok_s = np.real(s_p) > 0
    om = np.where(ok_a, s_p / _safe(a_p, ok_a), 0.0)
    gp = np.where(ok_s, sg_p / _safe(s_p, ok_s), 0.0)
    share = np.where(ok_s[:, None, :], sj / _safe(s_p, ok_s)[:, None, :], zero)
    return a_p, om, gp, share


# --- (c) transfer coefficients -------------------------------------------------------------------
def transfer_coefficients(ap, om, gp, dr):
    """`calculate_deltas` in tdts.c, §C.2.5 (c). Returns D (6, N, B), l_t (N, B), T' (N, B)."""
    dt = np.result_type(ap, om, gp, float)
    pos = np.real(ap) > 0
    ap_s = _safe(ap, pos)
    lt = 1.0 / (ap_s * (1.0 - om * gp))
    ft = np.exp(-dr * ap_s * (1.0 - om * gp))
    fa = np.exp(-dr * ap_s * (1.0 - om))
    f = MU1 * (lt / dr - ft / (1.0 - ft))
    D = MU1 * om * fa * np.sqrt(lt / (3.0 * dr))
    D1 = D * np.exp(-3.7 * (lt / dr) ** 0.75)
    D2 = D * np.exp(-3.7 * lt / dr)
    d0 = ft * (1.0 - MU1) + (fa - ft) * (0.5 - f) - D1
    d1 = 0.5 * (fa - ft) * (1.0 - f) - D2
    d2 = MU1 * ft + (fa - ft) * f + 0.5 * D1
    d3 = 0.25 * (fa - ft) * f + 0.5 * D2
    d4 = 0.5 * D1
    d5 = d3
    neg0 = np.real(d0) < 0
    d2 = np.where(neg0, d2 - 0.5 * d0, d2)
    d4 = np.where(neg0, d4 - 0.5 * d0, d4)
    d0 = np.where(neg0, 0.0, d0)
    neg1 = np.real(d1) < 0
    d3 = np.where(neg1, d3 - 0.5 * d1, d3)
    d5 = np.where(neg1, d5 - 0.5 * d1, d5)
    d1 = np.where(neg1, 0.0, d1)
    clear = (1.0 - MU1, 0.0, MU1, 0.0, 0.0, 0.0)
    Dk = np.stack([np.where(pos, dk, cv) for dk, cv in zip((d0, d1, d2, d3, d4, d5), clear)]).astype(dt)
    lt = np.where(pos, lt, 1e6 * dr)
    # T'_i: one-way transmittance of the scaled medium, bin mean
    od = np.where(pos, ap, 0.0) * dr
    cum = np.concatenate([np.zeros_like(od[:1]), np.cumsum(od, axis=0)[:-1]], axis=0)
    Tp = np.exp(-cum) * np.where(pos, -np.expm1(-_safe(od, pos)) / _safe(od, pos), 1.0)
    return Dk, lt, Tp


def step_streams(Im, Ip, D, lo, hi):
    """One application of the stencil (d) on padded arrays (N+2, B) for bins lo..hi (inclusive).

    Returns the new values for padded rows lo+1..hi+1 (I-, I+).
    """
    s = slice(lo + 1, hi + 2)
    up = slice(lo, hi + 1)          # i - 1
    dn = slice(lo + 2, hi + 3)      # i + 1
    D0, D1, D2, D3, D4, D5 = D
    im = D0[s] * Im[s] + D1[s] * Ip[s] + D2[dn] * Im[dn] + D3[up] * Ip[up] + D4[up] * Im[up] + D5[dn] * Ip[dn]
    ip = D0[s] * Ip[s] + D1[s] * Im[s] + D2[up] * Ip[up] + D3[dn] * Im[dn] + D4[dn] * Ip[dn] + D5[up] * Im[up]
    return im, ip


# --- (e) lateral variance --------------------------------------------------------------------------
def expansion(I, V, s2d, lt, dr, d_factor=D_FACTOR):
    """Xi of §C.2.5 (e), with the fast inversion of nu = (V/I - s_d^2) / l_t^2 as coded."""
    posI = np.real(I) > 0
    Is = _safe(I, posI)
    nu = (V / Is - s2d) / (lt * lt)
    grow = posI & (np.real(V) > np.real(I * s2d))
    nus = _safe(nu, grow)
    sq = np.sqrt(1.5 * nus)
    n = np.where(np.real(nus) < 1.0, sq * (1.0 + 0.2 * sq), 0.75 * nus + 1.0 - 0.24 / (nus * nus))
    n = np.where(grow, n, 0.0)
    xi = d_factor * lt * (dr + lt * (np.exp(-n - dr / lt) - np.exp(-n)))
    return np.where(posI, xi, 0.0)


# --- driver ------------------------------------------------------------------------------------------
def wide_angle(r, dr, ap, om, gp, share, sp, sm, s2d, rho_tr, rho_fov, n_out, receiver="top_hat",
               d_factor=D_FACTOR, c_exp=False):
    """Time stepping (d)-(f) and the class partition H1 (g).

    Parameters
    ----------
    r : (N, B) range of bin centres from the instrument
    dr : (B,) regular spacing
    ap, om, gp : (N, B) scaled medium
    share : (N, J, B) H1 class shares of the scaled scattering
    sp, sm : (N, B) sources S+, S-
    s2d : (N, B) source width
    n_out : int >= N
    c_exp : bool
        Use 1 - exp(-x) as the C code does instead of -expm1(-x) (compat).

    Returns
    -------
    beta_w : (n_out, J, B) wide-angle apparent backscatter per class (before H2)
    """
    N, J, B = share.shape
    dt = np.result_type(ap, om, gp, sp, sm, s2d, float)
    Dk, lt, Tp = transfer_coefficients(ap, om, gp, dr)
    pad = lambda a: np.concatenate([np.zeros_like(a[:1]), a, np.zeros_like(a[:1])], axis=0)
    D = np.stack([pad(Dk[k]) for k in range(6)])
    Im, Ip, Vm, Vp = (np.zeros((N + 2, B), dtype=dt) for _ in range(4))
    Imn, Ipn, Vmn, Vpn = (np.zeros((N + 2, B), dtype=dt) for _ in range(4))
    wgt = Tp * om * ap * ONE_OVER_FOUR_PI                      # (N, B)
    fwd = np.real(gp) >= 1.0 / (3.0 * MU1)
    cm = np.where(fwd, 2.0, 1.0 + 3.0 * gp * MU1)
    cp = np.where(fwd, 0.0, 1.0 - 3.0 * gp * MU1)
    rf2 = rho_fov ** 2 * r * r
    rt2 = rho_tr ** 2 * r * r
    one_minus_exp = (lambda x: 1.0 - np.exp(-x)) if c_exp else (lambda x: -np.expm1(-x))
    out = np.zeros((n_out, J, B), dtype=dt)
    nt = 2 * n_out
    idx = np.arange(N)
    for it in range(nt):
        hi = min(it, N - 1, nt - it - 1)
        if hi >= 0:
            s = slice(1, hi + 2)
            im, ip = step_streams(Im, Ip, D, 0, hi)
            vm, vp = step_streams(Vm, Vp, D, 0, hi)
            sl = slice(0, hi + 1)
            vm = vm + Im[s] * expansion(Im[s], Vm[s], s2d[sl], lt[sl], dr, d_factor)
            vp = vp + Ip[s] * expansion(Ip[s], Vp[s], s2d[sl], lt[sl], dr, d_factor)
            Imn[s], Ipn[s], Vmn[s], Vpn[s] = im, ip, vm, vp
            # (f) tally: bin i -> output k = floor((it + i) / 2)
            ok = (np.real(vm) > 0) & (np.real(vp) > 0)
            if receiver == "top_hat":
                lim = (((it - idx[sl])[:, None] * dr[None, :]) ** 2) + rt2[sl]
                R2 = np.where(np.real(lim) < np.real(rf2[sl]), lim, rf2[sl])
                fm = one_minus_exp(R2 * im / _safe(vm, ok))
                fp = one_minus_exp(R2 * ip / _safe(vp, ok))
            else:
                fm = 1.0 / (1.0 + _safe(vm, ok) / (_safe(im, ok) * rf2[sl]))
                fp = 1.0 / (1.0 + _safe(vp, ok) / (_safe(ip, ok) * rf2[sl]))
            inc = wgt[sl] * (im * cm[sl] * fm + ip * cp[sl] * fp)
            inc = np.where(ok, inc, 0.0)
            _tally(out, inc, share[sl], it, n_out)
        if it < N:
            Imn[it + 1] = Imn[it + 1] + sm[it]
            Ipn[it + 1] = Ipn[it + 1] + sp[it]
            Vmn[it + 1] = Vmn[it + 1] + sm[it] * s2d[it]
            Vpn[it + 1] = Vpn[it + 1] + sp[it] * s2d[it]
        Im, Imn = Imn, Im
        Ip, Ipn = Ipn, Ip
        Vm, Vmn = Vmn, Vm
        Vp, Vpn = Vpn, Vp
    return out


def _tally(out, inc, share, it, n_out):
    """Add inc (n, B) x share (n, J, B) into out[k] with k = (it + i) // 2, as contiguous pairs."""
    s = it % 2
    c = inc[:, None, :] * share
    if s:
        c = np.concatenate([np.zeros_like(c[:1]), c], axis=0)
    if c.shape[0] % 2:
        c = np.concatenate([c, np.zeros_like(c[:1])], axis=0)
    pairs = c.reshape((c.shape[0] // 2, 2) + c.shape[1:]).sum(axis=1)
    k0 = it // 2
    n = min(pairs.shape[0], n_out - k0)
    if n > 0:
        out[k0:k0 + n] += pairs[:n]


# --- coarsening (§C.3.5) -------------------------------------------------------------------------------
def coarsen_classes(k, r, alpha, omega, g, q_j):
    """Conservative k-fold coarsening of class arrays for the wide-angle grid (`grid.coarsen`).

    alpha mean; scattering alpha*omega mean; alpha*omega*g mean; alpha_j Theta_j^2 mean. Pads the end
    with clear bins. Returns r_c (Nc, B), alpha_c, omega_c, g_c, q_c (Nc, J, B).
    """
    r_c = _grid.coarsen_centres(r, k)
    a_c = _grid.coarsen(alpha, k)
    s_c = _grid.coarsen(alpha * omega, k)
    sg_c = _grid.coarsen(alpha * omega * g, k)
    q_c = _grid.coarsen(q_j, k)
    ok_a = np.real(a_c) > 0
    ok_s = np.real(s_c) > 0
    om_c = np.where(ok_a, s_c / _safe(a_c, ok_a), 1.0)
    g_c = np.where(ok_s, sg_c / _safe(s_c, ok_s), 0.0)
    return r_c, a_c, om_c, g_c, q_c


def expand(beta_c, k, n_out):
    """Piecewise-constant expansion of a coarse (Nc_out, J, B) density back to n_out fine bins."""
    return np.repeat(beta_c, k, axis=0)[:n_out]
