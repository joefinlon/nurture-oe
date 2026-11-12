"""Solver `pvc`: small-angle multiple scattering (Hogan 2008), appendix C §C.2.4 steps 1-6.

The recursion over range (steps 1-4 and the coangle variance of step 5) is one Python loop over
bins, vectorized over the batch axis B. It does not depend on the backscatter classes, so the class
factors F_j (step 5) and the bin means (step 6) are evaluated afterwards on whole arrays.

Branches use real parts only and both sides are computed with safe denominators, so every
function accepts float64 or complex128 (§C.2.8).
"""
import numpy as np

from .constants import COMPAT_ONE_THIRD, FORWARD_FRACTION_THRESHOLD, LOGMEAN_TOL


def _safe(x, cond):
    return np.where(cond, x, 1.0)


def propagate(dr, r_half, alpha_c, alpha_n, q, tu, rho_tr, rho_fov, kappa, receiver="top_hat",
              native=True, return_state=False):
    """Steps 1-5 of §C.2.4 at every half-level.

    Parameters
    ----------
    dr : (N, B) bin widths
    r_half : (N+1, B) half-level ranges, r_half[0] = r_{-1/2}
    alpha_c, alpha_n, q : (N, B) from Eq. (C.1)
    tu : (N+1, B) two-way unscattered transmittance at half-levels
    rho_tr, rho_fov, kappa : float
    receiver : {"top_hat", "gaussian"}
    native : bool
        False reproduces `fast_small_angle.c` (Dep4: correlation uses zeta^2; 1/3 written 0.333333).
    return_state : bool
        Also return the per-half-level moments (test A6).

    Returns
    -------
    Mt : (N+1, B)  M~ = T_u M, Eq. (C.3) without the division by T_u; Mt[0] = 0
    gamma2 : (N+1, B)  coangle variance of all forward-scattered photons, after the clip; 0 where invalid
    valid : (N+1, B) bool  T_a > 1e-12 T_r (some forward scattering in view)
    state : dict, only if return_state
    """
    N, B = dr.shape
    dt = np.result_type(dr, alpha_c, alpha_n, q, float)
    third = 1.0 / 3.0 if native else COMPAT_ONE_THIRD
    rt2, rf2 = rho_tr ** 2, rho_fov ** 2

    R0 = r_half[0]
    s2 = (R0 * R0 * rt2).astype(dt)
    z2 = np.full(B, rt2, dtype=dt)
    c = (R0 * rt2).astype(dt)
    T = np.ones(B, dtype=dt)
    s2r, z2r, cr, Tr = s2.copy(), z2.copy(), c.copy(), T.copy()

    Mt = np.zeros((N + 1, B), dtype=dt)
    gamma2 = np.zeros((N + 1, B), dtype=dt)
    valid = np.zeros((N + 1, B), dtype=bool)
    if return_state:
        keys = ("s2", "T", "Ta", "Tb", "s2a", "s2b", "s2u", "clipped", "Tr_pre", "s2r_pre")
        state = {k: np.zeros((N + 1, B), dtype=bool if k == "clipped" else dt) for k in keys}

    for i in range(N):
        d = dr[i]
        d2 = d * d
        R = r_half[i + 1]
        s2max = rf2 * R * R
        s2u = rt2 * R * R
        cu = rt2 * R
        tui = tu[i + 1]
        qi = q[i]
        att = np.exp(-(alpha_c[i] + 2.0 * alpha_n[i]) * d)

        # step 1: propagate the combined and the reduced distributions (H08 Eqs. 16-19)
        s2 = s2 + z2 * d2 + 2.0 * c * d + qi * d2 * d * third
        c = c + z2 * d + 0.5 * qi * d2
        z2 = z2 + qi * d
        T = T * att
        s2r = s2r + z2r * d2 + 2.0 * cr * d + qi * d2 * d * third
        cr = cr + z2r * d + 0.5 * qi * d2
        z2r = z2r + qi * d
        Tr = Tr * att

        # step 2: first Gaussian (forward-scattered photons kept in view)
        Ta = Tr - tui
        Ta = np.where(np.real(Ta) < 0, 0.0 * Ta, Ta)
        ok = np.real(Ta) > FORWARD_FRACTION_THRESHOLD * np.real(Tr)
        Tas = _safe(Ta, ok)
        s2a = (Tr * s2r - tui * s2u) / Tas
        z2a = (Tr * z2r - tui * rt2) / Tas
        ca = (Tr * cr - tui * cu) / Tas
        Tb = T - Tr                                   # before the clip, as coded (Dep10)
        if return_state:
            state["Tr_pre"][i + 1], state["s2r_pre"][i + 1] = Tr, s2r

        # step 3: field-of-view clip (H08 Eqs. 20-24)
        clip = ok & (np.real(s2a) > s2max)
        s2a_s = _safe(s2a, clip)
        z2a_s = _safe(z2a, clip)
        f = s2max / s2a_s
        rho2 = ca * ca / (s2a_s * z2a_s)
        Ta_c = Ta * f
        ca_c = ca * f
        z2a_c = z2a * (rho2 * f + 1.0 - rho2)
        Tr_c = tui + Ta_c
        Trs = _safe(Tr_c, clip)
        s2r_c = (Ta_c * s2max + tui * s2u) / Trs
        z2r_c = (Ta_c * z2a_c + tui * rt2) / Trs
        cr_c = (Ta_c * ca_c + tui * cu) / Trs
        Ta = np.where(clip, Ta_c, Ta)
        ca = np.where(clip, ca_c, ca)
        z2a = np.where(clip, z2a_c, z2a)
        s2a = np.where(clip, s2max + 0.0 * s2a, s2a)
        Tr = np.where(clip, Tr_c, Tr)
        s2r = np.where(clip, s2r_c, s2r)
        z2r = np.where(clip, z2r_c, z2r)
        cr = np.where(clip, cr_c, cr)

        # step 5: coangle variance of all forward-scattered photons (H08 Eqs. 32-33)
        Tg = T - tui
        Tgs = _safe(Tg, ok)
        s2g = (T * s2 - tui * s2u) / Tgs
        z2g = (T * z2 - tui * rt2) / Tgs
        cg = (T * c - tui * cu) / Tgs
        g2 = z2g + s2g / (R * R) - 2.0 * cg / R
        gclip = ok & (np.real(s2g) > s2max)
        s2gs = _safe(s2g, gclip)
        csg = cg - s2g / R
        den = s2gs * _safe(g2 if native else z2g, gclip)           # Dep4
        rho2g = csg * csg / den
        g2 = np.where(gclip, g2 * (rho2g * s2max / s2gs + 1.0 - rho2g), g2)

        # step 4: second Gaussian and multiple-scattering factor (H08 Eqs. 27-28)
        has_b = ok & (np.real(Tb) > 0)
        Tbs = _safe(Tb, has_b)
        s2b = (T * s2 - Tr * s2r) / Tbs
        s2b = np.where(np.real(s2b) <= np.real(s2a), s2a, s2b)
        s2as = _safe(s2a, ok)
        s2bs = _safe(s2b, has_b)
        if receiver == "top_hat":
            ma = -np.expm1(-s2max / s2as) * Ta
            mb = -np.expm1(-s2max / s2bs) * Tb
        else:
            ma = Ta / (1.0 + s2as / s2max)
            mb = Tb / (1.0 + s2bs / s2max)
        m = kappa * (ma + np.where(has_b, mb, 0.0))
        Mt[i + 1] = np.where(ok, m, 0.0)
        gamma2[i + 1] = np.where(ok, g2, 0.0)
        valid[i + 1] = ok
        if return_state:
            state["s2"][i + 1], state["T"][i + 1] = s2, T
            state["Ta"][i + 1], state["Tb"][i + 1] = np.where(ok, Ta, 0.0), np.where(has_b, Tb, 0.0)
            state["s2a"][i + 1], state["s2b"][i + 1] = s2a, s2b
            state["s2u"][i + 1], state["clipped"][i + 1] = s2u, clip

    if return_state:
        return Mt, gamma2, valid, state
    return Mt, gamma2, valid


def log_mean(beta, y_minus, y_plus):
    """Bin mean, numerical form of Eq. (C.4): beta (Y- - Y+) / ln(Y- / Y+).

    beta Y- where |Y- - Y+| <= 1e-12 Y-, 0 where Y- = 0 or Y+ = 0. When Y- and Y+ are close
    the logarithm is taken as -log1p(-(Y- - Y+)/Y-) with the same difference as the numerator, so
    the ratio stays accurate; when they are far apart (Y+ < Y-/2) it is ln(Y-/Y+).
    """
    d = y_minus - y_plus
    pos = (np.real(y_minus) > 0) & (np.real(y_plus) > 0)
    ym = _safe(y_minus, pos)
    yp = _safe(y_plus, pos)
    near = np.abs(np.real(d)) <= LOGMEAN_TOL * np.real(ym)
    far = np.real(yp) < 0.5 * np.real(ym)
    close = ~far & ~near & pos
    lg = np.where(far, np.log(ym / yp), -np.log1p(np.where(close, -d / ym, 0.0)))
    val = d / _safe(lg, ~near & pos)
    out = np.where(near, y_minus, val)
    return beta * np.where(pos, out, 0.0)


def bin_means(beta, tu, Mt, F_minus, F_plus):
    """Step 6 for every class: Y-_{j,i} = T_{u,i-1/2} + F-_{j,i} M~_{i-1/2}, Y+ likewise at i+1/2.

    beta, F_minus, F_plus : (N, J, B); tu, Mt : (N+1, B). Returns (N, J, B).
    """
    ym = tu[:-1, None, :] + F_minus * Mt[:-1, None, :]
    yp = tu[1:, None, :] + F_plus * Mt[1:, None, :]
    return log_mean(beta, ym, yp)


def bin_means_direct(beta, tu, Mt, F_minus, F_plus, tau):
    """Eq. (C.4) in its direct form (H08 Eqs. A1-A4), for comparison only (test N3). May be NaN."""
    with np.errstate(all="ignore"):
        M = Mt / tu
        xm = 1.0 + F_minus * M[:-1, None, :]
        xp = 1.0 + F_plus * M[1:, None, :]
        t = tau[:, None, :]
        return beta * tu[:-1, None, :] * (xm - xp * np.exp(-t)) / (t - np.log(xp / xm))
