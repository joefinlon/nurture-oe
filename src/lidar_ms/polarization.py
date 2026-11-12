"""Depolarization under multiple scattering: closures C1 and C2 and the selection rule
(appendix C §C.2.6, §C.3.4) [DECISION D26].

The fast solvers are scalar; these closures split each class's apparent backscatter into co- and
cross-polarized parts. Arrays are bin-major: (N, J, B) per class, (N, B) per bin, (N+1, B) at
half-levels. Complex-safe.
"""
import numpy as np

from .anisotropy import half_level_factors
from .constants import C1_DELTA_ACC_MAX, C1_MIN_FOOTPRINT, flag


def single_split(beta_hat, delta):
    """Single-scatter split: beta_perp = beta_hat delta / (1 + delta)."""
    return beta_hat * delta / (1.0 + delta)


def closure_c1(beta_liq, beta_ss_liq, dr):
    """Closure C1 (Hu-type), Eqs. (C.5)-(C.6), for the summed liquid classes.

    Parameters
    ----------
    beta_liq, beta_ss_liq : (N, B) apparent and single-scatter backscatter of the liquid classes
    dr : (N, B)

    Returns
    -------
    perp : (N, B)  cross-polarized liquid backscatter (may be negative: never clipped)
    a_s : (N, B)   accumulated single-scatter fraction I_S / I_T at r_{i+1/2} (real; NaN where I_T = 0)
    """
    it = np.cumsum(beta_liq * dr, axis=0)
    is_ = np.cumsum(beta_ss_liq * dr, axis=0)
    prod = it * is_
    pos = np.real(prod) > 0
    iperp = 0.5 * (it - np.where(pos, np.sqrt(np.where(pos, prod, 1.0)), 0.0))
    iperp_prev = np.concatenate([np.zeros_like(iperp[:1]), iperp[:-1]], axis=0)
    perp = (iperp - iperp_prev) / dr
    itr = np.real(it)
    with np.errstate(divide="ignore", invalid="ignore"):
        a_s = np.where(itr > 0, np.real(is_) / np.where(itr > 0, itr, 1.0), np.nan)
    return perp, a_s


def delta_acc_from_as(a_s):
    """Hu (2007) closed form inverted: delta_acc = (1 - sqrt(A_s)) / (1 + sqrt(A_s))."""
    s = np.sqrt(a_s)
    return (1.0 - s) / (1.0 + s)


def closure_c2_ratio(tu, Mt, F_minus, F_plus, Fperp_minus, Fperp_plus, c0):
    """Closure C2, Eq. (C.8): beta_perp_j / beta_hat_j = (X_perp-/X- + X_perp+/X+) / 2.

    With X = 1 + F M and X_perp = c_perp(0) + F_perp M, written with M~ = T_u M so that T_u
    underflow is harmless: X_perp / X = (T_u c0 + F_perp M~) / (T_u + F M~).
    """
    def ratio(t, m, f, fp):
        num = t[:, None, :] * c0 + fp * m[:, None, :]
        den = t[:, None, :] + f * m[:, None, :]
        ok = np.real(den) > 0
        return np.where(ok, num / np.where(ok, den, 1.0), c0)
    return 0.5 * (ratio(tu[:-1], Mt[:-1], F_minus, Fperp_minus) + ratio(tu[1:], Mt[1:], F_plus, Fperp_plus))


def lobe_above_liquid(alpha, lobe, kinds):
    """Per profile: does a non-liquid lobe-forming class lie above or inside the liquid? (N, J, B) -> (B,)."""
    liq = np.array([k == "liquid" for k in kinds])
    present = np.real(alpha) > 0
    liq_any = (present[:, liq]).any(axis=1)                          # (N, B)
    N = alpha.shape[0]
    idx = np.arange(N)[:, None]
    deepest = np.where(liq_any, idx, -1).max(axis=0)                 # (B,)
    other = (present & lobe & ~liq[None, :, None]).any(axis=1)        # (N, B)
    return (other & (idx <= deepest[None, :])).any(axis=0), liq_any.any(axis=0), liq_any


def tables_complete(alpha, lobe, kinds, tables):
    """Per profile: every present lobe-forming particulate class has a near-backscatter table."""
    missing = np.array([t is None for t in tables]) & np.array([k != "molecular" for k in kinds])
    need = (np.real(alpha) > 0) & lobe & missing[None, :, None]
    return ~need.any(axis=(0, 1))


def polarize(beta_ms, beta_ss, delta, kinds, descriptors, config, *, dr, r, alpha, lobe, tu, Mt, gamma2,
             theta2, F_minus, F_plus, tables, near_index, rho_fov):
    """Apply the configured closure per profile.

    Returns
    -------
    par, perp : (N, J, B)
    closure_used : (B,) object array of "hu", "near_backscatter" or "none"
    flags : (N, B) uint16  bits 0-3 (and 5 for tables used by C2)
    diag : dict or None  {"closure": (B,), "beta_perp": (N, J, B)} from the closure not selected
    """
    N, J, B = beta_ms.shape
    liq = np.array([k == "liquid" for k in kinds])
    above, has_liq, liq_bins = lobe_above_liquid(alpha, lobe, kinds)
    complete = tables_complete(alpha, lobe, kinds, tables)
    have_any_table = any(t is not None for t in tables)

    # selection (§C.3.4)
    mode = config.closure
    if mode == "auto":
        choice = np.where(has_liq & ~above, "hu",
                          np.where(has_liq & above, np.where(complete & have_any_table, "near_backscatter", "hu"),
                                   np.where(complete & have_any_table, "near_backscatter", "none")))
    else:
        choice = np.full(B, mode)
    choice = choice.astype(object)

    need_c1 = bool(np.any(choice == "hu")) or config.diagnose_closures
    need_c2 = bool(np.any(choice == "near_backscatter")) or config.diagnose_closures

    # candidate splits
    perp_none = single_split(beta_ms, delta)

    c1_perp = perp_none.copy()
    c1_flags = np.zeros((N, B), np.uint16)
    if liq.any() and need_c1:
        bl = beta_ms[:, liq].sum(axis=1)
        bs = beta_ss[:, liq].sum(axis=1)
        pl, a_s = closure_c1(bl, bs, dr)
        share = beta_ms[:, liq] / np.where(np.real(bl) != 0, bl, 1.0)[:, None, :]
        c1_perp[:, liq] = share * pl[:, None, :]
        dacc = delta_acc_from_as(a_s)
        with np.errstate(invalid="ignore"):
            c1_flags |= np.where(liq_bins & (dacc > C1_DELTA_ACC_MAX), flag("C1_DELTA_RANGE"), 0).astype(np.uint16)
        c1_flags |= np.where(liq_bins & (2.0 * rho_fov * r < C1_MIN_FOOTPRINT), flag("C1_SMALL_FOOTPRINT"),
                             0).astype(np.uint16)
        c1_flags |= np.where(liq_bins & (np.real(pl) < 0), flag("C1_NEGATIVE"), 0).astype(np.uint16)
        c1_flags |= np.where(liq_bins & above[None, :], flag("C1_UNDER_OTHER_LOBE"), 0).astype(np.uint16)

    c2_perp = perp_none.copy()
    if have_any_table and need_c2:
        c0 = delta / (1.0 + delta)
        Fpm = F_minus * c0
        Fpp = F_plus * c0
        for j, tab in enumerate(tables):
            if tab is None or kinds[j] == "molecular":
                continue
            tab.check_delta(delta[:, j], near_index[:, j] if near_index is not None else 0, name=str(j))
            fm, fp, _ = half_level_factors(gamma2, theta2[:, j:j + 1], descriptors[j:j + 1], kinds[j:j + 1],
                                           [tab], None if near_index is None else near_index[:, j:j + 1],
                                           component="perp")
            Fpm[:, j], Fpp[:, j] = fm[:, 0], fp[:, 0]
        c2_perp = closure_c2_ratio(tu, Mt, F_minus, F_plus, Fpm, Fpp, c0) * beta_ms

    pick = {"hu": c1_perp, "near_backscatter": c2_perp, "none": perp_none}
    perp = np.empty_like(beta_ms)
    flags = np.zeros((N, B), np.uint16)
    for name, arr in pick.items():
        m = choice == name
        if m.any():
            perp[..., m] = arr[..., m]
            if name == "hu":
                flags[:, m] |= c1_flags[:, m]
    diag = None
    if config.diagnose_closures:
        alt = np.where(choice == "hu", "near_backscatter", "hu").astype(object)
        alt_perp = np.where((alt == "hu")[None, None, :], c1_perp, c2_perp)
        diag = {"closure": alt, "beta_perp": alt_perp}
    return beta_ms - perp, perp, choice, flags, diag
