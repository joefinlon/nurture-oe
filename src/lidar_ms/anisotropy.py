"""Near-backscatter anisotropy factors F_j (appendix C §C.2.4 step 5; closure C2, §C.2.6).

F_j(gamma2) = int_0^inf p_j(sqrt(gamma2 t)) exp(-t) dt (H08 Eq. 34), where p_j(gamma) is the phase
function at pi - gamma normalized by its value at pi. Closed forms: `droplet` (H08 Eq. 36),
`yang` (H08 Eq. 37), `isotropic` (1). Tables: exact integral of the piecewise-linear-in-gamma^2
interpolant, Eq. (C.9). All functions are complex-safe (branches on real parts only).
"""
import numpy as np
from scipy.special import erfcx

from .constants import (DROPLET_U0, DROPLET_U1, DROPLET_U2, DROPLET_V1_SQ, DROPLET_V2_SQ, YANG_GAMMA0,
                        YANG_W)


def factor_droplet(gamma2, theta2):
    """H08 Eq. 36: U0 + U1 / (1 + gamma2 V1^2 / Theta^2) + U2 / (1 + gamma2 V2^2 / Theta^2)."""
    u = gamma2 / theta2
    return DROPLET_U0 + DROPLET_U1 / (1.0 + u * DROPLET_V1_SQ) + DROPLET_U2 / (1.0 + u * DROPLET_V2_SQ)


def factor_yang(gamma2):
    """H08 Eq. 37: 1 - W + W [1 - g sqrt(pi) erfcx(g)], g = sqrt(gamma2) / (2 gamma0); 1 at gamma2 = 0."""
    pos = np.real(gamma2) > 0
    g = np.sqrt(np.where(pos, gamma2, 1.0)) / (2.0 * YANG_GAMMA0)
    f = 1.0 - YANG_W + YANG_W * (1.0 - g * np.sqrt(np.pi) * erfcx(g))
    return np.where(pos, f, 1.0)


def p_droplet(gamma_over_theta):
    """Normalized near-180 phase function behind Eq. 36 (H08 Eq. 29), argument gamma / Theta."""
    x2 = np.asarray(gamma_over_theta) ** 2
    return DROPLET_U0 + DROPLET_U1 * np.exp(-DROPLET_V1_SQ * x2) + DROPLET_U2 * np.exp(-DROPLET_V2_SQ * x2)


def p_yang(gamma):
    """Normalized near-180 phase function behind Eq. 37: 1 - W + W exp(-gamma / gamma0)."""
    return 1.0 - YANG_W + YANG_W * np.exp(-np.asarray(gamma) / YANG_GAMMA0)


def table_factor(gamma2, x, p, scale2=1.0):
    """Eq. (C.9): exact integral of a tabulated p, linear in gamma^2, constant beyond the last node.

    Parameters
    ----------
    gamma2 : array (...)
        Coangle variance.
    x : (G,)
        Node abscissae, x_0 = 0 (gamma, or gamma / Theta).
    p : (..., G)
        Normalized function at the nodes, broadcast against gamma2[..., None].
    scale2 : array (...) or float
        Theta^2 for a table in gamma / Theta, 1 for a table in gamma.

    Returns
    -------
    F : (...)
    tail : (...) real
        exp(-t_{G-1}), the weight beyond the last node (flag TABLE_EXTRAPOLATED).
    """
    gamma2 = np.asarray(gamma2)
    pos = np.real(gamma2) > 0
    s = np.where(pos, scale2 / np.where(pos, gamma2, 1.0), 0.0)[..., None]   # t = x^2 * s
    x2 = np.asarray(x, dtype=float) ** 2
    t = x2 * s
    dt = np.diff(x2) * s
    dp = np.diff(p, axis=-1)
    small = np.real(dt) < 1e-300
    w = np.exp(-t[..., :-1]) * np.where(small, 1.0, -np.expm1(-dt) / np.where(small, 1.0, dt))
    F = p[..., 0] + np.sum(dp * w, axis=-1)
    F = np.where(pos, F, p[..., 0])
    tail = np.where(pos, np.exp(-np.real(t[..., -1])), 0.0)
    return F, tail


def class_factor(descriptor, kind, gamma2, theta2, table=None, near_index=None, component="total"):
    """F for one class at given coangle variances.

    Parameters
    ----------
    descriptor : str
    kind : str
        Molecular classes return 1 whatever the descriptor.
    gamma2, theta2 : arrays broadcast together
    table : NearBackscatter, needed for "table" and for the polarized components
    near_index : int array like gamma2, row of the table
    component : {"total", "par", "perp"}
        "par"/"perp" give F^par, F^perp of closure C2 (they need a table).

    Returns
    -------
    F, tail : arrays; tail is 0 for closed forms.
    """
    shape = np.broadcast(gamma2, theta2).shape
    if kind == "molecular" or (descriptor == "isotropic" and component == "total"):
        one = np.ones(shape, dtype=np.result_type(gamma2, theta2, float))
        return one, np.zeros(shape)
    if component == "total" and descriptor == "droplet":
        return factor_droplet(gamma2, theta2), np.zeros(shape)
    if component == "total" and descriptor == "yang":
        return factor_yang(np.broadcast_to(gamma2, shape)), np.zeros(shape)
    if table is None:
        raise ValueError(f"descriptor {descriptor!r} / component {component!r} needs a NearBackscatter table")
    rows = {"total": table.c_par + table.c_perp, "par": table.c_par, "perp": table.c_perp}[component]
    idx = np.zeros(shape, dtype=int) if near_index is None else np.clip(near_index, 0, table.n_types - 1)
    p = rows[np.broadcast_to(idx, shape)]
    scale2 = theta2 if table.coord == "gamma_over_theta" else 1.0
    return table_factor(np.broadcast_to(gamma2, shape), table.x, p,
                        np.broadcast_to(scale2, shape) if table.coord == "gamma_over_theta" else 1.0)


def half_level_factors(gamma2, theta2, descriptors, kinds, tables=None, near_index=None, component="total"):
    """F^-_{j,i} = F_j(gamma2_{i-1/2}) and F^+_{j,i} = F_j(gamma2_{i+1/2}), both with bin-i properties.

    Parameters
    ----------
    gamma2 : (N+1, B) coangle variance at half-levels
    theta2 : (N, J, B)
    tables : sequence of length J of NearBackscatter or None
    near_index : (N, J, B) int

    Returns
    -------
    F_minus, F_plus : (N, J, B)
    tail : (N, J, B) real, max over both half-levels
    """
    N, J, B = theta2.shape
    tables = tables if tables is not None else [None] * J
    dtype = np.result_type(gamma2, theta2, float)
    Fm = np.empty((N, J, B), dtype=dtype)
    Fp = np.empty((N, J, B), dtype=dtype)
    tail = np.zeros((N, J, B))
    for j in range(J):
        th = np.where(np.isfinite(np.real(theta2[:, j])), theta2[:, j], 1.0)
        ni = None if near_index is None else near_index[:, j]
        fm, tm = class_factor(descriptors[j], kinds[j], gamma2[:-1], th, tables[j], ni, component)
        fp, tp = class_factor(descriptors[j], kinds[j], gamma2[1:], th, tables[j], ni, component)
        Fm[:, j], Fp[:, j] = fm, fp
        tail[:, j] = np.maximum(tm, tp)
    return Fm, Fp, tail
