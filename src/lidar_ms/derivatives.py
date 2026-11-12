"""Directional derivatives by complex step (appendix C §C.2.8) [DECISION D27].

d beta_hat / dp . t = Im beta_hat(p + i h t) / h, h = 1e-30. No subtraction occurs, so the result is
exact to rounding. Tangents are a batch axis, like profiles: the kernels run once on P x K complex
columns.
"""
from __future__ import annotations

import numpy as np

from .constants import COMPLEX_STEP
from .forward import _FIELDS, assemble, check_inputs, solve, to_bin_major
from .schema import OpticalProfile, OpticalTangent, ReturnsTangent, SolverConfig

_OUT = ("beta_par", "beta_perp")


def _check(optics: OpticalProfile, tangents: OpticalTangent, config: SolverConfig):
    if not isinstance(tangents, OpticalTangent):
        raise TypeError("tangents must be an OpticalTangent")
    if config.solver == "mc":
        raise ValueError("jvp is not available for solver 'mc'")
    if optics.is_complex:
        raise ValueError("jvp needs real optics")
    P, K, N, J = tangents.shape
    if (P, N, J) != optics.shape:
        raise ValueError(f"tangents {tangents.shape} do not match optics {optics.shape} as (P, K, N, J)")
    return P, K


def _run(grid_, instrument, optics, tangents, config, near, step):
    """Run the kernels on x + step * t for all P x K directions; step may be complex."""
    r = check_inputs(grid_, instrument, optics, config)
    P, K, N, J = tangents.shape
    opt = {}
    for nm in _FIELDS:
        x = np.asarray(getattr(optics, nm))[:, None]
        opt[nm] = to_bin_major(x + step * getattr(tangents, nm))
    rB = np.repeat(r, K, axis=0).T.copy()
    ni = to_bin_major(np.repeat(optics.near_index[:, None], K, axis=1))
    out = solve(rB, opt, optics.kinds, optics.descriptors, instrument, config, near=near, near_index=ni,
                names=optics.names)
    return out, r


def _to_pk(a, P, K):
    """(N_out, J, P*K) -> (P, K, N_out, J)."""
    a = np.moveaxis(a, -1, 0)
    return a.reshape((P, K) + a.shape[1:])


def jvp(grid_, instrument, optics, tangents, config=SolverConfig(), *, near=None):
    """Returns and their directional derivatives along K tangents per profile (not for solver 'mc').

    Returns
    -------
    (Returns, ReturnsTangent)
        Real parts as `Returns` (identical to `simulate`), imaginary parts / h as `ReturnsTangent`.
    """
    P, K = _check(optics, tangents, config)
    h = COMPLEX_STEP
    out, r = _run(grid_, instrument, optics, tangents, config, near, 1j * h)
    first = {}
    for k, v in out.items():
        if isinstance(v, np.ndarray) and v.ndim >= 1 and v.shape[-1] == P * K:
            first[k] = np.moveaxis(v, -1, 0).reshape((P, K) + v.shape[:-1])[:, 0]
            first[k] = np.moveaxis(first[k], 0, -1)
        else:
            first[k] = v
    first = {k: (np.real(v) if np.iscomplexobj(v) else v) for k, v in first.items()}
    returns = assemble(first, P, r, optics, instrument, config, real=True)
    tan = ReturnsTangent(d_beta_par=np.imag(_to_pk(out["beta_par"], P, K)) / h,
                         d_beta_perp=np.imag(_to_pk(out["beta_perp"], P, K)) / h)
    return returns, tan


def central_difference(grid_, instrument, optics, tangents, config=SolverConfig(), *, near=None, rel_step=1e-6):
    """Central differences along each tangent; step eps_pk = rel_step * max|x| / max|t| over the
    entries the tangent touches. Returns d_beta_par, d_beta_perp of shape (P, K, N_out, J)."""
    P, K = _check(optics, tangents, config)
    eps = np.empty((P, K))
    for p in range(P):
        for k in range(K):
            tmax, xmax = 0.0, 0.0
            for nm in _FIELDS:
                t = np.abs(getattr(tangents, nm)[p, k])
                if t.max() > 0:
                    tmax = max(tmax, t.max())
                    xmax = max(xmax, np.abs(np.real(getattr(optics, nm)[p]))[t > 0].max())
            eps[p, k] = rel_step * (xmax if xmax > 0 else 1.0) / (tmax if tmax > 0 else 1.0)
    step = eps[:, :, None, None]
    scaled = OpticalTangent(**{nm: getattr(tangents, nm) * step for nm in _FIELDS})
    plus, _ = _run(grid_, instrument, optics, scaled, config, near, 1.0)
    minus, _ = _run(grid_, instrument, optics, scaled, config, near, -1.0)
    res = {}
    for nm in _OUT:
        res[nm] = (_to_pk(plus[nm], P, K) - _to_pk(minus[nm], P, K)) / (2.0 * step)
    return res["beta_par"], res["beta_perp"]


def check_jvp(grid_, instrument, optics, tangents, config=SolverConfig(), *, near=None, rel_step=1e-6):
    """Largest relative mismatch between `jvp` and central differences.

    Normwise per block: for each profile, tangent and output (par, perp), max|jvp - fd| / max|fd|.
    (Element-wise ratios are dominated by the central difference's roundoff, ~1e-16 |F| / step, in
    elements whose derivative is many orders below the block maximum.)
    """
    _, tan = jvp(grid_, instrument, optics, tangents, config, near=near)
    fd_par, fd_perp = central_difference(grid_, instrument, optics, tangents, config, near=near, rel_step=rel_step)
    worst = 0.0
    for a, b in ((tan.d_beta_par, fd_par), (tan.d_beta_perp, fd_perp)):
        scale = np.abs(b).max(axis=(2, 3))
        err = np.abs(a - b).max(axis=(2, 3))
        ok = scale > 0
        if ok.any():
            worst = max(worst, float(np.max(err[ok] / scale[ok])))
        if (~ok).any() and np.abs(a).max(axis=(2, 3))[~ok].max() > 0:
            worst = np.inf                      # fd is zero but jvp is not
    return worst
