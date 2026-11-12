"""`simulate()`: solver dispatch and assembly of `Returns` (appendix C §C.3.2).

The public layout is (P, N, J). Kernels work bin-major: (N, J, B) and (N, B), with B the batch
axis (profiles, or profiles x tangents inside `jvp`) (§C.3.5).
"""
from __future__ import annotations

import hashlib
from dataclasses import asdict

import numpy as np

from . import __version__, compat, grid, mixing
from .constants import EXTINGUISHED_TU, flag
from .schema import Instrument, OpticalProfile, RangeGrid, Returns, SolverConfig
from .single import single

_FIELDS = ("alpha", "omega", "g", "a_eq", "beta", "delta")


# --- layout ------------------------------------------------------------------------------------
def to_bin_major(a):
    """(P, N, J) -> (N, J, P); (P, K, N, J) -> (N, J, P*K)."""
    a = np.asarray(a)
    if a.ndim == 4:
        P, K, N, J = a.shape
        a = a.reshape(P * K, N, J)
    return np.ascontiguousarray(np.moveaxis(a, 0, -1))


def from_bin_major(a, P, K=None):
    """(N, J, B) -> (P, N, J) or (P, K, N, J); (N, B) -> (P, N) or (P, K, N)."""
    a = np.moveaxis(a, -1, 0)
    if K is not None:
        a = a.reshape((P, K) + a.shape[1:])
    return a


def check_inputs(grid_: RangeGrid, instrument: Instrument, optics: OpticalProfile, config: SolverConfig):
    if not isinstance(grid_, RangeGrid) or not isinstance(optics, OpticalProfile) \
            or not isinstance(instrument, Instrument) or not isinstance(config, SolverConfig):
        raise TypeError("expected RangeGrid, Instrument, OpticalProfile, SolverConfig")
    P, N, _ = optics.shape
    r = grid_.r
    if r.shape[1] != N:
        raise ValueError(f"grid has {r.shape[1]} bins, optics {N}")
    if r.shape[0] == 1 and P > 1:
        r = np.broadcast_to(r, (P, N))
    elif r.shape[0] != P:
        raise ValueError(f"grid has {r.shape[0]} profiles, optics {P}")
    return np.asarray(r)


def checksum(r, optics: OpticalProfile) -> str:
    h = hashlib.sha256(np.ascontiguousarray(r).tobytes())
    for nm in _FIELDS:
        h.update(np.ascontiguousarray(np.real(getattr(optics, nm))).tobytes())
    h.update(repr((optics.kinds, optics.descriptors)).encode())
    return h.hexdigest()[:16]


# --- core --------------------------------------------------------------------------------------
def solve(r, opt, kinds, descriptors, instrument, config, near=None, near_index=None):
    """Run the configured solver on bin-major arrays.

    Parameters
    ----------
    r : (N, B) bin centres
    opt : dict of (N, J, B) arrays: alpha, omega, g, a_eq, beta, delta
    kinds, descriptors : (J,) tuples

    Returns
    -------
    dict with beta_par, beta_perp, beta_wide (N_out, J, B); beta_ss (N, J, B); ms_factor, gamma2
    (N+1, B); a_s, flags (N, B); closure_used (B,) object.
    """
    N, J, B = opt["alpha"].shape
    is_compat = config.compat is not None
    dr = grid.bin_widths(r.T).T
    alpha = opt["alpha"]
    alpha_tot = alpha.sum(axis=1)
    beta = opt["beta"]

    if config.solver == "single":
        if is_compat:
            beta_ss, tu, tau = compat.single_compat(alpha_tot, beta, dr)
        else:
            beta_ss, tu, tau = single(alpha_tot, beta, dr)
        beta_ms = beta_ss
        M = np.zeros((N + 1, B), dtype=beta_ss.dtype)
        gamma2 = np.zeros((N + 1, B))
    else:
        raise NotImplementedError(f"solver {config.solver!r} is not built yet")

    n_out = N if config.n_out is None else max(int(config.n_out), N)
    beta_wide = np.zeros((n_out, J, B), dtype=beta_ms.dtype)

    par, perp, closure_used, pol_flags = _polarize(beta_ms, beta_ss, opt["delta"], kinds, config)
    flags = pol_flags | np.where(np.real(tu[:-1]) < EXTINGUISHED_TU, flag("EXTINGUISHED"), 0).astype(np.uint16)

    return dict(beta_par=_pad(par, n_out), beta_perp=_pad(perp, n_out), beta_ss=beta_ss,
                beta_wide=beta_wide, ms_factor=M, gamma2=gamma2,
                a_s=accumulated_ss_fraction(beta_ms, beta_ss, dr, kinds), flags=flags.astype(np.uint16),
                closure_used=closure_used)


def _pad(a, n_out):
    if a.shape[0] == n_out:
        return a
    return np.concatenate([a, np.zeros((n_out - a.shape[0],) + a.shape[1:], dtype=a.dtype)], axis=0)


def _polarize(beta_ms, beta_ss, delta, kinds, config):
    """Single-scatter split per class: beta_perp = beta delta / (1 + delta). Closures arrive at M6."""
    perp = beta_ms * delta / (1.0 + delta)
    B = beta_ms.shape[-1]
    return beta_ms - perp, perp, np.full(B, "none", dtype=object), np.zeros(beta_ms.shape[::2], np.uint16)


def accumulated_ss_fraction(beta_ms, beta_ss, dr, kinds):
    """A_s(r) = I_S / I_T of the liquid classes from the first bin to the bottom of each bin; NaN where I_T = 0."""
    liq = np.array([k == "liquid" for k in kinds])
    if not liq.any():
        return np.full(dr.shape, np.nan)
    it = np.cumsum(np.real(beta_ms[:, liq]).sum(axis=1) * dr, axis=0)
    is_ = np.cumsum(np.real(beta_ss[:, liq]).sum(axis=1) * dr, axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(it > 0, is_ / np.where(it > 0, it, 1.0), np.nan)


# --- public ------------------------------------------------------------------------------------
def simulate(grid_, instrument, optics, config=SolverConfig(), *, near=None, phase=None) -> Returns:
    """Apparent backscatter per class and polarization for P profiles.

    Parameters
    ----------
    grid_ : RangeGrid, (P, N) or (1, N) broadcast over profiles
    instrument : Instrument
    optics : OpticalProfile, (P, N, J)
    config : SolverConfig
    near : dict {class name or index: NearBackscatter}, optional
        Needed for descriptor "table" and closure C2.
    phase : PhaseTables, optional (solver "mc" only)
    """
    r = check_inputs(grid_, instrument, optics, config)
    P = optics.shape[0]
    opt = {nm: to_bin_major(getattr(optics, nm)) for nm in _FIELDS}
    out = solve(r.T.copy(), opt, optics.kinds, optics.descriptors, instrument, config, near=near,
                near_index=to_bin_major(optics.near_index))
    return assemble(out, P, r, optics, instrument, config)


def assemble(out, P, r, optics, instrument, config, real=False) -> Returns:
    f = (lambda a: np.real(a)) if real else (lambda a: a)
    n_out = out["beta_par"].shape[0]
    r_out = r
    if n_out > r.shape[1]:
        d = r[:, -1:] - r[:, -2:-1]
        r_out = np.concatenate([r, r[:, -1:] + d * np.arange(1, n_out - r.shape[1] + 1)], axis=1)
    meta = dict(solver=config.solver, config=_config_dict(config), version=__version__,
                compat=config.compat or "native", checksum=checksum(r, optics),
                instrument=asdict(instrument))
    return Returns(beta_par=f(from_bin_major(out["beta_par"], P)),
                   beta_perp=f(from_bin_major(out["beta_perp"], P)),
                   beta_ss=f(from_bin_major(out["beta_ss"], P)),
                   beta_wide=f(from_bin_major(out["beta_wide"], P)),
                   ms_factor=np.real(from_bin_major(out["ms_factor"], P)),
                   gamma2=np.real(from_bin_major(out["gamma2"], P)),
                   a_s=from_bin_major(out["a_s"], P), closure_used=np.asarray(out["closure_used"]),
                   flags=from_bin_major(out["flags"], P), kinds=optics.kinds, names=optics.names,
                   r=r_out, meta=meta)


def _config_dict(config):
    d = asdict(config)
    d["mc"] = dict(d["mc"])
    return d
