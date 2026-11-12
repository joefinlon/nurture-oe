"""Frozen data schemas of `lidar_ms` (appendix C §C.3.3) [DECISION D29].

All quantities are SI. Arrays are `float64`, or `complex128` inside `jvp`; constraints are tested on
real parts. Invalid inputs raise `ValueError`; nothing is silently fixed (§C.3.0 rule 6).

Dimensions: P profiles, N range bins, N_out output bins, J classes, K tangents per profile,
T table types, G coangle nodes, A scattering angles.
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass, field, fields
from typing import Literal, Optional

import numpy as np

from . import grid as _grid
from .constants import DELTA_TABLE_MISMATCH, THETA_MAX_DEFAULT

KINDS = ("liquid", "ice", "other", "molecular")
DESCRIPTORS = ("droplet", "yang", "isotropic", "table")
SOLVERS = ("single", "pvc", "pvc+tdts", "mc")
CLOSURES = ("auto", "hu", "near_backscatter", "none")
COMPAT_MODES = (None, "multiscatter-1.2.11")


def _real(x):
    return np.real(np.asarray(x))


def _num(x):
    """float64 or complex128 copy of x."""
    x = np.asarray(x)
    return np.array(x, dtype=complex if np.iscomplexobj(x) else float)


def _set(obj, name, value):
    object.__setattr__(obj, name, value)


# --------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class Instrument:
    """Lidar geometry (§C.2.1). Angles in radians, lengths in metres.

    Parameters
    ----------
    wavelength : float
    rho_tr : float
        Transmitter 1/e half-width; energy ~ exp(-theta^2 / rho_tr^2).
    rho_fov : float
        Receiver half-angle (top hat) or 1/e half-width (Gaussian receiver).
    receiver : {"top_hat", "gaussian"}
    name : str
        Provenance only.
    """

    wavelength: float
    rho_tr: float
    rho_fov: float
    receiver: Literal["top_hat", "gaussian"] = "top_hat"
    name: str = ""

    def __post_init__(self):
        for nm in ("wavelength", "rho_tr", "rho_fov"):
            v = getattr(self, nm)
            if not (np.isscalar(v) and np.isfinite(v) and v > 0):
                raise ValueError(f"Instrument.{nm} must be a positive finite number, got {v!r}")
            _set(self, nm, float(v))
        if self.receiver not in ("top_hat", "gaussian"):
            raise ValueError(f"Instrument.receiver must be 'top_hat' or 'gaussian', got {self.receiver!r}")

    @classmethod
    def from_full_angles(cls, wavelength, divergence_full_e2, fov_full, **kw):
        """Build from the quantities instrument papers quote.

        `divergence_full_e2` is the full-angle divergence at 1/e^2 of energy, so the 1/e half-width is
        divergence_full_e2 / (2 sqrt 2); `fov_full` is the full-angle field of view.
        HALO [GAP G25]: ``from_full_angles(532e-9, 0.8e-3, 1.0e-3)``. There is no F-Lidar preset [GAP G6].
        """
        return cls(wavelength, 0.5 * divergence_full_e2 / np.sqrt(2.0), 0.5 * fov_full, **kw)

    @property
    def kappa(self) -> float:
        """Calibration factor that makes single scattering return beta exp(-2 tau) (§C.2.1)."""
        if self.receiver == "top_hat":
            return 1.0 / (1.0 - np.exp(-(self.rho_fov / self.rho_tr) ** 2))
        return 1.0 + (self.rho_tr / self.rho_fov) ** 2


# --------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class RangeGrid:
    """Bin centres along the beam.

    Parameters
    ----------
    r : array_like, shape (P, N) or (N,)
        Positive, strictly increasing along N, N >= 3. A 1-D input is one profile.
    z : array_like, shape (P, N) or (N,), optional
        Altitude of bin centres, provenance only.
    """

    r: np.ndarray
    z: Optional[np.ndarray] = None

    def __post_init__(self):
        r = np.array(self.r, dtype=float)
        if r.ndim == 1:
            r = r[None, :]
        if r.ndim != 2 or r.shape[1] < 3:
            raise ValueError(f"RangeGrid.r must have shape (P, N) with N >= 3, got {r.shape}")
        if not np.all(np.isfinite(r)) or np.any(r <= 0):
            raise ValueError("RangeGrid.r must be positive and finite")
        if np.any(np.diff(r, axis=1) <= 0):
            raise ValueError("RangeGrid.r must be strictly increasing along range")
        _set(self, "r", r)
        if self.z is not None:
            z = np.array(self.z, dtype=float)
            if z.shape not in (r.shape, r.shape[1:]):
                raise ValueError(f"RangeGrid.z must have shape {r.shape} or {r.shape[1:]}, got {z.shape}")
            _set(self, "z", z)

    @property
    def n_profiles(self) -> int:
        return self.r.shape[0]

    @property
    def n_bins(self) -> int:
        return self.r.shape[1]

    @property
    def dr(self) -> np.ndarray:
        """Bin widths, (P, N) (§C.2.1)."""
        return _grid.bin_widths(self.r)

    @property
    def r_half(self) -> np.ndarray:
        """Half-level ranges, (P, N+1) (§C.2.1)."""
        return _grid.half_levels(self.r)

    def is_regular(self, tol=0.05) -> np.ndarray:
        """Per profile: spacing constant to within `tol`."""
        return _grid.is_regular(self.r, tol)


# --------------------------------------------------------------------------------------------
_OPT_FIELDS = ("alpha", "omega", "g", "a_eq", "beta", "delta")


@dataclass(frozen=True)
class OpticalProfile:
    """Per-bin, per-class optical properties, shape (P, N, J) (§C.2.2).

    Parameters
    ----------
    kinds : sequence of {"liquid", "ice", "other", "molecular"}
        Exactly one "molecular".
    descriptors : sequence of {"droplet", "yang", "isotropic", "table"}
        Near-180 deg form per class.
    alpha : extinction coefficient [m-1] (molecular: scattering + gas absorption), >= 0
    omega : single-scattering albedo, in [0, 1]
    g : asymmetry factor, in [-1, 1]
    a_eq : equivalent-area radius [m], pi a^2 = <G>; > 0 where alpha > 0; ignored for molecular
    beta : backscatter coefficient [m-1 sr-1], >= 0
    delta : single-scatter linear depolarization ratio, in [0, 1]
    near_index : int array, optional
        Row of the class's near-180 table; default 0.
    names : sequence of str, optional
        Free labels; default: the kinds.
    """

    kinds: tuple
    descriptors: tuple
    alpha: np.ndarray
    omega: np.ndarray
    g: np.ndarray
    a_eq: np.ndarray
    beta: np.ndarray
    delta: np.ndarray
    near_index: Optional[np.ndarray] = None
    names: tuple = field(default=())

    def __post_init__(self):
        kinds = tuple(self.kinds)
        desc = tuple(self.descriptors)
        J = len(kinds)
        bad = [k for k in kinds if k not in KINDS]
        if bad:
            raise ValueError(f"unknown kinds {bad}; allowed {KINDS}")
        if kinds.count("molecular") != 1:
            raise ValueError("OpticalProfile needs exactly one 'molecular' class")
        if len(desc) != J or any(d not in DESCRIPTORS for d in desc):
            raise ValueError(f"need one descriptor per class from {DESCRIPTORS}, got {desc}")
        names = tuple(self.names) if self.names else kinds
        if len(names) != J or len(set(names)) != J:
            raise ValueError("names must be unique, one per class")
        _set(self, "kinds", kinds)
        _set(self, "descriptors", desc)
        _set(self, "names", names)

        arrs = {nm: _num(getattr(self, nm)) for nm in _OPT_FIELDS}
        shape = arrs["alpha"].shape
        if len(shape) != 3 or shape[2] != J or shape[1] < 1:
            raise ValueError(f"arrays must have shape (P, N, J={J}), alpha has {shape}")
        for nm, a in arrs.items():
            if a.shape != shape:
                raise ValueError(f"{nm} has shape {a.shape}, expected {shape}")
            if not np.all(np.isfinite(a)):
                raise ValueError(f"{nm} must be finite")
        al, om, gg, ae, be, de = (_real(arrs[nm]) for nm in _OPT_FIELDS)
        if np.any(al < 0):
            raise ValueError("alpha must be >= 0")
        if np.any((om < 0) | (om > 1)):
            raise ValueError("omega must lie in [0, 1]")
        if np.any((gg < -1) | (gg > 1)):
            raise ValueError("g must lie in [-1, 1]")
        if np.any(be < 0):
            raise ValueError("beta must be >= 0")
        if np.any((de < 0) | (de > 1)):
            raise ValueError("delta must lie in [0, 1]")
        part = np.array([k != "molecular" for k in kinds])
        if np.any((al > 0) & (ae <= 0) & part):
            raise ValueError("a_eq must be > 0 wherever a particulate class has alpha > 0")
        for nm, a in arrs.items():
            _set(self, nm, a)
        if self.near_index is None:
            ni = np.zeros(shape, dtype=np.int64)
        else:
            ni = np.asarray(self.near_index)
            if ni.shape != shape or not np.issubdtype(ni.dtype, np.integer) or np.any(ni < 0):
                raise ValueError(f"near_index must be a non-negative integer array of shape {shape}")
            ni = ni.astype(np.int64)
        _set(self, "near_index", ni)

    @property
    def shape(self):
        return self.alpha.shape

    @property
    def is_complex(self) -> bool:
        return any(np.iscomplexobj(getattr(self, nm)) for nm in _OPT_FIELDS)

    @property
    def mol_index(self) -> int:
        return self.kinds.index("molecular")

    def replace(self, **kw) -> "OpticalProfile":
        d = {f.name: getattr(self, f.name) for f in fields(self)}
        d.update(kw)
        return OpticalProfile(**d)


# --------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class NearBackscatter:
    """Near-180 deg co-/cross-polarized functions of one class (C.7) [DECISION D30].

    Parameters
    ----------
    coord : {"gamma", "gamma_over_theta"}
        Abscissa: coangle gamma [rad] (ice) or gamma / Theta_j (droplets).
    x : array (G,)
        Starts at 0, strictly increasing.
    c_par, c_perp : array (T, G) or (G,)
        c_par[:, 0] + c_perp[:, 0] == 1.
    """

    coord: Literal["gamma", "gamma_over_theta"]
    x: np.ndarray
    c_par: np.ndarray
    c_perp: np.ndarray

    def __post_init__(self):
        if self.coord not in ("gamma", "gamma_over_theta"):
            raise ValueError("coord must be 'gamma' or 'gamma_over_theta'")
        x = np.array(self.x, dtype=float)
        if x.ndim != 1 or x.size < 2 or x[0] != 0 or np.any(np.diff(x) <= 0):
            raise ValueError("x must be 1-D, start at 0 and be strictly increasing")
        cp, cx = (np.array(a, dtype=float) for a in (self.c_par, self.c_perp))
        cp, cx = (a[None, :] if a.ndim == 1 else a for a in (cp, cx))
        if cp.ndim != 2 or cp.shape != cx.shape or cp.shape[1] != x.size:
            raise ValueError(f"c_par and c_perp must have shape (T, {x.size})")
        if not (np.all(np.isfinite(cp)) and np.all(np.isfinite(cx))):
            raise ValueError("c_par and c_perp must be finite")
        if np.any(np.abs(cp[:, 0] + cx[:, 0] - 1.0) > 1e-6):
            raise ValueError("c_par[:, 0] + c_perp[:, 0] must equal 1")
        _set(self, "x", x)
        _set(self, "c_par", cp)
        _set(self, "c_perp", cx)

    @property
    def n_types(self) -> int:
        return self.c_par.shape[0]

    def check_delta(self, delta, near_index, name=""):
        """Warn if the single-scatter split delta/(1+delta) differs from c_perp(0) by > 0.02 (§C.3.3)."""
        d = _real(delta)
        c0 = self.c_perp[np.clip(near_index, 0, self.n_types - 1), 0]
        if np.any(np.abs(d / (1 + d) - c0) > DELTA_TABLE_MISMATCH):
            warnings.warn(f"class {name!r}: delta/(1+delta) differs from the near-backscatter table's "
                          f"c_perp(0) by more than {DELTA_TABLE_MISMATCH}; delta from OpticalProfile is used",
                          stacklevel=3)


# --------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class PhaseTables:
    """Full phase matrices for solver `mc` (one profile). theta (A,) from 0 to pi; P_ij (T, A).

    The lobe-resolution constraint (first spacing <= Theta_min / 10) needs the optics and is checked
    by the Monte Carlo set-up.
    """

    theta: np.ndarray
    P11: np.ndarray
    P12: np.ndarray
    P22: np.ndarray
    P33: np.ndarray
    P34: np.ndarray
    P44: np.ndarray
    type_index: np.ndarray

    def __post_init__(self):
        th = np.array(self.theta, dtype=float)
        if th.ndim != 1 or th.size < 3 or th[0] != 0 or not np.isclose(th[-1], np.pi) \
                or np.any(np.diff(th) <= 0):
            raise ValueError("theta must be 1-D, strictly increasing, from 0 to pi")
        _set(self, "theta", th)
        shape = None
        for nm in ("P11", "P12", "P22", "P33", "P34", "P44"):
            a = np.array(getattr(self, nm), dtype=float)
            a = a[None, :] if a.ndim == 1 else a
            if a.ndim != 2 or a.shape[1] != th.size or (shape is not None and a.shape != shape):
                raise ValueError(f"{nm} must have shape (T, {th.size}) like P11")
            shape = a.shape
            _set(self, nm, a)
        norm = np.trapezoid(self.P11 * np.sin(th), th, axis=1)
        if np.any(np.abs(norm - 2.0) > 2e-3):
            raise ValueError(f"int P11 sin(theta) dtheta must equal 2 to 1e-3, got {norm}")
        ti = np.asarray(self.type_index)
        if ti.ndim != 2 or not np.issubdtype(ti.dtype, np.integer) or np.any(ti < 0) or np.any(ti >= shape[0]):
            raise ValueError("type_index must be an (N, J) integer array indexing the T types")
        _set(self, "type_index", ti.astype(np.int64))


# --------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class McConfig:
    """Monte Carlo settings (§C.2.7)."""

    n_photons: int = 10_000_000
    n_batches: int = 20
    chunk: int = 200_000
    epsilon_ddis: float = 0.1
    w_roulette: float = 1e-3
    max_order: int = 500
    polarized: bool = True
    seed: int = 0
    workers: int = 1

    def __post_init__(self):
        for nm in ("n_photons", "chunk", "max_order", "workers"):
            if int(getattr(self, nm)) < 1:
                raise ValueError(f"McConfig.{nm} must be >= 1")
        if int(self.n_batches) < 2:
            raise ValueError("McConfig.n_batches must be >= 2")
        if not 0 <= self.epsilon_ddis < 1:
            raise ValueError("McConfig.epsilon_ddis must lie in [0, 1)")
        if not 0 < self.w_roulette < 1:
            raise ValueError("McConfig.w_roulette must lie in (0, 1)")


@dataclass(frozen=True)
class SolverConfig:
    """Solver selection and switches (§C.3.3)."""

    solver: str = "pvc"
    closure: str = "auto"
    diagnose_closures: bool = False
    theta_max: float = THETA_MAX_DEFAULT
    n_out: Optional[int] = None
    wide_coarsen: int = 1
    delta_wide: float = 1.0
    compat: Optional[str] = None
    mc: McConfig = field(default_factory=McConfig)

    def __post_init__(self):
        if self.solver not in SOLVERS:
            raise ValueError(f"solver must be one of {SOLVERS}")
        if self.closure not in CLOSURES:
            raise ValueError(f"closure must be one of {CLOSURES}")
        if not self.theta_max > 0:
            raise ValueError("theta_max must be > 0")
        if self.n_out is not None and int(self.n_out) < 1:
            raise ValueError("n_out must be None or >= 1")
        if int(self.wide_coarsen) < 1 or int(self.wide_coarsen) != self.wide_coarsen:
            raise ValueError("wide_coarsen must be an integer >= 1")
        if not self.delta_wide >= 0:
            raise ValueError("delta_wide must be >= 0")
        if self.compat not in COMPAT_MODES:
            raise ValueError(f"compat must be one of {COMPAT_MODES}")
        if not isinstance(self.mc, McConfig):
            raise ValueError("mc must be a McConfig")


# --------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class OpticalTangent:
    """Directional inputs, shape (P, K, N, J) float64: row k is d(field)/d x_k (§C.3.3).

    Fields left as None are zero.
    """

    alpha: Optional[np.ndarray] = None
    omega: Optional[np.ndarray] = None
    g: Optional[np.ndarray] = None
    a_eq: Optional[np.ndarray] = None
    beta: Optional[np.ndarray] = None
    delta: Optional[np.ndarray] = None

    def __post_init__(self):
        given = {nm: getattr(self, nm) for nm in _OPT_FIELDS if getattr(self, nm) is not None}
        if not given:
            raise ValueError("OpticalTangent needs at least one field")
        shape = None
        for nm, a in given.items():
            a = np.array(a, dtype=float)
            if a.ndim != 4 or (shape is not None and a.shape != shape):
                raise ValueError(f"tangent field {nm} must have shape (P, K, N, J) like the others")
            if not np.all(np.isfinite(a)):
                raise ValueError(f"tangent field {nm} must be finite")
            shape = a.shape
            _set(self, nm, a)
        for nm in _OPT_FIELDS:
            if getattr(self, nm) is None:
                _set(self, nm, np.zeros(shape))

    @property
    def shape(self):
        return self.alpha.shape


# --------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class Returns:
    """Apparent backscatter by class and polarization (§C.3.3) [DECISION D31].

    `beta_par`, `beta_perp`, `beta_wide`: (P, N_out, J); `beta_ss`: (P, N, J); `ms_factor` (M) and
    `gamma2`: (P, N+1) at half-levels; `a_s`: (P, N); `closure_used`: (P,) str; `flags`: (P, N) uint16.
    `beta_order`, `stderr_par`, `stderr_perp` are for solver `mc` only.
    """

    beta_par: np.ndarray
    beta_perp: np.ndarray
    beta_ss: np.ndarray
    beta_wide: np.ndarray
    ms_factor: np.ndarray
    gamma2: np.ndarray
    a_s: np.ndarray
    closure_used: np.ndarray
    flags: np.ndarray
    kinds: tuple
    names: tuple
    r: np.ndarray
    meta: dict = field(default_factory=dict)
    beta_order: Optional[np.ndarray] = None
    stderr_par: Optional[np.ndarray] = None
    stderr_perp: Optional[np.ndarray] = None

    def __post_init__(self):
        P, n_out, J = np.shape(self.beta_par)
        if np.shape(self.beta_perp) != (P, n_out, J) or np.shape(self.beta_wide) != (P, n_out, J):
            raise ValueError("beta_par, beta_perp and beta_wide must share shape (P, N_out, J)")
        N = np.shape(self.beta_ss)[1]
        if np.shape(self.beta_ss) != (P, N, J) or N > n_out:
            raise ValueError("beta_ss must have shape (P, N, J) with N <= N_out")
        for nm, shp in (("ms_factor", (P, N + 1)), ("gamma2", (P, N + 1)), ("a_s", (P, N)),
                        ("flags", (P, N)), ("closure_used", (P,))):
            if np.shape(getattr(self, nm)) != shp:
                raise ValueError(f"{nm} must have shape {shp}, got {np.shape(getattr(self, nm))}")
        if len(self.kinds) != J or len(self.names) != J:
            raise ValueError("kinds and names need one entry per class")

    @property
    def mol_index(self) -> int:
        return self.kinds.index("molecular")

    def _pad(self, a):
        n_out = self.beta_par.shape[1]
        if a.shape[1] == n_out:
            return a
        pad = np.zeros((a.shape[0], n_out - a.shape[1]) + a.shape[2:], dtype=a.dtype)
        return np.concatenate([a, pad], axis=1)

    @property
    def total(self) -> np.ndarray:
        """Sum over class and polarization, (P, N_out)."""
        return (self.beta_par + self.beta_perp).sum(axis=-1)

    @property
    def total_by_class(self) -> np.ndarray:
        """Sum over polarization, (P, N_out, J)."""
        return self.beta_par + self.beta_perp

    @property
    def small_angle_ms(self) -> np.ndarray:
        """Small-angle multiple scattering: total - single - wide-angle, (P, N_out)."""
        return self.total - self._pad(self.beta_ss).sum(-1) - self.beta_wide.sum(-1)

    @property
    def delta_volume(self) -> np.ndarray:
        with np.errstate(divide="ignore", invalid="ignore"):
            return self.beta_perp.sum(-1) / self.beta_par.sum(-1)

    @property
    def delta_particulate(self) -> np.ndarray:
        m = np.array([k != "molecular" for k in self.kinds])
        with np.errstate(divide="ignore", invalid="ignore"):
            return self.beta_perp[..., m].sum(-1) / self.beta_par[..., m].sum(-1)

    @property
    def mol_par(self) -> np.ndarray:
        """HSRL molecular co-polarized return, (P, N_out)."""
        return self.beta_par[..., self.mol_index]


@dataclass(frozen=True)
class ReturnsTangent:
    """Directional derivatives of `Returns`: (P, K, N_out, J)."""

    d_beta_par: np.ndarray
    d_beta_perp: np.ndarray

    def __post_init__(self):
        if np.ndim(self.d_beta_par) != 4 or np.shape(self.d_beta_par) != np.shape(self.d_beta_perp):
            raise ValueError("d_beta_par and d_beta_perp must share shape (P, K, N_out, J)")

    @property
    def d_total(self) -> np.ndarray:
        return (self.d_beta_par + self.d_beta_perp).sum(axis=-1)
