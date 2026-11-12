"""Schema sketch for lidar_ms (Phase 3). Standard library + NumPy only."""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Literal, Optional
import numpy as np

Kind = Literal["liquid", "ice", "other", "molecular"]
Descriptor = Literal["droplet", "yang", "isotropic", "table"]


def _real(x):
    return np.real(np.asarray(x))


@dataclass(frozen=True)
class Instrument:
    """Lidar geometry. Angles in radians, lengths in metres."""
    wavelength: float
    rho_tr: float                      # transmitter 1/e half-width
    rho_fov: float                     # receiver half-angle
    receiver: Literal["top_hat", "gaussian"] = "top_hat"
    name: str = ""

    def __post_init__(self):
        if not (self.wavelength > 0 and self.rho_tr > 0 and self.rho_fov > 0):
            raise ValueError("wavelength, rho_tr and rho_fov must be positive")

    @classmethod
    def from_full_angles(cls, wavelength, divergence_full_e2, fov_full, **kw):
        """Build from a 1/e^2 full-angle divergence and a full-angle field of view."""
        return cls(wavelength, 0.5 * divergence_full_e2 / np.sqrt(2.0), 0.5 * fov_full, **kw)

    @property
    def kappa(self) -> float:
        if self.receiver == "top_hat":
            return 1.0 / (1.0 - np.exp(-(self.rho_fov / self.rho_tr) ** 2))
        return 1.0 + (self.rho_tr / self.rho_fov) ** 2


@dataclass(frozen=True)
class RangeGrid:
    """Bin centres along the beam, shape (P, N), increasing away from the instrument."""
    r: np.ndarray
    z: Optional[np.ndarray] = None     # altitude of bin centres, for provenance only

    def __post_init__(self):
        r = np.asarray(self.r, dtype=float)
        if r.ndim != 2 or r.shape[1] < 3:
            raise ValueError("r must have shape (P, N) with N >= 3")
        if np.any(np.diff(r, axis=1) <= 0) or np.any(r[:, 0] <= 0):
            raise ValueError("r must be positive and strictly increasing along range")

    @property
    def dr(self) -> np.ndarray:
        r = np.asarray(self.r, dtype=float)
        d = np.empty_like(r)
        d[:, 1:-1] = 0.5 * (r[:, 2:] - r[:, :-2])
        d[:, 0] = r[:, 1] - r[:, 0]
        d[:, -1] = r[:, -1] - r[:, -2]
        return d


@dataclass(frozen=True)
class OpticalProfile:
    """Per-bin, per-class optical properties, shape (P, N, J). float64 or complex128."""
    kinds: tuple[Kind, ...]
    descriptors: tuple[Descriptor, ...]
    alpha: np.ndarray                  # extinction, m-1
    omega: np.ndarray                  # single-scattering albedo
    g: np.ndarray                      # asymmetry factor
    a_eq: np.ndarray                   # equivalent-area radius, m (ignored for molecular)
    beta: np.ndarray                   # backscatter coefficient, m-1 sr-1
    delta: np.ndarray                  # single-scatter linear depolarization ratio
    near_index: Optional[np.ndarray] = None   # int (P, N, J): row of the class's near-180 table
    names: tuple[str, ...] = field(default=())

    def __post_init__(self):
        J = len(self.kinds)
        if self.kinds.count("molecular") != 1 or len(self.descriptors) != J:
            raise ValueError("need exactly one molecular class and one descriptor per class")
        shape = np.shape(self.alpha)
        if len(shape) != 3 or shape[2] != J:
            raise ValueError("arrays must have shape (P, N, J)")
        for nm in ("omega", "g", "a_eq", "beta", "delta"):
            if np.shape(getattr(self, nm)) != shape:
                raise ValueError(f"{nm} has shape {np.shape(getattr(self, nm))}, expected {shape}")
        if np.any(_real(self.alpha) < 0) or np.any(_real(self.beta) < 0):
            raise ValueError("alpha and beta must be non-negative")
        if np.any((_real(self.omega) < 0) | (_real(self.omega) > 1)):
            raise ValueError("omega must lie in [0, 1]")
        if np.any((_real(self.delta) < 0) | (_real(self.delta) > 1)):
            raise ValueError("delta must lie in [0, 1]")


if __name__ == "__main__":
    halo = Instrument.from_full_angles(532e-9, 0.8e-3, 1.0e-3, name="HALO 532 nm")
    print("HALO rho_tr = %.4f mrad, rho_fov = %.2f mrad, kappa = %.4f" % (1e3 * halo.rho_tr, 1e3 * halo.rho_fov, halo.kappa))
    P, N = 2, 120
    grid = RangeGrid(r=9000.0 + 15.0 * np.arange(N)[None, :] * np.ones((P, 1)))
    shape = (P, N, 3)
    opt = OpticalProfile(kinds=("liquid", "ice", "molecular"), descriptors=("droplet", "isotropic", "isotropic"),
                         alpha=np.zeros(shape), omega=np.ones(shape), g=np.zeros(shape), a_eq=np.full(shape, 8e-6),
                         beta=np.zeros(shape), delta=np.zeros(shape))
    print("grid dr[0, :3] =", grid.dr[0, :3], "| optics shape", np.shape(opt.alpha))
    try:
        OpticalProfile(kinds=("liquid", "ice"), descriptors=("droplet", "isotropic"), alpha=np.zeros((P, N, 2)),
                       omega=np.ones((P, N, 2)), g=np.zeros((P, N, 2)), a_eq=np.ones((P, N, 2)),
                       beta=np.zeros((P, N, 2)), delta=np.zeros((P, N, 2)))
    except ValueError as e:
        print("validation works:", e)
