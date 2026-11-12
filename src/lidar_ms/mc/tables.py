"""Precomputed tables for solver `mc` (appendix C §C.2.7, §C.3.5).

- Phase-matrix checks and evaluation on the common angle grid of `PhaseTables` (P_ij linear in theta).
- Inverse CDFs of P11 per type on 2^16 uniform probability nodes. Theta is tabulated rather than
  mu = cos(theta): it keeps full precision inside narrow forward lobes. A draw is one gather and a
  linear interpolation, with no search.
- The detector importance function p_ddis (normalized envelope of P11 over the types in use) and its
  inverse CDF.
- The vertical optical depth tau_z(z) of the layered medium, with the floor alpha >= 1e-12 m-1 so that
  it is strictly increasing and invertible by interpolation in both directions.
"""
from typing import NamedTuple

import numpy as np

from ..schema import PhaseTables

N_NODES = 2 ** 16
REFINE = 8
ALPHA_FLOOR = 1e-12          # m-1
LOBE_STEPS = 10              # first angle spacing <= narrowest forward lobe / 10 (§C.3.3)


# --- phase matrices ----------------------------------------------------------------------------
def elements(phase: PhaseTables):
    """(6, T, A) array of P11, P12, P22, P33, P34, P44."""
    return np.stack([phase.P11, phase.P12, phase.P22, phase.P33, phase.P34, phase.P44])


def check_lobe_resolution(theta, alpha, a_eq, kinds, wavelength, theta_max):
    """Raise if the first angle spacing exceeds one tenth of the narrowest forward lobe present.

    alpha, a_eq : (N, J). Lobe width Theta_j = lambda / (pi a_j) for lobe-forming particulate classes.
    """
    part = np.array([k != "molecular" for k in kinds])
    present = (np.real(alpha) > 0) & part[None, :] & (np.real(a_eq) > 0)
    if not present.any():
        return
    lobes = wavelength / (np.pi * np.real(a_eq[present]))
    lobes = lobes[lobes <= theta_max]
    if lobes.size and theta[1] > lobes.min() / LOBE_STEPS:
        raise ValueError(f"PhaseTables.theta: first spacing {theta[1]:.3g} rad exceeds one tenth of the narrowest "
                         f"forward lobe ({lobes.min():.3g} rad) (§C.3.3)")


def refine(theta, m=REFINE):
    """Insert m - 1 equally spaced points into every interval of theta."""
    t = theta[:-1, None] + np.diff(theta)[:, None] * (np.arange(m) / m)[None, :]
    return np.concatenate([t.ravel(), theta[-1:]])


def cumulative_sphere_integral(theta, f):
    """C(theta_k) = int_0^theta_k f sin, f linear in theta between nodes (closed form per interval)."""
    t0, t1 = theta[:-1], theta[1:]
    h = t1 - t0
    i0 = np.cos(t0) - np.cos(t1)
    i1 = (np.sin(t1) - np.sin(t0) - h * np.cos(t1)) / h
    seg = f[..., :-1] * i0 + (f[..., 1:] - f[..., :-1]) * i1
    return np.concatenate([np.zeros(f.shape[:-1] + (1,)), np.cumsum(seg, axis=-1)], axis=-1)


def inverse_cdf(theta, p11, n_nodes=N_NODES):
    """theta_of_u (..., n_nodes) for u = k / (n_nodes - 1); p11 (..., A) linear in theta on `theta`."""
    tf = refine(theta)
    p = np.stack([np.interp(tf, theta, row) for row in np.atleast_2d(p11)])
    cdf = cumulative_sphere_integral(tf, p)
    cdf = cdf / cdf[:, -1:]
    u = np.linspace(0.0, 1.0, n_nodes)
    out = np.stack([np.interp(u, c, tf) for c in cdf])
    return out.reshape(np.shape(p11)[:-1] + (n_nodes,))


def sample_theta(theta_of_u, rows, xi):
    """Draw theta from the inverse-CDF table: rows (n,) int (type), xi (n,) uniform on [0, 1)."""
    n_nodes = theta_of_u.shape[-1]
    pos = xi * (n_nodes - 1)
    i = np.minimum(pos.astype(np.int64), n_nodes - 2)
    f = pos - i
    if theta_of_u.ndim == 1:
        return theta_of_u[i] * (1 - f) + theta_of_u[i + 1] * f
    return theta_of_u[rows, i] * (1 - f) + theta_of_u[rows, i + 1] * f


def interp_index(theta_grid, th):
    """Index and weight for linear interpolation in theta: (i, f) with value = P[i] (1 - f) + P[i + 1] f."""
    i = np.clip(np.searchsorted(theta_grid, th, side="right") - 1, 0, theta_grid.size - 2)
    f = (th - theta_grid[i]) / (theta_grid[i + 1] - theta_grid[i])
    return i, np.clip(f, 0.0, 1.0)


def lookup(P, types, i, f):
    """P (..., T, A) evaluated at (i, f) for per-packet types (n,): returns (..., n)."""
    return P[..., types, i] * (1 - f) + P[..., types, i + 1] * f


def p11_at(theta_grid, p11, th):
    """Interpolate a single P11 row (A,) at th."""
    return np.interp(th, theta_grid, p11)


def ddis_envelope(theta, P11, types_used):
    """Normalized envelope max P11 over the types in use (int p_ddis dOmega = 4 pi), and its inverse CDF."""
    env = P11[np.asarray(sorted(set(int(t) for t in types_used)))].max(axis=0)
    env = env * 2.0 / cumulative_sphere_integral(theta, env)[-1]
    return env, inverse_cdf(theta, env)


def backscatter_from_tables(alpha, omega, phase: PhaseTables):
    """beta_j = alpha_j omega_j P11_j(pi) / (4 pi), (N, J), from the tables (for consistency checks)."""
    p180 = phase.P11[phase.type_index, -1]
    return alpha * omega * p180 / (4 * np.pi)


# --- layered medium ----------------------------------------------------------------------------
class Medium(NamedTuple):
    """Layered medium along the beam axis z (range); layer 0 is the gap from the instrument to the
    first half-level, layers 1..N are the range bins."""

    z_nodes: np.ndarray      # (N + 2,) 0, r_half_0, ..., r_half_N
    tau_nodes: np.ndarray    # (N + 2,) tau_z at z_nodes
    alpha: np.ndarray        # (N + 1,) total extinction per layer (floored)
    sigma: np.ndarray        # (N + 1, J) scattering coefficient alpha_j omega_j per class
    types: np.ndarray        # (N + 1, J) phase-table type per class (layer 0: copy of layer 1)
    r_half: np.ndarray       # (N + 1,) output bin edges
    z_max: float


def build_medium(r_half, alpha_j, omega_j, type_index):
    """Medium from bin edges (N + 1,), alpha_j, omega_j (N, J) and type_index (N, J)."""
    r_half = np.asarray(r_half, dtype=float)
    N, J = alpha_j.shape
    a_tot = np.maximum(alpha_j.sum(axis=1), ALPHA_FLOOR)
    alpha = np.concatenate([[ALPHA_FLOOR], a_tot])
    sigma = np.concatenate([np.zeros((1, J)), alpha_j * omega_j], axis=0)
    z_nodes = np.concatenate([[0.0], r_half])
    tau_nodes = np.concatenate([[0.0], np.cumsum(alpha * np.diff(z_nodes))])
    types = np.concatenate([type_index[:1], type_index], axis=0).astype(np.int64)
    return Medium(z_nodes, tau_nodes, alpha, sigma, types, r_half, float(r_half[-1]))


def tau_of_z(med: Medium, z):
    return np.interp(z, med.z_nodes, med.tau_nodes)


def z_of_tau(med: Medium, tau):
    return np.interp(tau, med.tau_nodes, med.z_nodes)


def layer_of_z(med: Medium, z):
    """Layer index 0..N of height z (clipped)."""
    return np.clip(np.searchsorted(med.z_nodes, z, side="right") - 1, 0, med.alpha.size - 1)
