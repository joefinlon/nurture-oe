"""Solver `single`: Beer-Lambert bin means, Eq. (C.2) (appendix C §C.2.3)."""

import numpy as np


def two_way_depth(alpha_tot, dr):
    """tau_i = 2 alpha_i Delta r_i, (N, B)."""
    return 2.0 * alpha_tot * dr


def unscattered_transmittance(tau):
    """Two-way unscattered transmittance at half-levels, (N+1, B); T_{u,-1/2} = 1."""
    cum = np.concatenate([np.zeros_like(tau[:1]), np.cumsum(tau, axis=0)], axis=0)
    return np.exp(-cum)


def bin_mean_factor(tau):
    """(1 - exp(-tau)) / tau, -> 1 as tau -> 0. Complex-safe."""
    pos = np.real(tau) > 0
    safe = np.where(pos, tau, 1.0)
    return np.where(pos, -np.expm1(-safe) / safe, 1.0)


def single(alpha_tot, beta, dr):
    """Single-scatter apparent backscatter per class, Eq. (C.2).

    Parameters
    ----------
    alpha_tot : (N, B)  total extinction (all classes, incl. molecular)
    beta : (N, J, B)    backscatter coefficient per class
    dr : (N, B)         bin widths

    Returns
    -------
    beta_ss : (N, J, B)
    tu : (N+1, B)  two-way unscattered transmittance at half-levels
    tau : (N, B)   two-way bin optical depth
    """
    tau = two_way_depth(alpha_tot, dr)
    tu = unscattered_transmittance(tau)
    fac = tu[:-1] * bin_mean_factor(tau)
    return beta * fac[:, None, :], tu, tau
