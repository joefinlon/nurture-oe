"""Class -> bin quantities (appendix C §C.2.2). Pure functions on bin-major arrays.

Internal layout: per-class arrays have shape (N, J, B), bin arrays (N, B); B is the batch axis
(profiles x tangents). Branches use real parts only (complex-step safe, §C.2.8).
"""

import numpy as np


def lobe_variance(a_eq, wavelength):
    """Theta_j^2 = (lambda / (pi a_j))^2, the squared forward-lobe width (H06). Inf where a_eq <= 0."""
    pos = np.real(a_eq) > 0
    safe = np.where(pos, a_eq, 1.0)
    th2 = (wavelength / (np.pi * safe)) ** 2
    return np.where(pos, th2, np.inf)


def lobe_mask(theta2, particulate, theta_max, all_particulate=False):
    """Lobe-forming classes: particulate and Theta_j <= theta_max (§C.2.2).

    Parameters
    ----------
    theta2 : (N, J, B)
    particulate : (J,) bool
        False for the molecular class.
    all_particulate : bool
        `compat` (multiscatter 1.2.11 has no lobe test in `fast_small_angle.c`): every particulate
        class with finite Theta is lobe-forming.
    """
    part = np.asarray(particulate)[None, :, None]
    finite = np.isfinite(np.real(theta2))
    if all_particulate:
        return part & finite
    return part & finite & (np.real(theta2) <= theta_max ** 2)


def pvc_bin(alpha, theta2, lobe):
    """Bin quantities for `pvc`, Eq. (C.1).

    Returns
    -------
    alpha_c : (N, B)  extinction of lobe-forming classes
    alpha_n : (N, B)  extinction of all other classes, incl. molecular
    q : (N, B)        sum over lobe-forming classes of alpha_j Theta_j^2
    """
    zero = np.zeros((), dtype=alpha.dtype)
    alpha_c = np.where(lobe, alpha, zero).sum(axis=1)
    alpha_n = np.where(lobe, zero, alpha).sum(axis=1)
    q = np.where(lobe, alpha * np.where(lobe, theta2, 0.0), zero).sum(axis=1)
    return alpha_c, alpha_n, q
