"""Range-bin geometry (appendix C §C.2.1). Pure functions on arrays whose last axis is range."""

import numpy as np

from .constants import RANGE_SPACING_TOL


def bin_widths(r):
    """Bin widths: Delta r_i = |r_{i+1} - r_{i-1}| / 2, one-sided at the ends (§C.2.1).

    Parameters
    ----------
    r : array_like, shape (..., N), N >= 2
        Bin centres along the beam.
    """
    r = np.asarray(r, dtype=float)
    d = np.empty_like(r)
    d[..., 1:-1] = 0.5 * np.abs(r[..., 2:] - r[..., :-2])
    d[..., 0] = np.abs(r[..., 1] - r[..., 0])
    d[..., -1] = np.abs(r[..., -1] - r[..., -2])
    return d


def half_levels(r):
    """Half-level ranges, shape (..., N+1): r_{-1/2} = r_0 - |r_1 - r_0|/2, r_{i+1/2} = r_i + Delta r_i / 2."""
    r = np.asarray(r, dtype=float)
    out = np.empty(r.shape[:-1] + (r.shape[-1] + 1,))
    out[..., 0] = r[..., 0] - 0.5 * np.abs(r[..., 1] - r[..., 0])
    out[..., 1:] = r + 0.5 * bin_widths(r)
    return out


def is_regular(r, tol=RANGE_SPACING_TOL):
    """True where the spacing is monotonic and constant to within `tol` (as multiscatter's test).

    Returns a boolean array of shape r.shape[:-1].
    """
    d = np.diff(np.asarray(r, dtype=float), axis=-1)
    dmax, dmin = d.max(axis=-1), d.min(axis=-1)
    monotonic = dmax * dmin > 0
    amax, amin = np.abs(dmax), np.abs(dmin)
    return monotonic & (amax <= amin * (1 + tol)) & (amin <= amax * (1 + tol))


def coarsen(a, k, axis=0, fill=0.0):
    """Mean over groups of `k` consecutive bins along `axis`, padding the end with `fill` (§C.3.5).

    Conservative for densities: the coarse value times k Delta r equals the sum over the group.
    """
    a = np.moveaxis(np.asarray(a), axis, 0)
    n = a.shape[0]
    nc = -(-n // k)
    if nc * k > n:
        a = np.concatenate([a, np.full((nc * k - n,) + a.shape[1:], fill, dtype=a.dtype)], axis=0)
    out = a.reshape((nc, k) + a.shape[1:]).mean(axis=1)
    return np.moveaxis(out, 0, axis)


def coarsen_centres(r, k):
    """Bin centres of the k-fold coarsened regular grid along axis 0 (end padded by extrapolation)."""
    r = np.asarray(r)
    n = r.shape[0]
    nc = -(-n // k)
    if nc * k > n:
        d = r[1] - r[0]
        r = np.concatenate([r, r[-1:] + d * np.arange(1, nc * k - n + 1).reshape((-1,) + (1,) * (r.ndim - 1))], axis=0)
    return r.reshape((nc, k) + r.shape[1:]).mean(axis=1)
