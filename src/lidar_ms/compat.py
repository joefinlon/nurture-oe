"""Everything that knows about `multiscatter` 1.2.11 (appendix C §C.2.9, §C.4.1).

- Reader for the ASCII input format (`.in`) and conversion to `lidar_ms` schema objects.
- `compat="multiscatter-1.2.11"` variants of the kernels where the C code departs from §C.2
  (Dep1-Dep5, Dep8, Dep11 and the C code's arithmetic forms needed for regression).
- An optional oracle runner for a compiled `multiscatter` executable (not needed by the tests).

The C source is in `legacy/multiscatter-1.2.11/` (reference only, Apache-2.0).
"""
from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .constants import COMPAT_AIR_BSCAT_MAX_WAVELENGTH, RAYLEIGH_BSCAT_EXT
from .schema import Instrument, OpticalProfile, RangeGrid

COMPAT = "multiscatter-1.2.11"

# Classes built from one multiscatter particulate: extinction split by droplet_fraction and
# pristine_ice_fraction (anisotropic.c); the remainder is isotropic near 180 deg.
MS_KINDS = ("liquid", "ice", "other", "molecular")
MS_DESCRIPTORS = ("droplet", "yang", "isotropic", "isotropic")
MS_NAMES = ("droplet", "pristine_ice", "isotropic", "air")


@dataclass(frozen=True)
class MsCase:
    """A multiscatter input converted to `lidar_ms` objects."""

    grid: RangeGrid
    instrument: Instrument
    optics: OpticalProfile
    altitude: float
    z: np.ndarray                 # gate heights as in the file
    has_air: bool                 # >= 5 columns on every row
    has_wide_angle: bool          # >= 8 columns (ssa, g, ssa_air) on every row
    has_near_180: bool            # 10 columns (droplet and pristine-ice fractions) on every row
    comments: tuple


def _parse(text: str):
    comments, rows = [], []
    for line in text.splitlines():
        s = line.strip()
        if not s:
            continue
        (comments if s.startswith("#") else rows).append(s if s.startswith("#") else s.split())
    head = rows[0]
    n = int(head[0])
    body = rows[1:1 + n]
    if len(body) != n:
        raise ValueError(f"header says {n} gates, found {len(body)} data rows")
    vals = np.full((n, 10), np.nan)
    ncols = np.empty(n, dtype=int)
    for i, row in enumerate(body):
        v = [float(x) for x in row[:10]]
        if len(v) < 4:
            raise ValueError(f"garbled input at range gate {i + 1}: only {len(v)} values")
        vals[i, :len(v)] = v
        ncols[i] = len(v)
    hdr = dict(n=n, wavelength=float(head[1]), altitude=float(head[2]), rho_tr=float(head[3]),
               rho_fov=[float(x) for x in head[4:]])
    return hdr, vals, ncols, tuple(comments)


def read_ms_input(source, *, air_backscatter="compat", delta_air=0.0) -> MsCase:
    """Read a multiscatter 1.2.11 ASCII input (version 2.x header) as a one-profile case.

    Parameters
    ----------
    source : path or str
        File path, or the file content if it contains a newline.
    air_backscatter : {"compat", "rayleigh"}
        "compat": beta_air = 3/(8 pi) alpha_air, set to zero above 1 um as the C code does (Dep1).
        "rayleigh": 3/(8 pi) alpha_air at every wavelength.
    delta_air : float
        Single-scatter depolarization assigned to the molecular class (not in the file).

    Notes
    -----
    Column defaults follow `multiscatter_ascii.c`: fewer than 10 values -> isotropic near 180 deg;
    fewer than 8 -> no wide-angle properties (omega = 1, g = 0 are placeholders); fewer than 5 ->
    no air. Range from the instrument is r = |altitude - z|. Only the first field of view is used.
    """
    text = str(source) if "\n" in str(source) else Path(source).read_text()
    hdr, v, ncols, comments = _parse(text)
    n, lam, alt = hdr["n"], hdr["wavelength"], hdr["altitude"]
    z, ext, rad, S = v[:, 0], v[:, 1], v[:, 2], v[:, 3]
    has_air, has_wide, has_near = bool(np.all(ncols >= 5)), bool(np.all(ncols >= 8)), bool(np.all(ncols >= 10))
    ext_air = v[:, 4] if has_air else np.zeros(n)
    ssa = v[:, 5] if has_wide else np.ones(n)
    g = v[:, 6] if has_wide else np.zeros(n)
    ssa_air = v[:, 7] if has_wide else np.ones(n)
    fd = v[:, 8] if has_near else np.zeros(n)
    fi = v[:, 9] if has_near else np.zeros(n)

    if air_backscatter == "compat":
        k_air = RAYLEIGH_BSCAT_EXT if lam <= COMPAT_AIR_BSCAT_MAX_WAVELENGTH else 0.0
    elif air_backscatter == "rayleigh":
        k_air = RAYLEIGH_BSCAT_EXT
    else:
        raise ValueError("air_backscatter must be 'compat' or 'rayleigh'")

    frac = np.stack([fd, fi, 1.0 - fd - fi], axis=1)                  # (n, 3)
    alpha = np.concatenate([ext[:, None] * frac, ext_air[:, None]], axis=1)
    beta = np.concatenate([(ext / S)[:, None] * frac, (ext_air * k_air)[:, None]], axis=1)
    omega = np.stack([ssa, ssa, ssa, ssa_air], axis=1)
    gg = np.stack([g, g, g, np.zeros(n)], axis=1)
    a_eq = np.stack([rad, rad, rad, np.zeros(n)], axis=1)
    delta = np.zeros((n, 4))
    delta[:, 3] = delta_air
    optics = OpticalProfile(kinds=MS_KINDS, descriptors=MS_DESCRIPTORS, names=MS_NAMES,
                            alpha=alpha[None], omega=omega[None], g=gg[None], a_eq=a_eq[None],
                            beta=beta[None], delta=delta[None])
    r = np.abs(alt - z)
    instrument = Instrument(lam, hdr["rho_tr"], hdr["rho_fov"][0], name="multiscatter input")
    return MsCase(grid=RangeGrid(r[None], z=z), instrument=instrument, optics=optics, altitude=alt,
                  z=z, has_air=has_air, has_wide_angle=has_wide, has_near_180=has_near, comments=comments)


# --- compat kernels ----------------------------------------------------------------------------
def single_compat(alpha_tot, beta, dr):
    """`singlescatter.c` arithmetic: factor = (T_{i-1/2} - T_{i+1/2}) / tau_i, or T_{i-1/2} if tau_i = 0.

    Algebraically equal to Eq. (C.2); kept because the subtraction's rounding (relative error
    ~ eps / tau in near-vacuum bins) is part of the frozen references.
    """
    tau = 2.0 * alpha_tot * dr
    cum = np.concatenate([np.zeros_like(tau[:1]), np.cumsum(tau, axis=0)], axis=0)
    tu = np.exp(-cum)
    pos = np.real(tau) > 0
    fac = np.where(pos, (tu[:-1] - tu[1:]) / np.where(pos, tau, 1.0), tu[:-1])
    return beta * fac[:, None, :], tu, tau


# --- oracle --------------------------------------------------------------------------------------
def run_multiscatter(executable, input_path, algorithms=("fast", "none"), extra=()):
    """Run a compiled multiscatter executable; returns the numeric output, (n, ncol).

    Optional: a debugging aid only. The regression tests use the frozen `.ref` files.
    """
    cmd = [str(executable), "-quiet", "-algorithms", *algorithms, *extra, str(input_path)]
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=120, check=True)
    rows = [ln.split() for ln in p.stdout.splitlines() if ln.strip() and not ln.startswith("#")]
    return np.array(rows, dtype=float)


# --- compat pvc pieces (Dep2, Dep3; anisotropic.c, fast_small_angle.c) ------------------------
def _ms_erfcx(x):
    """exp(x^2) erfc(x) for x >= 0 from the rational approximation `ms_erfc` in anisotropic.c."""
    t = 1.0 / (1.0 + 0.5 * x)
    poly = (-1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 + t * (-0.18628806 + t * (
        0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))))
    return t * np.exp(poly)


def factor_yang_compat(gamma2):
    """H08 Eq. 37 with the C code's erfc approximation (relative error < 1.2e-7)."""
    from .constants import YANG_GAMMA0, YANG_W
    pos = np.real(gamma2) > 0
    g = np.sqrt(np.where(pos, gamma2, 1.0)) * 0.5 / YANG_GAMMA0
    f = (1.0 - YANG_W) + YANG_W * (1.0 - g * np.sqrt(np.pi) * _ms_erfcx(g))
    return np.where(pos, f, 1.0)


def factors_compat(gamma2, theta2, alpha, beta, kinds, descriptors):
    """`ms_anisotropic_factor`: one factor per bin and half-level, applied to every class.

    The particulate factor is the extinction-weighted mean over particulate classes (Dep3); air is
    isotropic and the result is weighted by backscatter between particulate and air, then applied to
    the particulate *and* the air return (Dep2).

    Returns F_minus, F_plus : (N, J, B), identical along J.
    """
    from .anisotropy import factor_droplet
    part = np.array([k != "molecular" for k in kinds])
    mol = ~part
    a_p = alpha[:, part].sum(axis=1)
    b_p = beta[:, part].sum(axis=1)
    b_air = beta[:, mol].sum(axis=1)
    has_p = np.real(a_p) > 0

    def at(g2):
        fpart = np.zeros_like(a_p * g2)
        for j in np.flatnonzero(part):
            th = np.where(np.isfinite(np.real(theta2[:, j])), theta2[:, j], 1.0)
            if descriptors[j] == "droplet":
                fj = factor_droplet(g2, th)
            elif descriptors[j] == "yang":
                fj = factor_yang_compat(g2)
            elif descriptors[j] == "isotropic":
                fj = 1.0
            else:
                raise ValueError("compat mode supports the droplet, yang and isotropic descriptors only")
            fpart = fpart + alpha[:, j] / np.where(has_p, a_p, 1.0) * fj
        fpart = np.where(has_p, fpart, 1.0)
        den = b_p + b_air
        okd = np.real(den) > 0
        return np.where(okd, (b_p * fpart + b_air) / np.where(okd, den, 1.0), 1.0)

    J = alpha.shape[1]
    fm = at(gamma2[:-1])
    fp = at(gamma2[1:])
    return (np.repeat(fm[:, None, :], J, axis=1), np.repeat(fp[:, None, :], J, axis=1))


def bin_means_compat(beta, tu, Mt, F_minus, F_plus, valid, tau):
    """`fast_small_angle.c` step 6: where no forward-scattered light is in view at r_{i+1/2},
    the C code returns T_{u,i-1/2} (1 - exp(-tau)) / tau (zero if tau = 0), ignoring M_{i-1/2}.
    Elsewhere the logarithmic-mean form of Eq. (C.4), algebraically equal to the C code's expression.
    """
    from .pvc import bin_means
    lm = bin_means(beta, tu, Mt, F_minus, F_plus)
    pos = np.real(tau) > 0
    ss = np.where(pos, tu[:-1] * (1.0 - np.exp(-tau)) / np.where(pos, tau, 1.0), 0.0)
    return np.where(valid[1:, None, :], lm, beta * ss[:, None, :])
