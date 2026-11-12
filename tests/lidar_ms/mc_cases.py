"""Cases shared by the Monte Carlo tests (X5-X9, P1-P4)."""
import numpy as np

from lidar_ms import Instrument, McConfig, OpticalProfile, PhaseTables, RangeGrid, SolverConfig
from lidar_ms.compat import read_ms_input
from lidar_ms.mc import fixtures as fx
from lidar_ms.mc import tables as tb
from msref import REF_DIR

HALO = Instrument.from_full_angles(532e-9, 0.8e-3, 1.0e-3)


def mc_config(n, polarized=False, seed=0, epsilon=0.1, chunk=50_000, batches=20, **kw):
    return SolverConfig(solver="mc", mc=McConfig(n_photons=n, n_batches=batches, chunk=chunk, polarized=polarized,
                                                 seed=seed, epsilon_ddis=epsilon, **kw))


def phase(theta, mats, type_index):
    return PhaseTables(theta=theta, type_index=np.asarray(type_index), **fx.stack(*mats))


def consistent_beta(optics, ph):
    """Optics whose beta equals alpha omega P11(pi) / 4 pi of the tables."""
    return optics.replace(beta=tb.backscatter_from_tables(optics.alpha[0], optics.omega[0], ph)[None])


def truncate(c, n_bins):
    """First n_bins of a one-profile multiscatter case: (grid, optics)."""
    o = c.optics
    kw = {nm: getattr(o, nm)[:, :n_bins] for nm in ("alpha", "omega", "g", "a_eq", "beta", "delta", "near_index")}
    return RangeGrid(c.grid.r[:, :n_bins]), o.replace(**kw)


def lobe_case(name, n_extra=2, delta_pi=0.0, instrument=None, iso_fraction=None):
    """A multiscatter input with its particulate classes on the Gaussian-lobe fixture (Theta = lambda/(pi a)),
    air on Rayleigh, beta from the tables, grid ending n_extra bins below the last particulate bin.

    iso_fraction None: half of the extinction in the lobe, half scattered isotropically (omega = 1).
    iso_fraction f: half in the lobe, f isotropic, the rest absorbed (omega = 1/2 + f): light that leaves
    the lobe does not come back, the configuration in which A3 is the exact wide-field limit."""
    c = read_ms_input(REF_DIR / f"{name}.in")
    part = c.optics.alpha[0, :, :3].sum(-1) > 0
    last = int(np.flatnonzero(part)[-1])
    grid, opt = truncate(c, last + 1 + n_extra)
    a = float(c.optics.a_eq[0, part, 0][0])
    lobe = c.instrument.wavelength / (np.pi * a)
    th = fx.angle_grid(lobe / 20)
    frac = 0.5 if iso_fraction is None else 0.5 / (0.5 + iso_fraction)
    if iso_fraction is not None:
        om = np.array(opt.omega)
        om[..., :3] = 0.5 + iso_fraction
        opt = opt.replace(omega=om)
    ph = phase(th, [fx.gaussian_lobe(th, lobe, lobe_fraction=frac, delta_pi=delta_pi), fx.rayleigh(th)],
               np.tile([0, 0, 0, 1], (grid.n_bins, 1)))
    return grid, instrument or c.instrument, consistent_beta(opt, ph), ph, lobe


def tau_mid(grid, optics):
    a = optics.alpha[0].sum(-1)
    dr = grid.dr[0]
    return np.cumsum(a * dr) - 0.5 * a * dr
