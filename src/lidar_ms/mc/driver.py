"""Solver `mc` driver: batches, seeds, statistics and assembly of `Returns` (appendix C §C.2.7, §C.3.5).

Each batch has its own `SeedSequence.spawn` child, so batches are independent and reproducible on any
number of workers. Within a batch, packets run in chunks to extinction, one generation (scattering
event) at a time, with dead packets compacted away each generation. Standard errors come from the
batch means.
"""
from __future__ import annotations

import warnings
from concurrent.futures import ProcessPoolExecutor
from typing import NamedTuple

import numpy as np

from .. import grid
from ..schema import Instrument, McConfig, PhaseTables
from . import tables as tb
from . import transport as tr
from .estimator import N_ORDERS, local_estimate, tally

BETA_TABLE_TOL = 1e-3          # relative mismatch of beta and alpha omega P11(pi) / 4 pi that warns


class Tables(NamedTuple):
    theta: np.ndarray            # (A,)
    P: np.ndarray                # (6, T, A)
    theta_of_u: np.ndarray       # (T, U)
    ddis_p: np.ndarray           # (A,)
    ddis_theta_of_u: np.ndarray  # (U,)


def build_tables(phase: PhaseTables, types_used) -> Tables:
    P = tb.elements(phase)
    env, env_icdf = tb.ddis_envelope(phase.theta, phase.P11, types_used)
    return Tables(phase.theta, P, tb.inverse_cdf(phase.theta, phase.P11), env, env_icdf)


def run_chunk(n, med, tabs, inst, mc: McConfig, rng, acc, peak, pol_angle=0.0, record=None):
    """Transport n packets to extinction, adding their local estimates to acc (N, J, 3, 2) and peak (N,).

    `record`, if a list, receives one dict per generation (order, bin, r_app, dist, layer, contribution
    total): a test hook, off in production.
    """
    pol = mc.polarized
    p = tr.launch(n, inst.rho_tr, rng, polarized=pol, pol_angle=pol_angle)
    while p["W"].size:
        p = tr.move(p, med, rng)
        if not p["W"].size:
            break
        layer = tb.layer_of_z(med, p["x"][2])
        alive = med.sigma[layer].sum(axis=1) > 0
        if not alive.all():
            p = tr.compact(p, alive)
            layer = layer[alive]
            if not p["W"].size:
                break
        d, bin_, co, cross, r_app, dist = local_estimate(p, med, tabs, layer, inst.rho_fov, inst.receiver,
                                                         polarized=pol, pol_angle=pol_angle)
        tally(acc, peak, bin_, p["order"], co, cross)
        if record is not None:
            record.append(dict(order=p["order"].copy(), bin=bin_, r_app=r_app, dist=dist, layer=layer,
                               total=co.sum(axis=0) + (0 if cross is None else cross.sum(axis=0))))
        p, _ = tr.redirect(p, med, tabs, layer, d, mc.epsilon_ddis, rng, polarized=pol)
        p = tr.roulette(p, mc.w_roulette, rng)
        p = tr.compact(p, (p["order"] < mc.max_order) & (p["W"] > 0))


def run_batch(args):
    """One batch: (n_photons, seed sequence, med, tabs, inst, mc, pol_angle) -> (acc, peak)."""
    n, seed, med, tabs, inst, mc, pol_angle = args
    rng = np.random.Generator(np.random.PCG64(seed))
    N, J = med.sigma.shape[0] - 1, med.sigma.shape[1]
    acc = np.zeros((N, J, N_ORDERS, 2))
    peak = np.zeros(N)
    done = 0
    while done < n:
        m = min(mc.chunk, n - done)
        run_chunk(m, med, tabs, inst, mc, rng, acc, peak, pol_angle=pol_angle)
        done += m
    return acc, peak


def batch_sizes(n_photons, n_batches):
    base, extra = divmod(int(n_photons), int(n_batches))
    return [base + (1 if b < extra else 0) for b in range(int(n_batches))]


def run_profile(r, alpha, omega, a_eq, kinds, phase: PhaseTables, instrument: Instrument, mc: McConfig,
                theta_max, beta=None, pol_angle=0.0):
    """Monte Carlo for one profile.

    Parameters
    ----------
    r : (N,) bin centres; alpha, omega, a_eq, beta : (N, J)
    pol_angle : float
        Transmitted polarization angle from x-hat [rad]; co- and cross-polarized are referred to it.
        A test hook (P3); 0 in production.

    Returns
    -------
    dict with
      beta : (N, J, 3, 2) apparent backscatter by class, order (1, 2, >=3), polarization (co, cross)
      batches : (n_batches, N, J, 3, 2) the same per batch
      se : (N, J, 3, 2) batch-mean standard errors
      spike : (N,) largest single contribution / bin total
    """
    r = np.asarray(r, dtype=float)
    N, J = alpha.shape
    tb.check_lobe_resolution(phase.theta, alpha, a_eq, kinds, instrument.wavelength, theta_max)
    if phase.type_index.shape != (N, J):
        raise ValueError(f"PhaseTables.type_index has shape {phase.type_index.shape}, expected {(N, J)}")
    if beta is not None:
        bt = tb.backscatter_from_tables(alpha, omega, phase)
        bad = (alpha > 0) & (np.abs(beta - bt) > BETA_TABLE_TOL * np.maximum(bt, 1e-300))
        if bad.any():
            warnings.warn("OpticalProfile.beta differs from alpha omega P11(pi) / 4 pi of the PhaseTables by more "
                          f"than {BETA_TABLE_TOL:g} in {int(bad.sum())} bin-class entries; the Monte Carlo uses the "
                          "tables", stacklevel=3)
    r_half = grid.half_levels(r[None])[0]
    dr = grid.bin_widths(r[None])[0]
    med = tb.build_medium(r_half, alpha, omega, phase.type_index)
    used = phase.type_index[(alpha * omega) > 0]
    tabs = build_tables(phase, used if used.size else phase.type_index.ravel())
    seeds = np.random.SeedSequence(mc.seed).spawn(mc.n_batches)
    sizes = batch_sizes(mc.n_photons, mc.n_batches)
    args = [(n, s, med, tabs, instrument, mc, pol_angle) for n, s in zip(sizes, seeds)]
    if mc.workers > 1:
        with ProcessPoolExecutor(max_workers=mc.workers) as ex:
            res = list(ex.map(run_batch, args))
    else:
        res = [run_batch(a) for a in args]
    scale = instrument.kappa / dr[:, None, None, None]
    batches = np.stack([acc * scale / n for (acc, _), n in zip(res, sizes)])
    w = np.asarray(sizes, float) / np.sum(sizes)
    beta = np.tensordot(w, batches, axes=1)
    se = batches.std(axis=0, ddof=1) / np.sqrt(len(sizes))
    total = sum(acc.sum(axis=(1, 2, 3)) for acc, _ in res)
    peak = np.max([pk for _, pk in res], axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        spike = np.where(total > 0, peak / total, 0.0)
    return dict(beta=beta, batches=batches, se=se, spike=spike)


def simulate_mc(r, optics, instrument, config, phase):
    """Solver `mc` for P profiles sharing one PhaseTables; returns the dict `forward.assemble` expects
    plus beta_order and standard errors, bin-major (N, J, P)."""
    if phase is None:
        raise ValueError("solver 'mc' needs PhaseTables (argument `phase`)")
    if not isinstance(phase, PhaseTables):
        raise TypeError("phase must be a PhaseTables")
    if config.n_out is not None and int(config.n_out) > r.shape[1]:
        raise ValueError("solver 'mc' tallies on the input grid only; n_out must be None or <= N")
    if optics.is_complex:
        raise ValueError("solver 'mc' takes real inputs (no jvp)")
    mc = config.mc
    P, N, J = optics.shape
    outs = [run_profile(r[b], optics.alpha[b], optics.omega[b], optics.a_eq[b], optics.kinds, phase, instrument,
                        mc, config.theta_max, beta=optics.beta[b]) for b in range(P)]
    beta = np.stack([o["beta"] for o in outs], axis=-1)          # (N, J, 3, 2, P)
    se = np.stack([np.std(o["batches"].sum(axis=3), axis=0, ddof=1) / np.sqrt(mc.n_batches) for o in outs],
                  axis=-1)                                        # (N, J, 2, P)
    spike = np.stack([o["spike"] for o in outs], axis=-1)        # (N, P)
    return dict(beta_par=beta[:, :, :, 0].sum(axis=2), beta_perp=beta[:, :, :, 1].sum(axis=2),
                beta_ss=beta[:, :, 0].sum(axis=2), beta_order=beta.sum(axis=3),
                stderr_par=se[:, :, 0], stderr_perp=se[:, :, 1], spike=spike)
