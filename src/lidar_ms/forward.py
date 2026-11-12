"""`simulate()`: solver dispatch and assembly of `Returns` (appendix C §C.3.2).

The public layout is (P, N, J). Kernels work bin-major: (N, J, B) and (N, B), with B the batch
axis (profiles, or profiles x tangents inside `jvp`) (§C.3.5).
"""
from __future__ import annotations

import hashlib
from dataclasses import asdict, replace

import numpy as np

from . import __version__, anisotropy, compat, grid, mixing, polarization, pvc, tdts
from .constants import EXTINGUISHED_TU, LOBE_MIXTURE_RATIO, TABLE_TAIL_WEIGHT, flag
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
def solve(r, opt, kinds, descriptors, instrument, config, near=None, near_index=None, names=None):
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
    tables = near_tables(near, kinds, names) if near is not None or "table" in descriptors else [None] * J
    for j, d in enumerate(descriptors):
        if d == "table" and tables[j] is None:
            raise ValueError(f"class {j} uses descriptor 'table' but no NearBackscatter was given")
    dr = grid.bin_widths(r.T).T
    alpha = opt["alpha"]
    alpha_tot = alpha.sum(axis=1)
    beta = opt["beta"]

    if is_compat:
        beta_ss, tu, tau = compat.single_compat(alpha_tot, beta, dr)
    else:
        beta_ss, tu, tau = single(alpha_tot, beta, dr)
    theta2 = mixing.lobe_variance(opt["a_eq"], instrument.wavelength)
    part = np.array([k != "molecular" for k in kinds])
    lobe = mixing.lobe_mask(theta2, part, config.theta_max, all_particulate=is_compat)
    flags = np.where(np.real(tu[:-1]) < EXTINGUISHED_TU, flag("EXTINGUISHED"), 0).astype(np.uint16)

    if config.solver == "single":
        beta_ms = beta_ss
        Mt = np.zeros((N + 1, B), dtype=beta_ss.dtype)
        gamma2 = np.zeros((N + 1, B), dtype=beta_ss.dtype)
        Fm = Fp = np.ones((N, J, B))
    elif config.solver in ("pvc", "pvc+tdts"):
        alpha_c, alpha_n, q = mixing.pvc_bin(alpha, theta2, lobe)
        r_half = grid.half_levels(r.T).T
        Mt, gamma2, valid = pvc.propagate(dr, r_half, alpha_c, alpha_n, q, tu, instrument.rho_tr,
                                          instrument.rho_fov, instrument.kappa, instrument.receiver,
                                          native=not is_compat)
        if is_compat:
            Fm, Fp = compat.factors_compat(gamma2, theta2, alpha, beta, kinds, descriptors)
            beta_ms = compat.bin_means_compat(beta, tu, Mt, Fm, Fp, valid, tau)
            tail = np.zeros(Fm.shape)
        else:
            Fm, Fp, tail = anisotropy.half_level_factors(gamma2, theta2, descriptors, kinds, tables,
                                                         near_index)
            beta_ms = pvc.bin_means(beta, tu, Mt, Fm, Fp)
            # (C.4) reduces to (C.2) where M = 0 at both half-levels: use the single-scatter value
            no_ms = (np.real(Mt[:-1]) == 0) & (np.real(Mt[1:]) == 0)
            beta_ms = np.where(no_ms[:, None, :], beta_ss, beta_ms)
        flags |= pvc_flags(alpha, theta2, lobe, tail, descriptors, Mt)
    else:
        raise NotImplementedError(f"solver {config.solver!r} is not built yet")

    n_out = N if config.n_out is None else max(int(config.n_out), N)
    if config.solver == "pvc+tdts":
        beta_wide = wide_angle_part(r, alpha, opt["omega"], opt["g"], theta2, lobe, kinds, instrument, config,
                                    n_out)
    else:
        beta_wide = np.zeros((n_out, J, B), dtype=beta_ms.dtype)

    par, perp, closure_used, pol_flags, diag = polarization.polarize(
        beta_ms, beta_ss, opt["delta"], kinds, descriptors, config, dr=dr, r=r, alpha=alpha, lobe=lobe, tu=tu,
        Mt=Mt, gamma2=gamma2, theta2=theta2, F_minus=Fm, F_plus=Fp, tables=tables, near_index=near_index,
        rho_fov=instrument.rho_fov)
    flags |= pol_flags
    if config.solver == "pvc+tdts":
        dw = config.delta_wide
        par = _pad(par, n_out) + beta_wide / (1.0 + dw)
        perp = _pad(perp, n_out) + beta_wide * dw / (1.0 + dw)
        has_w = (np.real(beta_wide[:N]) > 0).any(axis=1)
        flags |= np.where(has_w, flag("WIDE_HEURISTIC"), 0).astype(np.uint16)
    with np.errstate(divide="ignore", invalid="ignore"):
        M = np.where(np.real(tu) > 0, Mt / np.where(np.real(tu) > 0, tu, 1.0), np.nan)

    return dict(beta_par=_pad(par, n_out), beta_perp=_pad(perp, n_out), beta_ss=beta_ss,
                beta_wide=beta_wide, ms_factor=M, gamma2=gamma2,
                a_s=accumulated_ss_fraction(beta_ms, beta_ss, dr, kinds), flags=flags.astype(np.uint16),
                closure_used=closure_used, diagnostic=diag)


def wide_angle_part(r, alpha, omega, g, theta2, lobe, kinds, instrument, config, n_out):
    """beta_w per class, (n_out, J, B), §C.2.5 with heuristic H1; H2 is applied by the caller.

    Native: regular grid required (Dep8), sources start at the first bin with particulate extinction,
    exact source width (Dep5), class-general sums (Dep11), optional coarsening `wide_coarsen`.
    compat: `wide_angle.c` from gate 0, branch inferred from the lobe test, variance.c width.
    """
    N, J, B = alpha.shape
    is_compat = config.compat is not None
    k = int(config.wide_coarsen)
    if is_compat:
        if k != 1:
            raise ValueError("wide_coarsen must be 1 in compat mode")
        dr = np.abs(r[2] - r[1])
        part = np.array([kk != "molecular" for kk in kinds])
        lobe_test = mixing.lobe_mask(theta2, part, config.theta_max)       # pvc compat treats all as lobe
        lobe_w = compat.tdts_branch_lobe(alpha, lobe_test, kinds)
        lidar = lobe_w.any(axis=(0, 1))
        sp, sm = compat.tdts_sources_compat(alpha, omega, g, kinds, lidar, dr, instrument.kappa)
        s2d = compat.source_width_compat(r, alpha, theta2, kinds, lidar, instrument.rho_tr)
        ap, om, gp, share = tdts.scaled_medium(alpha, omega, g, lobe_w)
        return tdts.wide_angle(r, dr, ap, om, gp, share, sp, sm, s2d, instrument.rho_tr, instrument.rho_fov,
                               n_out, instrument.receiver, d_factor=tdts.COMPAT_D_FACTOR, c_exp=True)
    if not np.all(grid.is_regular(r.T)):
        raise ValueError("pvc+tdts needs a regular range grid (spacing constant to 5 %); regrid first (Dep8)")
    part = np.array([kk != "molecular" for kk in kinds])
    _check_lobe_albedo(alpha, omega, g, lobe)
    q_j = np.where(lobe, alpha * np.where(lobe, theta2, 0.0), 0.0)
    if k > 1:
        r_w, a_w, om_w, g_w, q_w = tdts.coarsen_classes(k, r, alpha, omega, g, q_j)
        th2_w = np.where(np.real(a_w) > 0, q_w / np.where(np.real(a_w) > 0, a_w, 1.0), np.inf)
        lobe_w = lobe.any(axis=0)[None] & (np.real(a_w) > 0) & (np.real(th2_w) <= config.theta_max ** 2)
        lobe_w = lobe_w | (lobe.any(axis=0)[None] & (np.real(a_w) == 0))
        q_w = np.where(lobe_w, q_w, 0.0)
        n_out_w = -(-n_out // k)
    else:
        r_w, a_w, om_w, g_w, q_w, lobe_w, n_out_w = r, alpha, omega, g, q_j, lobe, n_out
    dr = (r_w[-1] - r_w[0]) / (r_w.shape[0] - 1)
    cloudy = (np.real(a_w[:, part]) > 0).any(axis=1)
    first = np.where(cloudy.any(axis=0), np.argmax(cloudy, axis=0), r_w.shape[0])
    sp, sm = tdts.source_terms(a_w, om_w, g_w, lobe_w, dr, instrument.kappa, first_bin=first)
    s2d = tdts.source_width(r_w, dr, q_w.sum(axis=1), instrument.rho_tr)
    ap, om, gp, share = tdts.scaled_medium(a_w, om_w, g_w, lobe_w)
    bw = tdts.wide_angle(r_w, dr, ap, om, gp, share, sp, sm, s2d, instrument.rho_tr, instrument.rho_fov,
                         n_out_w, instrument.receiver)
    return tdts.expand(bw, k, n_out) if k > 1 else bw


def _check_lobe_albedo(alpha, omega, g, lobe):
    bad = lobe & (np.real(alpha) > 0) & (np.real(omega * g) < 0.5)
    if np.any(bad):
        raise ValueError("a lobe-forming class must satisfy omega * g >= 1/2 where alpha > 0 (§C.2.2)")


def near_tables(near, kinds, names):
    """Normalize `near` ({class name or index: NearBackscatter}) to a list of length J."""
    J = len(kinds)
    names = tuple(names) if names else tuple(kinds)
    tables = [None] * J
    for key, tab in (near or {}).items():
        j = key if isinstance(key, (int, np.integer)) else names.index(key)
        if not 0 <= j < J:
            raise ValueError(f"near-backscatter key {key!r} does not name a class")
        tables[j] = tab
    return tables


def pvc_flags(alpha, theta2, lobe, tail, descriptors, Mt):
    """Flag bits 4 (LOBE_MIXTURE) and 5 (TABLE_EXTRAPOLATED), (N, B)."""
    present = lobe & (np.real(alpha) > 0)
    th = np.where(present, np.sqrt(np.real(np.where(present, theta2, 1.0))), np.nan)
    with np.errstate(invalid="ignore"):
        ratio = np.nanmax(np.where(present, th, -np.inf), axis=1) / np.nanmin(np.where(present, th, np.inf), axis=1)
    mix = (present.sum(axis=1) >= 2) & (ratio > LOBE_MIXTURE_RATIO)
    is_tab = np.array([d == "table" for d in descriptors])[None, :, None]
    ms = (np.real(Mt[:-1]) > 0) | (np.real(Mt[1:]) > 0)
    ext = (is_tab & (np.real(alpha) > 0) & (tail > TABLE_TAIL_WEIGHT)).any(axis=1) & ms
    return (np.where(mix, flag("LOBE_MIXTURE"), 0) | np.where(ext, flag("TABLE_EXTRAPOLATED"), 0)).astype(np.uint16)


def _pad(a, n_out):
    if a.shape[0] == n_out:
        return a
    return np.concatenate([a, np.zeros((n_out - a.shape[0],) + a.shape[1:], dtype=a.dtype)], axis=0)


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
    if config.solver == "mc":
        return simulate_mc(r, instrument, optics, config, phase, P)
    opt = {nm: to_bin_major(getattr(optics, nm)) for nm in _FIELDS}
    out = solve(r.T.copy(), opt, optics.kinds, optics.descriptors, instrument, config, near=near,
                near_index=to_bin_major(optics.near_index), names=optics.names)
    return assemble(out, P, r, optics, instrument, config)


def simulate_mc(r, instrument, optics, config, phase, P) -> Returns:
    """Solver `mc` (§C.2.7): Monte Carlo tallies assembled as `Returns`.

    beta_par/beta_perp are co-/cross-polarized tallies (scalar mode: all of I in beta_par, beta_perp = 0,
    closure_used "mc_scalar"); beta_ss is the order-1 tally; beta_order holds orders 1, 2, >= 3;
    stderr_* are batch-mean standard errors; ms_factor and gamma2 are NaN (no PVC moments);
    beta_wide is zero (the Monte Carlo does not separate wide-angle light). meta carries the largest
    single-contribution fraction of any bin with one-way tau <= 3 (spike diagnostic).
    """
    from .mc.driver import simulate_mc as _run
    out = _run(r, optics, instrument, config, phase)
    N, J = optics.shape[1:]
    dr = grid.bin_widths(r).T
    alpha = to_bin_major(optics.alpha)
    _, tu, _ = single(alpha.sum(axis=1), to_bin_major(optics.beta), dr)
    flags = np.where(tu[:-1] < EXTINGUISHED_TU, flag("EXTINGUISHED"), 0).astype(np.uint16)
    total = out["beta_par"] + out["beta_perp"]
    tau_mid = np.cumsum(alpha.sum(axis=1) * dr, axis=0) - 0.5 * alpha.sum(axis=1) * dr
    spike = np.where(tau_mid <= 3.0, out["spike"], 0.0)
    nan = np.full((N + 1, P), np.nan)
    res = dict(beta_par=out["beta_par"], beta_perp=out["beta_perp"], beta_ss=out["beta_ss"],
               beta_wide=np.zeros_like(out["beta_par"]), ms_factor=nan, gamma2=nan.copy(),
               a_s=accumulated_ss_fraction(total, out["beta_ss"], dr, optics.kinds), flags=flags,
               closure_used=np.array(["mc" if config.mc.polarized else "mc_scalar"] * P, dtype=object))
    ret = assemble(res, P, r, optics, instrument, config)
    meta = dict(ret.meta, mc_spike_max=float(spike.max()) if spike.size else 0.0)
    return replace(ret, meta=meta, beta_order=from_bin_major(out["beta_order"], P),
                   stderr_par=from_bin_major(out["stderr_par"], P), stderr_perp=from_bin_major(out["stderr_perp"], P))


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
                   r=r_out, meta=meta, diagnostic=_diag(out.get("diagnostic"), P, f))


def _diag(d, P, f):
    if d is None:
        return None
    return {"closure": np.asarray(d["closure"]), "beta_perp": f(from_bin_major(d["beta_perp"], P))}


def _config_dict(config):
    d = asdict(config)
    d["mc"] = dict(d["mc"])
    return d
