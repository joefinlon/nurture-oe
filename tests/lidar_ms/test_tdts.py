"""Solver `pvc+tdts`: analytical A7-A9, regression R6-R7, derivatives J5, flag bit 7 (appendix C §C.4.3)."""
import numpy as np
import pytest

from lidar_ms import (Instrument, OpticalProfile, OpticalTangent, RangeGrid, SolverConfig, check_jvp,
                      simulate)
from lidar_ms import tdts
from lidar_ms.compat import read_ms_input
from lidar_ms.constants import FLAG_BITS
from msref import REF_DIR, REF_PAIRS, load_ref

COMPAT = SolverConfig(solver="pvc+tdts", compat="multiscatter-1.2.11")
NATIVE = SolverConfig(solver="pvc+tdts")


def case(name):
    return read_ms_input(REF_DIR / f"{name}.in")


def _padded(D):
    pad = lambda a: np.concatenate([np.zeros_like(a[:1]), a, np.zeros_like(a[:1])])
    return np.stack([pad(D[k]) for k in range(6)])


def _uniform(N, ap, om, gp, dr):
    D, lt, _ = tdts.transfer_coefficients(np.full((N, 1), ap), np.full((N, 1), om), np.full((N, 1), gp),
                                          np.array([dr]))
    return _padded(D), lt


# --- A7 ------------------------------------------------------------------------------------------------
def test_a7_coefficients_sum_and_sign():
    ap = np.logspace(-4, 0, 40)[:, None]
    for om in (0.5, 0.8, 0.95, 0.99, 1.0):
        for gp in (0.0, 0.3, 0.6, 0.85):
            for dr in (1.0, 15.0, 200.0):
                D, _, _ = tdts.transfer_coefficients(ap, np.full_like(ap, om), np.full_like(ap, gp), np.array([dr]))
                fa = np.exp(-dr * ap * (1 - om))
                np.testing.assert_allclose(D.sum(axis=0), fa, rtol=1e-12)
                assert np.all(D >= 0)


def test_a7_conservative_energy_budget():
    N, dr = 60, 15.0
    rng = np.random.default_rng(0)
    ap = rng.uniform(1e-3, 0.1, (N, 1))
    D, _, _ = tdts.transfer_coefficients(ap, np.ones((N, 1)), rng.uniform(0, 0.8, (N, 1)), np.array([dr]))
    Dp = _padded(D)
    Im = np.zeros((N + 2, 1))
    Ip = np.zeros((N + 2, 1))
    injected = left = 0.0
    for it in range(4 * N):
        lo_m = (Dp[2, 1] + Dp[3, 1]) * Im[1] + (Dp[5, 1] + Dp[4, 1]) * Ip[1]           # bin 0 -> above the grid
        lo_p = (Dp[4, N] + Dp[5, N]) * Im[N] + (Dp[3, N] + Dp[2, N]) * Ip[N]           # bin N-1 -> below
        im, ip = tdts.step_streams(Im, Ip, Dp, 0, N - 1)
        Im[1:N + 1], Ip[1:N + 1] = im, ip
        left += float(lo_m[0] + lo_p[0])
        if it < N:
            s = rng.uniform(0, 1, 2)
            Im[it + 1] += s[0]
            Ip[it + 1] += s[1]
            injected += s.sum()
        total = float(Im.sum() + Ip.sum())
        assert total + left == pytest.approx(injected, rel=1e-12)


# --- A8 ------------------------------------------------------------------------------------------------
def test_a8_mirror_symmetry():
    N = 201
    c = N // 2
    Dp, _ = _uniform(N, 0.05, 0.97, 0.0, 15.0)
    Im = np.zeros((N + 2, 1))
    Ip = np.zeros((N + 2, 1))
    Im[c + 1] = Ip[c + 1] = 1.0
    for _ in range(c - 1):
        im, ip = tdts.step_streams(Im, Ip, Dp, 0, N - 1)
        Im[1:N + 1], Ip[1:N + 1] = im, ip
        a, b = Ip[1:N + 1, 0], Im[1:N + 1, 0][::-1]
        np.testing.assert_allclose(a, b, rtol=1e-12, atol=1e-12 * a.max())


# --- A9 ------------------------------------------------------------------------------------------------
@pytest.mark.parametrize("dr,lt", [(1.0, 20.0), (1.0, 2.0), (5.0, 1.0)])
def test_a9_lateral_spreading(dr, lt):
    """Thick uniform medium, uniformly filled (so V/I stays uniform): the coded fast inversion follows
    (4/3) l_t^2 (n + exp(-n) - 1) from n = 0.1 to 30.

    A single pulse in the discrete scheme falls up to 18 % below the law at n < 1, because the
    stencil mixes bins whose V/I differ (recorded in docs/findings_for_pk.md)."""
    steps = int(30 * lt / dr)
    N = 2 * steps + 41
    c = N // 2
    Dp, ltv = _uniform(N, 1 / lt, 1.0, 0.0, dr)
    Im = np.zeros((N + 2, 1))
    Im[1:-1] = 0.5
    Ip = Im.copy()
    Vm = np.zeros_like(Im)
    Vp = np.zeros_like(Im)
    s2d = np.zeros((N, 1))
    s = slice(1, N + 1)
    for it in range(1, steps + 1):
        im, ip = tdts.step_streams(Im, Ip, Dp, 0, N - 1)
        vm, vp = tdts.step_streams(Vm, Vp, Dp, 0, N - 1)
        vm = vm + Im[s] * tdts.expansion(Im[s], Vm[s], s2d, ltv, dr)
        vp = vp + Ip[s] * tdts.expansion(Ip[s], Vp[s], s2d, ltv, dr)
        Im[s], Ip[s], Vm[s], Vp[s] = im, ip, vm, vp
        n = it * dr / lt
        if n >= 0.1:
            law = 4 / 3 * lt ** 2 * (n + np.exp(-n) - 1)
            assert Vm[c + 1, 0] / Im[c + 1, 0] == pytest.approx(law, rel=0.03)


# --- R6, R7 --------------------------------------------------------------------------------------------
def _compare(got, ref, tol):
    ok = np.isfinite(ref) & (np.abs(ref) > 1e-12 * np.nanmax(np.abs(ref)))
    assert np.all(np.isfinite(got))
    np.testing.assert_allclose(got[ok], ref[ok], rtol=tol, atol=0)


R6 = [(nm, t) for nm, t in REF_PAIRS if t in ("tdts_only", "tdts_only_nolobe")]
R7 = [(nm, t) for nm, t in REF_PAIRS if t in ("pvc_tdts", "single_tdts_nolobe")]


def test_r6_r7_counts():
    assert len([t for _, t in R6 if t == "tdts_only"]) == 5
    assert len([t for _, t in R6 if t == "tdts_only_nolobe"]) == 3


@pytest.mark.parametrize("name,tag", R6)
def test_r6_tdts_alone(name, tag):
    c = case(name)
    out = simulate(c.grid, c.instrument, c.optics, COMPAT)
    _compare(out.beta_wide[0].sum(-1), load_ref(name, tag).data[:, 4], 1e-6)


@pytest.mark.parametrize("name,tag", R7)
def test_r7_pvc_plus_tdts(name, tag):
    """`pvc_tdts`: the whole return. `single_tdts_nolobe`: single scattering plus the no-lobe
    wide-angle part (in Returns: beta_ss + beta_wide)."""
    c = case(name)
    out = simulate(c.grid, c.instrument, c.optics, COMPAT)
    got = out.total[0] if tag == "pvc_tdts" else out.beta_ss[0].sum(-1) + out.beta_wide[0].sum(-1)
    _compare(got, load_ref(name, tag).data[:, 4], 1e-6)


# --- J5 ------------------------------------------------------------------------------------------------
J5 = [("hogan2008fig8", "alpha"), ("hogan2008fig8", "a_eq"), ("hogan2008fig8", "beta"),
      ("halo_ice_over_liquid", "alpha"), ("halo_ice_over_liquid", "a_eq"), ("halo_ice_over_liquid", "beta"),
      ("i3rc_lidar_exp5", "alpha"), ("i3rc_lidar_exp5", "beta")]
# Not included: i3rc_lidar_exp5 with a_eq tangents. They move the source width across the branch point
# nu -> 0+ of the inversion (e), where n ~ sqrt(1.5 nu) has unbounded slope (docs/findings_for_pk.md).


@pytest.mark.parametrize("name,field", J5)
def test_j5_tdts_vs_central_difference(name, field):
    c = case(name)
    P, N, J = c.optics.shape
    rng = np.random.default_rng(1)
    t = OpticalTangent(**{field: np.real(getattr(c.optics, field))[:, None] * rng.standard_normal((P, 3, N, J))})
    err = check_jvp(c.grid, c.instrument, c.optics, t, NATIVE)
    assert err <= 1e-5, err


# --- flag bit 7, coarsening, input checks ------------------------------------------------------------
def test_k6_wide_heuristic_flag_alone():
    c = case("halo_liquid_30")
    out = simulate(c.grid, c.instrument, c.optics, NATIVE)
    assert int(np.bitwise_or.reduce(out.flags[0])) == 1 << FLAG_BITS["WIDE_HEURISTIC"]
    assert not simulate(c.grid, c.instrument, c.optics, SolverConfig(solver="pvc")).flags.any()


def test_h2_split_and_class_shares():
    c = case("hogan2008fig8")
    out = simulate(c.grid, c.instrument, c.optics, SolverConfig(solver="pvc+tdts", delta_wide=0.5))
    pvc = simulate(c.grid, c.instrument, c.optics, SolverConfig(solver="pvc"))
    w = out.beta_wide[0]
    np.testing.assert_allclose(out.beta_perp[0] - pvc.beta_perp[0], w * 0.5 / 1.5, rtol=1e-12, atol=1e-30)
    np.testing.assert_allclose(out.total[0], pvc.total[0] + w.sum(-1), rtol=1e-12)
    assert w[:, 3].sum() > 0 and w[:, 2].sum() > 0          # both air and the particulate class get a share


@pytest.mark.parametrize("k", [2, 4])
def test_wide_coarsen_runs(k):
    c = case("halo_ice_over_liquid")
    full = simulate(c.grid, c.instrument, c.optics, NATIVE).beta_wide[0].sum(-1)
    coarse = simulate(c.grid, c.instrument, c.optics, SolverConfig(solver="pvc+tdts", wide_coarsen=k))
    w = coarse.beta_wide[0].sum(-1)
    assert w.shape == full.shape and np.all(np.isfinite(w)) and w.sum() > 0


def test_n_out_extends_returns():
    c = case("hogan2008fig8")
    out = simulate(c.grid, c.instrument, c.optics, SolverConfig(solver="pvc+tdts", n_out=60))
    assert out.beta_par.shape[1] == 60 and out.r.shape[1] == 60 and out.beta_wide[0, 41:].sum() > 0


def test_native_tdts_input_checks():
    c = case("halo_liquid_30")
    r = c.grid.r.copy()
    r[0, 60:] += 5.0                                       # 15 m -> 20 m step somewhere: irregular
    with pytest.raises(ValueError):
        simulate(RangeGrid(r), c.instrument, c.optics, NATIVE)
    g = np.array(c.optics.g)
    g[..., 0] = 0.3                                        # omega g < 1/2 for a lobe-forming class
    with pytest.raises(ValueError):
        simulate(c.grid, c.instrument, c.optics.replace(g=g), NATIVE)
    with pytest.raises(ValueError):
        simulate(c.grid, c.instrument, c.optics, SolverConfig(solver="pvc+tdts", compat="multiscatter-1.2.11",
                                                              wide_coarsen=2))
