"""Derivative tests J1-J3 (appendix C §C.4.3). J4 (closures) is in test_polarization.py."""
import numpy as np
import pytest

from lidar_ms import OpticalTangent, SolverConfig, check_jvp, jvp, simulate
from lidar_ms.compat import read_ms_input
from msref import REF_DIR

PVC = SolverConfig(solver="pvc")


def case(name):
    return read_ms_input(REF_DIR / f"{name}.in")


def random_tangents(optics, K, fields=("alpha", "a_eq", "beta", "delta"), seed=0):
    """Relative random perturbations of every class, only where the field is non-zero (delta: everywhere)."""
    rng = np.random.default_rng(seed)
    P, N, J = optics.shape
    out = {}
    for nm in fields:
        x = np.real(getattr(optics, nm))[:, None]
        if nm == "delta":
            t = 0.1 * rng.standard_normal((P, K, N, J))
        else:
            t = x * rng.standard_normal((P, K, N, J))
        out[nm] = t
    return OpticalTangent(**out)


@pytest.mark.parametrize("cfg", [SolverConfig(solver="single"), PVC,
                                 SolverConfig(solver="pvc", compat="multiscatter-1.2.11")])
def test_j1_real_part_equals_simulate(cfg):
    c = case("halo_ice_over_liquid")
    ret, tan = jvp(c.grid, c.instrument, c.optics, random_tangents(c.optics, 3), cfg)
    ref = simulate(c.grid, c.instrument, c.optics, cfg)
    # relative to each output's maximum: the real part carries the complex step's h^2 term (~1e-68
    # here), which no element-wise relative test passes where the true value is 0
    scale = np.abs(ref.total_by_class).max()
    np.testing.assert_allclose(ret.beta_par, ref.beta_par, rtol=1e-14, atol=1e-14 * scale)
    np.testing.assert_allclose(ret.beta_perp, ref.beta_perp, rtol=1e-14, atol=1e-14 * scale)
    np.testing.assert_allclose(ret.ms_factor, ref.ms_factor, rtol=1e-14, atol=0)
    assert tan.d_beta_par.shape == (1, 3) + ref.beta_par.shape[1:]


@pytest.mark.parametrize("name", ["halo_liquid_30", "halo_ice_over_liquid"])
@pytest.mark.parametrize("field", ["alpha", "a_eq", "beta", "delta"])
def test_j2_pvc_vs_central_difference(name, field):
    c = case(name)
    t = random_tangents(c.optics, 4, fields=(field,), seed=hash(field) % 1000)
    err = check_jvp(c.grid, c.instrument, c.optics, t, PVC)
    assert err <= 1e-6, err


def test_j2_all_fields_at_once():
    c = case("halo_ice_over_liquid")
    err = check_jvp(c.grid, c.instrument, c.optics, random_tangents(c.optics, 3, seed=7), PVC)
    assert err <= 1e-6, err


def test_j3_linearity():
    c = case("halo_ice_over_liquid")
    t = random_tangents(c.optics, 2, seed=3)
    both = OpticalTangent(**{nm: np.concatenate([getattr(t, nm), getattr(t, nm).sum(1, keepdims=True)], axis=1)
                             for nm in ("alpha", "omega", "g", "a_eq", "beta", "delta")})
    _, tan = jvp(c.grid, c.instrument, c.optics, both, PVC)
    for d in (tan.d_beta_par, tan.d_beta_perp):
        np.testing.assert_allclose(d[:, 2], d[:, 0] + d[:, 1], rtol=1e-10, atol=1e-10 * np.abs(d).max())


def test_jvp_rejects_mc_and_bad_shapes():
    c = case("halo_liquid_30")
    with pytest.raises(ValueError):
        jvp(c.grid, c.instrument, c.optics, random_tangents(c.optics, 2), SolverConfig(solver="mc"))
    with pytest.raises(ValueError):
        jvp(c.grid, c.instrument, c.optics, OpticalTangent(alpha=np.zeros((1, 2, 5, 4))), PVC)
