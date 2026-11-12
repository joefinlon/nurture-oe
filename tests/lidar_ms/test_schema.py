"""Structure tests ST1 (schemas) and ST2 (grid) (appendix C §C.4.3)."""
import numpy as np
import pytest

from lidar_ms import grid
from lidar_ms.schema import (Instrument, McConfig, NearBackscatter, OpticalProfile, OpticalTangent,
                             PhaseTables, RangeGrid, Returns, SolverConfig)

P, N, J = 2, 6, 3
SHAPE = (P, N, J)


def _optics(**kw):
    d = dict(kinds=("liquid", "ice", "molecular"), descriptors=("droplet", "yang", "isotropic"),
             alpha=np.full(SHAPE, 1e-3), omega=np.ones(SHAPE), g=np.full(SHAPE, 0.85),
             a_eq=np.full(SHAPE, 8e-6), beta=np.full(SHAPE, 5e-5), delta=np.zeros(SHAPE))
    d.update(kw)
    return OpticalProfile(**d)


def _with(name, idx, value):
    a = np.array(getattr(_optics(), name))
    a[idx] = value
    return {name: a}


# --- ST1 ---------------------------------------------------------------------------------------
def test_st1_halo_instrument():
    halo = Instrument.from_full_angles(532e-9, 0.8e-3, 1.0e-3)
    assert halo.rho_tr == pytest.approx(0.28284e-3, rel=1e-5)
    assert halo.rho_fov == pytest.approx(0.50e-3, rel=1e-5)
    assert halo.kappa == pytest.approx(1.04596, rel=1e-5)


@pytest.mark.parametrize("kw", [dict(wavelength=0), dict(rho_tr=-1e-3), dict(rho_fov=0.0),
                                dict(rho_fov=np.nan), dict(receiver="cone")])
def test_st1_instrument_raises(kw):
    d = dict(wavelength=532e-9, rho_tr=3e-4, rho_fov=5e-4)
    d.update(kw)
    with pytest.raises(ValueError):
        Instrument(**d)


@pytest.mark.parametrize("r", [np.arange(2.0) + 1, -np.arange(5.0), np.array([1.0, 3.0, 2.0, 4.0]),
                               np.array([1.0, 2.0, 2.0, 3.0]), np.ones((2, 2, 3)), np.array([1, np.nan, 3.0])])
def test_st1_rangegrid_raises(r):
    with pytest.raises(ValueError):
        RangeGrid(r)


def test_st1_rangegrid_z_shape():
    RangeGrid(np.arange(1.0, 5.0), z=np.arange(4.0))
    with pytest.raises(ValueError):
        RangeGrid(np.arange(1.0, 5.0), z=np.arange(3.0))


@pytest.mark.parametrize("kw", [
    dict(kinds=("liquid", "ice", "ice")),                       # no molecular
    dict(kinds=("molecular", "ice", "molecular")),              # two molecular
    dict(kinds=("liquid", "snow", "molecular")),                # unknown kind
    dict(descriptors=("droplet", "yang")),                      # wrong count
    dict(descriptors=("droplet", "mie", "isotropic")),          # unknown descriptor
    dict(names=("a", "a", "b")),                                # duplicate names
    dict(alpha=np.zeros((P, N, 2))),                            # wrong J
    dict(omega=np.ones((P, N + 1, J))),                         # shape mismatch
    dict(near_index=np.full(SHAPE, -1)),                        # negative index
    dict(near_index=np.zeros(SHAPE)),                           # not integer
])
def test_st1_optics_structure_raises(kw):
    with pytest.raises(ValueError):
        _optics(**kw)


@pytest.mark.parametrize("name,value", [("alpha", -1e-6), ("omega", 1.1), ("omega", -0.1), ("g", 1.5),
                                        ("g", -1.01), ("beta", -1e-9), ("delta", 1.2), ("delta", -0.1),
                                        ("a_eq", 0.0), ("alpha", np.inf)])
def test_st1_optics_value_raises(name, value):
    with pytest.raises(ValueError):
        _optics(**_with(name, (0, 2, 0), value))


def test_st1_a_eq_ignored_for_molecular():
    _optics(**_with("a_eq", (0, 2, 2), 0.0))


def test_st1_complex_accepted():
    alpha = np.full(SHAPE, 1e-3) + 1j * 1e-30
    opt = _optics(alpha=alpha)
    assert opt.alpha.dtype == np.complex128 and opt.is_complex
    assert opt.beta.dtype == np.float64
    with pytest.raises(ValueError):
        _optics(alpha=alpha * 0 - 1e-3 + 1j)        # constraint tested on real part


def test_st1_near_backscatter():
    x = np.linspace(0, 1, 5)
    NearBackscatter("gamma", x, np.full(5, 0.6), np.full(5, 0.4))
    for kw in (dict(coord="theta"), dict(x=x + 0.1), dict(x=x[::-1]), dict(c_par=np.full(5, 0.7)),
               dict(c_perp=np.full((2, 5), 0.4)), dict(c_par=np.full(4, 0.6))):
        d = dict(coord="gamma", x=x, c_par=np.full(5, 0.6), c_perp=np.full(5, 0.4))
        d.update(kw)
        with pytest.raises(ValueError):
            NearBackscatter(**d)


def test_st1_near_backscatter_delta_warning():
    nb = NearBackscatter("gamma", np.linspace(0, 1, 5), np.full(5, 0.9), np.full(5, 0.1))
    with pytest.warns(UserWarning):
        nb.check_delta(np.array([0.5]), np.array([0]))


def test_st1_phase_tables():
    th = np.linspace(0, np.pi, 181)
    iso = np.ones((1, th.size))
    ok = dict(theta=th, P11=iso, P12=0 * iso, P22=iso, P33=iso, P34=0 * iso, P44=iso,
              type_index=np.zeros((N, J), int))
    PhaseTables(**ok)
    for kw in (dict(theta=th[:-1], P11=iso[:, :-1]), dict(P11=2 * iso), dict(P22=np.ones((2, th.size))),
               dict(type_index=np.ones((N, J), int)), dict(theta=th[::-1])):
        with pytest.raises(ValueError):
            PhaseTables(**{**ok, **kw})


@pytest.mark.parametrize("kw", [dict(solver="fast"), dict(closure="hu2"), dict(theta_max=0), dict(n_out=0),
                                dict(wide_coarsen=0), dict(wide_coarsen=1.5), dict(delta_wide=-1),
                                dict(compat="1.2.10"), dict(mc=None)])
def test_st1_solver_config_raises(kw):
    with pytest.raises(ValueError):
        SolverConfig(**kw)


@pytest.mark.parametrize("kw", [dict(n_photons=0), dict(n_batches=1), dict(epsilon_ddis=1.0),
                                dict(w_roulette=0), dict(workers=0)])
def test_st1_mc_config_raises(kw):
    with pytest.raises(ValueError):
        McConfig(**kw)


def test_st1_tangent():
    t = OpticalTangent(alpha=np.ones((P, 4, N, J)))
    assert t.beta.shape == (P, 4, N, J) and not t.beta.any()
    for kw in (dict(), dict(alpha=np.ones((P, N, J))), dict(alpha=np.ones((P, 4, N, J)), beta=np.ones((P, 3, N, J))),
               dict(alpha=np.full((P, 4, N, J), np.nan))):
        with pytest.raises(ValueError):
            OpticalTangent(**kw)


def test_st1_returns_shapes():
    z = np.zeros
    kw = dict(beta_par=z((P, N, J)), beta_perp=z((P, N, J)), beta_ss=z((P, N, J)), beta_wide=z((P, N, J)),
              ms_factor=z((P, N + 1)), gamma2=z((P, N + 1)), a_s=z((P, N)), closure_used=np.array(["none"] * P),
              flags=z((P, N), np.uint16), kinds=("liquid", "ice", "molecular"), names=("l", "i", "m"),
              r=z((P, N)))
    Returns(**kw)
    with pytest.raises(ValueError):
        Returns(**{**kw, "beta_perp": z((P, N + 1, J))})
    with pytest.raises(ValueError):
        Returns(**{**kw, "flags": z((P, N - 1))})


# --- ST2 ---------------------------------------------------------------------------------------
def test_st2_irregular_grid():
    r = np.array([100.0, 110.0, 130.0, 160.0, 200.0])
    dr = grid.bin_widths(r)
    np.testing.assert_allclose(dr, [10.0, 15.0, 25.0, 35.0, 40.0], rtol=1e-14)
    rh = grid.half_levels(r)
    np.testing.assert_allclose(rh, [95.0, 105.0, 117.5, 142.5, 177.5, 220.0], rtol=1e-14)
    g = RangeGrid(np.stack([r, r + 1000]))
    np.testing.assert_allclose(g.dr[1], dr, rtol=1e-14)
    np.testing.assert_allclose(g.r_half[1], rh + 1000, rtol=1e-14)


def test_st2_is_regular():
    reg = 1000 + 15.0 * np.arange(10)
    assert grid.is_regular(reg)
    jit = reg.copy()
    jit[5] += 0.3                      # spacings 15.3, 14.7: within 5 %
    assert grid.is_regular(jit)
    bad = reg.copy()
    bad[5] += 1.0                      # spacings 16, 14: 14 % apart
    assert not grid.is_regular(bad)
    assert not grid.is_regular(np.array([1.0, 2.0, 1.5, 2.5]))
    np.testing.assert_array_equal(RangeGrid(np.stack([reg, bad])).is_regular(), [True, False])


def test_coarsen_conserves():
    a = np.arange(10.0)[:, None] * np.ones((1, 3))
    c = grid.coarsen(a, 4)
    assert c.shape == (3, 3)
    np.testing.assert_allclose(c.sum(0) * 4, a.sum(0))
    np.testing.assert_allclose(grid.coarsen_centres(np.arange(10.0), 4), [1.5, 5.5, 9.5])
