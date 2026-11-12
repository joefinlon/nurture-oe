"""Structure test ST4: Returns -> netCDF -> Returns (appendix C §C.4.3)."""
import numpy as np
import pytest

from lidar_ms import SolverConfig, simulate
from lidar_ms.compat import read_ms_input
from lidar_ms.constants import FLAG_MASKS
from lidar_ms.io import (dataset_to_inputs, inputs_to_dataset, read_returns, returns_to_dataset,
                         write_returns)
from msref import REF_DIR


def _returns(cfg):
    c = read_ms_input(REF_DIR / "halo_ice_over_liquid.in")
    return c, simulate(c.grid, c.instrument, c.optics, cfg)


def _equal(a, b):
    for nm in ("beta_par", "beta_perp", "beta_wide", "beta_ss", "ms_factor", "gamma2", "a_s", "flags", "r"):
        np.testing.assert_array_equal(getattr(a, nm), getattr(b, nm), err_msg=nm)
        assert getattr(a, nm).dtype == getattr(b, nm).dtype, nm
    assert list(a.closure_used) == list(b.closure_used)
    assert a.kinds == b.kinds and a.names == b.names
    assert a.meta == b.meta


@pytest.mark.parametrize("cfg", [SolverConfig(solver="pvc"), SolverConfig(solver="pvc+tdts", n_out=320,
                                                                          diagnose_closures=True)])
def test_st4_round_trip(tmp_path, cfg):
    _, ret = _returns(cfg)
    path = tmp_path / "returns.nc"
    write_returns(ret, path)
    back = read_returns(path)
    _equal(ret, back)
    if ret.diagnostic is not None:
        np.testing.assert_array_equal(ret.diagnostic["beta_perp"], back.diagnostic["beta_perp"])
        assert list(ret.diagnostic["closure"]) == list(back.diagnostic["closure"])


def test_st4_units_and_cf_flags(tmp_path):
    _, ret = _returns(SolverConfig(solver="pvc+tdts"))
    path = tmp_path / "returns.nc"
    write_returns(ret, path)
    import xarray as xr
    with xr.open_dataset(path) as ds:
        for nm in ds.data_vars:
            assert "long_name" in ds[nm].attrs, nm
            if nm not in ("closure_used", "flags", "diagnostic_closure"):
                assert "units" in ds[nm].attrs, nm
        assert ds["r"].attrs["units"] == "m"
        np.testing.assert_array_equal(ds["flags"].attrs["flag_masks"], FLAG_MASKS)
        assert ds["flags"].attrs["flag_meanings"].split()[6] == "extinguished"
        assert ds.attrs["solver"] == "pvc+tdts" and "config_json" in ds.attrs


def test_inputs_round_trip(tmp_path):
    c, _ = _returns(SolverConfig())
    ds = inputs_to_dataset(c.grid, c.instrument, c.optics)
    ds.to_netcdf(tmp_path / "in.nc")
    import xarray as xr
    with xr.open_dataset(tmp_path / "in.nc") as d2:
        g, i, o = dataset_to_inputs(d2.load())
    np.testing.assert_array_equal(g.r, c.grid.r)
    assert i == c.instrument
    for nm in ("alpha", "omega", "g", "a_eq", "beta", "delta", "near_index"):
        np.testing.assert_array_equal(getattr(o, nm), getattr(c.optics, nm))
    assert (o.kinds, o.descriptors, o.names) == (c.optics.kinds, c.optics.descriptors, c.optics.names)
