"""Reader for the HALO Aerosol + Water Vapor L2 product (HSRL/DIAL curtains).

See ``docs/data_map.md`` §1.1. The reader returns native (time, z) grids, native units
and NaN fills, with every attribute attached verbatim (including names with trailing
spaces). It does no screening, regridding, averaging or unit conversion.
"""

from pathlib import Path

import h5py
import numpy as np
import xarray as xr

from ._hdf5 import apply_cf_fill, check_units, decode_attrs, dims_of, parse_time_units, time_slice, to_datetime64

GROUPS = ("DataProducts", "DataUncertainty", "Nav_Data", "State", "UserInput")
_ROOT_COORDS = ("time", "z", "lat", "lon")


def read_halo_aerowv(path, groups=GROUPS, variables=None, time_range=None):
    """Read a HALO Aerosol/WV L2 granule into an ``xarray.Dataset``.

    Parameters
    ----------
    path : str or Path
        HDF5 granule (``nurture-HALO_G3_<date>_R0_L<n>.h5``).
    groups : sequence of str
        Groups to read. Variable names are kept verbatim (e.g. ``"532_ext"``); the
        source path is stored in each variable's ``hdf5_path`` attribute.
    variables : sequence of str, optional
        Subset of variable names (bare names, e.g. ``["532_total_attn_bsc", "532_dep"]``).
    time_range : (start, end), optional
        Inclusive UTC bounds (anything ``np.datetime64`` accepts). Only that slice is read.

    Returns
    -------
    xarray.Dataset
        Dims ``time`` (datetime64; raw seconds in coordinate ``time_native``) and ``z``
        (m, ascending on 2026-02-17). Datasets without dimension scales are mapped to
        dims by length; ambiguous lengths raise.

    Warns
    -----
    UserWarning
        When a legacy ``"Units "`` attribute disagrees with CF ``units``.
    """
    path = Path(path)
    with h5py.File(path, "r") as f:
        t_attrs = decode_attrs(f["time"])
        epoch, unit = parse_time_units(t_attrs["units"])
        t_native = f["time"][:]
        t = to_datetime64(t_native, epoch, unit)
        sl = time_slice(t, time_range)
        sizes = {"time": f["time"].shape[0], "z": f["z"].shape[0]}

        coords = {
            "time": ("time", t[sl], {k: v for k, v in t_attrs.items() if k != "units"}),
            "time_native": ("time", t_native[sl], t_attrs),
            "z": ("z", f["z"][:], decode_attrs(f["z"])),
            "lat": ("time", f["lat"][sl], decode_attrs(f["lat"])),
            "lon": ("time", f["lon"][sl], decode_attrs(f["lon"])),
        }
        data_vars = {}
        for g in groups:
            for name, ds in f[g].items():
                if not isinstance(ds, h5py.Dataset):
                    continue
                if variables is not None and name not in variables:
                    continue
                if name in data_vars or name in _ROOT_COORDS:
                    raise ValueError(f"variable name collision: {name!r} ({ds.name})")
                attrs = decode_attrs(ds)
                attrs["hdf5_path"] = ds.name
                check_units(attrs, ds.name)
                if ds.shape == (1,) and name == "State_Type":
                    data_vars[name] = ((), ds[0], attrs)
                    continue
                dims = dims_of(ds, sizes)
                idx = tuple(sl if d == "time" else slice(None) for d in dims)
                data = apply_cf_fill(ds[idx], attrs)
                if "_FillValue" in attrs:
                    attrs["source_FillValue"] = attrs.pop("_FillValue")
                data_vars[name] = (dims, data, attrs)
        gattrs = decode_attrs(f)
    gattrs["source_file"] = path.name
    return xr.Dataset(data_vars, coords=coords, attrs=gattrs)
