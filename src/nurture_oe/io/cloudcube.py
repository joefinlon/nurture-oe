"""Reader for CloudCube L2 (Ka/W/G) granules.

See ``docs/data_map.md`` §1.3. Native (time, range) arrays; ``range`` is slant range from
the aircraft (negative = above it). ``zm_*`` is NOT gas-corrected and ``katten_*`` is one-way
specific attenuation (JPL, 2026-10-08); the reader changes neither.
"""

import re
import warnings
from pathlib import Path

import netCDF4
import numpy as np
import xarray as xr

from ._hdf5 import parse_time_units, time_slice, to_datetime64

BANDS = ("ka", "w", "g")
_PROFILE = ("zm", "katten")
_RAY = ("s0",)
_NAV = ("alt", "lat", "lon", "roll", "pitch", "heading", "incidence_angle", "landmask")
PK_ZMIN_1KM_DBZ = {"ka": -19.26, "w": -32.64, "g": -37.58}  # appendix A §A.4 (v1.10 single values)
# v1.11 (appendix A §A.4): Ka is per granule on 02-17 (-16.32 for 15:13-18:59, -19.26 for 19:00-19:26)
PK_ZMIN_1KM_DBZ_RECORDED = {"ka": (-16.32, -19.26), "w": (-32.64,), "g": (-37.58,)}


def parse_sensitivity(text):
    """Parse ``'-16.32 dBZ @ 1 km'`` -> -16.32. Raises if the format differs."""
    m = re.fullmatch(r"\s*(-?\d+(?:\.\d+)?)\s*dBZ\s*@\s*1\s*km\s*", text)
    if not m:
        raise ValueError(f"unexpected sensitivity attribute {text!r}")
    return float(m.group(1))


def read_cloudcube(path, bands=BANDS, time_range=None, profiles=_PROFILE):
    """Read a CloudCube L2 granule into an ``xarray.Dataset``.

    Parameters
    ----------
    path : str or Path
    bands : sequence of {"ka", "w", "g"}
    time_range : (start, end), optional
        Inclusive UTC bounds. Only that slice is read (the arrays are contiguous, so a
        30 min slice of all three bands is ~0.3 GB per profile variable family).
    profiles : sequence of {"zm", "katten"}
        Which (time, range) families to read.

    Returns
    -------
    xarray.Dataset
        Dims ``time`` (datetime64, interval mid-time; raw seconds in ``time_native``),
        ``range`` (m) and ``nxyz`` (3). Global attributes are kept verbatim; the 1 km
        sensitivities are also parsed to floats as ``zmin_1km_dbz_<band>``.

    Warns
    -----
    UserWarning
        When a granule's declared 1 km sensitivity is none of the values recorded in PK
        (appendix A §A.4, v1.11). The file value governs (CLAUDE.md rule 3).
    """
    path = Path(path)
    with netCDF4.Dataset(path) as nc:
        nc.set_auto_mask(False)
        tv = nc["time"]
        epoch, unit = parse_time_units(tv.units)
        t_native = tv[:]
        t = to_datetime64(t_native, epoch, unit)
        sl = time_slice(t, time_range)

        def attrs_of(v):
            return {k: v.getncattr(k) for k in v.ncattrs()}

        def read(name):
            v = nc[name]
            a = attrs_of(v)
            data = v[sl] if v.dimensions and v.dimensions[0] == "time" else v[:]
            fv = a.get("_FillValue")
            if data.dtype.kind == "f":
                if fv is not None and np.isfinite(fv):
                    data = np.where(data == fv, np.nan, data)
                # netCDF default fill for variables without _FillValue
                data = np.where(data >= 9.9e36, np.nan, data)
            if "_FillValue" in a:
                a["source_FillValue"] = a.pop("_FillValue")
            return v.dimensions, data, a

        coords = {
            "time": ("time", t[sl], {k: v for k, v in attrs_of(tv).items() if k != "units"}),
            "time_native": ("time", t_native[sl], attrs_of(tv)),
            "range": read("range"),
        }
        data_vars = {n: read(n) for n in _NAV}
        data_vars["look_vector_ecef"] = read("look_vector_ecef")
        for b in bands:
            for fam in profiles:
                data_vars[f"{fam}_{b}"] = read(f"{fam}_{b}")
            for fam in _RAY:
                data_vars[f"{fam}_{b}"] = read(f"{fam}_{b}")
        gattrs = {k: nc.getncattr(k) for k in nc.ncattrs()}

    for b in BANDS:
        key = f"radar_sensitivity_{b}_dbz"
        if key in gattrs:
            val = parse_sensitivity(gattrs[key])
            gattrs[f"zmin_1km_dbz_{b}"] = val
            if not np.isclose(val, PK_ZMIN_1KM_DBZ_RECORDED[b]).any():
                warnings.warn(f"{path.name}: {key} = {gattrs[key]!r} is not among the PK values "
                              f"(appendix A §A.4: {PK_ZMIN_1KM_DBZ_RECORDED[b]} dBZ @ 1 km); "
                              f"the file value governs", stacklevel=2)
    gattrs["source_file"] = path.name
    return xr.Dataset(data_vars, coords=coords, attrs=gattrs)
