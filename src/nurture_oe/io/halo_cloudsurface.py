"""Reader for the HALO Cloud & Surface L2 product.

See ``docs/data_map.md`` §1.2. Datasets stored as (N, 1) are squeezed to (N,); that is the
only reshaping. Variable names (some with spaces) are kept verbatim. No class remapping:
``CloudMask/Type`` is returned as stored (float, NaN where missing).
"""

from pathlib import Path

import h5py
import xarray as xr

from ._hdf5 import decode_attrs, flight_date_from_name, time_slice, to_datetime64


def read_halo_cloudsurface(path, variables=None, time_range=None):
    """Read a HALO Cloud & Surface L2 granule into an ``xarray.Dataset``.

    Parameters
    ----------
    path : str or Path
        HDF5 granule (``NURTURE-HALO-CloudAndSurface_G3_<date>_R0_L<n>.h5``).
    variables : sequence of str, optional
        Subset by bare name (e.g. ``["Type", "Cloud Top Height"]``). Nav is always read.
    time_range : (start, end), optional
        Inclusive UTC bounds; only that slice is read.

    Returns
    -------
    xarray.Dataset
        Dims ``profile`` (0.5 s; coordinate ``time`` as datetime64, from ``Nav/gps_time``
        in hours UTC since 00 UTC of the flight date in the file name; raw hours in
        ``time_native``), ``z`` (``CloudMask/Altitude``, 15 m), ``wcp_dz``
        (``WaterCloudProfile/Altitude``, relative to cloud top) and ``layer`` (10).
        ``000_Readme`` and ``000_revision_date`` go into the global attributes.

    Notes
    -----
    The flight date is not stored in the file; it is taken from the file name.
    """
    path = Path(path)
    with h5py.File(path, "r") as f:
        hours = f["Nav/gps_time"][:, 0]
        t = to_datetime64(hours, flight_date_from_name(path), "hours")
        sl = time_slice(t, time_range)
        n_prof = f["Nav/gps_time"].shape[0]
        z = f["DataProducts/CloudMask/Altitude"][0, :]
        dz = f["DataProducts/WaterCloudProfile/Altitude"][:, 0]

        coords = {
            "time": ("profile", t[sl], {"description": "from Nav/gps_time (hours UTC) + flight date"}),
            "time_native": ("profile", hours[sl], decode_attrs(f["Nav/gps_time"])),
            "z": ("z", z, decode_attrs(f["DataProducts/CloudMask/Altitude"])),
            "wcp_dz": ("wcp_dz", dz, decode_attrs(f["DataProducts/WaterCloudProfile/Altitude"])),
        }
        data_vars = {}

        def visit(name, ds):
            if not isinstance(ds, h5py.Dataset) or name.startswith("000_"):
                return
            if name in ("Nav/gps_time", "DataProducts/CloudMask/Altitude",
                        "DataProducts/WaterCloudProfile/Altitude"):
                return
            bare = name.rsplit("/", 1)[-1]
            if variables is not None and bare not in variables and not name.startswith("Nav/"):
                return
            if bare in data_vars:
                raise ValueError(f"variable name collision: {bare!r} ({name})")
            if ds.shape[0] != n_prof:
                raise ValueError(f"{name}: leading axis {ds.shape[0]} != {n_prof} profiles")
            attrs = decode_attrs(ds)
            attrs["hdf5_path"] = "/" + name
            data = ds[sl]
            if ds.shape[1] == 1:
                dims, data = ("profile",), data[:, 0]
            elif ds.shape[1] == z.size:
                dims = ("profile", "z")
            elif ds.shape[1] == dz.size:
                dims = ("profile", "wcp_dz")
            elif ds.shape[1] == 10:
                dims = ("profile", "layer")
            else:
                raise ValueError(f"{name}: unmapped second axis {ds.shape}")
            data_vars[bare] = (dims, data, attrs)

        f.visititems(visit)
        gattrs = {
            "000_Readme": [s for s in (x.decode() if isinstance(x, bytes) else str(x)
                                       for x in f["000_Readme"][0, :])],
            "000_revision_date": f["000_revision_date"][0, 0].decode()
            if isinstance(f["000_revision_date"][0, 0], bytes) else str(f["000_revision_date"][0, 0]),
            "source_file": path.name,
        }
    return xr.Dataset(data_vars, coords=coords, attrs=gattrs)
