"""Reader for the HALO 0.5 s "cloud" profile file (``<YYYYMMDD>_cloud_data.mat``).

Delivered 2026-10-23 for 2026-02-17 (``Cloud_Plots/`` of ``nurture_example``); see
``docs/data_map.md`` §1.6 and ``docs/reports/cloud_0p5s.md``. MATLAB v7.3 (HDF5): one struct
``data`` with groups ``nav``, ``cloudphase``, ``cloud``, ``column`` (+ ``column/OD``),
``cprofile``, ``sprofile``, ``surface``, ``reference`` and ``QA``.

This interim file declares no units, descriptions or fill values [GAP]. The reader returns
native values only. Where this repo has evidence for a meaning or unit, it attaches
``units_inferred`` / ``meaning_inferred`` / ``evidence`` (table ``INFERRED``), never ``units``.
A CF-compliant version (.mat or .h5) is announced: file layout lives in ``LAYOUT``, and any
attribute the file declares governs and is checked against ``INFERRED`` (CLAUDE.md rule 3).
"""

import warnings
from pathlib import Path

import h5py
import numpy as np
import xarray as xr

from ._hdf5 import apply_cf_fill, decode_attrs, flight_date_from_name, time_slice, to_datetime64

ROOT = "data"
TIME = "nav/gps_time"          # hours UTC since 00 UTC of the flight date (same as C&S Nav/gps_time)
GROUPS = ("nav", "cloudphase", "cloud", "column", "cprofile", "sprofile", "surface", "reference", "QA")

# Second-axis coordinates per group: dim name -> dataset path below ROOT (or a fixed length).
LAYOUT = {
    "cloudphase": {"z": "cloudphase/alt"},
    "cprofile": {"wcp_dz": "cprofile/alt"},
    "sprofile": {"sfc_dz": "sprofile/alt"},
    "column": {"layer": 10, "od_coord": "column/OD/OD_coord", "od_level": "column/OD/OD_levels"},
}

_MATLAB_ATTRS = {"MATLAB_class", "MATLAB_fields", "H5PATH", "MATLAB_empty"}

_EV_BSC = ("10 s window means ([t-5 s, t+5 s)) equal HALO Aerosol/WV {ref} (km-1 sr-1) x 1e-3: median ratio "
           "{ratio}, 2026-02-17 15:30-16:00 (scripts/check_halo_cloud.py); PI figure colorbar 'beta_p (Mm-1 sr-1)'")
INFERRED = {
    "nav/gps_time": {"units_inferred": "hours since 00 UTC of the flight date",
                     "evidence": "equals C&S Nav/gps_time except 37 profiles (<= 0.072 s)"},
    "nav/gps_alt": {"units_inferred": "m", "meaning_inferred": "aircraft altitude, MSL [DECISION D23]",
                    "evidence": "identical to C&S Nav/gps_alt"},
    "cloudphase/alt": {"units_inferred": "m", "meaning_inferred": "altitude grid, MSL",
                       "evidence": "identical to HALO Aerosol/WV z and C&S CloudMask/Altitude"},
    "cloudphase/bsc": {"units_inferred": "m-1 sr-1",
                       "meaning_inferred": "532 nm particulate backscatter coefficient beta_p, 0.5 s, unscreened",
                       "evidence": _EV_BSC.format(ref="532_bsc", ratio="1.01")},
    "cloudphase/bsc1064": {"meaning_inferred": "1064 nm backscatter, 0.5 s, unscreened; calibration and units unknown",
                           "evidence": "10 s window means are 0.55 x HALO 1064_bsc x 1e-3 (median; 5-95 % 0.10-0.78), "
                                       "so not the Aerosol/WV 1064 product [GAP]"},
    "cloudphase/dep": {"units_inferred": "1",
                       "meaning_inferred": "532 nm particulate depolarization ratio delta_p, 0.5 s, unscreened",
                       "evidence": "beta_p-weighted 10 s means match HALO 532_aer_dep (median difference <= 0.002 "
                                   "in every 532_bsr bin), not volume 532_dep (+0.15 at bsr < 2); PI figure label"},
    "cloudphase/dep1064": {"units_inferred": "1",
                           "meaning_inferred": "1064 nm depolarization ratio, 0.5 s, unscreened; particulate or volume "
                                               "undetermined",
                           "evidence": "10 s means differ from both HALO 1064_aer_dep and 1064_dep by 0.006-0.04"},
    "cloudphase/cp": {"flag_values_inferred": "0 1 2 3 4 5",
                      "flag_meanings_inferred": "WATER MIX ICE SPEC AER NLM",
                      "evidence": "identical to C&S CloudMask/Type on every pixel; labels from the colorbar of the "
                                  "PI's 20260217_F1_cloudphase.fig; NLM undefined [GAP G17]"},
    "cloud/IAB": {"meaning_inferred": "integrated attenuated backscatter, 532 nm",
                  "evidence": "identical to C&S Cloud/IAB532"},
}


def _relname(path, group):
    return path[len(group) + 1:].replace("/", "_")


def _coord_dims(f, group):
    """{length: (dim, values or None)} for the second-axis coordinates of ``group``."""
    out = {}
    for dim, src in LAYOUT.get(group, {}).items():
        if isinstance(src, int):
            out[src] = (dim, None, None)
        else:
            out[f[f"{ROOT}/{src}"].size] = (dim, np.ravel(f[f"{ROOT}/{src}"][...]), src)
    return out


def _attrs(ds, relpath):
    attrs = {k: v for k, v in decode_attrs(ds).items() if k not in _MATLAB_ATTRS}
    attrs["hdf5_path"] = ds.name
    if "MATLAB_class" in ds.attrs:
        attrs["matlab_class"] = ds.attrs["MATLAB_class"].decode()
    inf = INFERRED.get(relpath, {})
    if "units" in attrs and "units_inferred" in inf and str(attrs["units"]).strip() != inf["units_inferred"]:
        warnings.warn(f"{relpath}: file declares units={attrs['units']!r}, repo inferred "
                      f"{inf['units_inferred']!r}; the file governs", stacklevel=3)
    attrs.update({k: v for k, v in inf.items() if k.replace("_inferred", "") not in attrs})
    return attrs


def read_halo_cloud(path, group="cloudphase", variables=None, time_range=None):
    """Read one group of the HALO 0.5 s cloud profile file into an ``xarray.Dataset``.

    Parameters
    ----------
    path : str or Path
        ``<YYYYMMDD>_cloud_data.mat`` (MATLAB v7.3 / HDF5).
    group : str
        One of ``GROUPS``. ``nav`` (``gps_alt``, ``gps_lat``, ``gps_lon``) is always added.
    variables : sequence of str, optional
        Subset by name within the group (path below the group, ``/`` -> ``_``, e.g.
        ``"OD_ext"`` for ``column/OD/ext``).
    time_range : (start, end), optional
        Inclusive UTC bounds; only that slice of profiles is read.

    Returns
    -------
    xarray.Dataset
        Dim ``profile`` (0.5 s; coordinate ``time`` from ``nav/gps_time`` hours + the flight
        date in the file name, raw hours in ``time_native``). Second axes per ``LAYOUT``:
        ``z`` (``cloudphase``), ``wcp_dz`` (``cprofile``), ``sfc_dz`` (``sprofile``),
        ``layer`` / ``od_coord`` / ``od_level`` (``column``); (N, 1) datasets are squeezed.
        Datasets without a profile axis (e.g. ``QA/block``) get dims ``<name>_d0``, ``<name>_d1``.
        Values are native (no units conversion, no screening, NaN as stored).
    """
    if group not in GROUPS:
        raise ValueError(f"unknown group {group!r}; expected one of {GROUPS}")
    path = Path(path)
    with h5py.File(path, "r") as f:
        hours = f[f"{ROOT}/{TIME}"][:, 0]
        t = to_datetime64(hours, flight_date_from_name(path), "hours")
        if not np.all(np.diff(t) > np.timedelta64(0, "ns")):
            raise ValueError(f"{path.name}: {TIME} is not strictly increasing")
        sl = time_slice(t, time_range)
        n_prof = hours.size
        cdims = _coord_dims(f, group)
        coords = {"time": ("profile", t[sl], {"description": f"from {TIME} (hours UTC) + flight date in file name"}),
                  "time_native": ("profile", hours[sl], _attrs(f[f"{ROOT}/{TIME}"], TIME))}
        for dim, vals, src in cdims.values():
            if vals is not None:
                coords[dim] = (dim, vals, _attrs(f[f"{ROOT}/{src}"], src))
        coord_srcs = {src for _, _, src in cdims.values()}
        data_vars = {}

        def add(relpath, ds, rel_to):
            name = _relname(relpath, rel_to)
            if name in data_vars or name in coords:
                raise ValueError(f"variable name collision: {name!r} ({relpath})")
            attrs = _attrs(ds, relpath)
            if ds.ndim == 2 and ds.shape[0] == n_prof:
                data = apply_cf_fill(ds[sl], attrs)
                if ds.shape[1] == 1:
                    dims, data = ("profile",), data[:, 0]
                elif ds.shape[1] in cdims:
                    dims = ("profile", cdims[ds.shape[1]][0])
                else:
                    raise ValueError(f"{relpath}: unmapped second axis {ds.shape}; update LAYOUT")
            else:
                data = apply_cf_fill(ds[...], attrs)
                dims = tuple(f"{name}_d{i}" for i in range(ds.ndim))
            data_vars[name] = (dims, data, attrs)

        def visit(name, obj):
            relpath = f"{group}/{name}"
            if not isinstance(obj, h5py.Dataset) or relpath in coord_srcs or relpath == TIME:
                return
            if variables is not None and _relname(relpath, group) not in variables:
                return
            add(relpath, obj, group)

        f[f"{ROOT}/{group}"].visititems(visit)
        if group != "nav":
            for k in ("gps_alt", "gps_lat", "gps_lon"):
                add(f"nav/{k}", f[f"{ROOT}/nav/{k}"], "nav")
    attrs = {"source_file": path.name,
             "note": "interim MATLAB file without metadata; *_inferred attributes are this repo's inferences [GAP]"}
    return xr.Dataset(data_vars, coords=coords, attrs=attrs)
