"""xarray / netCDF round trip of inputs and `Returns`, with provenance (appendix C §C.3.3).

This is the only `lidar_ms` module that imports xarray (§C.3.0 rule 1).

Dimensions: ``profile``, ``range`` (N_out output bins), ``bin`` (N input bins), ``half_level`` (N+1),
``cls``, ``order``. Every variable carries ``units`` and ``long_name``; ``flags`` carries CF
``flag_masks`` and ``flag_meanings``; `Returns.meta` is stored as global attributes (JSON for the
nested parts).
"""
from __future__ import annotations

import json

import numpy as np
import xarray as xr

from . import __version__
from .constants import FLAG_MASKS, FLAG_MEANINGS
from .schema import Instrument, OpticalProfile, RangeGrid, Returns

_BSC = "m-1 sr-1"
_RET_VARS = {
    "beta_par": (("profile", "range", "cls"), _BSC, "apparent backscatter, co-polarized"),
    "beta_perp": (("profile", "range", "cls"), _BSC, "apparent backscatter, cross-polarized"),
    "beta_wide": (("profile", "range", "cls"), _BSC, "wide-angle part of the apparent backscatter"),
    "beta_ss": (("profile", "bin", "cls"), _BSC, "single-scatter apparent backscatter (both polarizations)"),
    "ms_factor": (("profile", "half_level"), "1", "small-angle multiple-scattering factor M at half-levels"),
    "gamma2": (("profile", "half_level"), "rad2", "coangle variance of forward-scattered photons at half-levels"),
    "a_s": (("profile", "bin"), "1", "accumulated single-scatter fraction of the liquid classes"),
}
_OPTIONAL = {
    "beta_order": (("profile", "range", "cls", "order"), _BSC, "apparent backscatter by scattering order (1, 2, >=3)"),
    "stderr_par": (("profile", "range", "cls"), _BSC, "batch-mean standard error, co-polarized"),
    "stderr_perp": (("profile", "range", "cls"), _BSC, "batch-mean standard error, cross-polarized"),
}


def returns_to_dataset(ret: Returns) -> xr.Dataset:
    """`Returns` -> `xarray.Dataset`."""
    P, n_out, J = ret.beta_par.shape
    N = ret.beta_ss.shape[1]
    data = {}
    for nm, (dims, units, ln) in {**_RET_VARS, **_OPTIONAL}.items():
        val = getattr(ret, nm)
        if val is None:
            continue
        data[nm] = xr.Variable(dims, np.asarray(val), {"units": units, "long_name": ln})
    data["closure_used"] = xr.Variable(("profile",), np.asarray(ret.closure_used, dtype=str),
                                       {"long_name": "depolarization closure used (hu, near_backscatter, none)"})
    data["flags"] = xr.Variable(("profile", "bin"), np.asarray(ret.flags, dtype=np.uint16), {
        "long_name": "lidar_ms validity flags",
        "flag_masks": np.array(FLAG_MASKS, dtype=np.uint16),
        "flag_meanings": " ".join(m.lower() for m in FLAG_MEANINGS),
    })
    if ret.diagnostic is not None:
        data["diagnostic_closure"] = xr.Variable(("profile",), np.asarray(ret.diagnostic["closure"], dtype=str),
                                                 {"long_name": "closure not selected (diagnose_closures)"})
        data["diagnostic_beta_perp"] = xr.Variable(("profile", "bin", "cls"), np.asarray(ret.diagnostic["beta_perp"]),
                                                   {"units": _BSC,
                                                    "long_name": "cross-polarized backscatter from the closure not selected"})
    coords = {
        "r": xr.Variable(("profile", "range"), np.asarray(ret.r, dtype=float),
                         {"units": "m", "long_name": "range from the instrument to the bin centre"}),
        "cls": xr.Variable(("cls",), np.asarray(ret.names, dtype=str), {"long_name": "backscatterer class"}),
        "kind": xr.Variable(("cls",), np.asarray(ret.kinds, dtype=str), {"long_name": "class kind"}),
    }
    attrs = _meta_attrs(ret.meta)
    attrs.update(n_bins=np.int64(N), Conventions="CF-1.8")
    return xr.Dataset(data, coords=coords, attrs=attrs)


def dataset_to_returns(ds: xr.Dataset) -> Returns:
    """`xarray.Dataset` (from `returns_to_dataset`) -> `Returns`."""
    kw = {nm: ds[nm].values for nm in _RET_VARS}
    for nm in _OPTIONAL:
        kw[nm] = ds[nm].values if nm in ds else None
    diag = None
    if "diagnostic_closure" in ds:
        diag = {"closure": ds["diagnostic_closure"].values.astype(object),
                "beta_perp": ds["diagnostic_beta_perp"].values}
    return Returns(closure_used=ds["closure_used"].values.astype(object), flags=ds["flags"].values.astype(np.uint16),
                   kinds=tuple(str(k) for k in ds["kind"].values), names=tuple(str(c) for c in ds["cls"].values),
                   r=ds["r"].values, meta=_attrs_meta(ds.attrs), diagnostic=diag, **kw)


def write_returns(ret: Returns, path) -> None:
    """Write `Returns` to netCDF-4."""
    returns_to_dataset(ret).to_netcdf(path, engine="netcdf4", format="NETCDF4")


def read_returns(path) -> Returns:
    """Read `Returns` written by `write_returns`."""
    with xr.open_dataset(path, engine="netcdf4") as ds:
        return dataset_to_returns(ds.load())


def _meta_attrs(meta: dict) -> dict:
    attrs = {"lidar_ms_version": __version__}
    for k, v in meta.items():
        if isinstance(v, (str, int, float)) and not isinstance(v, bool):
            attrs[k] = v
        else:
            attrs[f"{k}_json"] = json.dumps(v, sort_keys=True, default=_jsonable)
    return attrs


def _attrs_meta(attrs: dict) -> dict:
    meta = {}
    for k, v in attrs.items():
        if k in ("lidar_ms_version", "n_bins", "Conventions"):
            continue
        if k.endswith("_json"):
            meta[k[:-5]] = json.loads(v)
        else:
            meta[k] = v.item() if isinstance(v, np.generic) else v
    return meta


def _jsonable(x):
    if isinstance(x, np.generic):
        return x.item()
    if isinstance(x, np.ndarray):
        return x.tolist()
    raise TypeError(f"not JSON serializable: {type(x)}")


# --- inputs ------------------------------------------------------------------------------------------
def inputs_to_dataset(grid: RangeGrid, instrument: Instrument, optics: OpticalProfile) -> xr.Dataset:
    """Grid, instrument and optics -> `xarray.Dataset` (real optics only)."""
    if optics.is_complex:
        raise ValueError("complex optics cannot be written")
    dims = ("profile", "bin", "cls")
    spec = {"alpha": ("m-1", "extinction coefficient"), "omega": ("1", "single-scattering albedo"),
            "g": ("1", "asymmetry factor"), "a_eq": ("m", "equivalent-area radius"),
            "beta": (_BSC, "backscatter coefficient"), "delta": ("1", "single-scatter linear depolarization ratio")}
    data = {nm: xr.Variable(dims, np.asarray(getattr(optics, nm)), {"units": u, "long_name": ln})
            for nm, (u, ln) in spec.items()}
    data["near_index"] = xr.Variable(dims, optics.near_index, {"long_name": "row of the class's near-180 table"})
    data["r"] = xr.Variable(("grid_profile", "bin"), grid.r, {"units": "m", "long_name": "range to bin centre"})
    if grid.z is not None:
        zd = ("grid_profile", "bin") if grid.z.ndim == 2 else ("bin",)
        data["z"] = xr.Variable(zd, grid.z, {"units": "m", "long_name": "altitude of bin centre (provenance)"})
    coords = {"cls": ("cls", np.asarray(optics.names, dtype=str)),
              "kind": ("cls", np.asarray(optics.kinds, dtype=str)),
              "descriptor": ("cls", np.asarray(optics.descriptors, dtype=str))}
    attrs = {"wavelength": instrument.wavelength, "rho_tr": instrument.rho_tr, "rho_fov": instrument.rho_fov,
             "receiver": instrument.receiver, "instrument_name": instrument.name,
             "lidar_ms_version": __version__, "Conventions": "CF-1.8"}
    return xr.Dataset(data, coords=coords, attrs=attrs)


def dataset_to_inputs(ds: xr.Dataset):
    """`xarray.Dataset` (from `inputs_to_dataset`) -> (RangeGrid, Instrument, OpticalProfile)."""
    grid = RangeGrid(ds["r"].values, z=ds["z"].values if "z" in ds else None)
    inst = Instrument(float(ds.attrs["wavelength"]), float(ds.attrs["rho_tr"]), float(ds.attrs["rho_fov"]),
                      receiver=str(ds.attrs["receiver"]), name=str(ds.attrs["instrument_name"]))
    opt = OpticalProfile(kinds=tuple(str(k) for k in ds["kind"].values),
                         descriptors=tuple(str(d) for d in ds["descriptor"].values),
                         names=tuple(str(c) for c in ds["cls"].values),
                         near_index=ds["near_index"].values.astype(np.int64),
                         **{nm: ds[nm].values for nm in ("alpha", "omega", "g", "a_eq", "beta", "delta")})
    return grid, inst, opt
