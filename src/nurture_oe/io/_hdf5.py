"""Shared helpers for HDF5-based readers."""

import re
import warnings

import h5py
import numpy as np

_HIDDEN = {"DIMENSION_LIST", "REFERENCE_LIST", "CLASS", "NAME", "_Netcdf4Dimid", "_Netcdf4Coordinates"}


def decode(v):
    """Decode an HDF5 attribute value into a plain Python / NumPy value."""
    if isinstance(v, h5py.Empty):
        return ""
    if isinstance(v, bytes):
        return v.decode("utf-8", "replace")
    if isinstance(v, np.ndarray):
        if v.dtype.kind in "SO":
            out = [decode(x) for x in v.ravel()]
            return out[0] if len(out) == 1 else out
        return v.item() if v.size == 1 else v
    if isinstance(v, np.generic):
        return v.item()
    return v


def decode_attrs(obj):
    """All attributes of an HDF5 object, decoded, names kept verbatim (incl. trailing spaces)."""
    return {k: decode(v) for k, v in obj.attrs.items() if k not in _HIDDEN}


def _norm_units(u):
    u = str(u).strip().lower().replace("^", "").replace(".", "").replace(" ", "")
    return "1" if u in ("1", "ratio", "none", "unitless", "") else u


def check_units(attrs, name):
    """Warn if a legacy ``"Units "`` attribute disagrees with CF ``units`` (notation-normalized).

    The CF ``units`` attribute governs; both are kept. Returns True if they disagree.
    """
    legacy = attrs.get("Units ", attrs.get("Units"))
    cf = attrs.get("units")
    if legacy is None or cf is None:
        return False
    if _norm_units(legacy) != _norm_units(cf):
        warnings.warn(f"{name}: 'Units ' = {legacy!r} disagrees with CF units = {cf!r}; "
                      f"using CF units", stacklevel=3)
        return True
    return False


def dims_of(ds, sizes):
    """Dimension names of HDF5 dataset ``ds``.

    Uses attached dimension scales when present; otherwise maps each axis by its length
    to ``sizes`` ({name: length}). Raises if a length is ambiguous or unknown.
    """
    if "DIMENSION_LIST" in ds.attrs:
        names = []
        for refs in ds.attrs["DIMENSION_LIST"]:
            names.append(ds.file[refs[0]].name.rsplit("/", 1)[-1])
        return tuple(names)
    by_len = {}
    for name, n in sizes.items():
        by_len.setdefault(n, []).append(name)
    out = []
    for n in ds.shape:
        cands = by_len.get(n, [])
        if len(cands) != 1:
            raise ValueError(f"{ds.name}: cannot map axis of length {n} to a dimension "
                             f"(candidates {cands}; sizes {sizes})")
        out.append(cands[0])
    if len(set(out)) != len(out):
        raise ValueError(f"{ds.name}: repeated dimension {out}")
    return tuple(out)


def apply_cf_fill(data, attrs):
    """Replace a non-NaN CF ``_FillValue`` / ``missing_value`` with NaN (float data only)."""
    if data.dtype.kind != "f":
        return data
    for key in ("_FillValue", "missing_value"):
        fv = attrs.get(key)
        if fv is not None and np.isfinite(fv):
            data = np.where(data == fv, np.nan, data)
    return data


def parse_time_units(units):
    """Epoch (datetime64[ns]) and unit from ``'<unit> since YYYY-MM-DD HH:MM:SS'``."""
    m = re.match(r"\s*(\w+)\s+since\s+(\d{4}-\d{2}-\d{2})(?:[ T](\d{2}:\d{2}:\d{2}))?", units)
    if not m:
        raise ValueError(f"unparseable time units {units!r}")
    epoch = np.datetime64(f"{m.group(2)}T{m.group(3) or '00:00:00'}", "ns")
    return epoch, m.group(1)


def to_datetime64(values, epoch, unit):
    """Convert numeric offsets to datetime64[ns]; NaN -> NaT."""
    scale = {"seconds": 1e9, "second": 1e9, "s": 1e9, "hours": 3600e9, "hour": 3600e9}[unit]
    v = np.asarray(values, dtype=float)
    out = np.full(v.shape, np.datetime64("NaT", "ns"))
    ok = np.isfinite(v)
    out[ok] = epoch + np.round(v[ok] * scale).astype("timedelta64[ns]")
    return out


def time_slice(t, time_range):
    """Contiguous index slice of sorted datetime64 ``t`` within [start, end] (inclusive)."""
    if time_range is None:
        return slice(None)
    start, end = (np.datetime64(x, "ns") if x is not None else None for x in time_range)
    i0 = 0 if start is None else int(np.searchsorted(t, start, side="left"))
    i1 = len(t) if end is None else int(np.searchsorted(t, end, side="right"))
    return slice(i0, i1)
