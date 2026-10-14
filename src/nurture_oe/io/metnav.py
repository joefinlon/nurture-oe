"""Reader for the G-III MetNav ICARTT (FFI 1001) file.

See ``docs/data_map.md`` §1.4. Missing (-99999) becomes NaN. LLOD/ULOD flags (-8888/-7777)
also become NaN in the data, but each variable gets boolean ``<name>_llod`` / ``<name>_ulod``
companions so the reason is not lost. Scale factors are applied as declared.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


def _split(line):
    return [x.strip() for x in line.split(",")]


def parse_icartt_header(path):
    """Parse an ICARTT FFI 1001 header.

    Returns
    -------
    dict
        ``nhead``, ``date`` (datetime64 of the data date), ``independent`` (name, units),
        ``variables`` (list of dicts: name, units, description, scale, missing),
        ``special_comments``, ``normal_comments`` (dict of KEY: value plus ``raw`` list).
    """
    with open(path, encoding="utf-8", errors="replace") as f:
        lines = [f.readline().rstrip("\n")]
        nhead, ffi = (int(x) for x in _split(lines[0])[:2])
        if ffi != 1001:
            raise ValueError(f"only FFI 1001 supported, got {ffi}")
        lines += [f.readline().rstrip("\n") for _ in range(nhead - 1)]
    d = _split(lines[6])
    date = np.datetime64(f"{int(d[0]):04d}-{int(d[1]):02d}-{int(d[2]):02d}", "ns")
    indep = _split(lines[8])
    nvar = int(lines[9])
    scales = [float(x) for x in _split(lines[10])]
    missing = [float(x) for x in _split(lines[11])]
    variables = []
    for i in range(nvar):
        parts = _split(lines[12 + i])
        variables.append({"name": parts[0], "units": parts[1] if len(parts) > 1 else "",
                          "description": ", ".join(parts[2:]), "scale": scales[i], "missing": missing[i]})
    k = 12 + nvar
    nspecial = int(lines[k])
    special = lines[k + 1:k + 1 + nspecial]
    k = k + 1 + nspecial
    nnormal = int(lines[k])
    normal_raw = lines[k + 1:k + 1 + nnormal]
    normal = {}
    for ln in normal_raw[:-1]:
        if ":" in ln:
            key, val = ln.split(":", 1)
            normal[key.strip()] = val.strip()
    return {"nhead": nhead, "date": date, "independent": (indep[0], indep[1] if len(indep) > 1 else ""),
            "variables": variables, "special_comments": special, "normal_comments": normal,
            "column_header": normal_raw[-1], "header_lines": lines}


def _flag(text):
    try:
        return float(str(text).split(";")[0])
    except ValueError:
        return None


def read_metnav(path):
    """Read a MetNav ICARTT file into an ``xarray.Dataset`` on dim ``time``.

    Returns
    -------
    xarray.Dataset
        ``time`` is datetime64 (data date + ``Time_Start`` seconds; raw seconds in
        ``time_native``). Each variable carries ``units`` and ``description`` as declared.
    """
    path = Path(path)
    h = parse_icartt_header(path)
    names = [h["independent"][0]] + [v["name"] for v in h["variables"]]
    cols = _split(h["column_header"])
    if cols != names:
        raise ValueError(f"column header {cols} does not match declared variables {names}")
    df = pd.read_csv(path, skiprows=h["nhead"], header=None, names=names)
    llod = _flag(h["normal_comments"].get("LLOD_FLAG"))
    ulod = _flag(h["normal_comments"].get("ULOD_FLAG"))
    t_native = df[names[0]].to_numpy(float)
    time = h["date"] + np.round(t_native * 1e9).astype("timedelta64[ns]")
    data_vars = {}
    for v in h["variables"]:
        raw = df[v["name"]].to_numpy(float)
        attrs = {"units": v["units"], "description": v["description"], "scale_factor_applied": v["scale"],
                 "missing_value_source": v["missing"]}
        is_llod = raw == llod if llod is not None else np.zeros(raw.shape, bool)
        is_ulod = raw == ulod if ulod is not None else np.zeros(raw.shape, bool)
        bad = (raw == v["missing"]) | is_llod | is_ulod
        data_vars[v["name"]] = ("time", np.where(bad, np.nan, raw * v["scale"]), attrs)
        if is_llod.any():
            data_vars[f"{v['name']}_llod"] = ("time", is_llod, {"description": "below LLOD flag"})
        if is_ulod.any():
            data_vars[f"{v['name']}_ulod"] = ("time", is_ulod, {"description": "above ULOD flag"})
    gattrs = {"source_file": path.name, **h["normal_comments"],
              "special_comments": h["special_comments"], "pi_line": h["header_lines"][1]}
    coords = {"time": ("time", time), "time_native": ("time", t_native,
                                                       {"units": f"{h['independent'][1]} since data date 00 UTC"})}
    return xr.Dataset(data_vars, coords=coords, attrs=gattrs)
