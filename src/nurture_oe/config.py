"""Data locations from ``config/paths.yaml``."""

import re
from pathlib import Path

import numpy as np
import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
PATHS_YAML = REPO_ROOT / "config" / "paths.yaml"


def load_paths(path=PATHS_YAML):
    """Read ``config/paths.yaml`` and expand ``~``.

    Returns
    -------
    dict of str -> pathlib.Path
    """
    with open(path) as f:
        raw = yaml.safe_load(f)
    return {k: Path(v).expanduser() for k, v in raw.items()}


def cache_dir():
    """Return ``cache_dir``, creating it if needed. Never inside a data store."""
    d = load_paths()["cache_dir"]
    d.mkdir(parents=True, exist_ok=True)
    return d


def find_granule(pattern, key="nurture_example"):
    """Return the single file in data store ``key`` matching glob ``pattern``.

    Raises
    ------
    FileNotFoundError
        If zero or more than one file matches.
    """
    root = load_paths()[key]
    hits = sorted(root.glob(pattern))
    if len(hits) != 1:
        raise FileNotFoundError(f"{len(hits)} files match {pattern!r} in {root}: {hits}")
    return hits[0]


EXAMPLE_PATTERNS = {
    "halo_aerowv": "nurture-HALO_G3_20260217_*.h5",
    "halo_cloudsurface": "NURTURE-HALO-CloudAndSurface_G3_20260217_*.h5",
    "metnav": "NURTURE-MetNav_G3_20260217_*.ict",
    "cloudcube_1513": "NURTURE-CLOUDCUBE_G3_20260217_*_S20260217a151300_*.nc",
    "cloudcube_1521": "NURTURE-CLOUDCUBE_G3_20260217_*_S20260217a152158_*.nc",
    "cloudcube_1900": "NURTURE-CLOUDCUBE_G3_20260217_*_S20260217a190000_*.nc",
    "halo_cloud": "Cloud_Plots/*_cloud_data.mat",
}
CLOUDCUBE_GLOB = "NURTURE-CLOUDCUBE_G3_*_KaWG.nc"


def example_file(name):
    """Path of a 2026-02-17 example granule by short name (see ``EXAMPLE_PATTERNS``)."""
    return find_granule(EXAMPLE_PATTERNS[name])


def cloudcube_span(path):
    """(start, end) datetime64[ns] of a CloudCube granule from its ``_S<date>a<HHMMSS>_E...`` name."""
    m = re.search(r"_S(\d{8})a(\d{6})_E(\d{8})a(\d{6})", Path(path).name)
    if not m:
        raise ValueError(f"no _S..._E... time stamps in {Path(path).name!r}")

    def t(d, hms):
        return np.datetime64(f"{d[:4]}-{d[4:6]}-{d[6:]}T{hms[:2]}:{hms[2:4]}:{hms[4:]}", "ns")
    return t(m.group(1), m.group(2)), t(m.group(3), m.group(4))


def cloudcube_granules(start, end, key="nurture_example"):
    """CloudCube granules in data store ``key`` that overlap ``[start, end]``, sorted by start.

    Returns
    -------
    list of (Path, start, end)
    """
    t0, t1 = np.datetime64(start, "ns"), np.datetime64(end, "ns")
    out = []
    for p in sorted(load_paths()[key].glob(CLOUDCUBE_GLOB)):
        s, e = cloudcube_span(p)
        if s <= t1 and e >= t0:
            out.append((p, s, e))
    return sorted(out, key=lambda x: x[1])
