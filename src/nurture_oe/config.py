"""Data locations from ``config/paths.yaml``."""

from pathlib import Path

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
}


def example_file(name):
    """Path of a 2026-02-17 example granule by short name (see ``EXAMPLE_PATTERNS``)."""
    return find_granule(EXAMPLE_PATTERNS[name])
