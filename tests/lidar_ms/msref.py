"""Raw loaders for the frozen `multiscatter` 1.2.11 inputs (`.in`) and references (`.ref`).

These read the raw text only; turning an input into `lidar_ms` schema objects is
`lidar_ms.compat`'s job. The files under `reference/` are never written (appendix C §C.4.0).
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


REF_DIR = Path(__file__).parent / "reference"
FROZEN_DIR = Path(__file__).parent / "frozen"

IN_NAMES = tuple(sorted(p.stem for p in REF_DIR.glob("*.in")))
REF_PAIRS = tuple(sorted(tuple(p.stem.split("__")) for p in REF_DIR.glob("*.ref")))


@dataclass(frozen=True)
class MsInput:
    """Raw content of a multiscatter input file (version 2.x header)."""

    name: str
    comments: tuple[str, ...]
    n: int
    wavelength: float
    altitude: float
    rho_tr: float
    rho_fov: tuple[float, ...]
    data: np.ndarray        # (n, 10), NaN where a row has fewer columns
    ncols: np.ndarray       # (n,) number of values on each row


@dataclass(frozen=True)
class MsRef:
    """Raw content of a multiscatter reference output."""

    name: str
    tag: str
    options: tuple[str, ...]
    data: np.ndarray        # (n, 5) or (n, 6 with -hsrl)


def _numeric_lines(path: Path):
    comments, rows = [], []
    for line in path.read_text().splitlines():
        s = line.strip()
        if not s:
            continue
        if s.startswith("#"):
            comments.append(s)
        else:
            rows.append(s.split())
    return comments, rows


def load_input(name: str) -> MsInput:
    comments, rows = _numeric_lines(REF_DIR / f"{name}.in")
    head = rows[0]
    n = int(head[0])
    body = rows[1:1 + n]
    if len(body) != n:
        raise ValueError(f"{name}.in: header says {n} gates, found {len(body)}")
    data = np.full((n, 10), np.nan)
    ncols = np.zeros(n, dtype=int)
    for i, row in enumerate(body):
        vals = [float(v) for v in row[:10]]
        data[i, :len(vals)] = vals
        ncols[i] = len(vals)
    return MsInput(name=name, comments=tuple(comments), n=n, wavelength=float(head[1]),
                   altitude=float(head[2]), rho_tr=float(head[3]),
                   rho_fov=tuple(float(v) for v in head[4:]), data=data, ncols=ncols)


def load_ref(name: str, tag: str) -> MsRef:
    comments, rows = _numeric_lines(REF_DIR / f"{name}__{tag}.ref")
    opts = comments[0].split("options:", 1)[1].split() if comments and "options:" in comments[0] else []
    return MsRef(name=name, tag=tag, options=tuple(opts), data=np.array(rows, dtype=float))


def ref_tags(name: str) -> tuple[str, ...]:
    return tuple(t for nm, t in REF_PAIRS if nm == name)
