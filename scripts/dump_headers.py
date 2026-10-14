"""Dump the full header of every 2026-02-17 example granule to docs/reports/headers/.

Task (a) step 1. Read-only. Attribute names are printed with repr() so trailing spaces show.
"""

import sys
from pathlib import Path

import h5py
import numpy as np

from nurture_oe.config import REPO_ROOT, EXAMPLE_PATTERNS, example_file

OUT = REPO_ROOT / "docs" / "reports" / "headers"


def _val(v, maxlen=400):
    if isinstance(v, bytes):
        v = v.decode("utf-8", "replace")
    if isinstance(v, np.ndarray):
        if v.dtype.kind in "SO":
            v = [x.decode("utf-8", "replace") if isinstance(x, bytes) else x for x in v.ravel()]
            v = v[0] if len(v) == 1 else v
        elif v.size == 1:
            v = v.item()
    s = repr(v)
    return s if len(s) <= maxlen else s[:maxlen] + f"... [{len(s)} chars]"


def _dims(ds):
    if "DIMENSION_LIST" not in ds.attrs:
        return None
    names = []
    for refs in ds.attrs["DIMENSION_LIST"]:
        try:
            names.append(ds.file[refs[0]].name)
        except Exception:
            names.append("?")
    return names


SKIP_ATTRS = {"DIMENSION_LIST", "REFERENCE_LIST", "CLASS", "NAME", "_Netcdf4Dimid", "_Netcdf4Coordinates"}


def dump_hdf5(path):
    lines = [f"# Header dump: {path.name}", f"# size: {path.stat().st_size:,} bytes", ""]
    with h5py.File(path, "r") as f:
        lines.append("## Global attributes")
        for k, v in f.attrs.items():
            if k not in SKIP_ATTRS:
                lines.append(f"  {k!r}: {_val(v)}")
        lines.append("")

        def visit(name, obj):
            if isinstance(obj, h5py.Group):
                lines.append(f"## Group /{name}")
                for k, v in obj.attrs.items():
                    lines.append(f"    {k!r}: {_val(v)}")
                return
            hidden = {k: _val(v) for k, v in obj.attrs.items() if k in ("CLASS", "NAME")}
            lines.append(f"- /{name}  shape={obj.shape} dtype={obj.dtype} chunks={obj.chunks} "
                         f"compression={obj.compression} fillvalue={obj.fillvalue!r}")
            dims = _dims(obj)
            if dims:
                lines.append(f"    dims: {dims}")
            if hidden.get("CLASS") == "'DIMENSION_SCALE'":
                lines.append(f"    dimension scale NAME: {hidden.get('NAME')}")
            for k, v in obj.attrs.items():
                if k not in SKIP_ATTRS:
                    lines.append(f"    {k!r}: {_val(v)}")
            if obj.ndim <= 1 and obj.size <= 4 and obj.dtype.kind not in "V":
                lines.append(f"    value: {_val(obj[()])}")

        f.visititems(visit)
    return "\n".join(lines) + "\n"


def dump_icartt(path):
    with open(path, encoding="utf-8", errors="replace") as f:
        first = f.readline()
        nhead = int(first.split(",")[0])
        head = [first.rstrip("\n")] + [f.readline().rstrip("\n") for _ in range(nhead - 1)]
        nrows = sum(1 for _ in f)
    return (f"# Header dump: {path.name}\n# size: {path.stat().st_size:,} bytes\n"
            f"# header lines: {nhead}; data rows: {nrows}\n\n" + "\n".join(head) + "\n")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for name in EXAMPLE_PATTERNS:
        path = example_file(name)
        text = dump_icartt(path) if path.suffix == ".ict" else dump_hdf5(path)
        out = OUT / f"{path.stem}.txt"
        out.write_text(text)
        print(f"{name}: {out.relative_to(REPO_ROOT)} ({text.count(chr(10))} lines)")


if __name__ == "__main__":
    sys.exit(main())
