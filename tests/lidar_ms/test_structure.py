"""Structure tests ST0 and ST3 (appendix C §C.4.3)."""
import subprocess
import sys
import textwrap

import numpy as np
import pytest

from msref import IN_NAMES, REF_PAIRS, load_input, load_ref


def _modules_after(stmt: str) -> set[str]:
    code = textwrap.dedent(f"""
        import sys
        {stmt}
        print("\\n".join(sorted(sys.modules)))
    """)
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    return set(out.stdout.split())


def _top(mods):
    return {m.split(".")[0] for m in mods}


def test_st0_import_boundary():
    """ST0: `import lidar_ms` loads no nurture_oe, numba, jax or torch; xarray only via lidar_ms.io.

    Run in a fresh interpreter: the repository's root conftest imports nurture_oe.
    """
    mods = _top(_modules_after("import lidar_ms"))
    for banned in ("nurture_oe", "numba", "jax", "torch", "xarray"):
        assert banned not in mods, banned


def test_st0_io_is_only_xarray_entry():
    """ST0: every lidar_ms submodule except io imports without xarray."""
    import pkgutil

    import lidar_ms
    names = [m.name for m in pkgutil.walk_packages(lidar_ms.__path__, "lidar_ms.")]
    others = [n for n in names if n != "lidar_ms.io"]
    if others:
        mods = _top(_modules_after("; ".join(f"import {n}" for n in others)))
        for banned in ("nurture_oe", "numba", "jax", "torch", "xarray"):
            assert banned not in mods, banned
    if "lidar_ms.io" in names:
        mods = _top(_modules_after("import lidar_ms.io"))
        assert "nurture_oe" not in mods


def test_st3_counts():
    assert len(IN_NAMES) == 15
    assert len(REF_PAIRS) == 50
    assert {nm for nm, _ in REF_PAIRS} <= set(IN_NAMES)


@pytest.mark.parametrize("name", IN_NAMES)
def test_st3_inputs(name):
    inp = load_input(name)
    assert inp.data.shape == (inp.n, 10)
    assert inp.wavelength > 0 and inp.rho_tr > 0 and len(inp.rho_fov) == 1
    assert np.all(inp.ncols >= 4)
    assert np.all(np.isfinite(inp.data[:, :4]))


@pytest.mark.parametrize("name,tag", REF_PAIRS)
def test_st3_refs(name, tag):
    inp = load_input(name)
    ref = load_ref(name, tag)
    ncol = 6 if "-hsrl" in ref.options else 5
    assert ref.data.shape == (inp.n, ncol)
    assert ref.options[0] == "-algorithms"
    np.testing.assert_array_equal(ref.data[:, 0], np.arange(1, inp.n + 1))
    np.testing.assert_allclose(ref.data[:, 1], inp.data[:, 0], rtol=1e-6)
