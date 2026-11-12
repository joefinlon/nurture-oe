"""Fixtures for the lidar_ms tests; loaders live in `msref.py`."""
import pytest

from msref import load_input, load_ref


@pytest.fixture(scope="session")
def ms_input():
    """Factory: name -> MsInput."""
    return load_input


@pytest.fixture(scope="session")
def ms_ref():
    """Factory: (name, tag) -> MsRef."""
    return load_ref
