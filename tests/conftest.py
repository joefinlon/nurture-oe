import warnings

import pytest

from nurture_oe.config import example_file
from nurture_oe.io.halo_aerowv import read_halo_aerowv
from nurture_oe.io.halo_cloudsurface import read_halo_cloudsurface


@pytest.fixture(scope="session")
def aerowv():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return read_halo_aerowv(example_file("halo_aerowv"))


@pytest.fixture(scope="session")
def cs_type():
    return read_halo_cloudsurface(example_file("halo_cloudsurface"), variables=["Type", "altSurface"])
