"""Reader tests on the 2026-02-17 example flight (numbers: appendix A §A.2.4, §A.10.5)."""

import numpy as np
import pytest

from nurture_oe.config import example_file, load_paths
from nurture_oe.io.cloudcube import parse_sensitivity, read_cloudcube
from nurture_oe.io.halo_aerowv import read_halo_aerowv
from nurture_oe.io.halo_cloudsurface import read_halo_cloudsurface
from nurture_oe.io.metnav import read_metnav


def test_paths_resolve():
    p = load_paths()
    assert p["nurture_example"].is_dir()
    for name in ("halo_aerowv", "halo_cloudsurface", "metnav", "cloudcube_1521"):
        assert example_file(name).is_file()


def test_aerowv_grid_and_state(aerowv):
    assert aerowv.sizes == {"time": 1689, "z": 845}
    assert np.all(np.diff(aerowv.z.values) > 0)
    assert aerowv.z.values[0] == pytest.approx(-300.0)
    assert float(aerowv["State_Type"]) == 2.0
    assert aerowv["Pressure"].attrs["units"] == "atm"
    assert aerowv["532_total_attn_bsc"].dims == ("time", "z")
    assert aerowv["Temperature"].dims == ("time", "z")  # mapped by length (no dimension scales)


def test_aerowv_attrs_verbatim(aerowv):
    a = aerowv["532_total_attn_bsc"].attrs
    assert "Units " in a and a["units"] == "km-1.sr-1"
    assert a["hdf5_path"] == "/DataProducts/532_total_attn_bsc"
    assert aerowv["z"].attrs["positive"] == "up "


def test_aerowv_units_conflict_warns():
    with pytest.warns(UserWarning, match="Angstrom_Spherical"):
        read_halo_aerowv(example_file("halo_aerowv"), groups=("DataProducts",), variables=["Angstrom_Spherical"])


def test_aerowv_time_range():
    ds = read_halo_aerowv(example_file("halo_aerowv"), groups=("DataProducts",), variables=["532_dep"],
                          time_range=("2026-02-17T15:30", "2026-02-17T16:00"))
    assert ds.time.values.min() >= np.datetime64("2026-02-17T15:30")
    assert ds.time.values.max() <= np.datetime64("2026-02-17T16:00")
    assert 170 <= ds.sizes["time"] <= 190


def test_cloudsurface(aerowv, cs_type):
    assert cs_type.sizes["profile"] == 33765
    assert np.array_equal(cs_type.z.values, aerowv.z.values)
    assert cs_type["altSurface"].dims == ("profile",)
    assert len(cs_type.attrs["000_Readme"]) == 18
    assert np.nanmedian(np.diff(cs_type.time_native.values * 3600)) == pytest.approx(0.5, abs=0.01)


def test_cloudsurface_wcp_coordinate():
    ds = read_halo_cloudsurface(example_file("halo_cloudsurface"), variables=["ExtinctionProfile532"],
                                time_range=("2026-02-17T15:30", "2026-02-17T15:31"))
    dz = ds.wcp_dz.values
    assert dz.size == 501
    assert np.median(np.diff(dz)) == pytest.approx(1.249, abs=0.002)
    assert dz.min() == pytest.approx(-374.7, abs=0.5)


@pytest.mark.parametrize("text,val", [("-16.32 dBZ @ 1 km", -16.32), ("-19.26 dBZ @ 1 km", -19.26)])
def test_parse_sensitivity(text, val):
    assert parse_sensitivity(text) == val


def test_cloudcube_slice_and_sensitivity():
    seg = ("2026-02-17T15:30:00", "2026-02-17T15:30:30")
    with pytest.warns(UserWarning, match="radar_sensitivity_ka_dbz"):
        cc = read_cloudcube(example_file("cloudcube_1521"), time_range=seg)
    assert cc.attrs["zmin_1km_dbz_ka"] == -16.32
    assert cc.sizes["range"] == 751
    assert np.median(np.diff(cc.range.values)) == pytest.approx(19.986, abs=0.001)
    assert cc["zm_w"].dims == ("time", "range")
    assert cc["look_vector_ecef"].dims == ("time", "nxyz")
    assert 290 <= cc.sizes["time"] <= 300
    assert np.isfinite(cc["alt"].values).all()


def test_cloudcube_1900_matches_pk_ka():
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        cc = read_cloudcube(example_file("cloudcube_1900"), profiles=(), time_range=("2026-02-17T19:00", "2026-02-17T19:00:05"))
    assert cc.attrs["zmin_1km_dbz_ka"] == -19.26


def test_metnav():
    m = read_metnav(example_file("metnav"))
    assert m.sizes["time"] == 17733
    assert np.all(np.diff(m.time_native.values) == 1.0)
    assert m["Pitch_Angle"].attrs["units"].startswith("deg")
    assert np.isnan(m["Track_Angle"].values).any()   # -99999 -> NaN
    assert np.nanmin(m["Track_Angle"].values) > -1000
