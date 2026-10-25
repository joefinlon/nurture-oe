"""HALO 0.5 s cloud file: reader and footprint co-registration (2026-10-25).

Reader tests check content against the other HALO products rather than the interim file's
layout, so they should survive the announced CF-compliant version.
"""

import warnings

import numpy as np
import pytest
import xarray as xr

from nurture_oe import coreg
from nurture_oe.config import example_file
from nurture_oe.io.halo_aerowv import read_halo_aerowv
from nurture_oe.io.halo_cloud import read_halo_cloud

SEG = ("2026-02-17T15:40:00", "2026-02-17T15:45:00")


@pytest.fixture(scope="module")
def cloud():
    return read_halo_cloud(example_file("halo_cloud"), time_range=SEG)


def test_grid_and_time(cloud):
    halo = read_halo_aerowv(example_file("halo_aerowv"), groups=("DataProducts",), variables=["532_dep"],
                            time_range=SEG)
    assert np.array_equal(cloud.z.values, halo.z.values)
    t = cloud.time.values
    assert t[0] >= np.datetime64(SEG[0]) and t[-1] <= np.datetime64(SEG[1])
    assert np.median(np.diff(t).astype(float)) / 1e9 == pytest.approx(0.5007, abs=1e-3)


def test_native_values_no_declared_units(cloud):
    assert "units" not in cloud.bsc.attrs                    # the interim file declares none
    assert cloud.bsc.attrs["units_inferred"] == "m-1 sr-1"
    assert np.nanmin(cloud.bsc.values) < 0                     # unscreened, native noise kept


def test_phase_equals_cloudsurface_type(cloud, cs_type):
    typ = cs_type.Type.sel(profile=(cs_type.time >= cloud.time[0]) & (cs_type.time <= cloud.time[-1]))
    assert np.array_equal(cloud.cp.values, typ.values, equal_nan=True)


def test_column_group_names_and_dims():
    ds = read_halo_cloud(example_file("halo_cloud"), group="column", time_range=(SEG[0], "2026-02-17T15:40:05"))
    assert ds.OD_ext.dims == ("profile", "od_coord")
    assert ds.ext.dims == ("profile",)
    assert ds.cths.dims == ("profile", "layer")
    with pytest.raises(ValueError):
        read_halo_cloud(example_file("halo_cloud"), group="nope")


def test_bsc_matches_aerowv_532_bsc(cloud):
    """10 s window means of bsc reproduce 532_bsc (km-1 sr-1) x 1e-3."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        halo = read_halo_aerowv(example_file("halo_aerowv"), groups=("DataProducts",), variables=["532_bsc"],
                                time_range=SEG)
    halo = halo.isel(time=slice(1, -1))
    w = coreg.halo_windows(halo.time.values)
    b = cloud.bsc.values
    s, _ = coreg.window_sums(np.nan_to_num(b), cloud.time.values, w)
    n, _ = coreg.window_sums(np.isfinite(b).astype(float), cloud.time.values, w)
    with np.errstate(invalid="ignore", divide="ignore"):
        r = (s / n) / halo["532_bsc"].values
    ok = np.isfinite(r) & (halo["532_bsc"].values > 1e-2)
    assert np.median(r[ok]) * 1e3 == pytest.approx(1.0, abs=0.05)


# --- footprint co-registration (synthetic) ---------------------------------------------

def _moving_cc(tilt_deg, n=200, alt=10000.0, v=240.0, lat0=54.0, lon=-70.0):
    a, e2 = 6378137.0, 6.69437999014e-3
    t_s = np.arange(n) * 0.1
    phi = np.deg2rad(lat0)
    m_rad = a * (1 - e2) / (1 - e2 * np.sin(phi) ** 2) ** 1.5
    lat = lat0 + np.rad2deg(v * t_s / (m_rad + alt))            # due north at v (at flight altitude)
    up = coreg.local_up_ecef(lat, np.full(n, lon))
    la, lo = np.deg2rad(lat), np.deg2rad(lon)
    north = np.stack([-np.sin(la) * np.cos(lo), -np.sin(la) * np.sin(lo), np.cos(la)], axis=-1)
    th = np.deg2rad(tilt_deg)
    lv = -up * np.cos(th) + north * np.sin(th)
    rng = np.arange(0.0, 12000.0, 20.0)
    zm = np.zeros((n, rng.size))
    t = np.datetime64("2026-02-17T15:30:00", "ns") + (t_s * 1e9).astype("timedelta64[ns]")
    return xr.Dataset({"zm_w": (("time", "range"), zm), "alt": ("time", np.full(n, alt)), "lat": ("time", lat),
                       "lon": ("time", np.full(n, lon)), "look_vector_ecef": (("time", "nxyz"), lv)},
                      coords={"time": t, "range": rng})


def test_footprint_time_forward_tilt():
    cc = _moving_cc(2.4)
    z = np.array([8000.0, 2000.0])
    alt = cc.alt.values
    cth = coreg.cloudcube_cos_offnadir(cc).values
    _, r = coreg.gate_index_on_grid(cc, z, alt, cth)
    along, cross, speed = coreg.footprint_offsets(cc, r)
    _, dt = coreg.footprint_time(cc, z, alt, along, speed, lidar_pointing="nadir")
    mid = slice(30, -30)                                        # away from smoothing edges
    expect = (alt[0] - z) * np.tan(np.deg2rad(2.4)) / 240.0     # h tan(theta) / v
    assert np.allclose(dt[mid], expect, rtol=2e-3)
    assert np.abs(cross[mid]).max() < 0.5
    _, dt_bf = coreg.footprint_time(cc, z, alt, along, speed)            # body-fixed lidar: same beam
    assert np.allclose(dt_bf[mid], 0.0)


def test_footprint_nadir_equals_time_mode():
    cc = _moving_cc(0.0)
    cc["zm_w"][:] = np.linspace(-10, 10, cc.sizes["time"])[:, None]
    lid_t = cc.time.values[::5][4:-4] + np.timedelta64(250, "ms")
    halo = xr.Dataset({"gps_alt": ("time", np.full(lid_t.size, 10000.0))}, coords={"time": lid_t, "z": np.array([3000.0, 6000.0])})
    w = coreg.subprofile_windows(lid_t)
    a = coreg.cloudcube_to_halo_grid(cc, halo, w, fields=("zm_w",), aircraft_alt=cc.alt, assign="time")
    b = coreg.cloudcube_to_halo_grid(cc, halo, w, fields=("zm_w",), aircraft_alt=cc.alt, assign="footprint")
    assert np.allclose(a.zm_w.values, b.zm_w.values, equal_nan=True)
    assert np.all(b.n_rays.values == a.n_rays.values[:, None])


def test_footprint_rejects_overlapping_windows():
    w = coreg.halo_windows(np.datetime64("2026-02-17T15:30:00", "ns") + np.array([0, 9], dtype="timedelta64[s]"))
    with pytest.raises(ValueError):
        coreg.window_index(np.array([0.0]), w)
