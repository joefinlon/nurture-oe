"""Co-registration of CloudCube and HALO Cloud & Surface onto the HALO Aerosol/WV grid.

Task (d); appendix A §A.8(4), §A.11(b). Target grid: HALO 10 s profiles x 15 m ``z``.

Conventions established from the 2026-02-17 data (``docs/reports/task_d_coreg.md``):

* A HALO 10 s profile averages the C&S 0.5 s sub-profiles in ``[t - 5 s, t + 5 s)``
  (empirical; best IAB532 agreement, 20 sub-profiles in >95 % of windows) [GAP G23].
* HALO altitudes (``z``, ``Nav_Data/gps_alt``) are MSL; CloudCube ``alt`` is above the
  WGS84 ellipsoid (~38 m lower here). Gate altitudes are therefore referenced to HALO's
  aircraft altitude by default.

Radar values are averaged in linear Z over each window. Because 20 m gates are mapped to
15 m bins, neighbouring output rows share gates: their errors are correlated [OPEN O29].
"""

import numpy as np
import xarray as xr

HALF_WINDOW = np.timedelta64(5, "s")


def halo_windows(halo_time, half_width=HALF_WINDOW):
    """Averaging window ``[t - half_width, t + half_width)`` of each HALO profile.

    Parameters
    ----------
    halo_time : array-like of datetime64
    half_width : numpy.timedelta64

    Returns
    -------
    xarray.Dataset
        ``t_lo``, ``t_hi`` on dim ``time``. Windows may overlap where HALO spacing < 10 s.
    """
    t = np.asarray(halo_time, dtype="datetime64[ns]")
    return xr.Dataset({"t_lo": ("time", t - half_width), "t_hi": ("time", t + half_width)},
                      coords={"time": t}, attrs={"half_width_s": half_width / np.timedelta64(1, "s"),
                                                 "rule": "[t - half_width, t + half_width)"})


def _bounds(sample_time, windows):
    t = np.asarray(sample_time, dtype="datetime64[ns]")
    i0 = np.searchsorted(t, windows.t_lo.values, side="left")
    i1 = np.searchsorted(t, windows.t_hi.values, side="left")
    return i0, i1


def window_sums(values, sample_time, windows):
    """Sum of ``values`` (n_samples, ...) over each window via cumulative sums.

    Returns
    -------
    sums : ndarray (n_windows, ...)
    counts : ndarray (n_windows,)  number of samples in each window
    """
    i0, i1 = _bounds(sample_time, windows)
    cs = np.concatenate([np.zeros((1,) + values.shape[1:]), np.cumsum(values, axis=0)])
    return cs[i1] - cs[i0], i1 - i0


def cloudmask_fractions(cs, windows, classes=(0, 1, 2, 3, 4, 5)):
    """Fraction of C&S sub-profiles in each window carrying each ``CloudMask/Type`` class.

    The denominator is all sub-profiles in the window (NaN counts as "no class").

    Returns
    -------
    xarray.Dataset
        ``frac_<k>`` (time, z) per class, ``n_sub`` (time).
    """
    typ = cs["Type"].values
    out = {}
    n = None
    for k in classes:
        s, n = window_sums((typ == k).astype(np.float64), cs["time"].values, windows)
        with np.errstate(invalid="ignore", divide="ignore"):
            out[f"frac_{k}"] = (("time", "z"), s / n[:, None])
    out["n_sub"] = ("time", n)
    return xr.Dataset(out, coords={"time": windows.time.values, "z": cs.z.values},
                      attrs={"source": cs.attrs.get("source_file", ""), "window_rule": windows.attrs["rule"]})


def phase_masks(fracs, threshold=0.5):
    """Liquid / ice / cloud masks as in appendix A §A.10.5.

    ``cloud``: classes 0-3 in >= ``threshold`` of sub-profiles; ``liquid``: class 0 in
    >= ``threshold``; ``ice``: classes 1-2 in >= ``threshold``.
    """
    cloud = (fracs.frac_0 + fracs.frac_1 + fracs.frac_2 + fracs.frac_3) >= threshold
    liquid = fracs.frac_0 >= threshold
    ice = (fracs.frac_1 + fracs.frac_2) >= threshold
    return xr.Dataset({"cloud": cloud, "liquid": liquid, "ice": ice}, attrs={"threshold": threshold})


def local_up_ecef(lat_deg, lon_deg):
    """Geodetic local-up unit vector in ECEF, shape (..., 3)."""
    lat, lon = np.deg2rad(lat_deg), np.deg2rad(lon_deg)
    return np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)], axis=-1)


def cloudcube_cos_offnadir(cc):
    """cos of the beam off-nadir angle, ``-look_vector_ecef . up`` (dim ``time``).

    The look vector is renormalized (98 rays of the 15:21 granule deviate from unit length
    by up to 4e-4, i.e. ~5 m of altitude at 12 km range). The flat-Earth form
    ``z = alt - r cos(theta)`` matches exact ECEF geolocation to 1.4 m (99th pct) at 12 km.
    """
    up = local_up_ecef(cc["lat"].values, cc["lon"].values)
    lv = cc["look_vector_ecef"].values
    lv = lv / np.linalg.norm(lv, axis=-1, keepdims=True)
    c = -(lv * up).sum(axis=-1)
    return xr.DataArray(c, dims="time", coords={"time": cc.time.values},
                        attrs={"long_name": "cosine of beam off-nadir angle from look_vector_ecef", "units": "1"})


def aircraft_alt_from_halo(cc, halo):
    """HALO ``Nav_Data/gps_alt`` (MSL) interpolated to CloudCube ray times."""
    th = halo["time"].values.astype("datetime64[ns]").astype(np.int64)
    tc = cc["time"].values.astype("datetime64[ns]").astype(np.int64)
    alt = np.interp(tc, th, halo["gps_alt"].values, left=np.nan, right=np.nan)
    return xr.DataArray(alt, dims="time", coords={"time": cc.time.values},
                        attrs={"long_name": "aircraft altitude from HALO Nav_Data/gps_alt", "units": "m",
                               "datum": "MSL (matches MetNav GPS_Altitude_MSL)"})


def gas_corrected_display(cc, band):
    """Display-only gas correction ``Z_gc(r) = zm + 2 * int_0^r katten dr'`` (dBZ).

    ``katten`` is one-way specific attenuation in dB/km (JPL, 2026-10-08), integrated
    outward from the aircraft with the trapezoid rule; gates at r <= 0 get no correction,
    and NaN ``katten`` contributes zero. ``zm`` itself enters y as delivered; this is for
    figures only (appendix A §A.4.1).
    """
    r_km = cc["range"].values / 1000.0
    k = np.nan_to_num(cc[f"katten_{band}"].values)
    seg = 0.5 * (k[:, 1:] + k[:, :-1]) * np.diff(r_km)[None, :]
    seg = np.where(r_km[None, :-1] >= 0, seg, 0.0)   # only segments starting at or beyond the aircraft
    one_way = np.concatenate([np.zeros((k.shape[0], 1)), np.cumsum(seg, axis=1)], axis=1)
    out = cc[f"zm_{band}"] + 2.0 * xr.DataArray(one_way, dims=("time", "range"))
    out.attrs = {"long_name": f"{band.upper()}-band reflectivity, gas-corrected only (display)",
                 "units": "dBZ", "note": "zm + 2*int katten dr from the aircraft; not for y"}
    return out


def gate_index_on_grid(cc, z, aircraft_alt, cos_theta, min_range=0.0):
    """Nearest-gate index for every (ray, HALO z bin).

    Slant range needed to reach altitude ``z`` is ``(aircraft_alt - z) / cos_theta``; the
    nearest gate on the uniform ``range`` grid is taken. Bins outside ``[min_range,
    range.max()]`` get -1.

    Returns
    -------
    idx : ndarray of int (n_rays, n_z)
    r_need : ndarray (n_rays, n_z), slant range (m) of each bin centre
    """
    rng = cc["range"].values
    dr = np.median(np.diff(rng))
    if not np.allclose(np.diff(rng), dr, rtol=1e-3):
        raise ValueError("range grid is not uniform")
    r_need = (np.asarray(aircraft_alt)[:, None] - np.asarray(z)[None, :]) / np.asarray(cos_theta)[:, None]
    idx = np.rint((r_need - rng[0]) / dr)
    bad = ~np.isfinite(idx) | (r_need < min_range) | (idx < 0) | (idx > rng.size - 1)
    idx = np.where(bad, -1, idx).astype(np.int64)
    return idx, r_need


def cloudcube_to_halo_grid(cc, halo, windows, fields=("zm_ka", "zm_w", "zm_g"), extra=None,
                           aircraft_alt="halo", min_range=0.0, min_detect_frac=0.5):
    """Window-average CloudCube reflectivities onto the HALO (time, z) grid in linear Z.

    Parameters
    ----------
    cc : xarray.Dataset
        From ``read_cloudcube`` (already sliced to the segment).
    halo : xarray.Dataset
        HALO Aerosol/WV with ``gps_alt`` (only needed when ``aircraft_alt="halo"``).
    windows : xarray.Dataset
        From ``halo_windows`` (restricted to the segment).
    fields : sequence of str
        dBZ variables of ``cc`` to aggregate.
    extra : dict of str -> DataArray (time, range), optional
        Additional dBZ fields (e.g. ``gas_corrected_display`` output) to aggregate.
    aircraft_alt : {"halo", "cloudcube"} or DataArray
        Reference for gate altitudes. ``"halo"`` (default) uses HALO MSL ``gps_alt`` so the
        result sits on HALO's own vertical datum; ``"cloudcube"`` uses ellipsoidal ``alt``.
    min_range : float
        Minimum slant range (m) used. Radar near-field limit not documented [GAP G1].
    min_detect_frac : float
        A pixel is valid if at least this fraction of its rays detect (finite value).

    Returns
    -------
    xarray.Dataset
        Dims (time, z) of the HALO grid. Per field: ``<f>`` (dBZ of the linear mean over
        detecting rays, NaN where detect fraction < ``min_detect_frac``), ``<f>_detect_frac``,
        and ``<f>_all`` (no detection threshold). Also ``n_rays`` (time), ``range_mean``
        (time, z), ``cos_theta_mean`` (time).
    """
    if isinstance(aircraft_alt, str):
        alt = aircraft_alt_from_halo(cc, halo) if aircraft_alt == "halo" else cc["alt"]
        alt_ref = aircraft_alt
    else:
        alt, alt_ref = aircraft_alt, "user"
    cth = cloudcube_cos_offnadir(cc)
    z = halo["z"].values
    idx, r_need = gate_index_on_grid(cc, z, alt.values, cth.values, min_range=min_range)
    valid = idx >= 0
    safe = np.where(valid, idx, 0)
    t = cc["time"].values

    n_rays = window_sums(np.ones((t.size, 1)), t, windows)[0][:, 0]
    out = {"n_rays": ("time", n_rays.astype(int), {"long_name": "CloudCube rays in window"})}
    rsum, _ = window_sums(np.where(valid, r_need, 0.0), t, windows)
    rcnt, _ = window_sums(valid.astype(float), t, windows)
    with np.errstate(invalid="ignore", divide="ignore"):
        out["range_mean"] = (("time", "z"), rsum / rcnt, {"long_name": "mean slant range of mapped gates", "units": "m"})
        out["cos_theta_mean"] = ("time", window_sums(cth.values[:, None], t, windows)[0][:, 0] / n_rays,
                                 {"long_name": "mean cosine of off-nadir angle", "units": "1"})

    items = [(f, cc[f]) for f in fields] + list((extra or {}).items())
    for name, da in items:
        v = np.take_along_axis(da.values, safe, axis=1)
        v = np.where(valid, v, np.nan)
        det = np.isfinite(v)
        lin = np.where(det, 10.0 ** (v / 10.0), 0.0)
        s, _ = window_sums(lin, t, windows)
        nd, _ = window_sums(det.astype(float), t, windows)
        with np.errstate(invalid="ignore", divide="ignore"):
            frac = nd / rcnt
            dbz_all = 10.0 * np.log10(s / nd)
        dbz = np.where(frac >= min_detect_frac, dbz_all, np.nan)
        attrs = dict(da.attrs)
        attrs.update({"aggregation": "10*log10(mean of linear Z over detecting rays in window)",
                      "min_detect_frac": min_detect_frac})
        out[name] = (("time", "z"), dbz, attrs)
        out[f"{name}_all"] = (("time", "z"), dbz_all, {**attrs, "min_detect_frac": 0.0})
        out[f"{name}_detect_frac"] = (("time", "z"), frac, {"long_name": f"fraction of rays with finite {name}"})

    ds = xr.Dataset(out, coords={"time": windows.time.values, "z": z})
    ds.attrs = {"source_file": cc.attrs.get("source_file", ""), "window_rule": windows.attrs["rule"],
                "aircraft_alt_reference": alt_ref, "min_range_m": min_range,
                "gate_assignment": "nearest gate in altitude per ray",
                **{k: v for k, v in cc.attrs.items() if k.startswith("zmin_1km_dbz_")}}
    return ds
