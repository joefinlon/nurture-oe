"""Co-registration of CloudCube and HALO Cloud & Surface onto HALO lidar grids.

Task (d); appendix A §A.8(4), §A.11(b). Target grids: HALO Aerosol/WV 10 s profiles x 15 m
``z`` (the y-row prototype), or the 0.5 s profiles of C&S / the HALO cloud file x the same ``z``.

Conventions established from the 2026-02-17 data (``docs/reports/task_d_coreg.md``):

* A HALO 10 s profile averages the C&S 0.5 s sub-profiles in ``[t - 5 s, t + 5 s)``
  (empirical; best IAB532 agreement, 20 sub-profiles in >95 % of windows) [GAP G23].
* HALO altitudes (``z``, ``Nav_Data/gps_alt``) are MSL; CloudCube ``alt`` is above the
  WGS84 ellipsoid (~38 m lower here). Gate altitudes are therefore referenced to HALO's
  aircraft altitude by default.

Radar values are averaged in linear Z over each window. Because 20 m gates are mapped to
15 m bins, neighbouring output rows share gates: their errors are correlated [OPEN O29].

Rays are assigned to windows by ray time (``assign="time"``) or by the ground position of
each gate (``assign="footprint"``). Radar gates lie 0.2-0.5 km ahead of the aircraft nadir
(beam ~3 deg off nadir with pitch). On 02-17 the lidar is processed as body-fixed too, so the
relative offset is about zero and footprint matching reduces to time matching; a ~1 s
empirical lag remains unexplained [GAP G27].
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


def subprofile_windows(sample_time, half_width=None):
    """Windows ``[t - h, t + h)`` around each 0.5 s lidar profile time (C&S or cloud file).

    ``h`` defaults to half the median profile spacing, so regularly spaced profiles tile
    time without gaps or overlap. That the 0.5 s time stamp is the centre of its interval
    is assumed by analogy with the HALO 10 s rule [GAP G23].
    """
    t = np.asarray(sample_time, dtype="datetime64[ns]")
    if half_width is None:
        half_width = (np.median(np.diff(t).astype(np.int64)) // 2).astype("timedelta64[ns]")
    w = halo_windows(t, half_width)
    w.attrs["rule"] = "[t - half_width, t + half_width) around each 0.5 s profile time (centre assumed)"
    return w


def _bounds(sample_time, windows):
    """Member index range ``[i0, i1)`` of each window in sorted ``sample_time``.

    This is also the window -> member map a 10 s row needs on a 0.5 s retrieval grid
    (``y = mean_k F(x_k)`` over its members); overlapping windows share members.
    """
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


def _ecef(lat_deg, lon_deg, h):
    a, e2 = 6378137.0, 6.69437999014e-3
    lat, lon = np.deg2rad(lat_deg), np.deg2rad(lon_deg)
    n = a / np.sqrt(1.0 - e2 * np.sin(lat) ** 2)
    return np.stack([(n + h) * np.cos(lat) * np.cos(lon), (n + h) * np.cos(lat) * np.sin(lon),
                     (n * (1.0 - e2) + h) * np.sin(lat)], axis=-1)


def track_unit_and_speed(cc, smooth_s=2.0):
    """Horizontal ground-track unit vector (ECEF) and ground speed of each CloudCube ray.

    From finite differences of the WGS84 position (``lat``, ``lon``, ellipsoidal ``alt``,
    which is the right height for ECEF), smoothed with a ``smooth_s`` running mean because
    CloudCube nav is 1 Hz MetNav interpolated to rays. Altitudes for gates still come from
    HALO [DECISION D23]; this is used for direction and speed only.

    Returns
    -------
    unit : ndarray (n_rays, 3)
    speed : ndarray (n_rays,), m s-1
    """
    t = cc["time"].values.astype("datetime64[ns]").astype(np.int64) / 1e9
    pos = _ecef(cc["lat"].values, cc["lon"].values, cc["alt"].values)
    vel = np.gradient(pos, t, axis=0)
    n = max(1, int(round(smooth_s / np.median(np.diff(t)))))
    if n > 1:
        k = np.ones(n) / n
        num = np.stack([np.convolve(vel[:, i], k, mode="same") for i in range(3)], axis=-1)
        vel = num / np.convolve(np.ones(len(t)), k, mode="same")[:, None]
    up = local_up_ecef(cc["lat"].values, cc["lon"].values)
    vh = vel - (vel * up).sum(axis=-1, keepdims=True) * up
    speed = np.linalg.norm(vh, axis=-1)
    return vh / speed[:, None], speed


def footprint_offsets(cc, r_need, smooth_s=2.0):
    """Along- and cross-track horizontal offset (m) of each mapped gate from the aircraft nadir.

    ``offset = r * (l . e)`` with ``l`` the unit look vector and ``e`` the horizontal track
    (along) or ``up x track`` (cross, positive to the left) unit vector. Positive along-track
    means the gate is ahead of the aircraft.

    Returns
    -------
    along, cross : ndarray (n_rays, n_z)
    speed : ndarray (n_rays,), ground speed in m s-1
    """
    unit, speed = track_unit_and_speed(cc, smooth_s)
    up = local_up_ecef(cc["lat"].values, cc["lon"].values)
    lv = cc["look_vector_ecef"].values
    lv = lv / np.linalg.norm(lv, axis=-1, keepdims=True)
    left = np.cross(up, unit)
    a = (lv * unit).sum(axis=-1)
    c = (lv * left).sum(axis=-1)
    return np.asarray(r_need) * a[:, None], np.asarray(r_need) * c[:, None], speed


def footprint_time(cc, z, aircraft_alt, along, speed, lidar_pointing="body_fixed", time_offset_s=0.0):
    """Time (int ns, n_rays x n_z) at which the lidar footprint reaches each radar gate.

    ``t_eff = t_ray + (along - along_lidar) / v_ground + time_offset_s``, where ``along_lidar``
    is the along-track offset of the lidar beam at altitude ``z``:

    * ``"body_fixed"`` (default): the lidar points along the body z-axis, like the radar look
      vector (``look_vector_ecef`` is the body z-axis, follow-up §A.8(5f)); HALO's
      ``surface/offset_angle`` equals cos(pitch) cos(roll) on every 02-17 profile, so HALO
      processing assumes the same. Then ``along_lidar = along`` and only ``time_offset_s`` remains.
    * ``"nadir"``: ``along_lidar = 0``.
    * a float: along-track tangent of the lidar pointing, ``along_lidar = (alt - z) * tan``.

    Lever arm and boresights are unknown [GAP G27]. An empirical radar-lidar lag of about 1 s on
    02-17 is unexplained (docs/reports/cloud_0p5s.md); ``time_offset_s`` exposes it, default 0.
    Cross-track offsets cannot be corrected on a curtain; they are only reported.

    Returns
    -------
    t_eff : ndarray (n_rays, n_z), float ns (NaN where undefined)
    dt : ndarray (n_rays, n_z), s
    """
    if lidar_pointing == "body_fixed":
        lid = along
    elif lidar_pointing == "nadir":
        lid = 0.0
    else:
        lid = (np.asarray(aircraft_alt)[:, None] - np.asarray(z)[None, :]) * float(lidar_pointing)
    dt = (along - lid) / np.asarray(speed)[:, None] + time_offset_s
    t = cc["time"].values.astype("datetime64[ns]").astype(np.int64)
    with np.errstate(invalid="ignore"):
        te = t[:, None] + np.rint(dt * 1e9)
    return te, dt


def window_index(t_eff_ns, windows):
    """Index of the (non-overlapping) window containing each time in ``t_eff_ns`` (int ns); -1 if none."""
    lo = windows.t_lo.values.astype("datetime64[ns]").astype(np.int64)
    hi = windows.t_hi.values.astype("datetime64[ns]").astype(np.int64)
    if np.any(lo[1:] < hi[:-1]):
        raise ValueError("footprint assignment needs non-overlapping windows (e.g. subprofile_windows)")
    te = np.asarray(t_eff_ns, dtype=float)
    ok = np.isfinite(te)
    k = np.full(te.shape, -1, dtype=np.int64)
    kk = np.searchsorted(lo, te[ok].astype(np.int64), side="right") - 1
    inside = (kk >= 0) & (te[ok].astype(np.int64) < hi[np.clip(kk, 0, None)])
    k[ok] = np.where(inside, kk, -1)
    return k


def _grid_sums(values, k, nw):
    """Sum ``values`` (n_rays, n_z) into (nw, n_z) by window index ``k`` (n_rays, n_z); k < 0 dropped."""
    nz = values.shape[1]
    j = np.broadcast_to(np.arange(nz), k.shape)
    ok = k >= 0
    flat = k[ok] * nz + j[ok]
    return np.bincount(flat, weights=values[ok], minlength=nw * nz).reshape(nw, nz)


def cloudcube_to_halo_grid(cc, halo, windows, fields=("zm_ka", "zm_w", "zm_g"), extra=None,
                           aircraft_alt="halo", min_range=0.0, min_detect_frac=0.5,
                           assign="time", lidar_pointing="body_fixed", time_offset_s=0.0, smooth_s=2.0):
    """Window-average CloudCube reflectivities onto a HALO (time, z) grid in linear Z.

    Parameters
    ----------
    cc : xarray.Dataset
        From ``read_cloudcube`` (already sliced to the segment).
    halo : xarray.Dataset
        HALO product with coordinates ``time`` and ``z`` and variable ``gps_alt`` (MSL):
        Aerosol/WV, or the cloud file / C&S nav for 0.5 s grids. ``gps_alt`` is only needed
        when ``aircraft_alt="halo"``.
    windows : xarray.Dataset
        From ``halo_windows`` or ``subprofile_windows`` (restricted to the segment).
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
    assign : {"time", "footprint"}
        ``"time"``: a ray belongs to every window containing its time (10 s default; windows
        may overlap). ``"footprint"``: each (ray, z) gate goes to the window containing the
        time the lidar footprint reaches it (``footprint_time``); windows must not overlap.
    lidar_pointing : {"body_fixed", "nadir"} or float
        Lidar pointing for ``"footprint"`` (see ``footprint_time``) [GAP G27].
    time_offset_s : float
        Added radar-to-lidar time for ``"footprint"`` (default 0; see ``footprint_time``).
    smooth_s : float
        Smoothing (s) of the ground-track velocity for ``"footprint"``.

    Returns
    -------
    xarray.Dataset
        Dims (time, z) of the target grid. Per field: ``<f>`` (dBZ of the linear mean over
        detecting rays, NaN where detect fraction < ``min_detect_frac``), ``<f>_detect_frac``,
        and ``<f>_all`` (no detection threshold). Also ``n_rays``, ``range_mean``,
        ``cos_theta_mean``: (time,) in time mode; (time, z) in footprint mode, which adds
        ``along_track_offset_mean``, ``cross_track_offset_mean`` and ``footprint_dt_mean``.
    """
    if assign not in ("time", "footprint"):
        raise ValueError(f"assign must be 'time' or 'footprint', not {assign!r}")
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
    nw = windows.sizes["time"]
    out = {}

    if assign == "time":
        def gsum(v):
            return window_sums(v, t, windows)[0]
        n_rays = gsum(np.ones((t.size, 1)))[:, 0]
        out["n_rays"] = ("time", n_rays.astype(int), {"long_name": "CloudCube rays in window"})
        rcnt = gsum(valid.astype(float))
        with np.errstate(invalid="ignore", divide="ignore"):
            out["cos_theta_mean"] = ("time", gsum(cth.values[:, None])[:, 0] / n_rays,
                                     {"long_name": "mean cosine of off-nadir angle", "units": "1"})
        assignment = "ray time in window"
    else:
        along, cross, speed = footprint_offsets(cc, r_need, smooth_s)
        te, dt = footprint_time(cc, z, alt.values, along, speed, lidar_pointing, time_offset_s)
        k = window_index(np.where(valid, te, np.nan), windows)
        valid = valid & (k >= 0)

        def gsum(v):
            return _grid_sums(np.broadcast_to(v, k.shape), k, nw)
        rcnt = gsum(valid.astype(float))
        out["n_rays"] = (("time", "z"), rcnt.astype(int), {"long_name": "CloudCube rays mapped to pixel (footprint)"})
        with np.errstate(invalid="ignore", divide="ignore"):
            for name, v, attrs in (
                    ("cos_theta_mean", np.broadcast_to(cth.values[:, None], k.shape),
                     {"long_name": "mean cosine of off-nadir angle", "units": "1"}),
                    ("along_track_offset_mean", along,
                     {"long_name": "mean along-track offset of gates from aircraft nadir (+ ahead)", "units": "m"}),
                    ("cross_track_offset_mean", cross,
                     {"long_name": "mean cross-track offset of gates from aircraft nadir (+ left)", "units": "m"}),
                    ("footprint_dt_mean", dt,
                     {"long_name": "mean time from ray to lidar footprint (t_eff - t_ray)", "units": "s"})):
                out[name] = (("time", "z"), gsum(np.where(valid, v, 0.0)) / rcnt, attrs)
        assignment = (f"footprint: lidar footprint time of each gate in window; lidar_pointing={lidar_pointing}, "
                      f"time_offset_s={time_offset_s} [GAP G27]; track smoothing {smooth_s} s")

    with np.errstate(invalid="ignore", divide="ignore"):
        out["range_mean"] = (("time", "z"), gsum(np.where(valid, r_need, 0.0)) / rcnt,
                             {"long_name": "mean slant range of mapped gates", "units": "m"})

    items = [(f, cc[f]) for f in fields] + list((extra or {}).items())
    for name, da in items:
        v = np.take_along_axis(da.values, safe, axis=1)
        v = np.where(valid, v, np.nan)
        det = np.isfinite(v)
        s = gsum(np.where(det, 10.0 ** (v / 10.0), 0.0))
        nd = gsum(det.astype(float))
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
                "gate_assignment": f"nearest gate in altitude per ray; window by {assignment}",
                **{k: v for k, v in cc.attrs.items() if k.startswith("zmin_1km_dbz_")}}
    return ds
