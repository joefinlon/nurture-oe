"""Follow-up §A.8(5a): geoid confirmation of [DECISION D23] on 2026-02-17.

For each CloudCube granule, compares HALO ``gps_alt`` (MSL) - CloudCube ``alt`` (WGS84
ellipsoid) with -N (EGM2008, EGM96) along track, and does the same with MetNav
``GPS_Altitude_MSL`` at 1 Hz. Writes ``docs/reports/followup_a_geoid.json`` and
``docs/reports/figures/followup_a_geoid.png``.
"""

import json
import sys
import warnings

import numpy as np

from nurture_oe.config import REPO_ROOT, example_file
from nurture_oe.geoid import geoid_undulation
from nurture_oe.io.cloudcube import read_cloudcube
from nurture_oe.io.halo_aerowv import read_halo_aerowv
from nurture_oe.io.metnav import read_metnav

GRANULES = ("cloudcube_1513", "cloudcube_1521", "cloudcube_1900")
PCT = [1, 5, 25, 50, 75, 95, 99]


def _ns(t):
    return np.asarray(t).astype("datetime64[ns]").astype(np.int64)


def stats(x):
    x = np.asarray(x)[np.isfinite(x)]
    return {"n": int(x.size), "pct_1_5_25_50_75_95_99": np.percentile(x, PCT).round(2).tolist(),
            "mean": round(float(x.mean()), 2), "std": round(float(x.std()), 2)}


def compare(t_ref, alt_msl, cc, vs=None):
    """Interpolate CloudCube nav to ``t_ref`` (inside the granule) and compare."""
    tc = _ns(cc.time.values)
    tr = _ns(t_ref)
    inside = (tr >= tc[0]) & (tr <= tc[-1]) & np.isfinite(alt_msl)
    tr = tr[inside]
    alt_cc = np.interp(tr, tc, cc["alt"].values)
    lat = np.interp(tr, tc, cc["lat"].values)
    lon = np.interp(tr, tc, cc["lon"].values)
    d = alt_msl[inside] - alt_cc
    out = {"msl_minus_cc_alt_m": stats(d)}
    res = {}
    for model in ("EGM2008", "EGM96"):
        n = geoid_undulation(lat, lon, model)
        res[model] = n
        out[f"minus_N_{model}_m"] = stats(-n)
        out[f"residual_{model}_m (msl - cc_alt + N)"] = stats(d + n)
        if vs is not None:
            level = np.abs(vs[inside]) < 0.5   # m/s, from d(GPS_Altitude_MSL)/dt
            out[f"residual_{model}_m_level_flight(|dz/dt|<0.5 m/s)"] = stats((d + n)[level])
    out["EGM2008_minus_EGM96_m"] = stats(res["EGM2008"] - res["EGM96"])
    return out, {"t": tr, "d": d, "N08": res["EGM2008"], "N96": res["EGM96"], "lat": lat, "lon": lon}


def d_level(halo, met, th, tm, ok, level):
    d = halo["gps_alt"].values[ok] - np.interp(th[ok], tm, met["GPS_Altitude_MSL"].values)
    return d[level]


def main():
    halo = read_halo_aerowv(example_file("halo_aerowv"), groups=("Nav_Data",), variables=["gps_alt"])
    met = read_metnav(example_file("metnav"))
    results = {"metnav_altitude_fields": [v for v in met.data_vars if "Alt" in v],
               "note": "MSL = ellipsoidal - N, so (MSL - CloudCube alt) should equal -N; residual = MSL - alt_CC + N"}
    series = {}
    vrate = np.gradient(met["GPS_Altitude_MSL"].values, met["time_native"].values)   # m/s
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for g in GRANULES:
            cc = read_cloudcube(example_file(g), profiles=())
            rh, sh = compare(halo.time.values, halo["gps_alt"].values, cc)
            rm, sm = compare(met.time.values, met["GPS_Altitude_MSL"].values, cc, vs=vrate)
            results[g] = {"file": cc.attrs["source_file"], "HALO_gps_alt_10s": rh, "MetNav_GPS_Altitude_MSL_1Hz": rm,
                          "lat_range": [round(float(np.nanmin(cc.lat)), 2), round(float(np.nanmax(cc.lat)), 2)],
                          "lon_range": [round(float(np.nanmin(cc.lon)), 2), round(float(np.nanmax(cc.lon)), 2)]}
            series[g] = (sh, sm)
    # HALO gps_alt vs MetNav MSL (both MSL; consistency of the reference)
    tm = _ns(met.time.values)
    th = _ns(halo.time.values)
    ok = (th >= tm[0]) & (th <= tm[-1])
    lag = {}
    for dt_s in (-1.0, -0.5, 0.0, 0.5, 1.0):
        d = halo["gps_alt"].values[ok] - np.interp(th[ok] + int(dt_s * 1e9), tm, met["GPS_Altitude_MSL"].values)
        lag[f"{dt_s:+.1f} s"] = stats(d)
    level = np.abs(np.interp(th[ok], tm, vrate)) < 0.5
    results["HALO_gps_alt_minus_MetNav_MSL_m"] = {"MetNav sampled at t + lag": lag,
                                                  "lag 0, level flight (|dz/dt|<0.5 m/s)": stats(d_level(halo, met, th, tm, ok, level))}
    out = REPO_ROOT / "docs" / "reports" / "followup_a_geoid.json"
    out.write_text(json.dumps(results, indent=1))
    print(json.dumps(results, indent=1))
    figure(series)


def figure(series):
    import matplotlib.dates as mdates
    import ultraplot as uplt
    fig, axs = uplt.subplots(nrows=2, ncols=3, refwidth=3.2, refaspect=1.4, sharex=False, sharey=False)
    for j, (g, (sh, sm)) in enumerate(series.items()):
        t = sm["t"].astype("datetime64[ns]")
        ax = axs[0, j]
        ax.plot(t, sm["d"], color="#3987e5", lw=1, label="MetNav MSL − CC alt")
        ax.plot(sh["t"].astype("datetime64[ns]"), sh["d"], ".", color="#0b0b0b", ms=2, label="HALO gps_alt − CC alt")
        ax.plot(t, -sm["N08"], color="#eb6834", lw=1.2, label="−N (EGM2008)")
        ax.plot(t, -sm["N96"], color="#1baf7a", lw=1, ls="--", label="−N (EGM96)")
        ax.format(title=g.replace("cloudcube_", "granule "), ylabel="m")
        axs[1, j].plot(t, sm["d"] + sm["N08"], color="#eb6834", lw=1, label="EGM2008")
        axs[1, j].plot(t, sm["d"] + sm["N96"], color="#1baf7a", lw=1, label="EGM96")
        axs[1, j].format(ylabel="residual (m)")
        for ax in (axs[0, j], axs[1, j]):
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[0, 0].legend(loc="ll", ncols=1, fontsize=7)
    axs[1, 0].legend(loc="ll", ncols=1, fontsize=7)
    axs.format(suptitle="MSL − CloudCube alt vs −N along track (2026-02-17); residual = MSL − alt + N",
               xlabel="UTC", xrotation=0)
    p = REPO_ROOT / "docs" / "reports" / "figures" / "followup_a_geoid.png"
    fig.savefig(p, dpi=150)
    print(f"wrote {p}")


if __name__ == "__main__":
    sys.exit(main())
