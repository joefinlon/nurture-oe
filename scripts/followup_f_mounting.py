"""Follow-up §A.8(5f): CloudCube antenna mounting, body-fixed vs stabilized [GAP G27].

For every CloudCube ray (all three 02-17 granules, subsampled to ~1 Hz):

* off-nadir angle theta from ``look_vector_ecef`` (geodetic local vertical);
* the body-fixed prediction theta_b = arccos(cos(pitch) cos(roll)) from MetNav attitude
  (CloudCube ``pitch``/``roll``/``heading`` are MetNav interpolated at lag 0, task (a));
* the look vector rotated into the aircraft body frame (x forward, y right, z down; heading,
  pitch up+, roll right+). A body-fixed antenna has a constant body-frame vector through turns;
  a stabilized one has a constant (near-vertical) theta instead.

Fits theta = a + b * theta_b and reports residuals overall, straight-and-level and in turns
(|roll| > 5 deg). Because ``look_vector_ecef`` turns out to be computed from the attitude, an
independent test compares the measured W surface-echo range with the height above the DEM:
body-fixed predicts r = h / cos(theta_b), stabilized r = h. Writes ``docs/reports/followup_f_mounting.json`` and
``figures/followup_f_mounting.png``.
"""

import json
import sys
import warnings

import numpy as np

from nurture_oe import coreg
from nurture_oe.config import REPO_ROOT, example_file
from nurture_oe.io.cloudcube import read_cloudcube
from nurture_oe.io.halo_aerowv import read_halo_aerowv

GRANULES = ("cloudcube_1513", "cloudcube_1521", "cloudcube_1900")
STEP = 10           # rays (~1 s)
TURN_DEG = 5.0
SFC_PEAK_MIN_DBZ = 10.0


def enu_basis(lat_deg, lon_deg):
    la, lo = np.deg2rad(lat_deg), np.deg2rad(lon_deg)
    e = np.stack([-np.sin(lo), np.cos(lo), np.zeros_like(lo)], -1)
    n = np.stack([-np.sin(la) * np.cos(lo), -np.sin(la) * np.sin(lo), np.cos(la)], -1)
    u = np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], -1)
    return e, n, u


def body_frame(lv, lat, lon, heading, pitch, roll):
    """Rotate ECEF unit vectors into the body frame (x fwd, y right, z down)."""
    e, n, u = enu_basis(lat, lon)
    ned = np.stack([(lv * n).sum(-1), (lv * e).sum(-1), -(lv * u).sum(-1)], -1)
    ps, th, ph = np.deg2rad(heading), np.deg2rad(pitch), np.deg2rad(roll)
    cps, sps, cth, sth, cph, sph = np.cos(ps), np.sin(ps), np.cos(th), np.sin(th), np.cos(ph), np.sin(ph)
    c = np.empty(lv.shape[:-1] + (3, 3))
    c[..., 0, :] = np.stack([cth * cps, cth * sps, -sth], -1)
    c[..., 1, :] = np.stack([sph * sth * cps - cph * sps, sph * sth * sps + cph * cps, sph * cth], -1)
    c[..., 2, :] = np.stack([cph * sth * cps + sph * sps, cph * sth * sps - sph * cps, cph * cth], -1)
    return np.einsum("...ij,...j->...i", c, ned)


def pct(x, q=(1, 5, 50, 95, 99)):
    x = np.asarray(x)[np.isfinite(x)]
    return np.percentile(x, q).round(3).tolist() if x.size else None


def main():
    warnings.simplefilter("ignore")
    halo = read_halo_aerowv(example_file("halo_aerowv"), groups=("Nav_Data", "UserInput"), variables=["gps_alt", "DEM_altitude"])
    th = halo.time.values.astype("datetime64[ns]").astype(np.int64)
    cols = {k: [] for k in ("t", "theta", "theta_b", "pitch", "roll", "b", "inc", "r_sfc", "pk_w", "h_agl")}
    for g in GRANULES:
        cc = read_cloudcube(example_file(g), bands=("w",), profiles=("zm",)).isel(time=slice(None, None, STEP))
        r = cc.range.values
        zw = np.where(r[None, :] < 3000, -np.inf, np.nan_to_num(cc["zm_w"].values, nan=-np.inf))
        j = np.argmax(zw, axis=1)
        cols["r_sfc"].append(r[j])
        cols["pk_w"].append(zw[np.arange(j.size), j])
        tc = cc.time.values.astype("datetime64[ns]").astype(np.int64)
        cols["h_agl"].append(np.interp(tc, th, halo["gps_alt"].values) - np.interp(tc, th, halo["DEM_altitude"].values))
        lv = cc["look_vector_ecef"].values
        lv = lv / np.linalg.norm(lv, axis=-1, keepdims=True)
        cth = coreg.cloudcube_cos_offnadir(cc).values
        cols["t"].append(cc.time.values)
        cols["theta"].append(np.degrees(np.arccos(np.clip(cth, -1, 1))))
        p, r = cc["pitch"].values, cc["roll"].values
        cols["theta_b"].append(np.degrees(np.arccos(np.cos(np.deg2rad(p)) * np.cos(np.deg2rad(r)))))
        cols["pitch"].append(p)
        cols["roll"].append(r)
        cols["inc"].append(cc["incidence_angle"].values)
        cols["b"].append(body_frame(lv, cc["lat"].values, cc["lon"].values, cc["heading"].values, p, r))
    d = {k: np.concatenate(v) for k, v in cols.items()}
    ok = np.isfinite(d["theta"]) & np.isfinite(d["theta_b"])
    turn = ok & (np.abs(d["roll"]) > TURN_DEG)
    level = ok & (np.abs(d["roll"]) < 1.0)
    A = np.stack([np.ones(ok.sum()), d["theta_b"][ok]], 1)
    coef, *_ = np.linalg.lstsq(A, d["theta"][ok], rcond=None)
    res = d["theta"] - (coef[0] + coef[1] * d["theta_b"])
    diff = d["theta"] - d["theta_b"]
    b = d["b"]
    fore = np.degrees(np.arctan2(b[:, 0], b[:, 2]))      # + = beam tilted forward
    side = np.degrees(np.arctan2(b[:, 1], b[:, 2]))      # + = beam tilted to starboard
    out = {
        "n_samples": int(ok.sum()), "n_turning(|roll|>5)": int(turn.sum()), "n_level(|roll|<1)": int(level.sum()),
        "roll_deg_pct_1_5_50_95_99": pct(d["roll"][ok]), "pitch_deg_pct_1_5_50_95_99": pct(d["pitch"][ok]),
        "theta_deg_pct": {"all": pct(d["theta"][ok]), "turning": pct(d["theta"][turn]), "level": pct(d["theta"][level])},
        "max_abs_roll_deg": round(float(np.nanmax(np.abs(d["roll"][ok]))), 2),
        "fit theta = a + b*theta_body": {"a_deg": round(float(coef[0]), 4), "b": round(float(coef[1]), 4),
                                        "residual_rms_deg_all": round(float(np.sqrt(np.nanmean(res[ok] ** 2))), 4),
                                        "residual_rms_deg_turning": round(float(np.sqrt(np.nanmean(res[turn] ** 2))), 4),
                                        "residual_pct_turning": pct(res[turn])},
        "theta_minus_theta_body_deg": {"all": pct(diff[ok]), "turning": pct(diff[turn]), "level": pct(diff[level])},
        "body_frame_look_vector": {
            "fore_aft_tilt_deg (+fwd)": {"all": pct(fore[ok]), "turning": pct(fore[turn]), "level": pct(fore[level])},
            "lateral_tilt_deg (+starboard)": {"all": pct(side[ok]), "turning": pct(side[turn]), "level": pct(side[level])},
            "angle_from_body_z_deg": {"all": pct(np.degrees(np.arccos(np.clip(b[ok, 2], -1, 1)))),
                                      "turning": pct(np.degrees(np.arccos(np.clip(b[turn, 2], -1, 1))))}},
        "theta_minus_incidence_angle_deg_pct": pct((d["theta"] - d["inc"])[ok]),
        "note": "look_vector_ecef equals the body z-axis from heading/pitch/roll: the geolocation assumes a body-fixed, "
                "exactly nadir beam, so the theta test above is circular. Independent test below.",
    }
    # independent test: measured surface range vs height above the DEM
    sfc = ok & np.isfinite(d["h_agl"]) & (d["pk_w"] > SFC_PEAK_MIN_DBZ) & (d["h_agl"] > 5000)
    with np.errstate(invalid="ignore"):
        th_eff = np.degrees(np.arccos(np.clip(d["h_agl"] / d["r_sfc"], -1, 1)))
    ratio = d["r_sfc"] / d["h_agl"]
    out["surface_range_test"] = {
        "selection": f"W surface peak > {SFC_PEAK_MIN_DBZ} dBZ beyond 3 km range; HALO gps_alt - DEM_altitude (nadir) > 5 km",
        "n_all": int(sfc.sum()), "n_turning": int((sfc & turn).sum()),
        "r_sfc/h_agl pct (turning)": pct(ratio[sfc & turn]),
        "1/cos(theta_body) pct (turning)": pct(1 / np.cos(np.deg2rad(d["theta_b"][sfc & turn]))),
        "theta_eff = arccos(h_agl/r_sfc) minus theta_body, deg": {"turning": pct((th_eff - d["theta_b"])[sfc & turn]),
                                                                 "level": pct((th_eff - d["theta_b"])[sfc & level])},
        "theta_eff deg (turning)": pct(th_eff[sfc & turn]),
        "theta_body deg (turning)": pct(d["theta_b"][sfc & turn]),
        "corr(theta_eff, theta_body) turning": round(float(np.corrcoef(th_eff[sfc & turn], d["theta_b"][sfc & turn])[0, 1]), 3)
        if (sfc & turn).sum() > 10 else None,
    }
    d["th_eff"], d["sfc"] = th_eff, sfc
    (REPO_ROOT / "docs" / "reports" / "followup_f_mounting.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))
    figure(d, ok, turn, coef, fore, side)


def figure(d, ok, turn, coef, fore, side):
    level = ok & ~turn
    import ultraplot as uplt
    fig, axs = uplt.subplots(ncols=3, refwidth=3.0, refaspect=1.1, share=False)
    ax = axs[0]
    ax.scatter(d["theta_b"][ok & ~turn], d["theta"][ok & ~turn], s=1, c="#86b6ef", label="|roll| ≤ 5°")
    ax.scatter(d["theta_b"][turn], d["theta"][turn], s=2, c="#b9471d", label="|roll| > 5°")
    x = np.array([0, np.nanmax(d["theta_b"][ok])])
    ax.plot(x, x, color="#0b0b0b", lw=0.8, label="1:1")
    ax.format(xlabel="arccos(cos pitch · cos roll) (°)", ylabel="off-nadir from look_vector_ecef (°)",
              title=f"fit: θ = {coef[0]:.3f} + {coef[1]:.4f} θ$_b$")
    ax.legend(loc="ul", ncols=1, fontsize=7)
    ax = axs[1]
    ax.scatter(d["roll"][ok], fore[ok], s=1, c="#3987e5", label="fore/aft")
    ax.scatter(d["roll"][ok], side[ok], s=1, c="#eb6834", label="lateral")
    ax.format(xlabel="MetNav roll (°)", ylabel="beam tilt in body frame (°)", title="look_vector_ecef in body axes")
    ax.legend(loc="ul", ncols=1, fontsize=7)
    ax = axs[2]
    m = d["sfc"]
    ax.scatter(d["theta_b"][m & ~turn], d["th_eff"][m & ~turn], s=1, c="#86b6ef", label="|roll| ≤ 5°")
    ax.scatter(d["theta_b"][m & turn], d["th_eff"][m & turn], s=3, c="#b9471d", label="|roll| > 5°")
    ax.plot(x, x, color="#0b0b0b", lw=0.8, label="body-fixed")
    ax.plot(x, [0, 0], color="#898781", lw=0.8, ls="--", label="stabilized (nadir)")
    ax.format(xlabel="arccos(cos pitch · cos roll) (°)", ylabel="arccos(h$_{AGL}$ / r$_{sfc}$) (°)",
              title="Independent: surface-echo range")
    ax.legend(loc="ul", ncols=1, fontsize=7)
    axs.format(suptitle="CloudCube beam vs aircraft attitude, 2026-02-17 (all granules, ~1 Hz)")
    p = REPO_ROOT / "docs" / "reports" / "figures" / "followup_f_mounting.png"
    fig.savefig(p, dpi=150)
    print(f"wrote {p}")


if __name__ == "__main__":
    sys.exit(main())
