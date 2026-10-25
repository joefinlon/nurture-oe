"""Four-panel STM figure (appendix A §A.8(4)) from the cached co-registration.

Run ``scripts/run_task_d.py`` first. Writes ``docs/reports/figures/stm_20260217_1530-1600.png``.
With ``--0p5s``, plots the 0.5 s x 15 m file of ``scripts/build_coreg_0p5s.py`` instead
(``..._0p5s.png``): native 0.5 s particulate beta_p and delta_p at 532 nm from the HALO cloud
file; radar by footprint with body-fixed beams (= ray time).
"""

import sys

import matplotlib.colors as mcolors
import matplotlib.dates as mdates
import numpy as np
import ultraplot as uplt
import xarray as xr

from nurture_oe.config import REPO_ROOT, cache_dir

# One-hue sequential ramps and a blue<->red diverging map with a gray midpoint
# (dataviz reference palette).
BLUE = ["#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]
ORANGE = ["#fbe0d4", "#f3a07e", "#eb6834", "#b9471d", "#7a2c10"]
AQUA = ["#d2f1e5", "#7fd4b3", "#1baf7a", "#13805a", "#0b5139"]
DIVERGING = ["#0d366b", "#3987e5", "#9ec5f4", "#f0efec", "#f3b3b2", "#e34948", "#a3201f"]
SURFACE, INK, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#898781", "#e1e0d9"


def ramp(colors, name):
    cm = mcolors.LinearSegmentedColormap.from_list(name, colors)
    cm.set_bad(SURFACE)
    return cm


RADAR_PANELS = [
    ("zgc_w", r"CloudCube $Z_W$ — gas-corrected only (display); masked <250 m above lidar surface", "dBZ",
     dict(cmap=ramp(AQUA, "aqua"), norm=mcolors.Normalize(-30, 20)), None),
    ("dwr_wg_gc", r"DWR$_{W-G}$ — gas-corrected only; contour: HALO liquid (phase class 0)", "dB",
     dict(cmap=ramp(DIVERGING, "div"), norm=mcolors.TwoSlopeNorm(vcenter=0, vmin=-6, vmax=12)), [-6, -3, 0, 3, 6, 9, 12]),
]
BSC_TICKS = [1e-4, 1e-3, 1e-2, 1e-1]
PANELS_10S = [
    ("532_total_attn_bsc", r"HALO $\beta'_{532}$", r"km$^{-1}$ sr$^{-1}$",
     dict(cmap=ramp(BLUE, "blue"), norm=mcolors.LogNorm(1e-4, 1e-1)), BSC_TICKS),
    ("532_dep", r"HALO $\delta_{532}$ (volume)", "1",
     dict(cmap=ramp(ORANGE, "orange"), norm=mcolors.Normalize(0, 0.5)), None),
] + RADAR_PANELS
PANELS_0P5S = [
    ("bsc532_km", r"HALO $\beta_{p,532}$ (particulate, 0.5 s; units inferred)", r"km$^{-1}$ sr$^{-1}$",
     dict(cmap=ramp(BLUE, "blue"), norm=mcolors.LogNorm(1e-4, 1e-1)), BSC_TICKS),
    ("dep532_shown", r"HALO $\delta_{p,532}$ (particulate, 0.5 s; shown where $\beta_p \geq 10^{-3}$ km$^{-1}$ sr$^{-1}$)", "1",
     dict(cmap=ramp(ORANGE, "orange"), norm=mcolors.Normalize(0, 0.5)), None),
] + RADAR_PANELS


def add_display_fields(ds, dep_bsc_min_km=1e-3):
    """0.5 s file: beta_p in km-1 sr-1 (inferred m-1 sr-1 x 1e3) and delta_p masked where beta_p is small (display)."""
    bsc_km = ds["bsc532"] * 1e3
    return ds.assign(bsc532_km=bsc_km, dep532_shown=ds["dep532"].where(bsc_km >= dep_bsc_min_km))


def main(path=None, ymax_km=6.0, radar_clearance_m=250.0, suffix="", grid_label="HALO 10 s × 15 m", panels=PANELS_10S):
    ds = xr.open_dataset(path or cache_dir() / f"coreg_20260217_1530-1600{suffix}.nc")
    if "bsc532" in ds:
        ds = add_display_fields(ds)
    t = ds["time"].values
    zkm = ds["z"].values / 1000.0
    liquid = ds["liquid"].values.astype(float)
    sfc = ds["surface_alt_lidar"].to_series().interpolate(limit_direction="both").to_numpy() / 1000.0
    # display-only: hide radar gates within radar_clearance_m of the lidar surface (surface clutter)
    near_sfc = (zkm[None, :] - sfc[:, None]) * 1000.0 < radar_clearance_m

    fig, axs = uplt.subplots(nrows=4, refwidth=7.5, refaspect=3.6, sharex=True, sharey=True,
                             facecolor=SURFACE)
    for ax, (var, title, units, kw, ticks) in zip(axs, panels):
        data = ds[var].values
        if var in ("zgc_w", "dwr_wg_gc"):
            data = np.where(near_sfc, np.nan, data)
        m = ax.pcolormesh(t, zkm, data.T, discrete=False, **kw)
        cb = ax.colorbar(m, loc="r", label=units, width=0.12)
        if ticks is not None:
            cb.set_ticks(ticks)
            if var in ("532_total_attn_bsc", "bsc532_km"):
                cb.set_ticklabels([f"$10^{{{int(np.log10(v))}}}$" for v in ticks])
        ax.plot(t, sfc, color=MUTED, lw=1.0)
        if var == "dwr_wg_gc":
            ax.contour(t, zkm, liquid.T, levels=[0.5], colors=INK, linewidths=1.0)
        ax.format(title=title, titleloc="l", facecolor=SURFACE, gridcolor=GRID)
    axs.format(ylim=(0, ymax_km), ylabel="Altitude (km MSL)", xlabel="UTC, 2026-02-17",
               suptitle=f"NURTURE G-III, HALO + CloudCube co-registered on the {grid_label} grid "
                        "(15:30–16:00 UTC)")
    axs.format(xrotation=0)
    axs[-1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    out = REPO_ROOT / "docs" / "reports" / "figures" / f"stm_20260217_1530-1600{suffix}.png"
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    print(f"wrote {out}")


if __name__ == "__main__":
    if "--0p5s" in sys.argv:
        sys.exit(main(path=cache_dir() / "coreg_cloud_20260217_1530-1600_0p5s.nc", suffix="_0p5s",
                      grid_label="HALO 0.5 s × 15 m", panels=PANELS_0P5S))
    sys.exit(main())
