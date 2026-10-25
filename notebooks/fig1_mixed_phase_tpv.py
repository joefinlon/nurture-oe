"""Conceptual schematic of a liquid-topped mixed-phase cloud near a TPV.

Modeled on the structure sampled by HALO (lidar) and CloudCube (radar) on
NURTURE RF17 (17 Feb 2026). Schematic only: no altitudes, no measured values.

Run ``python fig1_mixed_phase_tpv.py`` to write ``fig1_mixed_phase_tpv.svg``
(editable text) and a PNG preview.
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Circle, FancyArrowPatch, Rectangle

mpl.rcParams.update({
    "svg.fonttype": "none",
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 12,
})

COL = {
    "lidar": "#009E73",
    "both": "#E69F00",
    "radar": "#CC79A7",
    "heat": "#D55E00",
    "cool": "#0072B2",
    "ink": "#2F3A45",
    "streak": "#8A96A3",
    "surface": "#DCE3EA",
    "radar_line": "#882255",
    "ka": "#4B2C7A",
}
ALPHA = {"lidar": 0.30, "both": 0.45, "radar": 0.30}

X0, X1 = 24.0, 93.0
LIQ_PATCHES = [(24.0, 41.0, "dense"), (43.0, 47.5, "thin"), (50.0, 52.5, "thin"),
               (55.0, 60.0, "thin"), (62.5, 65.0, "thin"), (67.5, 70.5, "dense"),
               (73.0, 78.5, "thin")]
LIQ_DEPTH = 3.6


def _smooth(y, w):
    """Running mean with edge padding.

    Parameters
    ----------
    y : ndarray
        Profile to smooth.
    w : int
        Odd window length in samples.

    Returns
    -------
    ndarray
        Smoothed profile, same length as ``y``.
    """
    pad = w // 2
    return np.convolve(np.pad(y, pad, mode="edge"), np.ones(w) / w, mode="valid")


def deck_geometry(x):
    """Build the along-track cloud and instrument-coverage profiles.

    Parameters
    ----------
    x : ndarray
        Along-track coordinate (schematic units).

    Returns
    -------
    dict
        ``top0`` liquid-top height, ``liquid`` mask, ``ct`` cloud top,
        ``base`` hydrometeor lower boundary, ``lidar_loss`` height where the
        lidar signal is fully attenuated, ``radar_top`` radar echo top,
        ``ka_top`` Ka-band echo top (NaN where Ka does not detect).

    Notes
    -----
    Dense liquid stops the lidar near cloud top. Thin (low droplet number)
    liquid lets it through into the ice below, and gaps let it reach deep.
    """
    top0 = 44.0 + 0.6 * np.sin(2 * np.pi * (x - X0) / 16.0)
    top0 += np.interp(x, [X0, 80, 86, X1], [0.0, 0.0, -2.0, -5.5])
    dense = np.zeros_like(x, dtype=bool)
    thin = np.zeros_like(x, dtype=bool)
    for a, b, kind in LIQ_PATCHES:
        m = (x >= a) & (x <= b)
        if kind == "dense":
            dense |= m
        else:
            thin |= m
    liquid = dense | thin
    ct = _smooth(np.where(liquid, top0, top0 - 1.6), 25)

    bx = [X0, 31, 34, 46, 50, 55, 58, 68, 72, 80, 84, 88, X1]
    by = [13, 12.5, 4, 3.5, 9, 9, 4.5, 5, 12, 16, 30, 34, 40]
    base = _smooth(np.interp(x, bx, by), 61)
    base = np.minimum(base, ct)

    gap_depth = np.interp(x, [42.0, 48.75, 53.75, 61.25, 66.25, 71.75, 79.5],
                          [15.0, 12.5, 14.0, 11.5, 13.0, 12.0, 16.0])
    thin_depth = np.interp(x, [45.25, 57.5, 63.75, 75.75], [29.0, 26.0, 28.0, 27.0])
    ll = _smooth(np.where(thin, thin_depth, gap_depth), 121)
    ll = np.where(dense, top0 - 3.0, ll)
    ll = np.where(x > 80, -10.0, ll)
    lidar_loss = np.clip(_smooth(ll, 15), base, ct)

    rt = np.select([dense, thin], [top0 - 1.8, top0 - 2.0], ct - 1.0)
    rt = np.where(x > 82, -10.0, rt)
    radar_top = np.clip(_smooth(rt, 21), base, ct)

    kx = [X0, 30, 35, 44, 49, 56, 60.5, 66, 71, X1]
    ky = [-1, -1, 26, 27, -1, -1, 24.5, 25.5, -1, -1]
    ka = _smooth(np.interp(x, kx, ky), 81)
    ka_top = np.where(ka > base + 0.6, ka, np.nan)
    return dict(top0=top0, liquid=liquid, ct=ct, base=base,
                lidar_loss=lidar_loss, radar_top=radar_top, ka_top=ka_top)


def hexagon(ax, x, y, r, rot=0.0, z=4, **kw):
    """Draw a hexagonal plate outline."""
    t = np.deg2rad(rot) + np.arange(6) * np.pi / 3
    kw.setdefault("lw", 1.1)
    ax.add_patch(Polygon(np.c_[x + r * np.cos(t), y + r * np.sin(t)], closed=True,
                         fc="white", ec=COL["ink"], zorder=z, **kw))


def dendrite(ax, x, y, r, rot=0.0, z=4, alpha=1.0, lw=1.1):
    """Draw a six-armed stellar dendrite."""
    for k in range(6):
        t = np.deg2rad(rot) + k * np.pi / 3
        ax.plot([x, x + r * np.cos(t)], [y, y + r * np.sin(t)], color=COL["ink"],
                lw=lw, alpha=alpha, zorder=z, solid_capstyle="round")
        for f in (0.45, 0.72):
            bx, by = x + f * r * np.cos(t), y + f * r * np.sin(t)
            for s in (-1, 1):
                tb = t + s * np.pi / 3
                ax.plot([bx, bx + 0.28 * r * np.cos(tb)], [by, by + 0.28 * r * np.sin(tb)],
                        color=COL["ink"], lw=lw * 0.8, alpha=alpha, zorder=z,
                        solid_capstyle="round")


def rimed(ax, x, y, r, rng, rot=0.0, z=4):
    """Draw a hexagonal plate with frozen droplets on its surface."""
    hexagon(ax, x, y, r, rot, z=z)
    for _ in range(9):
        a, d = rng.uniform(0, 2 * np.pi), r * rng.uniform(0.35, 1.0)
        ax.add_patch(Circle((x + d * np.cos(a), y + d * np.sin(a)), 0.17 * r,
                            fc=COL["ink"], ec="none", zorder=z + 0.1))


def aggregate(ax, x, y, r, rng, z=4, alpha=1.0, lw=1.0):
    """Draw a snowflake aggregate as overlapping dendrites."""
    for dx, dy in ((0, 0), (0.75, 0.45), (-0.7, 0.4), (0.1, -0.75)):
        dendrite(ax, x + dx * r, y + dy * r, 0.62 * r, rot=rng.uniform(0, 60),
                 z=z, alpha=alpha, lw=lw)


def droplet(ax, x, y, r=0.22, z=4):
    """Draw a supercooled droplet."""
    ax.add_patch(Circle((x, y), r, fc=COL["ink"], ec="none", zorder=z))


def heat_marker(ax, x, y, sign, size, z=8):
    """Draw a heating (+) or cooling (-) badge.

    Parameters
    ----------
    sign : {'+', '-'}
        Heating or cooling.
    size : float
        Badge radius; encodes relative magnitude (qualitative).
    """
    c = COL["heat"] if sign == "+" else COL["cool"]
    ax.add_patch(Circle((x, y), size, fc=c, ec="white", lw=1.8, zorder=z))
    ax.text(x, y - 0.05 * size, "+" if sign == "+" else "−", color="white",
            ha="center", va="center", fontsize=11 * size, fontweight="bold", zorder=z + 1)


def wavy_arrow(ax, x, y0, y1, color, z=6):
    """Draw an upward wavy arrow (longwave emission)."""
    yy = np.linspace(y0, y1 - 0.8, 60)
    ax.plot(x + 0.35 * np.sin((yy - y0) * 2.4), yy, color=color, lw=1.5, zorder=z)
    ax.add_patch(FancyArrowPatch((x, y1 - 0.9), (x, y1), arrowstyle="-|>",
                                 mutation_scale=11, color=color, lw=1.5, zorder=z))


def aircraft(ax, x, y, s=1.0, z=7):
    """Draw a side-view jet silhouette pointing right."""
    parts = [
        [(-4.0, 0.2), (-3.4, 0.55), (3.0, 0.55), (4.0, 0.15), (3.6, -0.25), (-3.6, -0.35)],
        [(-3.9, 0.4), (-4.5, 2.1), (-3.6, 2.1), (-2.7, 0.5)],
        [(-0.4, 0.0), (1.3, 0.0), (-1.3, -1.3), (-2.0, -1.3)],
        [(-3.0, 0.7), (-1.6, 0.7), (-1.6, 1.15), (-3.0, 1.15)],
    ]
    for p in parts:
        ax.add_patch(Polygon([(x + s * px, y + s * py) for px, py in p], closed=True,
                             fc=COL["ink"], ec="none", zorder=z))


def draw_scene(ax, rng):
    """Draw the cloud cross-section, coverage zones, microphysics and heating."""
    x = np.linspace(X0, X1, 1381)
    g = deck_geometry(x)
    ct, base, ll, rt = g["ct"], g["base"], g["lidar_loss"], g["radar_top"]
    top0 = g["top0"]

    ax.add_patch(Rectangle((0, 0), 100, 2.5, fc=COL["surface"], ec="none", zorder=1))
    ax.text(99, 1.25, "snow / sea-ice surface", ha="right", va="center", fontsize=10,
            color="#5A6672", zorder=2)

    lid_lo = np.maximum(rt, ll)
    ax.fill_between(x, lid_lo, ct, where=ct > lid_lo + 1e-6, fc=COL["lidar"],
                    alpha=ALPHA["lidar"], lw=0, zorder=1)
    ax.fill_between(x, ll, rt, where=rt > ll + 1e-6, fc=COL["both"],
                    alpha=ALPHA["both"], lw=0, zorder=1)
    rad_hi = np.minimum(ll, rt)
    ax.fill_between(x, base, rad_hi, where=rad_hi > base + 1e-6, fc=COL["radar"],
                    alpha=ALPHA["radar"], lw=0, zorder=1)
    ax.plot(x, ct, color="#7D8894", lw=0.9, zorder=2)
    ax.plot(x, base, color="#7D8894", lw=0.9, zorder=2)

    for x0 in (35, 40, 45, 49, 62, 67):
        i0 = np.argmin(abs(x - x0))
        x1 = x0 - 5.0
        i1 = np.argmin(abs(x - x1))
        p0, p2 = np.array([x0, top0[i0] - 4.5]), np.array([x1, base[i1] + 0.6])
        p1 = np.array([x0 - 0.5, (p0[1] + p2[1]) / 2])
        t = np.linspace(0, 1, 50)[:, None]
        c = (1 - t) ** 2 * p0 + 2 * (1 - t) * t * p1 + t ** 2 * p2
        ax.plot(c[:, 0], c[:, 1], color=COL["streak"], lw=1.3, alpha=0.45, zorder=2)

    for a, b, kind in LIQ_PATCHES:
        m = (x >= a) & (x <= b)
        ax.fill_between(x[m], top0[m] - LIQ_DEPTH, ct[m], fc="none", ec=COL["ink"],
                        lw=0.9, ls="-" if kind == "dense" else (0, (3, 2)), zorder=3)
        n = int((b - a) * (1.6 if kind == "dense" else 0.55))
        xs = rng.uniform(a + 0.5, b - 0.5, n)
        for xd in xs:
            yt = np.interp(xd, x, top0)
            droplet(ax, xd, rng.uniform(yt - LIQ_DEPTH + 0.5, yt - 0.5))

    for xr in (27.5, 35.5, 51.2, 57.5):
        rimed(ax, xr, np.interp(xr, x, top0) - 1.85, 1.05, rng, rot=rng.uniform(0, 60))

    for xp, yp, kind in ((26.0, 36.0, "d"), (29.0, 38.2, "h"), (35.8, 33.4, "d")):
        if kind == "h":
            hexagon(ax, xp, yp, 0.85, rot=rng.uniform(0, 60))
        else:
            dendrite(ax, xp, yp, 1.05, rot=rng.uniform(0, 60))
    for xp, kind in ((48.7, "d"), (61.2, "d"), (66.2, "h"), (71.7, "d")):
        yp = np.interp(xp, x, ct) - 3.2
        if kind == "h":
            hexagon(ax, xp, yp, 0.85, rot=rng.uniform(0, 60))
        else:
            dendrite(ax, xp, yp, 0.95, rot=rng.uniform(0, 60))
    for xp, yp in ((83.5, 37.5), (87.5, 38.0), (85.5, 35.2)):
        hexagon(ax, xp, yp, 0.7, rot=rng.uniform(0, 60))

    for xa, ya in ((26.5, 22.0), (32.5, 25.0), (43.5, 20.0), (48.0, 25.5), (53.0, 19.0),
                   (59.5, 23.0), (66.0, 18.5), (72.0, 21.0), (78.0, 27.5)):
        aggregate(ax, xa, ya, 1.45, rng)
    for xa, ya in ((35.0, 10.0), (42.0, 7.5), (47.5, 12.0), (61.0, 9.5), (66.5, 7.0)):
        aggregate(ax, xa, ya, 2.0, rng)
    for xa, ya in ((25.3, 9.0), (30.6, 8.6), (75.0, 9.5), (78.5, 6.5)):
        aggregate(ax, xa, ya, 0.95, rng, alpha=0.45, lw=0.9)

    ka = g["ka_top"]
    ax.plot(x, ka, color=COL["ka"], lw=2.0, ls=(0, (1.2, 1.6)), zorder=5)
    ax.text(46.0, 28.2, "Ka-band echo top", color=COL["ka"], fontsize=10.5,
            fontweight="bold", ha="center", va="bottom", zorder=7,
            bbox=dict(fc="white", ec="none", alpha=0.8, pad=1.2))

    for xw in (25.8, 36.6):
        yt = np.interp(xw, x, top0)
        wavy_arrow(ax, xw, yt + 0.6, yt + 5.0, COL["cool"])
    heat_marker(ax, 29.6, 48.4, "-", 2.0)
    heat_marker(ax, 64.0, 46.6, "-", 1.05)
    heat_marker(ax, 32.3, 36.4, "+", 2.0)
    heat_marker(ax, 61.6, 31.0, "+", 1.3)
    heat_marker(ax, 28.0, 5.4, "-", 1.5)
    heat_marker(ax, 72.0, 5.5, "-", 1.1)

    xs0 = 39.2
    ys0 = np.interp(xs0, x, top0) - 2.6
    # ax.add_patch(FancyArrowPatch((xs0, ys0), (39.6, 14.5), connectionstyle="arc3,rad=0.12",
#                                  arrowstyle="-|>", mutation_scale=18, color=COL["ink"],
#                                  lw=2.0, zorder=9))
    #ax.add_patch(Circle((xs0, ys0), 0.5, fc=COL["ink"], ec="white", lw=0.8, zorder=10))
    #ax.plot([43.6, xs0], [52.0, ys0], color="#5A6672", lw=0.9, zorder=9)

    aircraft(ax, 33.0, 63.5, s=1.05)
    ax.plot([33.0, 33.0], [62.6, np.interp(33.0, x, ct) + 0.2], color=COL["lidar"], lw=2.2,
            zorder=6)
    ax.plot([34.0, 34.0], [62.6, np.interp(34.0, x, ct) + 0.2], color=COL["radar_line"],
            lw=2.2, zorder=6)
    ax.add_patch(FancyArrowPatch((39.0, 63.6), (52.0, 63.6), arrowstyle="-|>",
                                 mutation_scale=14, color=COL["ink"], lw=1.3,
                                 ls="--", zorder=6))
    ax.text(53.0, 64.3, "NASA G-III above cloud", ha="left", va="center", fontsize=11.5,
            fontweight="bold", color=COL["ink"])
    ax.text(53.0, 61.6, "HALO lidar and CloudCube radar",
            ha="left", va="center", fontsize=10.5, color=COL["ink"])

    ax.add_patch(FancyArrowPatch((84.0, 51.0), (99.0, 51.0), arrowstyle="-|>",
                                 mutation_scale=16, color="#5A6672", lw=1.6, zorder=6))
    ax.text(99.0, 52.6, "toward TPV core\n(drier air, fewer clouds)", ha="right",
            va="bottom", fontsize=10.5, color="#5A6672")
    ax.text(88.0, 41.0, "thin ice", ha="center", va="bottom", fontsize=10, color=COL["ink"],
            zorder=7)
    ax.text(71.5, 46.0, "broken liquid", ha="center", va="bottom", fontsize=10.5,
            color=COL["ink"], zorder=7)

    labels = [
        (48.5, "Supercooled liquid\nat cloud top", (24.6, 43.6)),
        (41.0, "Riming", (27.5, 42.75)),
        (33.0, "Vapor deposition", (26.0, 36.0)),
        (21.5, "Aggregation", (26.5, 22.0)),
        (8.0, "Sublimation", (25.3, 9.0)),
    ]
    for ly, txt, (tx, ty) in labels:
        ax.text(20.0, ly, txt, ha="right", va="center", fontsize=12, color=COL["ink"],
                linespacing=1.15)
        ax.plot([20.6, tx], [ly, ty], color="#5A6672", lw=0.8, zorder=6)
        ax.add_patch(Circle((tx, ty), 0.35, fc="#5A6672", ec="none", zorder=6))

    ax.text(0.6, 70.8, "Schematic based on RF-17 scene (17 Feb 2026). Not to scale.",
            ha="left", va="center", fontsize=12, color="#7D8894", style="italic")


def draw_key(ax, rng):
    """Draw the legend panel: instruments, coverage zones, hydrometeors, heating."""
    ink = COL["ink"]

    def head(y, s):
        ax.text(0.5, y, s, ha="left", va="center", fontsize=12, fontweight="bold", color=ink)

    def item(y, s, x=6.5, fs=10.5):
        ax.text(x, y, s, ha="left", va="center", fontsize=fs, color=ink, linespacing=1.15)

    head(69.5, "Instruments")
    ax.plot([0.8, 4.8], [66.0, 66.0], color=COL["lidar"], lw=3)
    item(66.0, "HALO: droplets + ice,\nmost sensitive to small particles")
    ax.plot([0.8, 4.8], [61.0, 61.0], color=COL["radar_line"], lw=3)
    item(61.0, "CloudCube (W, G, Ka):\nlarger ice and snow, weak\necho from small droplets")

    head(54.5, "What each instrument sees")
    for y, k, s in ((51.0, "lidar", "lidar only"), (48.0, "both", "lidar + radar"),
                    (45.0, "radar", "radar only")):
        ax.add_patch(Rectangle((0.8, y - 1.0), 4.0, 2.0, fc=COL[k], alpha=ALPHA[k], ec="none"))
        item(y, s)
    ax.plot([0.8, 4.8], [42.0, 42.0], color=COL["ka"], lw=2.0, ls=(0, (1.2, 1.6)))
    item(42.0, "Ka band: mostly larger snow")

    head(36.5, "Hydrometeors")
    droplet(ax, 1.8, 33.3, 0.3)
    droplet(ax, 3.0, 32.6, 0.3)
    droplet(ax, 3.7, 33.6, 0.3)
    item(33.1, "supercooled droplets")
    hexagon(ax, 1.9, 29.6, 0.85, rot=10)
    dendrite(ax, 4.2, 29.6, 1.05, rot=10)
    item(29.6, "ice crystals (plates, dendrites)")
    rimed(ax, 2.8, 26.0, 1.1, rng, rot=10)
    item(26.0, "rimed crystals")
    aggregate(ax, 2.8, 22.0, 1.45, rng)
    item(22.0, "snowflakes (aggregates)")

    head(16.5, "Heating and cooling")
    heat_marker(ax, 2.8, 13.0, "+", 1.15)
    item(13.0, "latent heating\n(condensation, deposition, riming)", fs=10)
    heat_marker(ax, 2.8, 8.6, "-", 1.15)
    item(8.6, "cooling: longwave emission at\ncloud top, sublimation below", fs=10)
    item(4.2, "Marker size = relative strength\n(qualitative)", x=0.8, fs=9.5)


def main(out="figures/fig1_mixed_phase_tpv"):
    """Render the schematic to SVG and PNG.

    Parameters
    ----------
    out : str
        Output file stem.
    """
    rng = np.random.default_rng(17)
    fig = plt.figure(figsize=(13.333, 7.5), facecolor="white")
    ax = fig.add_axes([0.005, 0.02, 0.745, 0.96])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 72.5)
    ax.set_aspect("equal")
    ax.axis("off")
    draw_scene(ax, rng)

    kx = fig.add_axes([0.765, 0.02, 0.23, 0.96])
    kx.set_xlim(0, 30.6)
    kx.set_ylim(0, 72.0)
    kx.set_aspect("equal")
    kx.axis("off")
    draw_key(kx, rng)

    fig.savefig(f"{out}.svg")
    fig.savefig(f"{out}.png", dpi=192)


if __name__ == "__main__":
    main()
