"""Block diagram of the NURTURE radar-lidar optimal estimation (OE) framework.

Written for an audience new to OE: plain-language box titles, with the OE
term in small italics. Layout borrows the organization of CAPTIVATE Fig. 1
(Mason et al. 2023): inputs on the left, the iterative retrieval inside one
bounding box, uncertainties carried with the measurements, and the
instrument simulators drawn inside the forward-model box. NURTURE 2026 inputs
are solid; future plug-in inputs (NURTURE 2027, FALCON) are dashed.

Run ``python fig2_oe_framework_diagram.py`` to write
``fig2_oe_framework_diagram.svg`` (editable text) and a PNG preview.
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle

mpl.rcParams.update({
    "svg.fonttype": "none",
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 11,
})

INK = "#2F3A45"
MUTED = "#5A6672"
TERM = "#7D8894"
LOOP = "#B97800"
STYLE = {
    "radar": dict(ec="#882255", fc="#F6E6EE"),
    "lidar": dict(ec="#009E73", fc="#E3F4EE"),
    "plugin": dict(ec=MUTED, fc="white", ls=(0, (4, 3))),
    "process": dict(ec=MUTED, fc="#F3F5F7"),
    "loop": dict(ec=LOOP, fc="white"),
    "sub": dict(ec=MUTED, fc="#F3F5F7"),
    "output": dict(ec=INK, fc="#EEF1F4"),
}


def box(ax, x0, y0, w, h, title, body=None, term=None, style="process", lw=1.4,
        title_fs=14.0, body_fs=12.0, body_dy=3.1):
    """Draw a rounded box with a bold title, body text and an OE term.

    Parameters
    ----------
    x0, y0, w, h : float
        Lower-left corner, width and height (data units).
    title : str
        Bold title at the top of the box.
    body : str, optional
        Body text, placed ``body_dy`` below the top edge.
    term : str, optional
        OE term, shown in small gray italics at the bottom.
    style : str
        Key into ``STYLE``.
    """
    s = dict(STYLE[style])
    ls = s.pop("ls", "-")
    ax.add_patch(FancyBboxPatch((x0, y0), w, h, boxstyle="round,pad=0,rounding_size=0.8",
                                lw=lw, ls=ls, zorder=2, **s))
    ax.text(x0 + w / 2, y0 + h - 1.05, title, ha="center", va="top", fontsize=title_fs,
            fontweight="bold", color=INK, zorder=3)
    if body:
        ax.text(x0 + w / 2, y0 + h - body_dy, body, ha="center", va="top",
                fontsize=body_fs, color=INK, linespacing=1.2, zorder=3)
    if term:
        ax.text(x0 + w / 2, y0 + 0.8, term, ha="center", va="bottom", fontsize=10,
                color=TERM, style="italic", zorder=3)


def arrow(ax, p0, p1, color=INK, lw=1.6, ls="-", rad=0.0, ms=14, z=4):
    """Draw a straight or gently curved arrow from ``p0`` to ``p1``."""
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=ms, color=color,
                                 lw=lw, ls=ls, connectionstyle=f"arc3,rad={rad}",
                                 shrinkA=0, shrinkB=0, zorder=z))


def elbow(ax, pts, color=INK, lw=1.6, ms=14, z=4):
    """Draw a polyline whose last segment ends in an arrowhead."""
    xs, ys = zip(*pts)
    ax.plot(xs[:-1], ys[:-1], color=color, lw=lw, zorder=z, solid_joinstyle="miter")
    arrow(ax, pts[-2], pts[-1], color=color, lw=lw, ms=ms, z=z)


def label(ax, x, y, s, fs=9.5, color=MUTED, ha="center", va="center", **kw):
    """Place a small annotation."""
    ax.text(x, y, s, ha=ha, va=va, fontsize=fs, color=color, zorder=5, **kw)


def draw_inputs(ax):
    """Measurements column (NURTURE 2026 solid, plug-ins dashed) and their routes."""
    ax.text(6.0, 54.6, "Measurements", ha="left", va="center", fontsize=16,
            fontweight="bold", color=INK)
    ax.text(5.0, 52.2, "NURTURE 2026 (G-III)", ha="left", va="center",
            fontsize=12, color=MUTED)
    ax.text(7.0, 50.7, "(observation vector)", ha="left", va="center", fontsize=10,
            color=TERM, style="italic")
    box(ax, 1.0, 41.4, 20.5, 7.0, "CloudCube radar", "G, W, Ka reflectivity", style="radar")
    box(ax, 1.0, 31.6, 20.5, 8.4, "HALO DIAL/HSRL",
        "backscatter, depolarization,\nwater vapor above cloud", style="lidar")
    ax.plot([21.5, 23.5], [44.9, 44.9], color=INK, lw=1.4, zorder=1)
    ax.plot([21.5, 23.5], [35.8, 35.8], color=INK, lw=1.4, zorder=1)
    ax.plot([23.5, 23.5], [35.8, 46.0], color=INK, lw=1.4, zorder=1)

    ax.text(1.8, 27.4, "Plug-in inputs, same framework", ha="left", va="center",
            fontsize=12, fontweight="bold", color=MUTED)
    for y0, t, b in ((17.6, "NURTURE 2027 (B-777)", "APR-3: Ku/Ka/W + Doppler\nVIPR: in-cloud water vapor\nCoSMIR/CoSSIR: microwave +\nsub-mm radiances"),
                     (7.8, "FALCON constellation", "radar with Doppler, lidar w/ DIAL,\nmicrowave radiometer")):
        box(ax, 1.0, y0, 20.5, 9.0, t, b, style="plugin")
        ax.plot([21.5, 23.5], [y0 + 3.7, y0 + 3.7], color=MUTED, lw=1.4, ls=(0, (3, 2)),
                zorder=1)
    ax.plot([23.5, 23.5], [12.3, 35.8], color=MUTED, lw=1.4, ls=(0, (3, 2)), zorder=1)

    arrow(ax, (23.5, 46.0), (25.5, 46.0), ms=10)
    elbow(ax, [(23.5, 35.8), (42.6, 35.8), (42.6, 12.6), (47.0, 12.6)])
    label(ax, 33.0, 36.6, "measurements +\nuncertainties", fs=10, color=INK, va="bottom")


def draw_start(ax):
    """Starting-point box: phase identification plus prior knowledge."""
    box(ax, 25.5, 40.0, 15.0, 12.0, "Starting point",
        "cloud phase from lidar\ndepolarization + typical\nproperties for that\n"
        "phase and temperature", "(prior)")
    arrow(ax, (40.5, 46.0), (47.0, 46.0))
    label(ax, 42.9, 46.6, "first\nguess", fs=9.5, va="bottom")


def draw_loop(ax):
    """Iterative retrieval: guess, simulate, compare inside one bounding box."""
    ax.add_patch(FancyBboxPatch((45.3, 5.0), 35.4, 50.6,
                                boxstyle="round,pad=0,rounding_size=1.4", fc="#FCEFD6",
                                alpha=0.55, ec="none", zorder=0))
    ax.add_patch(FancyBboxPatch((45.3, 5.0), 35.4, 50.6,
                                boxstyle="round,pad=0,rounding_size=1.4", fc="none",
                                ec=LOOP, lw=1.6, ls=(0, (5, 3)), zorder=1))
    label(ax, 63.0, 54.0, "Repeat until simulated ≈ measured", fs=12, color=LOOP,
          fontweight="bold")
    label(ax, 63.0, 52.2, "(iterative retrieval)", fs=9.5, color=TERM, style="italic")

    box(ax, 47.0, 41.0, 27.5, 9.4, "Current cloud guess",
        "ice amount and size, liquid water and\nwater vapor in each layer",
        "(state vector)", style="loop")
    arrow(ax, (60.75, 41.0), (60.75, 37.6), color=LOOP, lw=2.2, ms=16)

    box(ax, 47.0, 18.8, 27.5, 18.8, "Simulate the measurements", None, "(forward model)",
        style="loop")
    box(ax, 48.2, 25.4, 11.4, 9.0, "Radar",
        "ice: ARTS (DDA) or\nSSRGA, under\nevaluation\ndroplets: Mie",
        style="sub", lw=1.0, title_fs=10.5, body_fs=9.0, body_dy=2.6)
    box(ax, 61.9, 25.4, 11.4, 9.0, "Lidar",
        "ice: TAMUice2016\ndroplets: Mie\nincl. multiple\nscattering",
        style="sub", lw=1.0, title_fs=10.5, body_fs=9.0, body_dy=2.6)
    ax.add_patch(FancyArrowPatch((59.7, 29.9), (61.8, 29.9), arrowstyle="<|-|>",
                                 mutation_scale=9, color=MUTED, lw=1.1, zorder=4))
    label(ax, 60.75, 24.2, "same particle mass and shape at all wavelengths", fs=9)
    label(ax, 60.75, 22.3, "gas absorption: temp, pressure (MERRA-2) + water vapor",
          fs=9, color=INK)
    arrow(ax, (60.75, 18.8), (60.75, 16.4), color=LOOP, lw=2.2, ms=16)

    box(ax, 47.0, 7.0, 27.5, 9.4, "Compare and adjust",
        "weigh the disagreement by the uncertainties\nand the starting point; update the guess",
        "(cost function, solver)", style="loop")
    elbow(ax, [(74.5, 13.4), (77.6, 13.4), (77.6, 45.7), (74.5, 45.7)], color=LOOP, lw=2.2,
          ms=16)
    label(ax, 78.9, 29.5, "not yet: adjust the guess", fs=9.5, color=LOOP, rotation=90,
          fontweight="bold")


def draw_outputs(ax):
    """Outputs panel: planned retrieved variables, maturity, and diagnostics."""
    x0, y0, w, h = 83.0, 8.6, 16.2, 33.4
    ax.add_patch(FancyBboxPatch((x0, y0), w, h, boxstyle="round,pad=0,rounding_size=0.8",
                                lw=1.6, zorder=2, **STYLE["output"]))
    xt = x0 + 1.0
    ax.text(x0 + w / 2, y0 + h - 1.1, "Planned retrieved\nvariables", ha="center", va="top",
            fontsize=12, fontweight="bold", color=INK, zorder=3, linespacing=1.1)
    ax.text(x0 + w / 2, y0 + h - 5.2, "profiles along the flight track", ha="center",
            va="top", fontsize=9.5, color=MUTED, zorder=3)

    def bullet(y, s, filled=True):
        ax.add_patch(Circle((xt + 0.5, y), 0.42, fc=INK if filled else "white", ec=INK,
                            lw=1.1, zorder=3))
        ax.text(xt + 1.5, y, s, ha="left", va="center", fontsize=10.5, color=INK, zorder=3)

    bullet(33.6, "ice water content")
    bullet(31.4, "ice particle size")
    bullet(29.2, "cloud phase")
    bullet(26.4, "liquid water", filled=False)
    bullet(24.2, "riming (where Ka signal\nexists)", filled=False)
    ax.plot([x0 + 1.0, x0 + w - 1.0], [21.4, 21.4], color="#B9C2CB", lw=0.9, zorder=3)
    ax.text(xt, 20.4, "With every profile:", ha="left", va="top", fontsize=10.5,
            fontweight="bold", color=INK, zorder=3)
    ax.text(xt, 18.3, "• uncertainty\n• how much each\n   instrument contributed",
            ha="left", va="top", fontsize=10.5, color=INK, linespacing=1.25, zorder=3)
    ax.add_patch(Circle((xt + 0.5, 10.4), 0.42, fc=INK, ec=INK, lw=1.1, zorder=3))
    ax.text(xt + 1.3, 10.4, "established", ha="left", va="center", fontsize=9, color=MUTED,
            zorder=3)
    ax.add_patch(Circle((xt + 7.5, 10.4), 0.42, fc="white", ec=INK, lw=1.1, zorder=3))
    ax.text(xt + 8.3, 10.4, "developmental", ha="left", va="center", fontsize=9,
            color=MUTED, zorder=3)
    ax.text(x0 + w / 2, 12.6, "(posterior, averaging kernel)", ha="center", va="center",
            fontsize=9.5, color=TERM, style="italic", zorder=3)

    elbow(ax, [(74.5, 9.6), (79.8, 9.6), (79.8, 9.6), (83.0, 9.6)])
    label(ax, 77.4, 8.8, "agree", fs=9.5, color=INK, va="top")
    arrow(ax, (91.1, 8.6), (91.1, 5.4), color=MUTED)
    label(ax, 91.1, 4.8, "latent and radiative heating\nstudies (SQ 2.1)", fs=9.5, va="top")


def main(out="figures/fig2_oe_framework_diagram"):
    """Render the diagram to SVG and PNG.

    Parameters
    ----------
    out : str
        Output file stem.
    """
    fig = plt.figure(figsize=(13.333, 7.5), facecolor="white")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 56.25)
    ax.set_aspect("equal")
    ax.axis("off")
    draw_inputs(ax)
    draw_start(ax)
    draw_loop(ax)
    draw_outputs(ax)
    fig.savefig(f"{out}.svg")
    fig.savefig(f"{out}.png", dpi=150)


if __name__ == "__main__":
    main()
