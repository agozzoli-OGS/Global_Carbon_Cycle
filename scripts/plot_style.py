"""
================================================================================
 PAPER-GRADE PLOT STYLE TEMPLATE
 Domain: Oceanography · Biogeochemistry · Ecology · Geophysics
 Target: Journal papers + presentations
 Stack:  Python · Matplotlib · Seaborn · cmocean · (optional) Crameri scico
================================================================================

USAGE
-----
At the top of any plotting script:

    from plot_style import *          # imports everything below
    apply_style()                     # sets rcParams globally

Then use the helper functions and colormap aliases directly.

DEPENDENCIES
------------
    pip install matplotlib seaborn cmocean
    pip install scicomap            # optional, for Crameri maps
================================================================================
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import seaborn as sns

try:
    import cmocean
    import cmocean.cm as cmo
    HAS_CMOCEAN = True
except ImportError:
    HAS_CMOCEAN = False
    print("cmocean not found – cmocean colormaps unavailable.")

try:
    import scicomap as scm
    HAS_SCICO = False   # flag; set True if you install scicomap
except ImportError:
    HAS_SCICO = False


# ══════════════════════════════════════════════════════════════════════════════
# 1. GLOBAL rcPARAMS
# ══════════════════════════════════════════════════════════════════════════════

def apply_style(fontsize=12, spine_width=1.5, tick_width=1.2, tick_length=5):
    """
    Apply the base rcParams style globally.

    Parameters
    ----------
    fontsize    : base font size (12 for papers/presentations)
    spine_width : thickness of axis spines (1.5 = thick, crisp)
    tick_width  : tick mark linewidth
    tick_length : tick mark length in points
    """
    sns.set_theme(
        style="ticks",          # minimal: ticks + spines only, no box
        context="paper",        # scales everything; override below
        font_scale=fontsize / 10,
        rc={
            # ── Spines ──────────────────────────────────────────────────────
            "axes.spines.top":    True,
            "axes.spines.right":  True,
            "axes.linewidth":     spine_width,

            # ── Ticks ───────────────────────────────────────────────────────
            "xtick.direction":    "inout",   # bidirectional
            "ytick.direction":    "inout",
            "xtick.major.width":  tick_width,
            "ytick.major.width":  tick_width,
            "xtick.minor.width":  tick_width * 0.7,
            "ytick.minor.width":  tick_width * 0.7,
            "xtick.major.size":   tick_length,
            "ytick.major.size":   tick_length,
            "xtick.minor.size":   tick_length * 0.6,
            "ytick.minor.size":   tick_length * 0.6,
            "xtick.top":          True,      # ticks on all 4 sides
            "ytick.right":        True,

            # ── Font ────────────────────────────────────────────────────────
            "font.size":          fontsize,
            "axes.labelsize":     fontsize,
            "xtick.labelsize":    fontsize - 1,
            "ytick.labelsize":    fontsize - 1,
            "legend.fontsize":    fontsize - 1,
            "axes.titlesize":     fontsize,

            # ── Lines ───────────────────────────────────────────────────────
            "lines.linewidth":    2.0,       # thick lines throughout
            "lines.markersize":   6,

            # ── Grid ────────────────────────────────────────────────────────
            "axes.grid":          True,
            "grid.color":         "#CCCCCC",
            "grid.linestyle":     "--",
            "grid.linewidth":     0.5,
            "grid.alpha":         0.6,
            "axes.axisbelow":     True,      # grid behind data

            # ── Figure ──────────────────────────────────────────────────────
            "figure.dpi":         150,
            "savefig.dpi":        300,
            "savefig.bbox":       "tight",
            "savefig.pad_inches": 0.05,
            "figure.facecolor":   "white",
            "axes.facecolor":     "white",

            # ── Legend ──────────────────────────────────────────────────────
            "legend.frameon":     True,
            "legend.framealpha":  0.85,
            "legend.edgecolor":   "#AAAAAA",
        }
    )


# ══════════════════════════════════════════════════════════════════════════════
# 2. COLOR PALETTES  (qualitative – for line plots / categories)
# ══════════════════════════════════════════════════════════════════════════════

# Seaborn "deep" – your default, moderated hues, good contrast
PALETTE_DEEP    = sns.color_palette("deep")

# Seaborn "colorblind" – like deep but CVD-safe
PALETTE_CB      = sns.color_palette("colorblind")

# Nature Methods palette (Wong 2010) – 8 CVD-safe colors, publication staple
PALETTE_WONG    = [
    "#000000",  # black
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#009E73",  # bluish green
    "#F0E442",  # yellow  (use with care on white)
    "#0072B2",  # blue
    "#D55E00",  # vermilion
    "#CC79A7",  # reddish purple
]

# Active palette (change this to switch defaults)
PALETTE = PALETTE_DEEP


def set_line_palette(palette=None, n_colors=None):
    """Set the active line-plot color cycle. Pass a list or palette name."""
    if palette is None:
        palette = PALETTE
    if n_colors:
        palette = sns.color_palette(palette, n_colors)
    sns.set_palette(palette)


# ══════════════════════════════════════════════════════════════════════════════
# 3. COLORMAP ALIASES  (for 2D fields)
# ══════════════════════════════════════════════════════════════════════════════
#
# Convention (matches cmocean philosophy):
#   cmap_*         → sequential
#   dcmap_*        → diverging  (always center at 0 / midpoint)
#   ccmap_*        → cyclic     (phase, angle)
#
# ─────────────────────────────────────────────────────────────
# SEQUENTIAL – physical oceanography
# ─────────────────────────────────────────────────────────────
if HAS_CMOCEAN:
    cmap_temp     = cmo.thermal     # temperature (dark→warm yellow)
    cmap_sal      = cmo.haline      # salinity    (blue→yellow-green)
    cmap_dens     = cmo.dense       # density     (white→purple)
    cmap_oxy      = cmo.oxy         # oxygen      (sequential + diverging tail)
    cmap_depth    = cmo.deep        # bathymetry  (yellow→dark blue-purple)
    cmap_speed    = cmo.speed       # current speed (greenish-yellow)
    cmap_chl      = cmo.algae       # chlorophyll (white-green → dark green)
    cmap_turb     = cmo.turbid      # turbidity / TSM (light→dark brown)
    cmap_matter   = cmo.matter      # CDOM / DOM  (yellow → pink)
    cmap_rain     = cmo.rain        # precipitation
    cmap_solar    = cmo.solar       # irradiance  (dark brown → bright yellow)
    cmap_ice      = cmo.ice         # sea ice     (dark blue → white)
    cmap_amp      = cmo.amp         # amplitude / wave height (white → red)
    cmap_topo     = cmo.topo        # bathy+topo combined
    # ─────────────────────────────────────────────────────────
    # DIVERGING
    # ─────────────────────────────────────────────────────────
    dcmap_bal     = cmo.balance     # SSH / sea level anomaly (blue–white–red)
    dcmap_delta   = cmo.delta       # velocity anomaly        (blue–white–green)
    dcmap_curl    = cmo.curl        # vorticity               (teal–white–magenta)
    dcmap_diff    = cmo.diff        # generic diff            (blue–white–brown)
    dcmap_tarn    = cmo.tarn        # rain anomaly            (brown–white–green)
    # ─────────────────────────────────────────────────────────
    # CYCLIC
    # ─────────────────────────────────────────────────────────
    ccmap_phase   = cmo.phase       # tidal / wave phase

else:
    # Graceful fallback to matplotlib perceptually-uniform maps
    cmap_temp   = "inferno"
    cmap_sal    = "YlGnBu"
    cmap_dens   = "Purples"
    cmap_oxy    = "RdYlBu_r"
    cmap_depth  = "Blues_r"
    cmap_speed  = "YlGn"
    cmap_chl    = "Greens"
    dcmap_bal   = "RdBu_r"
    dcmap_delta = "BrBG"
    dcmap_curl  = "PiYG"
    dcmap_diff  = "RdBu_r"
    ccmap_phase = "hsv"


# ══════════════════════════════════════════════════════════════════════════════
# 4. FIGURE SIZE HELPERS
# ══════════════════════════════════════════════════════════════════════════════

# Common journal widths in inches
FIG_SINGLE  = (3.5, 2.8)    # single column
FIG_DOUBLE  = (7.0, 4.0)    # double column
FIG_FULL    = (10.0, 5.5)   # full page width
FIG_SQUARE  = (4.5, 4.5)    # square panel
FIG_TALL    = (7.0, 8.0)    # tall multi-panel

def figsize(cols=1, rows=1, aspect=1.4, col_width=3.5):
    """
    Compute a figure size for an (rows × cols) panel layout.

    Parameters
    ----------
    cols      : number of subplot columns
    rows      : number of subplot rows
    aspect    : width/height ratio per panel
    col_width : base width per column in inches

    Returns
    -------
    (width, height) tuple
    """
    w = cols * col_width
    h = rows * col_width / aspect
    return (w, h)


# ══════════════════════════════════════════════════════════════════════════════
# 5. contourf WRAPPER
# ══════════════════════════════════════════════════════════════════════════════

def pcontourf(ax, x, y, z, cmap, levels=20, vmin=None, vmax=None,
              extend="both", **kwargs):
    """
    Thin wrapper around ax.contourf with paper-grade defaults.

    Parameters
    ----------
    ax      : matplotlib Axes
    x, y    : coordinate arrays (1D or 2D)
    z       : data array (2D)
    cmap    : colormap (use aliases above)
    levels  : int (number of levels) or array of explicit levels
    vmin/vmax : colormap limits (required for diverging maps)
    extend  : "both" | "min" | "max" | "neither"

    Returns
    -------
    cf : the QuadContourSet (for colorbar)
    """
    if isinstance(levels, int):
        if vmin is not None and vmax is not None:
            levels = np.linspace(vmin, vmax, levels)
        # else let matplotlib auto-select
    cf = ax.contourf(x, y, z, levels=levels, cmap=cmap,
                     vmin=vmin, vmax=vmax, extend=extend, **kwargs)
    return cf


def add_colorbar(fig, cf, ax=None, label="", orientation="horizontal",
                 shrink=0.85, pad=0.08, aspect=30, **kwargs):
    """
    Add a formatted colorbar below (horizontal, default) or to the right.

    Parameters
    ----------
    fig         : Figure
    cf          : mappable (output of contourf / pcolormesh)
    ax          : Axes or list of Axes to steal space from
    label       : colorbar label string
    orientation : "horizontal" (bottom) or "vertical" (right)
    shrink      : fraction of axes size
    pad         : gap between axes and colorbar
    aspect      : aspect ratio of colorbar

    Returns
    -------
    cbar : Colorbar object
    """
    cbar = fig.colorbar(cf, ax=ax, orientation=orientation,
                        shrink=shrink, pad=pad, aspect=aspect, **kwargs)
    cbar.set_label(label, fontsize=mpl.rcParams["axes.labelsize"])
    cbar.ax.tick_params(labelsize=mpl.rcParams["xtick.labelsize"],
                        width=1.0, length=3)
    return cbar


# ══════════════════════════════════════════════════════════════════════════════
# 6. MULTI-PANEL HELPERS
# ══════════════════════════════════════════════════════════════════════════════

PANEL_LABELS = list("abcdefghijklmnopqrstuvwxyz")

def label_panels(axes, labels=None, x=0.03, y=0.96, bold=True,
                 fontsize=None, va="top", ha="left"):
    """
    Add bold panel labels (a, b, c …) inside each axes, top-left corner.

    Parameters
    ----------
    axes   : list or array of Axes
    labels : custom labels; defaults to lowercase a, b, c …
    x, y   : position in axes fraction
    bold   : whether to use fontweight="bold"
    """
    if fontsize is None:
        fontsize = mpl.rcParams["axes.labelsize"]
    if labels is None:
        labels = PANEL_LABELS
    axes_flat = np.array(axes).flatten()
    for ax, lbl in zip(axes_flat, labels):
        ax.text(x, y, lbl, transform=ax.transAxes,
                fontsize=fontsize,
                fontweight="bold" if bold else "normal",
                va=va, ha=ha)


def shared_colorbar(fig, axes_row, cf, label="", orientation="horizontal",
                    pad=0.08, shrink=0.9, aspect=35):
    """
    Add one colorbar shared across a row (or list) of axes.
    Ideal when all panels share the same colormap and limits.
    """
    return add_colorbar(fig, cf, ax=list(axes_row), label=label,
                        orientation=orientation, pad=pad,
                        shrink=shrink, aspect=aspect)


# ══════════════════════════════════════════════════════════════════════════════
# 7. AXIS UTILITIES
# ══════════════════════════════════════════════════════════════════════════════

def set_minor_ticks(ax, x_every=None, y_every=None):
    """Enable minor ticks. Pass spacing or use AutoMinorLocator."""
    if x_every:
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(x_every))
    else:
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    if y_every:
        ax.yaxis.set_minor_locator(ticker.MultipleLocator(y_every))
    else:
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())


def clean_axis(ax, grid=True):
    """
    Apply spine/tick cleanup: remove top-right spines if desired,
    ensure grid is behind data. Call after plotting.
    """
    if grid:
        ax.set_axisbelow(True)
        ax.grid(True, color="#CCCCCC", linestyle="--",
                linewidth=0.5, alpha=0.6)
    # Ensure spine width matches rcParams
    for spine in ax.spines.values():
        spine.set_linewidth(mpl.rcParams["axes.linewidth"])


def diverging_levels(vmax, n=21):
    """
    Generate symmetric levels for a diverging colormap centered at 0.

    Parameters
    ----------
    vmax : maximum absolute value
    n    : number of levels (odd number recommended)

    Returns
    -------
    levels : np.ndarray symmetric around 0
    """
    return np.linspace(-abs(vmax), abs(vmax), n)


# ══════════════════════════════════════════════════════════════════════════════
# 8. QUIVER / VECTOR OVERLAY
# ══════════════════════════════════════════════════════════════════════════════

def add_quiver(ax, x, y, u, v, stride=5, scale=None, color="k",
               alpha=0.7, width=0.003, **kwargs):
    """
    Overlay velocity vectors on a 2D field, subsampled by stride.

    Parameters
    ----------
    x, y    : coordinate arrays (2D)
    u, v    : velocity component arrays (2D)
    stride  : subsample every N points in both x and y
    scale   : arrow scale (larger = shorter arrows; None = auto)
    """
    xs = x[::stride, ::stride]
    ys = y[::stride, ::stride]
    us = u[::stride, ::stride]
    vs = v[::stride, ::stride]
    q = ax.quiver(xs, ys, us, vs, color=color, alpha=alpha,
                  scale=scale, width=width, **kwargs)
    return q


# ══════════════════════════════════════════════════════════════════════════════
# 9. SIGNIFICANCE HATCHING
# ══════════════════════════════════════════════════════════════════════════════

def add_significance_hatching(ax, x, y, mask, hatch="///", color="none",
                               edgecolor="0.3", alpha=0.4, linewidth=0.5):
    """
    Overlay hatching where mask == True (e.g., p < 0.05).

    Parameters
    ----------
    mask    : 2D boolean array (True = significant)
    hatch   : hatch pattern string ("///" | "..." | "xxx" etc.)
    color   : facecolor of hatched region ("none" = transparent)
    """
    # Draw a contourf with only one level, using hatching
    mpl.rcParams["hatch.linewidth"] = linewidth
    hatch_data = np.where(mask, 1.0, np.nan)
    cs = ax.contourf(x, y, hatch_data, levels=[0.5, 1.5],
                     colors=color, hatches=[hatch], alpha=alpha)
    for collection in cs.collections:
        collection.set_edgecolor(edgecolor)
    return cs


# ══════════════════════════════════════════════════════════════════════════════
# 10. SAVE HELPER
# ══════════════════════════════════════════════════════════════════════════════

def save(fig, path, dpi=300, formats=("pdf", "png"), **kwargs):
    """
    Save a figure in multiple formats. Default: PDF + PNG at 300 dpi.

    Parameters
    ----------
    path    : base path WITHOUT extension  (e.g. "figures/fig01")
    formats : tuple of extensions to save
    """
    for fmt in formats:
        fig.savefig(f"{path}.{fmt}", dpi=dpi, bbox_inches="tight",
                    pad_inches=0.05, **kwargs)
        print(f"Saved: {path}.{fmt}")


# ══════════════════════════════════════════════════════════════════════════════
# 11. QUICK DEMO  (run this file directly to verify setup)
# ══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    apply_style()
    set_line_palette()

    # ── Demo 1: line plot ────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=figsize(cols=2, rows=1))

    x = np.linspace(0, 4 * np.pi, 200)
    for i, (ax, label) in enumerate(zip(axes, ["Site A", "Site B"])):
        for j in range(4):
            ax.plot(x, np.sin(x + j * 0.5) * np.exp(-0.05 * j * x),
                    label=f"Depth {j*10} m")
        ax.set_xlabel("Time (days)")
        ax.set_ylabel("Temperature anomaly (°C)")
        ax.legend(loc="upper right")
        clean_axis(ax)
        set_minor_ticks(ax)

    label_panels(axes)
    fig.suptitle("Demo – line plot style", y=1.02)
    plt.tight_layout()

    # ── Demo 2: contourf ────────────────────────────────────────────────────
    fig2, axes2 = plt.subplots(1, 2, figsize=figsize(cols=2, rows=1))

    lon = np.linspace(-180, 180, 100)
    lat = np.linspace(-90,   90, 100)
    LON, LAT = np.meshgrid(lon, lat)
    Z_seq = np.sin(np.deg2rad(LON)) * np.cos(np.deg2rad(LAT))   # sequential
    Z_div = np.sin(np.deg2rad(LON * 2)) * np.cos(np.deg2rad(LAT))  # diverging

    cf1 = pcontourf(axes2[0], LON, LAT, Z_seq, cmap=cmap_temp,
                    levels=20, vmin=-1, vmax=1)
    add_colorbar(fig2, cf1, ax=axes2[0], label="Temperature (°C)")

    cf2 = pcontourf(axes2[1], LON, LAT, Z_div, cmap=dcmap_bal,
                    levels=diverging_levels(1.0, n=21), vmin=-1, vmax=1)
    add_colorbar(fig2, cf2, ax=axes2[1], label="SSH anomaly (m)")

    for ax in axes2:
        clean_axis(ax)
        ax.set_xlabel("Longitude (°)")
        ax.set_ylabel("Latitude (°)")

    label_panels(axes2)
    fig2.suptitle("Demo – contourf style", y=1.02)
    plt.tight_layout()
    plt.show()
