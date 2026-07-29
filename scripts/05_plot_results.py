"""
05_plot_results.py
==================
Produces all main result figures for Stage 1. Reads from the computed
outputs in data/ and output/ — no physics computation happens here.

Style follows plot_style.py guidelines: paper-grade rcParams, CVD-safe
colours, 4-sided spines with bidirectional inout ticks, bordered legends,
300 dpi output.

Figures produced
----------------
Primary results (global):
    fig01_flux_vs_co2.png
        Twin-axis: global net ocean CO₂ uptake J_net [Pg C yr⁻¹] (left) vs
        atmospheric CO₂ [ppm] (right). Monthly (faint) + annual mean (solid).
        The primary figure from the meeting.

    fig02_mean_flux_map.png
        Time-mean air-sea flux [mol C m⁻² yr⁻¹]. Blue = uptake, red =
        outgassing. Robinson projection.

    fig03_trend_significance_map.png
        Per-pixel Sen's slope [mol C m⁻² yr⁻¹ per decade] with Mann-Kendall
        significance hatching (/// where p > 0.05 = NOT significant).

    fig04_delta_pco2_map.png
        Time-mean ΔpCO₂ = pCO₂(ocean) − pCO₂(atm) [µatm].

    fig05_seasonal_cycle.png
        Global climatological seasonal cycle: J_net (bars) + atmospheric
        CO₂ (line) on twin axes.

Multi-year variability:
    fig06_loess_trend.png
        Raw monthly (faint) + LOESS smooth (bold) of global J_net.
        PELT breakpoint detection; piecewise linear trends per segment with
        annotated slopes [Pg C yr⁻²].

    fig07_sink_saturation.png
        Scatter + linear fit: annual J_net vs annual atmospheric CO₂. Shows
        whether the ocean uptake is growing proportionally with CO₂ (linear)
        or lagging (saturation). Includes Pearson r and slope annotation.

Regional analyses (Fay 2014 domains):
    fig08_fay_ts.png
        17-panel monthly timeseries of area-integrated flux per Fay domain
        [Pg C yr⁻¹] + global biome map.

    fig09_fay_clim.png
        17-panel climatological seasonal cycle ±1σ per Fay domain.

    fig10_fay_trends.png
        17-panel Sen's slope trend + Mann-Kendall significance per domain.
        Allows identification of which domains are driving the global trend.

    fig07_climate_regression_global.png
        Three-panel: (a) annual J_net timeseries overlaid with normalised
        ONI and SAM; (b) separate Spearman scatter J_net vs ONI and SAM;
        (c) multilinear OLS combination α·ONI + β·SAM scatter coloured by
        year. Following Bellacicco et al. (2025, Nat. Commun.).

Changelog
---------
    v1.0.0 — Initial implementation (5 figures, global only).
    v1.3.0 — Added statistical significance to trend map (Mann-Kendall p-value
              + significance hatching). Added LOESS + breakpoint figure.
              Added sink-saturation scatter. Fay domain figures moved to
              06_fay_analysis.py. 6 figures total.
    v1.3.1 — Added load_climate_indices() (ONI + SAM download + parsing).
              Added fig07_climate_index_regression() — global J_net vs
              ENSO/SAM analysis following Bellacicco et al. (2025).
              7 figures total.

Usage
-----
    python scripts/05_plot_results.py

    Requires data/flux_3d.nc, data/processed_surface.nc,
    output/global_flux.nc, data/co2_mm_gl.csv.
    ONI and SAM files are auto-downloaded to data/ on first run.
    Fay domain figures → run scripts/06_fay_analysis.py.
"""

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import seaborn as sns
from scipy import stats
from scipy.stats import theilslopes

sys.path.insert(0, str(Path(__file__).parent.parent))
import config as cfg

# ── Optional imports ───────────────────────────────────────────────────────────
try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    HAS_CARTOPY = True
except ImportError:
    HAS_CARTOPY = False
    print("[warn] cartopy not found — maps will use plain imshow.")

try:
    import cmocean.cm as cmo
    HAS_CMOCEAN = True
except ImportError:
    HAS_CMOCEAN = False
    print("[warn] cmocean not found — falling back to matplotlib colormaps.")

try:
    from statsmodels.nonparametric.smoothers_lowess import lowess
    HAS_STATSMODELS = True
except ImportError:
    HAS_STATSMODELS = False
    print("[warn] statsmodels not found — LOESS figure (fig06) will be skipped.")

try:
    import ruptures as rpt
    HAS_RUPTURES = True
except ImportError:
    HAS_RUPTURES = False
    print("[warn] ruptures not found — breakpoint detection disabled in fig06.")


# ===========================================================================
# STYLE
# ===========================================================================

sns.set_theme(
    style="ticks", context="paper", font_scale=1.2,
    rc={
        "axes.spines.top":    True,  "axes.spines.right":  True,
        "axes.linewidth":     1.5,
        "xtick.direction":    "inout", "ytick.direction":    "inout",
        "xtick.major.width":  1.2,    "ytick.major.width":  1.2,
        "xtick.minor.width":  0.84,   "ytick.minor.width":  0.84,
        "xtick.major.size":   5,      "ytick.major.size":   5,
        "xtick.minor.size":   3,      "ytick.minor.size":   3,
        "xtick.top":          True,   "ytick.right":        True,
        "font.size":          12,     "axes.labelsize":     12,
        "xtick.labelsize":    11,     "ytick.labelsize":    11,
        "legend.fontsize":    11,     "axes.titlesize":     12,
        "lines.linewidth":    2.0,    "lines.markersize":   6,
        "axes.grid":          True,   "grid.color":         "#CCCCCC",
        "grid.linestyle":     "--",   "grid.linewidth":     0.5,
        "grid.alpha":         0.6,    "axes.axisbelow":     True,
        "savefig.dpi":        300,    "savefig.bbox":       "tight",
        "savefig.pad_inches": 0.05,   "figure.facecolor":   "white",
        "axes.facecolor":     "white","legend.frameon":     True,
        "legend.framealpha":  0.85,   "legend.edgecolor":   "#AAAAAA",
    }
)

# ── Colormaps ──────────────────────────────────────────────────────────────────
CMAP_FLUX  = cmo.balance if HAS_CMOCEAN else "RdBu_r"
CMAP_TREND = cmo.balance if HAS_CMOCEAN else "RdBu_r"
CMAP_PCO2  = cmo.balance if HAS_CMOCEAN else "RdBu_r"

# ── Hardcoded colours ─────────────────────────────────────────────────────────
C_FLUX   = "#0558AD"   # deep blue — primary reconstruction
C_CO2    = "#CC2936"   # deep red  — atmospheric CO2
C_LOESS  = "#1A1A2E"   # near-black — LOESS smooth overlay

# LOESS segment colours (visually distinct, CVD-aware, avoids C_FLUX / C_CO2)
SEG_COLORS = [
    "#E69F00",  # warm amber
    "#2E8B57",  # sea green
    "#9400D3",  # vivid violet
    "#FF6F20",  # burnt orange
    "#1ABDE8",  # cerulean
    "#8B0000",  # dark red
    "#3D3D3D",  # charcoal
    "#B8860B",  # dark goldenrod
]

# ── Labels ────────────────────────────────────────────────────────────────────
REC_LABEL = "GLOBAL_MULTIYEAR_BGC_001_029"


# ===========================================================================
# SHARED HELPERS
# ===========================================================================

def _clean(ax, minor=True):
    """Apply spine width and optional minor ticks."""
    for spine in ax.spines.values():
        spine.set_linewidth(mpl.rcParams["axes.linewidth"])
    if minor:
        ax.xaxis.set_minor_locator(mticker.AutoMinorLocator())
        ax.yaxis.set_minor_locator(mticker.AutoMinorLocator())


def _label_panel(ax, lbl, x=0.015, y=0.97, fontsize=15):
    """Bold panel label (a, b, c …) in the top-left corner."""
    ax.text(x, y, lbl, transform=ax.transAxes,
            fontsize=fontsize, fontweight="bold", va="top", ha="left",
            zorder=10)


def compute_grid_cell_area(lat, lon):
    """Spherical grid-cell area [m²]."""
    R       = cfg.EARTH_RADIUS_M
    d_lat   = float(np.abs(lat.diff("latitude").mean()))
    d_lon   = float(np.abs(lon.diff("longitude").mean()))
    area_1d = (R**2) * np.deg2rad(d_lon) * np.deg2rad(d_lat) * np.cos(np.deg2rad(lat))
    area_2d = area_1d.expand_dims({"longitude": lon}).transpose("latitude", "longitude")
    return area_2d.rename("cell_area")


def mann_kendall_pvalue(y):
    """
    Two-sided Mann-Kendall trend test.

    Returns the p-value (float). Small p → statistically significant trend.
    Uses scipy.stats.kendalltau against a monotone integer sequence, which
    is mathematically equivalent to the standard MK test for a linear trend.

    Reference: Mann (1945); Kendall (1975); Sen (1968).
    """
    x = np.arange(len(y))
    ok = np.isfinite(y)
    if ok.sum() < 5:
        return np.nan
    _, p = stats.kendalltau(x[ok], y[ok])
    return p


def loess_smooth(x, y, frac=0.15):
    """LOESS smoother wrapper (requires statsmodels)."""
    return lowess(y, x, frac=frac, it=3, return_sorted=False)


# ===========================================================================
# DATA LOADING
# ===========================================================================

def load_all() -> dict:
    """
    Load all required output files and apply a global time cutoff
    1993-01-01 to 2024-01-01 (exclusive) to every time-indexed object.
    All downstream figures use this pre-cropped data automatically.
    """
    T_START = "1993-01-01"
    T_END   = "2023-12-31"   # xarray slice is inclusive on both ends

    ds_flux = xr.open_dataset(cfg.DATA_DIR / "flux_3d.nc")
    ds_glob = xr.open_dataset(cfg.OUT_DIR  / "global_flux.nc")
    ds_surf = xr.open_dataset(cfg.DATA_DIR / "processed_surface.nc")

    # Slice all time-indexed xarray objects to the shared window
    flux_3d     = ds_flux["fgco2"].sel(    time=slice(T_START, T_END))
    J_net       = ds_glob["J_net_PgC"].sel(time=slice(T_START, T_END))
    spco2_ocean = ds_surf["spco2_ocean"].sel(time=slice(T_START, T_END))
    spco2_atm   = ds_surf["spco2_atm"].sel( time=slice(T_START, T_END))

    # NOAA atmospheric CO₂ in ppm
    noaa = pd.read_csv(
        cfg.NOAA_CO2_FILE, comment="#",
        names=["year","month","decimal","average","average_unc","trend","trend_unc"],
    )
    noaa["average"] = pd.to_numeric(noaa["average"], errors="coerce")
    noaa = noaa[noaa["average"] > 0].copy()
    noaa["time"] = pd.to_datetime(
        noaa["year"].astype(str) + "-" + noaa["month"].astype(str).str.zfill(2)
    )
    noaa = noaa.set_index("time").sort_index()
    noaa_ppm = noaa["average"].loc[T_START:T_END]

    print(f"[load] Time window applied: {T_START} → {T_END}")
    print(f"       flux_3d : {flux_3d.time.values[0]} → {flux_3d.time.values[-1]}")
    print(f"       J_net   : {J_net.time.values[0]} → {J_net.time.values[-1]}")

    return {
        "flux_3d":      flux_3d,
        "J_net":        J_net,
        "spco2_ocean":  spco2_ocean,
        "spco2_atm":    spco2_atm,
        "ocean_mask":   ds_surf["ocean_mask"],
        "noaa_ppm":     noaa_ppm,
    }


def load_climate_indices() -> dict | None:
    """
    Load ONI and SAM monthly indices, returning annual-mean pandas Series
    aligned to the project time window, or None if the files are unavailable.

    On HPC systems with no outbound internet (e.g. Leonardo), the compute
    nodes cannot reach external URLs.  In that case this function prints
    clear download instructions and returns None so the rest of the pipeline
    continues unaffected.

    Manual download (run once from a login node with internet access):
    ------------------------------------------------------------------
    ONI:
        wget -O <DATA_DIR>/oni_monthly.txt \\
            https://psl.noaa.gov/data/correlation/oni.data

    SAM:
        wget -O <DATA_DIR>/sam_monthly.txt \\
            https://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/aao/monthly.aao.index.b79.current.ascii.table

    Sources
    -------
    ONI : NOAA PSL — ERSSTv5 Niño-3.4 3-month running mean anomaly.
          Format: year  Jan Feb … Dec   (space-separated, missing = -99.90)
    SAM : Marshall (2003) station-based index, NOAA CPC update.
          Format: year  Jan Feb … Dec   (space-separated)

    Scientific note
    ---------------
    Bellacicco et al. (2025) found that interannual variability in ocean
    carbon export correlates with a linear combination of ONI and SAM
    (Spearman r = 0.617, p = 0.001 for PIP export at depth). Here we test
    whether the same climate modes modulate the surface air-sea CO₂ flux.
    """
    import requests

    T_START_YR = 1993
    T_END_YR   = 2023

    oni_file = cfg.DATA_DIR / "oni_monthly.txt"
    sam_file = cfg.DATA_DIR / "sam_monthly.txt"

    # ── Download if absent (skipped gracefully on HPC without internet) ────────
    _DOWNLOAD_NOTE = (
        "\n[climate] *** File not found and download failed (no internet). ***\n"
        "[climate] On a login node with internet access, run:\n"
        f"[climate]   wget -O {oni_file} \\\n"
        "[climate]       https://psl.noaa.gov/data/correlation/oni.data\n"
        f"[climate]   wget -O {sam_file} \\\n"
        "[climate]       https://www.cpc.ncep.noaa.gov/products/precip/CWlink/"
        "daily_ao_index/aao/monthly.aao.index.b79.current.ascii.table\n"
        "[climate] Then re-run 05_plot_results.py — fig07 will be generated.\n"
    )

    for label, url, fpath in [
        ("ONI", "https://psl.noaa.gov/data/correlation/oni.data", oni_file),
        ("SAM",
         "https://www.cpc.ncep.noaa.gov/products/precip/CWlink/"
         "daily_ao_index/aao/monthly.aao.index.b79.current.ascii.table",
         sam_file),
    ]:
        if not fpath.exists():
            print(f"[download] {label} index ...")
            try:
                r = requests.get(url, timeout=30)
                r.raise_for_status()
                fpath.write_bytes(r.content)
                print(f"[ok] {fpath}")
            except Exception as e:
                print(f"[warn] Could not download {label}: {e}")
                print(_DOWNLOAD_NOTE)
                return None

    # ── Parse ONI ─────────────────────────────────────────────────────────────
    oni_rows = []
    with open(oni_file) as f:
        for line in f:
            parts = line.split()
            if len(parts) == 13:
                try:
                    yr = int(parts[0])
                    vals = [float(v) for v in parts[1:]]
                    for m, v in enumerate(vals, 1):
                        if v != -99.9:
                            oni_rows.append({"year": yr, "month": m, "oni": v})
                except ValueError:
                    pass
    if not oni_rows:
        print(f"[warn] ONI file {oni_file} parsed empty — check format.")
        return None
    oni_df = pd.DataFrame(oni_rows)
    oni_df["time"] = pd.to_datetime(
        oni_df["year"].astype(str) + "-" + oni_df["month"].astype(str).str.zfill(2)
    )
    oni_monthly = oni_df.set_index("time")["oni"].sort_index()

    # ── Parse SAM ─────────────────────────────────────────────────────────────
    sam_rows = []
    with open(sam_file) as f:
        for line in f:
            parts = line.split()
            if len(parts) == 13:
                try:
                    yr = int(parts[0])
                    vals = [float(v) for v in parts[1:]]
                    for m, v in enumerate(vals, 1):
                        sam_rows.append({"year": yr, "month": m, "sam": v})
                except ValueError:
                    pass
    if not sam_rows:
        print(f"[warn] SAM file {sam_file} parsed empty — check format.")
        return None
    sam_df = pd.DataFrame(sam_rows)
    sam_df["time"] = pd.to_datetime(
        sam_df["year"].astype(str) + "-" + sam_df["month"].astype(str).str.zfill(2)
    )
    sam_monthly = sam_df.set_index("time")["sam"].sort_index()

    # ── Annual means within project window ────────────────────────────────────
    t_range    = pd.date_range(f"{T_START_YR}", f"{T_END_YR}", freq="YE")
    oni_annual = oni_monthly.resample("YE").mean().reindex(t_range, method="nearest")
    sam_annual = sam_monthly.resample("YE").mean().reindex(t_range, method="nearest")

    print(f"[load] ONI annual : {int(oni_annual.index[0].year)} → {int(oni_annual.index[-1].year)}")
    print(f"[load] SAM annual : {int(sam_annual.index[0].year)} → {int(sam_annual.index[-1].year)}")

    return {
        "oni_monthly": oni_monthly,
        "sam_monthly": sam_monthly,
        "oni_annual":  oni_annual,
        "sam_annual":  sam_annual,
    }



# ===========================================================================
# FIG 01 — GLOBAL FLUX EVOLUTION: LOESS + BREAKPOINTS + ATMOSPHERIC CO₂
# ===========================================================================

def fig01_flux_vs_co2(data: dict) -> None:
    """
    Merged primary figure combining the net flux evolution with LOESS
    smoothing, piecewise trend analysis and atmospheric CO₂.

    Layout
    ------
    Left axis  : J_net [Pg C yr⁻¹]
        - Raw monthly values (faint fill + thin line)
        - LOESS smooth (frac=0.15, bold)
        - Piecewise trend lines for 3 segments defined by B1 and B4
        - Segment trend slope labels ABOVE each trend line
        - Breakpoint vertical lines at B1 and B4

    Right axis : Atmospheric CO₂ [ppm] (NOAA GML, annual mean, dashed red)
        - NOAA trend label BELOW the CO₂ line

    Breakpoint strategy
    -------------------
    PELT detects all breakpoints on the LOESS smooth. We keep only B1 (first)
    and B4 (fourth) — or the first and last if fewer than 4 exist — and use
    them to define three segments:
        Segment 1 : start → B1
        Segment 2 : B1 → B4
        Segment 3 : B4 → end
    Each segment gets an OLS trend line [Pg C yr⁻²] annotated above the line.
    """
    if not HAS_STATSMODELS:
        print("[skip] fig01 LOESS requires statsmodels — plotting simple twin-axis only.")
        _fig01_simple(data)
        return

    J_net = data["J_net"]
    noaa  = data["noaa_ppm"]

    t     = pd.to_datetime(J_net.time.values)
    v     = J_net.values.copy()
    x_num = np.arange(len(t), dtype=float)
    ok    = np.isfinite(v)

    # LOESS smooth
    sm = loess_smooth(x_num[ok], v[ok], frac=0.15)

    # ── Breakpoints ────────────────────────────────────────────────────────────
    t_ok = t[ok]

    if HAS_RUPTURES:
        signal = sm.reshape(-1, 1)
        algo   = rpt.Pelt(model="rbf").fit(signal)
        bkps   = algo.predict(pen=3)   # includes sentinel at end
        interior = bkps[:-1]           # actual breakpoints, drop sentinel

        print(f"\n[fig01] All PELT breakpoints detected (n={len(interior)}):")
        for k_bp, bp in enumerate(interior):
            print(f"  B{k_bp+1}: {t_ok[bp-1].strftime('%Y-%m')}")

        # Keep B1 (index 0) and B4 (index 3) if they exist,
        # otherwise fall back to first and last interior breakpoint
        if len(interior) >= 4:
            kept = [interior[0], interior[3]]
        elif len(interior) >= 2:
            kept = [interior[0], interior[-1]]
        elif len(interior) == 1:
            kept = [interior[0], interior[0]]
        else:
            kept = []

        if len(kept) == 2:
            b1_idx, b4_idx = kept
            b1_date = t_ok[b1_idx - 1]
            b4_date = t_ok[b4_idx - 1]
            print(f"  Using B1={b1_date.strftime('%Y-%m')}, "
                  f"B4={b4_date.strftime('%Y-%m')} for 3-segment analysis")
            segs = [0, b1_idx, b4_idx, len(sm)]
        else:
            print("  [warn] Fewer than 2 breakpoints — no piecewise trend drawn")
            segs = []
    else:
        print("[warn] ruptures not installed — breakpoints skipped in fig01.")
        segs = []

    # ── Figure ─────────────────────────────────────────────────────────────────
    fig, ax1 = plt.subplots(figsize=(15, 6))

    # Monthly flux — faint fill + thin line
    ax1.fill_between(t, v, 0, alpha=0.12, color=C_FLUX)
    ax1.plot(t, v, color=C_FLUX, lw=0.7, alpha=0.35, label="_nolegend_")

    # LOESS smooth
    ax1.plot(t_ok, sm, color=C_FLUX, lw=2.5,
             label=f"{REC_LABEL}  (LOESS, frac=0.15)")

    # ── Piecewise trend segments ────────────────────────────────────────────────
    seg_labels = ["Segment 1", "Segment 2", "Segment 3"]
    for seg_idx in range(len(segs) - 1):
        i0  = segs[seg_idx]
        i1  = segs[seg_idx + 1]
        col = SEG_COLORS[seg_idx % len(SEG_COLORS)]

        t_seg = t_ok[i0:i1]
        y_seg = sm[i0:i1]
        x_seg = np.arange(i0, i1, dtype=float)

        if len(t_seg) < 3:
            continue

        sl, ic, *_ = stats.linregress(x_seg, y_seg)
        sl_yr      = sl * 12     # per month → per year [Pg C yr⁻²]
        y_fit      = ic + sl * x_seg

        ax1.plot(t_seg, y_fit, color=col, lw=2.0, ls=":", alpha=0.95,
                 label=f"{seg_labels[seg_idx]}: {sl_yr:+.3f} Pg C yr⁻²")

        # Slope annotation ABOVE the trend line midpoint
        mid = len(x_seg) // 2
        ax1.annotate(
            f"{sl_yr:+.3f} Pg C yr⁻²",
            xy=(t_seg[mid], y_fit[mid]),
            xytext=(0, 14), textcoords="offset points",
            fontsize=10, fontweight="bold", color=col,
            ha="center", va="bottom",
            arrowprops=dict(arrowstyle="-", color=col, lw=1.6),
        )

    # ── Breakpoint vertical lines + top labels ─────────────────────────────────
    bp_tick_dates  = []
    bp_tick_labels = []
    if len(segs) > 0:
        for k, (idx, label) in enumerate(zip([segs[1], segs[2]],
                                              ["B1", "B4"])):
            if idx == 0 or idx >= len(t_ok):
                continue
            bp_date = t_ok[idx - 1]
            ax1.axvline(bp_date, color="#444444", lw=1.4, ls="-.", alpha=0.85,
                        zorder=3)
            bp_tick_dates.append(bp_date)
            bp_tick_labels.append(label)

    # B1 / B4 labels on top x-axis
    if bp_tick_dates:
        ax_top = ax1.twiny()
        ax_top.set_xlim(ax1.get_xlim())
        ax_top.set_xticks(bp_tick_dates)
        ax_top.set_xticklabels(bp_tick_labels, fontsize=10, fontweight="bold",
                               color="#444444")
        ax_top.tick_params(axis="x", length=6, width=1.2, direction="in",
                           color="#444444")
        ax_top.spines["top"].set_linewidth(0)

    ax1.axhline(0, color="black", lw=0.8, ls=":")
    ax1.set_xlabel("Year")
    ax1.set_ylabel("Net ocean CO₂ uptake  [Pg C yr⁻¹]\n(positive = ocean sink)",
                   color=C_FLUX)
    ax1.tick_params(axis="y", labelcolor=C_FLUX)
    ax1.set_ylim(-0.2, 3.0)

    # ── Right axis: atmospheric CO₂ ────────────────────────────────────────────
    noaa_annual = noaa.resample("YE").mean()

    ax2 = ax1.twinx()
    ax2.plot(noaa_annual.index, noaa_annual.values, color=C_CO2,
             lw=2.2, ls="--", label="Atmospheric CO₂  (NOAA GML)")
    ax2.set_ylabel("Atmospheric CO₂  [ppm]", color=C_CO2)
    ax2.tick_params(axis="y", labelcolor=C_CO2)

    # NOAA trend annotation BELOW the CO₂ line (mid-record)
    noaa_v   = noaa_annual.dropna()
    noaa_t   = noaa_v.index
    x_noaa   = np.arange(len(noaa_t), dtype=float)
    sl_n, ic_n, *_ = stats.linregress(x_noaa, noaa_v.values)
    sl_n_yr  = sl_n               # already per year [ppm yr⁻¹]
    mid_n    = len(noaa_t) // 2
    y_fit_n  = ic_n + sl_n * x_noaa
    ax2.annotate(
        f"+{sl_n_yr:.2f} ppm yr⁻¹",
        xy=(noaa_t[mid_n], y_fit_n[mid_n]),
        xytext=(0, -18), textcoords="offset points",
        fontsize=10, fontweight="bold", color=C_CO2,
        ha="center", va="top",
        arrowprops=dict(arrowstyle="-", color=C_CO2, lw=1.6),
    )

    # ── Legend ─────────────────────────────────────────────────────────────────
    lines1, labs1 = ax1.get_legend_handles_labels()
    lines2, labs2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labs1 + labs2,
               loc="upper left", fontsize=9,
               framealpha=0.88, edgecolor="#AAAAAA", ncol=2)

    ax1.set_title(
        "Global ocean CO₂ uptake — multi-year variability and trends\n"
        f"{REC_LABEL}  ·  Wanninkhof (2014) k  ·  Weiss (1974) K₀  "
        "·  LOESS frac=0.15  ·  B1/B4 breakpoints",
        fontsize=11,
    )
    _clean(ax1, minor=False)
    fig.tight_layout()

    out = cfg.FIG_DIR / "fig01_flux_vs_co2.png"
    fig.savefig(out, dpi=cfg.FIGURE_DPI)
    plt.close(fig)
    print(f"[fig] {out}")


def _fig01_simple(data: dict) -> None:
    """Fallback fig01 without LOESS (statsmodels not available)."""
    J_net = data["J_net"]
    noaa  = data["noaa_ppm"]
    t     = pd.to_datetime(J_net.time.values)
    J_da  = xr.DataArray(J_net.values, coords={"time": t}, dims="time")
    J_ann = J_da.resample(time="1YE").mean()

    fig, ax1 = plt.subplots(figsize=(15, 6))
    ax1.fill_between(t, J_net.values, 0, alpha=0.15, color=C_FLUX)
    ax1.plot(t, J_net.values, color=C_FLUX, lw=0.8, alpha=0.4)
    ax1.plot(pd.to_datetime(J_ann.time.values), J_ann.values,
             color=C_FLUX, lw=2.5,
             label="Net ocean uptake  $J_\\mathrm{net}$  (annual mean)")
    ax1.axhline(0, color="black", lw=0.8, ls=":")
    ax1.set_ylabel("Net ocean CO₂ uptake  [Pg C yr⁻¹]", color=C_FLUX)
    ax1.tick_params(axis="y", labelcolor=C_FLUX)
    ax1.set_ylim(-0.2, 3.0)
    ax2 = ax1.twinx()
    ax2.plot(noaa.resample("YE").mean().index,
             noaa.resample("YE").mean().values,
             color=C_CO2, lw=2.2, ls="--", label="Atmospheric CO₂ (NOAA GML)")
    ax2.set_ylabel("Atmospheric CO₂  [ppm]", color=C_CO2)
    ax2.tick_params(axis="y", labelcolor=C_CO2)
    lines1, l1 = ax1.get_legend_handles_labels()
    lines2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, l1 + l2, loc="upper left", fontsize=10)
    ax1.set_title("Global ocean CO₂ uptake vs. atmospheric CO₂")
    _clean(ax1, minor=False)
    fig.tight_layout()
    fig.savefig(cfg.FIG_DIR / "fig01_flux_vs_co2.png", dpi=cfg.FIGURE_DPI)
    plt.close(fig)
    print(f"[fig] {cfg.FIG_DIR}/fig01_flux_vs_co2.png")


# ===========================================================================
# FIG 02 — TIME-MEAN SPATIAL FLUX MAP
# ===========================================================================

def fig02_mean_flux_map(data: dict) -> None:
    """
    Time-mean air-sea CO₂ flux [mol C m⁻² yr⁻¹].
    Diverging colormap with white at zero: blue = uptake, red = outgassing.
    Fixed symmetric range ±2 mol C m⁻² yr⁻¹ to keep the colour scale
    comparable across runs and highlight the dominant signal.
    """
    flux_mean = data["flux_3d"].mean(dim="time").where(data["ocean_mask"] == 1)

    # White-centred diverging colormap: RdBu_r has blue=positive, red=negative
    # which matches our convention (blue=uptake=positive). If cmocean available
    # use balance (also white-centred), else fall back to RdBu_r.
    cmap_fig02 = cmo.balance if HAS_CMOCEAN else "RdBu_r"

    fig, ax = _map_axes()
    _plot_map(fig, ax, flux_mean, cmap_fig02, -2.5, 2.5,
              cbar_label="Air-sea CO₂ flux  [mol C m⁻² yr⁻¹]  (+ = uptake, − = outgassing)",
              title="Time-mean air-sea CO₂ flux  [mol C m⁻² yr⁻¹]")

    out = cfg.FIG_DIR / "fig02_mean_flux_map.png"
    fig.savefig(out, dpi=cfg.FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {out}")


# ===========================================================================
# FIG 03 — TREND MAP WITH MANN-KENDALL SIGNIFICANCE
# ===========================================================================

def fig03_trend_significance_map(data: dict) -> None:
    """
    Per-pixel Sen's slope of annual mean flux [mol C m⁻² yr⁻¹ per decade]
    with Mann-Kendall significance overlay.

    Hatching (///) indicates pixels where the trend is NOT statistically
    significant at α = 0.05 (two-sided Mann-Kendall p > 0.05).
    Unhatched pixels have a significant trend.

    Scientific rationale:
        Sen's slope is robust to outliers and non-normality (Theil-Sen
        estimator). Mann-Kendall is the standard non-parametric significance
        test for monotone trends and is consistent with Sen's slope.
        Serial autocorrelation in monthly data is avoided by working on
        annual means (reduces effective n but increases independence).
    """
    print("[compute] Per-pixel Sen's slope + Mann-Kendall significance ...")
    flux_annual = data["flux_3d"].resample(time="1YE").mean()
    n_yr   = len(flux_annual.time)
    times  = np.arange(n_yr, dtype=float)

    lat    = flux_annual.latitude.values
    lon    = flux_annual.longitude.values
    fv     = flux_annual.values          # (years, lat, lon)

    slope_arr = np.full((len(lat), len(lon)), np.nan)
    pval_arr  = np.full((len(lat), len(lon)), np.nan)

    for i in range(len(lat)):
        for j in range(len(lon)):
            y  = fv[:, i, j]
            ok = np.isfinite(y)
            if ok.sum() < 5:
                continue
            res = theilslopes(y[ok], times[ok])
            slope_arr[i, j] = res.slope * 10   # yr⁻¹ → decade⁻¹
            pval_arr[i, j]  = mann_kendall_pvalue(y)

    slope_da = xr.DataArray(
        slope_arr, coords={"latitude": lat, "longitude": lon},
        dims=["latitude", "longitude"],
        attrs={"long_name": "Sen's slope of annual flux",
               "units":     "mol C m-2 yr-1 per decade"},
    ).where(data["ocean_mask"] == 1)

    # Significance mask: True where NOT significant at α = 0.01 (stricter)
    # AND restricted to ocean pixels only — land pixels must never be hatched.
    ocean_2d   = data["ocean_mask"].values.astype(bool)
    insig_mask = ((pval_arr > 0.01) | ~np.isfinite(pval_arr)) & ocean_2d

    fig, ax = _map_axes()
    _plot_map(fig, ax, slope_da, CMAP_TREND, -0.5, 0.5,
              cbar_label="Trend in CO₂ flux  [mol C m⁻² yr⁻¹ per decade]\n"
                         "+ = increasing uptake  |  − = increasing outgassing",
              title="Per-pixel trend in air-sea CO₂ flux (Sen's slope, annual means)\n"
                    "Hatching (///) = NOT significant at α = 0.01 (Mann-Kendall)")

    # Overlay significance hatching + outline — ocean pixels only
    lon2d, lat2d = np.meshgrid(lon, lat)
    hatch_data   = np.where(insig_mask, 1.0, np.nan)

    # Hatching fill (///) — set colour via rcParams for matplotlib >= 3.8 compat
    with mpl.rc_context({"hatch.linewidth": 0.4, "hatch.color": "0.35"}):
        if HAS_CARTOPY:
            ax.contourf(lon2d, lat2d, hatch_data,
                        levels=[0.5, 1.5], colors="none",
                        hatches=["///"],
                        transform=ccrs.PlateCarree(), zorder=4)
        else:
            ax.contourf(lon2d, lat2d, hatch_data,
                        levels=[0.5, 1.5], colors="none",
                        hatches=["///"])

    # Outline the non-significant boundary with a darker gray contour line
    # Use a float array so contour has a clean single iso-level to trace
    outline_data = np.where(insig_mask, 1.0, 0.0)
    if HAS_CARTOPY:
        ax.contour(lon2d, lat2d, outline_data,
                   levels=[0.5], colors=["#555555"], linewidths=[0.6],
                   transform=ccrs.PlateCarree(), zorder=5)
    else:
        ax.contour(lon2d, lat2d, outline_data,
                   levels=[0.5], colors=["#555555"], linewidths=[0.6])

    out = cfg.FIG_DIR / "fig03_trend_significance_map.png"
    fig.savefig(out, dpi=cfg.FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {out}")


# ===========================================================================
# FIG 04 — ΔpCO₂ MAP
# ===========================================================================

def fig04_delta_pco2_map(data: dict) -> None:
    """
    Time-mean ΔpCO₂ = pCO₂(ocean) − pCO₂(atm) [µatm].
    Positive = ocean is a CO₂ source; negative = potential sink.
    """
    dpco2 = (data["spco2_ocean"] - data["spco2_atm"]).mean(dim="time")
    dpco2 = (dpco2 * 1e6).where(data["ocean_mask"] == 1)
    dpco2.attrs["units"] = "µatm"

    fig, ax = _map_axes()
    _plot_map(fig, ax, dpco2, CMAP_PCO2, -100.0, 100.0,
              cbar_label="ΔpCO₂ = pCO₂(ocean) − pCO₂(atm)  [µatm]  (+ = source, − = sink)",
              title="Time-mean ΔpCO₂ — thermodynamic driver of the air-sea flux")

    out = cfg.FIG_DIR / "fig04_delta_pco2_map.png"
    fig.savefig(out, dpi=cfg.FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {out}")


# ===========================================================================
# FIG 05 — GLOBAL SEASONAL CYCLE
# ===========================================================================

def fig05_seasonal_cycle(data: dict) -> None:
    """
    Climatological monthly mean J_net [Pg C yr⁻¹] (bars, left axis)
    and atmospheric CO₂ [ppm] (line, right axis).
    """
    J_net = data["J_net"]
    noaa  = data["noaa_ppm"]

    J_da   = xr.DataArray(J_net.values,
                           coords={"time": pd.to_datetime(J_net.time.values)},
                           dims="time")
    J_clim = J_da.groupby("time.month").mean()

    noaa.index = pd.to_datetime(noaa.index)
    noaa_clim  = noaa.groupby(noaa.index.month).mean()

    months      = np.arange(1, 13)
    month_labs  = ["Jan","Feb","Mar","Apr","May","Jun",
                   "Jul","Aug","Sep","Oct","Nov","Dec"]

    fig, ax1 = plt.subplots(figsize=(10, 5))
    ax1.bar(months, J_clim.values, color=C_FLUX, alpha=0.7,
            label="$J_\\mathrm{net}$ seasonal cycle")
    ax1.axhline(0, color="black", lw=0.8)
    ax1.set_xlabel("Month")
    ax1.set_ylabel("Global net uptake  [Pg C yr⁻¹]", color=C_FLUX)
    ax1.tick_params(axis="y", labelcolor=C_FLUX)
    ax1.set_xticks(months)
    ax1.set_xticklabels(month_labs)

    ax2 = ax1.twinx()
    ax2.plot(months, noaa_clim.values, color=C_CO2, lw=2.2, marker="o",
             label="Atm CO₂ (ppm)")
    ax2.set_ylabel("Atmospheric CO₂  [ppm]", color=C_CO2)
    ax2.tick_params(axis="y", labelcolor=C_CO2)

    lines1, labs1 = ax1.get_legend_handles_labels()
    lines2, labs2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labs1 + labs2,
               loc="lower right", fontsize=10, edgecolor="#AAAAAA")
    ax1.set_title("Climatological seasonal cycle — global ocean uptake vs. atmospheric CO₂")
    _clean(ax1)
    fig.tight_layout()

    out = cfg.FIG_DIR / "fig05_seasonal_cycle.png"
    fig.savefig(out, dpi=cfg.FIGURE_DPI)
    plt.close(fig)
    print(f"[fig] {out}")


# ===========================================================================
# FIG 06 — SINK SATURATION SCATTER  (renumbered from fig07)
# ===========================================================================

def fig06_sink_saturation(data: dict) -> None:
    """
    Scatter plot of annual J_net vs annual atmospheric CO₂.

    Scientific rationale:
        If the ocean sink scales linearly with rising atmospheric CO₂, the
        scatter should follow a straight line with positive slope (more CO₂
        → more uptake). If the sink is saturating, the relationship will
        flatten or bend over time. Colouring by year reveals temporal
        evolution. The regression slope β₁ [Pg C yr⁻¹ / ppm] quantifies
        ocean sensitivity to atmospheric CO₂ change.

        This implements the analysis requested in the meeting:
        "if [flux] is constant while [CO₂] is rising, saturation is occurring"
    """
    # Data already cropped 1993–2024 by load_all()
    J_net = data["J_net"]
    noaa  = data["noaa_ppm"]

    t_monthly   = pd.to_datetime(J_net.time.values)
    J_da        = xr.DataArray(J_net.values, coords={"time": t_monthly}, dims="time")
    J_annual    = J_da.resample(time="1YE").mean()
    t_annual    = pd.to_datetime(J_annual.time.values)
    years       = t_annual.year.values
    noaa_annual = noaa.resample("YE").mean()

    # Align NOAA annual means onto the J_net annual time axis
    co2_vals = noaa_annual.reindex(t_annual, method="nearest").values
    j_vals   = J_annual.values

    ok = np.isfinite(j_vals) & np.isfinite(co2_vals)
    x  = co2_vals[ok]
    y  = j_vals[ok]
    yr = years[ok]

    # Linear regression
    sl, ic, r, p, se = stats.linregress(x, y)
    x_fit = np.linspace(x.min(), x.max(), 100)
    y_fit = ic + sl * x_fit

    fig, ax = plt.subplots(figsize=(8, 6))

    # Scatter coloured by year
    sc = ax.scatter(x, y, c=yr, cmap="plasma", s=60, zorder=5,
                    edgecolors="white", linewidths=0.5)
    cbar = fig.colorbar(sc, ax=ax, pad=0.02, shrink=0.85)
    cbar.set_label("Year", fontsize=11)

    # Regression line
    ax.plot(x_fit, y_fit, color=C_FLUX, lw=2.2, ls="--", zorder=4,
            label=f"Linear fit\n$\\beta$ = {sl:.4f} Pg C yr⁻¹ ppm⁻¹\n"
                  f"$r$ = {r:.3f},  $p$ = {p:.3e}")

    ax.axhline(0, color="black", lw=0.7, ls=":")
    ax.set_xlabel("Atmospheric CO₂  [ppm]")
    ax.set_ylabel("Net ocean CO₂ uptake  $J_\\mathrm{net}$  [Pg C yr⁻¹]")
    ax.set_title(
        "Ocean sink sensitivity to atmospheric CO₂\n"
        "Positive slope = growing sink  |  Flat = saturation"
    )
    ax.legend(fontsize=10, loc="upper left", edgecolor="#AAAAAA")
    _clean(ax)
    fig.tight_layout()

    out = cfg.FIG_DIR / "fig06_sink_saturation.png"
    fig.savefig(out, dpi=cfg.FIGURE_DPI)
    plt.close(fig)
    print(f"[fig] {out}")


# ===========================================================================
# MAP HELPER UTILITIES
# ===========================================================================

# ===========================================================================
# FIG 07 — GLOBAL J_NET vs ONI + SAM CLIMATE INDICES
# ===========================================================================

def fig07_climate_index_regression(data: dict, idx: dict) -> None:
    """
    Three-panel figure linking global annual J_net to ENSO (ONI) and SAM.

    Panel a — Timeseries
        Annual J_net (blue, left axis) overlaid with normalised ONI (red)
        and normalised SAM (green) on the right axis. Reveals years where
        ENSO or SAM events coincide with anomalous ocean uptake.

    Panel b — Separate Spearman correlations
        Scatter of J_net vs ONI and J_net vs SAM, with Spearman ρ and
        p-value annotated.  Mirrors Bellacicco et al. (2025) Fig. 3 approach.

    Panel c — Multilinear combination
        J_net vs the OLS-fitted linear combination α·ONI + β·SAM.
        Spearman correlation and p-value of the combined predictor annotated.
        Points coloured by year.

    Scientific rationale (Bellacicco et al. 2025)
    -----------------------------------------------
    Interannual PIP export correlates with ONI (ρ = 0.57) and improves
    with a multilinear ONI+SAM combination (ρ = 0.617). Testing the same
    relationship for surface J_net reveals whether the ENSO/SAM teleconnection
    already manifests at the air-sea interface or only emerges at depth through
    the physical injection pump.
    """
    J_net = data["J_net"]
    t     = pd.to_datetime(J_net.time.values)
    J_da  = xr.DataArray(J_net.values, coords={"time": t}, dims="time")
    J_ann = J_da.resample(time="1YE").mean()
    t_ann = pd.to_datetime(J_ann.time.values)
    years = t_ann.year.values

    oni = idx["oni_annual"].reindex(t_ann, method="nearest").values
    sam = idx["sam_annual"].reindex(t_ann, method="nearest").values
    j   = J_ann.values

    # Drop years where any field is NaN
    ok  = np.isfinite(j) & np.isfinite(oni) & np.isfinite(sam)
    j, oni, sam, years_ok = j[ok], oni[ok], sam[ok], years[ok]
    t_ok = t_ann[ok]

    # Z-score normalise indices for plotting and multilinear fit
    oni_z = (oni - oni.mean()) / oni.std()
    sam_z = (sam - sam.mean()) / sam.std()

    # Spearman correlations — separate
    rho_oni, p_oni = stats.spearmanr(j, oni)
    rho_sam, p_sam = stats.spearmanr(j, sam)

    # OLS multilinear fit: j ~ alpha*oni_z + beta*sam_z + intercept
    X     = np.column_stack([oni_z, sam_z, np.ones(len(j))])
    coef, *_ = np.linalg.lstsq(X, j, rcond=None)
    j_pred   = X @ coef
    rho_ml, p_ml = stats.spearmanr(j, j_pred)

    C_ONI = "#CC2936"   # red
    C_SAM = "#2E8B57"   # green
    C_ML  = "#0558AD"   # blue

    fig = plt.figure(figsize=(18, 6))
    gs  = gridspec.GridSpec(1, 3, figure=fig, wspace=0.35,
                            left=0.06, right=0.97, top=0.88, bottom=0.13)

    # ── Panel a: timeseries ───────────────────────────────────────────────────
    ax_a  = fig.add_subplot(gs[0])
    ax_a2 = ax_a.twinx()

    ax_a.bar(t_ok, j, color=C_ML, alpha=0.55, width=300,
             label="$J_\\mathrm{net}$ (annual)")
    ax_a.axhline(0, color="black", lw=0.7, ls=":")
    ax_a.set_ylabel("$J_\\mathrm{net}$  [Pg C yr⁻¹]", color=C_ML)
    ax_a.tick_params(axis="y", labelcolor=C_ML)
    ax_a.set_ylim(-0.5, max(j) * 1.35)

    ax_a2.plot(t_ok, oni_z, color=C_ONI, lw=1.8, ls="-",  label="ONI (z-score)")
    ax_a2.plot(t_ok, sam_z, color=C_SAM, lw=1.8, ls="--", label="SAM (z-score)")
    ax_a2.axhline(0, color="gray", lw=0.5, ls=":")
    ax_a2.set_ylabel("Climate index (z-score)", color="gray")
    ax_a2.tick_params(axis="y", labelcolor="gray")

    lines_a, labs_a   = ax_a.get_legend_handles_labels()
    lines_a2, labs_a2 = ax_a2.get_legend_handles_labels()
    ax_a.legend(lines_a + lines_a2, labs_a + labs_a2,
                fontsize=9, loc="upper left", edgecolor="#AAAAAA")
    ax_a.set_title("(a)  Annual $J_\\mathrm{net}$ with climate indices",
                   fontsize=11, fontweight="bold")
    ax_a.tick_params(axis="x", rotation=30)
    _clean(ax_a, minor=False)

    # ── Panel b: separate scatter ─────────────────────────────────────────────
    ax_b = fig.add_subplot(gs[1])
    ax_b.scatter(oni, j, color=C_ONI, s=55, zorder=4, label=
                 f"vs ONI  ρ = {rho_oni:.2f}, p = {p_oni:.3f}",
                 edgecolors="white", linewidths=0.5)
    ax_b.scatter(sam, j, color=C_SAM, s=55, marker="^", zorder=4, label=
                 f"vs SAM  ρ = {rho_sam:.2f}, p = {p_sam:.3f}",
                 edgecolors="white", linewidths=0.5)

    for x_arr, col in [(oni, C_ONI), (sam, C_SAM)]:
        sl, ic, *_ = stats.linregress(x_arr, j)
        xs = np.linspace(x_arr.min(), x_arr.max(), 80)
        ax_b.plot(xs, ic + sl * xs, color=col, lw=1.6, ls="--", alpha=0.75)

    ax_b.axhline(0, color="black", lw=0.7, ls=":")
    ax_b.axvline(0, color="gray", lw=0.5, ls=":")
    ax_b.set_xlabel("Climate index value")
    ax_b.set_ylabel("Annual $J_\\mathrm{net}$  [Pg C yr⁻¹]")
    ax_b.legend(fontsize=9, edgecolor="#AAAAAA", loc="upper left")
    ax_b.set_title("(b)  Separate Spearman correlations",
                   fontsize=11, fontweight="bold")
    _clean(ax_b)

    # ── Panel c: multilinear scatter ──────────────────────────────────────────
    ax_c = fig.add_subplot(gs[2])
    sc = ax_c.scatter(j_pred, j, c=years_ok, cmap="plasma",
                      s=60, zorder=5,
                      edgecolors="white", linewidths=0.5)
    cbar = fig.colorbar(sc, ax=ax_c, pad=0.02, shrink=0.85)
    cbar.set_label("Year", fontsize=10)

    sl_ml, ic_ml, *_ = stats.linregress(j_pred, j)
    xs = np.linspace(j_pred.min(), j_pred.max(), 80)
    ax_c.plot(xs, ic_ml + sl_ml * xs, color=C_ML, lw=2.0, ls="--",
              label=f"ρ = {rho_ml:.2f}, p = {p_ml:.3f}")

    ax_c.set_xlabel(f"α·ONI + β·SAM  (α={coef[0]:.3f}, β={coef[1]:.3f})  [Pg C yr⁻¹]")
    ax_c.set_ylabel("Annual $J_\\mathrm{net}$  [Pg C yr⁻¹]")
    ax_c.legend(fontsize=10, edgecolor="#AAAAAA", loc="upper left")
    ax_c.set_title("(c)  Multilinear combination (ONI + SAM)",
                   fontsize=11, fontweight="bold")
    _clean(ax_c)

    fig.suptitle(
        "Global surface CO₂ flux — ENSO & SAM climate index regression\n"
        f"{REC_LABEL}  ·  following Bellacicco et al. (2025)",
        fontsize=12,
    )

    out = cfg.FIG_DIR / "fig07_climate_regression_global.png"
    fig.savefig(out, dpi=cfg.FIGURE_DPI)
    plt.close(fig)
    print(f"[fig] {out}")


def _map_axes():
    """Return (fig, ax) with Robinson projection if cartopy is available."""
    if HAS_CARTOPY:
        fig, ax = plt.subplots(
            figsize=(14, 7),
            subplot_kw={"projection": ccrs.Robinson()},
        )
    else:
        fig, ax = plt.subplots(figsize=(14, 6))
    return fig, ax


def _plot_map(fig, ax, da, cmap, vmin, vmax, cbar_label, title):
    """Shared map plotting logic for figs 02, 03, 04."""
    lat = da["latitude"].values
    lon = da["longitude"].values
    Z   = da.values
    levels = np.linspace(vmin, vmax, 21)

    if HAS_CARTOPY:
        cf = ax.contourf(lon, lat, Z, levels=levels, cmap=cmap,
                         extend="both", transform=ccrs.PlateCarree())
        ax.add_feature(cfeature.LAND,      facecolor="#d3d3d3", zorder=2)
        ax.add_feature(cfeature.COASTLINE, lw=0.4, zorder=3)
        ax.set_global()
    else:
        cf = ax.contourf(lon, lat, Z, levels=levels, cmap=cmap, extend="both")
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")

    cb = plt.colorbar(cf, ax=ax, orientation="horizontal",
                      shrink=0.6, pad=0.04, aspect=35)
    cb.set_label(cbar_label, fontsize=11)
    cb.ax.tick_params(labelsize=10, width=1.0, length=3)
    ax.set_title(title, fontsize=12)



# ===========================================================================
# MAIN
# ===========================================================================

def main():
    cfg.FIG_DIR.mkdir(parents=True, exist_ok=True)

    print("[load] Output files ...")
    data = load_all()

    print("[load] Climate indices (ONI, SAM) ...")
    idx = load_climate_indices()

    print("\n[figure 01] Global flux evolution — LOESS + breakpoints + CO₂ ...")
    fig01_flux_vs_co2(data)

    print("[figure 02] Time-mean flux map ...")
    fig02_mean_flux_map(data)

    print("[figure 03] Trend map + Mann-Kendall significance ...")
    fig03_trend_significance_map(data)

    print("[figure 04] ΔpCO₂ map ...")
    fig04_delta_pco2_map(data)

    print("[figure 05] Seasonal cycle ...")
    fig05_seasonal_cycle(data)

    print("[figure 06] Sink saturation scatter ...")
    fig06_sink_saturation(data)

    if idx is not None:
        print("[figure 07] Climate index regression (global) ...")
        fig07_climate_index_regression(data, idx)
    else:
        print("[skip] fig07 — climate index files not available (see [climate] messages above).")

    n_figs = 7 if idx is not None else 6
    print(f"\n[done] {n_figs} figures saved to {cfg.FIG_DIR}")
    print(f"       Fay domain figures (+ climate regression per domain)"
          f" → run scripts/06_fay_analysis.py\n")


if __name__ == "__main__":
    main()
