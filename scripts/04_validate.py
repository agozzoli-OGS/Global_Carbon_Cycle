"""
04_validate.py
==============
Validation of the GLOBAL_MULTIYEAR_BGC_001_029 surface CO2 flux reconstruction
against the CMEMS MULTIOBS SOCAT-NN observation-based product.

All data cut off at 2024-01-01. All figures follow plot_style.py guidelines.

Changelog:
    v1.2.0 — Full rewrite.
    v1.2.1 — Cutoff, plot_style, 2-row ts, LOESS breakpoints, Fay split,
              spectra ticks, CPS coherence fix, cmocean maps.
    v1.2.2 — MULTIOBS colour fixed to pure red (#CC0000);
              ts legends moved below plots (2-col, bordered), larger a/b labels;
              map titles multiline, units removed from title, tighter cbar padding;
              LOESS breakpoints: -. vertical lines darker, B1-Bn x-axis labels,
              breakpoint dates printed to stdout, unique hardcoded segment colours,
              same breakpoints applied to MULTIOBS, labels below/above for rec/obs,
              same colour per segment, thicker trend lines;
              spectra: rightmost tick = actual Nyquist (fs/2 → 2 months);
              CPS: thresholds red+thick, same rightmost period fix, sharex removed;
              Fay: tighter hspace/wspace, larger font or smaller subplots.

Usage:
    python scripts/04_validate.py
    Requires: data/Time_Varying_Biomes.cmems.nc
"""

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mticker
from matplotlib.colors import ListedColormap
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from scipy.signal import welch, csd
import xarray as xr

sys.path.insert(0, str(Path(__file__).parent.parent))
import config as cfg

# ── optional imports ──────────────────────────────────────────────────────────
try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    HAS_CARTOPY = True
except ImportError:
    HAS_CARTOPY = False
    print("[warn] cartopy not found.")

try:
    import cmocean.cm as cmo
    HAS_CMOCEAN = True
except ImportError:
    HAS_CMOCEAN = False
    print("[warn] cmocean not found — falling back to matplotlib colormaps.")

try:
    import ruptures as rpt
    HAS_RUPTURES = True
except ImportError:
    HAS_RUPTURES = False
    print("[warn] ruptures not found — breakpoint detection disabled.")

# ── apply style ───────────────────────────────────────────────────────────────
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

# ── colormaps ─────────────────────────────────────────────────────────────────
CMAP_RMSD = cmo.amp     if HAS_CMOCEAN else "Reds"
CMAP_BIAS = cmo.balance if HAS_CMOCEAN else "RdBu_r"

# ── hardcoded main colours ────────────────────────────────────────────────────
C_REC  = "#0558AD"   # reconstruction — deep blue
C_OBS  = "#CD2859"   # reference      — deep rose/red
C_RMSD = "#CC79A7"   # Wong reddish-purple
C_CORR = "#009E73"   # Wong bluish-green

# ── hardcoded segment colours (skip C_REC and C_OBS) ─────────────────────────
# Visually distinct palette, none overlapping with the two main curve colours
SEG_COLORS = [
    "#E69F00",   # warm amber
    "#2E8B57",   # sea green
    "#9400D3",   # vivid violet
    "#FF6F20",   # burnt orange
    "#1ABDE8",   # cerulean
    "#8B0000",   # dark red
    "#3D3D3D",   # charcoal
    "#B8860B",   # dark goldenrod
]

# ── labels ───────────────────────────────────────────────────────────────────
REC_LABEL = "GLOBAL_MULTIYEAR_BGC_001_029"
OBS_LABEL = "MULTIOBS (SOCAT-NN)"

# ── cutoff ───────────────────────────────────────────────────────────────────
CUTOFF = np.datetime64("2024-01-01")


# ===========================================================================
# HELPERS
# ===========================================================================

def _clean(ax, minor=True):
    for spine in ax.spines.values():
        spine.set_linewidth(mpl.rcParams["axes.linewidth"])
    if minor:
        ax.xaxis.set_minor_locator(mticker.AutoMinorLocator())
        ax.yaxis.set_minor_locator(mticker.AutoMinorLocator())


def _label_panel(ax, lbl, x=0.015, y=0.97, fontsize=16):
    ax.text(x, y, lbl, transform=ax.transAxes,
            fontsize=fontsize, fontweight="bold", va="top", ha="left")


def compute_grid_cell_area(lat, lon):
    R       = cfg.EARTH_RADIUS_M
    d_lat   = float(np.abs(lat.diff("latitude").mean()))
    d_lon   = float(np.abs(lon.diff("longitude").mean()))
    area_1d = (R**2) * np.deg2rad(d_lon) * np.deg2rad(d_lat) * np.cos(np.deg2rad(lat))
    area_2d = area_1d.expand_dims({"longitude": lon}).transpose("latitude", "longitude")
    return area_2d.rename("cell_area")


def loess_smooth(x, y, frac=0.15):
    from statsmodels.nonparametric.smoothers_lowess import lowess
    return lowess(y, x, frac=frac, it=3, return_sorted=False)


def _set_period_ticks(ax, min_period_yr):
    """
    Apply custom period ticks (30d, 90d, 180d, 1y, 5y, 10y) filtered to
    those >= min_period_yr. Rightmost tick = actual shortest resolvable period.
    """
    ticks_yr = np.array([30/365.25, 90/365.25, 180/365.25, 1.0, 5.0, 10.0])
    labels   = ["30 d", "90 d", "180 d", "1 yr", "5 yr", "10 yr"]
    # Keep only ticks within the data range
    valid = ticks_yr >= min_period_yr * 0.99
    ticks_yr = ticks_yr[valid]
    labels   = [l for l, v in zip(labels, valid) if v]
    # Add the actual Nyquist period as rightmost tick if not already covered
    if min_period_yr < ticks_yr[0] * 0.99:
        ticks_yr = np.concatenate([[min_period_yr], ticks_yr])
        # Format: days if < 1 yr
        d = min_period_yr * 365.25
        labels = [f"{d:.0f} d" if d < 365 else f"{min_period_yr:.1f} yr"] + labels
    ax.set_xticks(ticks_yr)
    ax.set_xticklabels(labels, fontsize=10)
    ax.xaxis.set_minor_locator(mticker.NullLocator())


# ===========================================================================
# LOAD DATA
# ===========================================================================

def load_data():
    ds_flux = xr.open_dataset(cfg.DATA_DIR / "flux_3d.nc",          chunks="auto")
    ds_surf = xr.open_dataset(cfg.DATA_DIR / "processed_surface.nc", chunks="auto")
    ds_glob = xr.open_dataset(cfg.OUT_DIR  / "global_flux.nc")

    if "fgco2_obs" not in ds_surf:
        raise FileNotFoundError("fgco2_obs not found in processed_surface.nc.")

    fgco2_rec = ds_flux["fgco2"]
    fgco2_obs = ds_surf["fgco2_obs"]
    J_rec     = ds_glob["J_net_PgC"]

    common_time = fgco2_rec.time[np.isin(fgco2_rec.time.values, fgco2_obs.time.values)]
    common_time = common_time.sel(time=common_time < CUTOFF)

    fgco2_rec = fgco2_rec.sel(time=common_time)
    fgco2_obs = fgco2_obs.sel(time=common_time)
    J_rec     = J_rec.sel(time=common_time)

    cell_area = compute_grid_cell_area(fgco2_obs["latitude"], fgco2_obs["longitude"])
    F_obs_mol = (fgco2_obs.where(ds_surf["ocean_mask"] == 1) * cell_area
                 ).sum(dim=["latitude", "longitude"], skipna=True)
    J_obs = (F_obs_mol * cfg.MOL_C_TO_PG).rename("J_obs_PgC")

    return fgco2_rec, fgco2_obs, J_rec, J_obs, ds_surf


# ===========================================================================
# SKILL METRICS
# ===========================================================================

def pixel_rmsd(rec, obs):
    rmsd = np.sqrt(((rec - obs)**2).mean(dim="time"))
    rmsd.attrs = {"long_name": "RMSD", "units": "mol C m-2 yr-1"}
    return rmsd

def pixel_bias(rec, obs):
    bias = (rec - obs).mean(dim="time")
    bias.attrs = {"long_name": "Bias (rec − obs)", "units": "mol C m-2 yr-1"}
    return bias

def scalar_metrics(rec, obs):
    r = rec.values.ravel(); o = obs.values.ravel()
    ok = np.isfinite(r) & np.isfinite(o); r, o = r[ok], o[ok]
    corr, pval = stats.pearsonr(r, o)
    return {
        "bias_mol_C_m2_yr": float(np.mean(r - o)),
        "rmsd_mol_C_m2_yr": float(np.sqrt(np.mean((r - o)**2))),
        "pearson_r":        float(corr),
        "pearson_pval":     float(pval),
        "std_rec":          float(np.std(r)),
        "std_obs":          float(np.std(o)),
    }

def rolling_rmsd_corr(J_rec, J_obs, window=12):
    r = pd.Series(J_rec.values, index=pd.to_datetime(J_rec.time.values))
    o = pd.Series(J_obs.values, index=pd.to_datetime(J_obs.time.values))
    aligned  = pd.DataFrame({"rec": r, "obs": o}).dropna()
    diff     = aligned["rec"] - aligned["obs"]
    roll_rmsd = diff.pow(2).rolling(window).mean().apply(np.sqrt)
    roll_corr = aligned["rec"].rolling(window).corr(aligned["obs"])
    return roll_rmsd, roll_corr


# ===========================================================================
# FIGURE 1 — TIME SERIES + ROLLING RMSD/r  (2-row)
# ===========================================================================

def plot_timeseries(J_rec, J_obs):
    roll_rmsd, roll_corr = rolling_rmsd_corr(J_rec, J_obs, window=12)

    fig = plt.figure(figsize=(14, 9))
    gs  = gridspec.GridSpec(2, 1, hspace=0.45, height_ratios=[1.6, 1])

    # ── Row 1: time series ────────────────────────────────────────────────────
    ax1 = fig.add_subplot(gs[0])
    t_rec = pd.to_datetime(J_rec.time.values)
    t_obs = pd.to_datetime(J_obs.time.values)

    ax1.plot(t_rec, J_rec.values, color=C_REC, lw=1.5, alpha=0.9, label=REC_LABEL)
    ax1.plot(t_obs, J_obs.values, color=C_OBS, lw=1.5, alpha=0.9, ls="--",
             label=f"{OBS_LABEL}  (reference)")

    ax1.set_ylabel("Net ocean CO₂ uptake  [Pg C yr⁻¹]")
    ax1.set_title("Global net ocean CO₂ uptake — reconstruction vs reference")
    # Legend below the axes, 2 columns, bordered
    ax1.legend(
        loc="upper center", bbox_to_anchor=(0.5, -0.07),
        ncol=2, frameon=True, framealpha=0.9,
        edgecolor="#555555", fontsize=11,
    )
    _clean(ax1)
    _label_panel(ax1, "a", fontsize=18)
    plt.setp(ax1.get_xticklabels(), visible=False)

    # ── Row 2: rolling RMSD + r ───────────────────────────────────────────────
    ax2  = fig.add_subplot(gs[1], sharex=ax1)
    ax2r = ax2.twinx()

    ax2.plot(roll_rmsd.index, roll_rmsd.values,
             color=C_RMSD, lw=1.8, label="Rolling RMSD (12-month)")
    ax2r.plot(roll_corr.index, roll_corr.values,
              color=C_CORR, lw=1.8, ls="-", label="Rolling Pearson r (12-month)")

    ax2.set_ylabel("RMSD  [Pg C yr⁻¹]", color=C_RMSD)
    ax2.tick_params(axis="y", labelcolor=C_RMSD)
    ax2r.set_ylabel("Pearson  r", color=C_CORR)
    ax2r.tick_params(axis="y", labelcolor=C_CORR)
    ax2r.set_ylim(-0.1, 1.05)
    ax2r.axhline(0, color="gray", lw=0.7, ls=":")
    ax2.set_xlabel("Year")
    ax2.set_title("Rolling annual skill metrics (12-month window)")

    h1, l1 = ax2.get_legend_handles_labels()
    h2, l2 = ax2r.get_legend_handles_labels()
    ax2.legend(
        h1 + h2, l1 + l2,
        loc="upper center", bbox_to_anchor=(0.5, -0.22),
        ncol=2, frameon=True, framealpha=0.9,
        edgecolor="#555555", fontsize=11,
    )
    _clean(ax2)
    _label_panel(ax2, "b", fontsize=18)

    fig.subplots_adjust(bottom=0.15, hspace=0.30)
    out = cfg.FIG_DIR / "fig_validation_ts.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"[fig] {out}")


# ===========================================================================
# FIGURE 2 — LOESS + BREAKPOINT TREND ANALYSIS
# ===========================================================================

def plot_loess(J_rec, J_obs, frac=0.15):
    """
    Raw monthly (faint) + LOESS smooth (bold).
    Breakpoints detected on reconstruction LOESS, then applied to both curves.
    Same segment index → same colour. Trend labels below for rec, above for obs.
    Vertical breakpoint lines are -. style, dark gray.
    Breakpoint dates printed to stdout.
    """
    fig, ax = plt.subplots(figsize=(15, 6))

    t_rec = pd.to_datetime(J_rec.time.values)
    t_obs = pd.to_datetime(J_obs.time.values)
    x_num = np.arange(len(t_rec), dtype=float)

    rec_v = J_rec.values.copy()
    obs_v = J_obs.values.copy()
    ok_r  = np.isfinite(rec_v)
    ok_o  = np.isfinite(obs_v)

    # Raw monthly — faint
    ax.plot(t_rec, rec_v, color=C_REC, lw=0.7, alpha=0.22)
    ax.plot(t_obs, obs_v, color=C_OBS, lw=0.7, alpha=0.22, ls="--")

    # LOESS smooth
    sm_rec = loess_smooth(x_num[ok_r], rec_v[ok_r], frac=frac)
    sm_obs = loess_smooth(x_num[ok_o], obs_v[ok_o], frac=frac)

    ax.plot(t_rec[ok_r], sm_rec, color=C_REC, lw=2.5,
            label=f"{REC_LABEL}  (LOESS)")
    ax.plot(t_obs[ok_o], sm_obs, color=C_OBS, lw=2.5, ls="--",
            label=f"{OBS_LABEL}  (LOESS)")

    # ── Breakpoints ───────────────────────────────────────────────────────────
    if not HAS_RUPTURES:
        print("[warn] ruptures not installed — breakpoint detection skipped.")
        ax.set_xlabel("Year")
        ax.set_ylabel("Net ocean CO₂ uptake  [Pg C yr⁻¹]")
        ax.set_title("Multi-year variability — LOESS smooth\nPale: raw monthly  |  Bold: LOESS")
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13),
                  ncol=2, frameon=True, edgecolor="#555555")
        _clean(ax)
        fig.savefig(cfg.FIG_DIR / "fig_validation_loess.png")
        plt.close(fig)
        return

    signal = sm_rec.reshape(-1, 1)
    algo   = rpt.Pelt(model="rbf").fit(signal)
    bkps   = algo.predict(pen=3)   # pen=3 → moderate sensitivity
    # bkps includes the final index (len(signal)) — that's the "end" sentinel
    segments = [0] + bkps          # e.g. [0, 48, 120, 372] where last = len

    print(f"\n[loess] Breakpoints detected (n={len(bkps)-1}):")
    bp_times = []
    for bp in bkps[:-1]:   # skip the last sentinel
        bp_t = t_rec[ok_r][bp - 1]
        bp_times.append(bp_t)
        print(f"  B{len(bp_times)}: {bp_t.strftime('%Y-%m')}")

    # x-tick positions for breakpoints (stored for labelling)
    bp_tick_dates  = []
    bp_tick_labels = []

    for seg_idx in range(len(segments) - 1):
        i0 = segments[seg_idx]
        i1 = segments[seg_idx + 1]
        col = SEG_COLORS[seg_idx % len(SEG_COLORS)]

        # ── Reconstruction segment ────────────────────────────────────────────
        t_seg_r  = t_rec[ok_r][i0:i1]
        y_seg_r  = sm_rec[i0:i1]
        x_seg    = np.arange(i0, i1, dtype=float)

        if len(t_seg_r) >= 3:
            sl_r, ic_r, *_ = stats.linregress(x_seg, y_seg_r)
            sl_r_yr = sl_r * 12
            y_fit_r = ic_r + sl_r * x_seg
            ax.plot(t_seg_r, y_fit_r, color=col, lw=2.2, ls=":", alpha=0.9)

            # Arrow originates from the midpoint of the TREND LINE (not LOESS)
            mid_idx = len(x_seg) // 2
            t_mid_r   = t_seg_r[mid_idx]
            y_trend_r = y_fit_r[mid_idx]   # point on the trend line
            ax.annotate(
                f"{sl_r_yr:+.3f} Pg C yr⁻²",
                xy=(t_mid_r, y_trend_r),
                xytext=(0, -22), textcoords="offset points",
                fontsize=10, fontweight="bold", color=col,
                ha="center", va="top",
                arrowprops=dict(arrowstyle="-", color=col, lw=1.8),
            )

        # ── MULTIOBS segment (same breakpoints, same colour) ─────────────────
        t_seg_o  = t_obs[ok_o][i0:min(i1, len(t_obs[ok_o]))]
        y_seg_o  = sm_obs[i0:min(i1, len(sm_obs))]
        x_seg_o  = np.arange(i0, i0 + len(t_seg_o), dtype=float)

        if len(t_seg_o) >= 3:
            sl_o, ic_o, *_ = stats.linregress(x_seg_o, y_seg_o)
            sl_o_yr = sl_o * 12
            y_fit_o = ic_o + sl_o * x_seg_o
            ax.plot(t_seg_o, y_fit_o, color=col, lw=2.2, ls=":", alpha=0.9)

            # Arrow originates from the midpoint of the TREND LINE
            mid_idx_o = len(x_seg_o) // 2
            t_mid_o   = t_seg_o[mid_idx_o]
            y_trend_o = y_fit_o[mid_idx_o]   # point on the trend line
            ax.annotate(
                f"{sl_o_yr:+.3f} Pg C yr⁻²",
                xy=(t_mid_o, y_trend_o),
                xytext=(0, 22), textcoords="offset points",
                fontsize=10, fontweight="bold", color=col,
                ha="center", va="bottom",
                arrowprops=dict(arrowstyle="-", color=col, lw=1.8),
            )

        # ── Breakpoint vertical line (all but last segment) ───────────────────
        if seg_idx < len(segments) - 2:
            bp_date = t_rec[ok_r][i1 - 1]
            ax.axvline(bp_date, color="#444444", lw=1.4, ls="-.", alpha=0.85)
            bp_tick_dates.append(bp_date)
            bp_tick_labels.append(f"B{seg_idx + 1}")

    # Add B1-Bn labels on the x-axis
    ax2_bp = ax.twiny()
    ax2_bp.set_xlim(ax.get_xlim())
    ax2_bp.set_xticks(bp_tick_dates)
    ax2_bp.set_xticklabels(bp_tick_labels, fontsize=9, fontweight="bold",
                            color="#444444")
    ax2_bp.tick_params(axis="x", length=6, width=1.2,
                       direction="in", color="#444444")
    ax2_bp.spines["top"].set_visible(True)
    ax2_bp.spines["top"].set_linewidth(0)   # hide spine, keep ticks

    ax.set_xlabel("Year")
    ax.set_ylabel("Net ocean CO₂ uptake  [Pg C yr⁻¹]")
    ax.set_title(
        "Multi-year variability — LOESS smooth with piecewise trend\n"
        "Pale: raw monthly  |  Bold: LOESS  |  Dotted: segment trend  "
        "|  B: breakpoints"
    )
    ax.legend(
        loc="upper center", bbox_to_anchor=(0.5, -0.12),
        ncol=2, frameon=True, edgecolor="#555555", fontsize=11,
    )
    _clean(ax, minor=False)

    out = cfg.FIG_DIR / "fig_validation_loess.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"[fig] {out}")


# ===========================================================================
# FIGURE 3 — VALIDATION MAP
# ===========================================================================

def plot_validation_map(rmsd, bias):
    rmsd_levels = np.linspace(0,                cfg.VAL_RMSD_MAX, 21)
    bias_levels = np.linspace(cfg.VAL_BIAS_MIN, cfg.VAL_BIAS_MAX, 21)

    fig = plt.figure(figsize=(11, 10))
    gs  = gridspec.GridSpec(2, 1, hspace=0.20,
                            left=0.03, right=0.93, top=0.93, bottom=0.04)
    cax_rmsd = fig.add_axes([0.945, 0.55, 0.015, 0.35])
    cax_bias = fig.add_axes([0.945, 0.08, 0.015, 0.35])

    panels = [
        (rmsd, rmsd_levels, CMAP_RMSD, cax_rmsd,
         f"RMSD\n{REC_LABEL}\nvs {OBS_LABEL}",
         "RMSD  [mol C m⁻² yr⁻¹]", "a"),
        (bias, bias_levels, CMAP_BIAS,  cax_bias,
         f"Bias (rec − obs)\n{REC_LABEL}\nvs {OBS_LABEL}",
         "Bias  [mol C m⁻² yr⁻¹]", "b"),
    ]

    for row, (da, levels, cmap, cax, title, cbar_label, lbl) in enumerate(panels):
        lat = da["latitude"].values
        lon = da["longitude"].values
        Z   = da.values

        if HAS_CARTOPY:
            ax = fig.add_subplot(gs[row], projection=ccrs.Robinson())
            cf = ax.contourf(lon, lat, Z, levels=levels, cmap=cmap,
                             extend="both", transform=ccrs.PlateCarree())
            ax.add_feature(cfeature.LAND,      facecolor="#d3d3d3", zorder=3)
            ax.add_feature(cfeature.COASTLINE, lw=0.3,              zorder=4)
            ax.set_global()
        else:
            ax = fig.add_subplot(gs[row])
            cf = ax.contourf(lon, lat, Z, levels=levels, cmap=cmap, extend="both")
            ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude")

        ax.set_title(title, pad=5, fontsize=10, linespacing=1.4)
        cb = plt.colorbar(cf, cax=cax, label=cbar_label)
        cb.ax.tick_params(labelsize=10, width=1.0, length=3)

    out = cfg.FIG_DIR / "fig_validation_map.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"[fig] {out}")


# ===========================================================================
# FAY DOMAIN HELPERS
# ===========================================================================

def load_fay_biomes():
    fay_file = cfg.DATA_DIR / "Time_Varying_Biomes.cmems.nc"
    if not fay_file.exists():
        raise FileNotFoundError(f"Fay biome file not found: {fay_file}")
    ds = xr.open_dataset(fay_file)
    mask = ds["MeanBiomes"].isel(year=0).values if "year" in ds["MeanBiomes"].dims \
           else ds["MeanBiomes"].values
    lon = ds["lon"].values[0, :] if ds["lon"].ndim == 2 else ds["lon"].values
    lat = ds["lat"].values[:, 0] if ds["lat"].ndim == 2 else ds["lat"].values
    return mask.astype(float), lon, lat


def domain_mean_timeseries(flux, biome_mask, domain_id):
    mask_2d = (biome_mask == domain_id).astype(float)
    mask_da = xr.DataArray(mask_2d, dims=["latitude", "longitude"],
                           coords={"latitude": flux["latitude"],
                                   "longitude": flux["longitude"]})
    cell_area = compute_grid_cell_area(flux["latitude"], flux["longitude"])
    weight    = cell_area * mask_da
    total_w   = float(weight.sum())
    if total_w == 0:
        return np.full(flux.sizes["time"], np.nan)
    return ((flux * weight).sum(dim=["latitude", "longitude"],
                                skipna=True) / total_w).values


def compute_climatology(ts, time):
    da = pd.Series(ts, index=time)
    return da.groupby(da.index.month).mean().values


def _draw_domain_inset(fig, ax, biome_mask, lon_fay, lat_fay, domain_id):
    pos     = ax.get_position()
    iw, ih  = 0.19 * pos.width, 0.26 * pos.height
    ix      = pos.x0 + 0.79 * pos.width
    iy      = pos.y0 + 0.70 * pos.height
    lon2d, lat2d = np.meshgrid(lon_fay, lat_fay)
    masked  = np.where(biome_mask == float(domain_id), 1, np.nan)
    if HAS_CARTOPY:
        axins = fig.add_axes([ix, iy, iw, ih], projection=ccrs.PlateCarree())
        axins.pcolormesh(lon2d, lat2d, masked,
                         cmap=ListedColormap(["darkred"]),
                         transform=ccrs.PlateCarree())
        axins.add_feature(cfeature.LAND, color="lightgray", zorder=2)
        axins.coastlines(linewidth=0.25)
        axins.set_global()
    else:
        axins = fig.add_axes([ix, iy, iw, ih])
        axins.contourf(lon2d, lat2d, masked, colors=["darkred"])
        axins.axis("off")


def _fay_global_biome_map(fig, subplot_spec, biome_mask, lon_fay, lat_fay, n=17):
    cmap_jet = plt.cm.jet
    bounds   = np.arange(0.5, n + 1.5, 1)
    norm     = mpl.colors.BoundaryNorm(bounds, cmap_jet.N)
    lon2d, lat2d = np.meshgrid(lon_fay, lat_fay)
    proj = ccrs.Robinson() if HAS_CARTOPY else None
    ax   = fig.add_subplot(subplot_spec, projection=proj)
    if HAS_CARTOPY:
        mm = ax.pcolormesh(lon2d, lat2d, biome_mask, cmap=cmap_jet,
                           norm=norm, transform=ccrs.PlateCarree())
        ax.add_feature(cfeature.LAND,      color="gray")
        ax.add_feature(cfeature.COASTLINE, linewidth=0.4)
        ax.set_global()
    else:
        mm = ax.contourf(lon_fay, lat_fay, biome_mask,
                         levels=bounds, cmap=cmap_jet, norm=norm)
    ax.set_title("Fay (2014) Mean Biomes", fontsize=9)
    plt.colorbar(mm, ax=ax, shrink=0.8, ticks=np.arange(1, n + 1),
                 orientation="horizontal", pad=0.06)
    return ax


# ===========================================================================
# FIGURE 4a — FAY DOMAIN TIME SERIES
# ===========================================================================

def plot_fay_ts(fgco2_rec, fgco2_obs):
    print("[load] Fay biome mask ...")
    biome_mask, lon_fay, lat_fay = load_fay_biomes()
    n_domains = 17

    # Per-domain colour for the reconstruction curve (jet, same as biome map)
    cmap_jet   = plt.cm.jet
    dom_colors = {i: cmap_jet((i - 1) / (n_domains - 1)) for i in range(1, n_domains + 1)}

    fgco2_rec_r = fgco2_rec.sel(latitude=lat_fay, longitude=lon_fay, method="nearest")
    fgco2_obs_r = fgco2_obs.sel(latitude=lat_fay, longitude=lon_fay, method="nearest")
    time_pd     = pd.to_datetime(fgco2_rec_r.time.values)

    fig = plt.figure(figsize=(20, 24))
    gs  = gridspec.GridSpec(6, 3, figure=fig,
                            hspace=0.32, wspace=0.25,
                            top=0.93, bottom=0.04,
                            left=0.06, right=0.97)

    domain_metrics = {}

    for i in range(1, n_domains + 1):
        print(f"  [fay-ts] domain {i:2d}/{n_domains} ...")
        row, col = divmod(i - 1, 3)
        ax = fig.add_subplot(gs[row, col])

        ts_rec = domain_mean_timeseries(fgco2_rec_r, biome_mask, i)
        ts_obs = domain_mean_timeseries(fgco2_obs_r, biome_mask, i)

        dc = dom_colors[i]
        ax.plot(time_pd, ts_rec, color=dc,   lw=1.4, alpha=0.95,
                label=REC_LABEL)
        ax.plot(time_pd, ts_obs, color=C_OBS, lw=1.4, alpha=0.95, ls="--",
                label=OBS_LABEL)
        ax.set_ylabel("mol C m⁻² yr⁻¹", fontsize=11)
        ax.tick_params(axis="x", labelsize=10, rotation=30)
        ax.tick_params(axis="y", labelsize=10)
        ax.set_title(f"Domain {i}", fontsize=12, pad=5, fontweight="bold")
        _clean(ax, minor=False)

        ok = np.isfinite(ts_rec) & np.isfinite(ts_obs)
        if ok.sum() > 5:
            r, _ = stats.pearsonr(ts_rec[ok], ts_obs[ok])
            domain_metrics[i] = {
                "bias":      float(np.mean(ts_rec[ok] - ts_obs[ok])),
                "rmsd":      float(np.sqrt(np.mean((ts_rec[ok] - ts_obs[ok])**2))),
                "pearson_r": float(r),
            }

    _fay_global_biome_map(fig, gs[5, 2], biome_mask, lon_fay, lat_fay)

    # Legend: use global colours for the two lines
    handles = [
        plt.Line2D([0],[0], color="gray",  lw=2.0,
                   label=f"{REC_LABEL}  (domain colour)"),
        plt.Line2D([0],[0], color=C_OBS, lw=2.0, ls="--", label=OBS_LABEL),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=2,
               fontsize=11, bbox_to_anchor=(0.5, 0.975),
               framealpha=0.9, edgecolor="#555555")
    fig.suptitle(
        f"Fay (2014) domain validation — monthly time series\n"
        f"{REC_LABEL} vs {OBS_LABEL}",
        fontsize=14, y=0.995,
    )

    out = cfg.FIG_DIR / "fig_fay_ts.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"[fig] {out}")
    return domain_metrics


# ===========================================================================
# FIGURE 4b — FAY DOMAIN CLIMATOLOGIES
# ===========================================================================

def plot_fay_clim(fgco2_rec, fgco2_obs):
    print("[load] Fay biome mask (clim) ...")
    biome_mask, lon_fay, lat_fay = load_fay_biomes()
    n_domains = 17

    cmap_jet   = plt.cm.jet
    dom_colors = {i: cmap_jet((i - 1) / (n_domains - 1)) for i in range(1, n_domains + 1)}

    fgco2_rec_r = fgco2_rec.sel(latitude=lat_fay, longitude=lon_fay, method="nearest")
    fgco2_obs_r = fgco2_obs.sel(latitude=lat_fay, longitude=lon_fay, method="nearest")
    time_pd     = pd.to_datetime(fgco2_rec_r.time.values)

    fig = plt.figure(figsize=(20, 24))
    gs  = gridspec.GridSpec(6, 3, figure=fig,
                            hspace=0.32, wspace=0.25,
                            top=0.93, bottom=0.04,
                            left=0.06, right=0.97)

    months_abbr = ["Jan","Feb","Mar","Apr","May","Jun",
                   "Jul","Aug","Sep","Oct","Nov","Dec"]
    m_idx = np.arange(1, 13)

    for i in range(1, n_domains + 1):
        print(f"  [fay-clim] domain {i:2d}/{n_domains} ...")
        row, col = divmod(i - 1, 3)
        ax = fig.add_subplot(gs[row, col])

        ts_rec   = domain_mean_timeseries(fgco2_rec_r, biome_mask, i)
        ts_obs   = domain_mean_timeseries(fgco2_obs_r, biome_mask, i)
        clim_rec = compute_climatology(ts_rec, time_pd)
        clim_obs = compute_climatology(ts_obs, time_pd)

        def _std(ts, time):
            s = pd.Series(ts, index=time)
            return s.groupby(s.index.month).std().values

        std_rec = _std(ts_rec, time_pd)
        std_obs = _std(ts_obs, time_pd)

        dc = dom_colors[i]
        ax.plot(m_idx, clim_rec, color=dc,    lw=2.0, marker="o", ms=5)
        ax.plot(m_idx, clim_obs, color=C_OBS, lw=2.0, marker="s", ms=5, ls="--")
        ax.fill_between(m_idx, clim_rec - std_rec, clim_rec + std_rec,
                        color=dc,    alpha=0.15)
        ax.fill_between(m_idx, clim_obs - std_obs, clim_obs + std_obs,
                        color=C_OBS, alpha=0.15)

        ax.set_xticks(m_idx)
        ax.set_xticklabels(months_abbr, fontsize=9, rotation=45)
        ax.set_ylabel("mol C m⁻² yr⁻¹", fontsize=11)
        ax.tick_params(axis="y", labelsize=10)
        ax.set_title(f"Domain {i}", fontsize=12, pad=5, fontweight="bold")
        _clean(ax, minor=False)

    _fay_global_biome_map(fig, gs[5, 2], biome_mask, lon_fay, lat_fay)

    handles = [
        plt.Line2D([0],[0], color="gray", lw=2.0, marker="o", ms=5,
                   label=f"{REC_LABEL}  (domain colour)"),
        plt.Line2D([0],[0], color=C_OBS, lw=2.0, marker="s", ms=5, ls="--",
                   label=OBS_LABEL),
        mpatches.Patch(color="gray", alpha=0.2, label="±1σ"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=3,
               fontsize=11, bbox_to_anchor=(0.5, 0.975),
               framealpha=0.9, edgecolor="#555555")
    fig.suptitle(
        f"Fay (2014) domain validation — climatological seasonal cycle\n"
        f"{REC_LABEL} vs {OBS_LABEL}",
        fontsize=14, y=0.995,
    )

    out = cfg.FIG_DIR / "fig_fay_clim.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"[fig] {out}")


# ===========================================================================
# FIGURE 5 — POWER SPECTRAL DENSITY
# ===========================================================================

def plot_spectra(J_rec, J_obs):
    fs  = 12.0   # samples per year
    rec = J_rec.values; obs = J_obs.values
    ok  = np.isfinite(rec) & np.isfinite(obs)
    rec = rec[ok] - np.mean(rec[ok])
    obs = obs[ok] - np.mean(obs[ok])

    nperseg = max(len(rec) // 3, 24)
    f_r, P_r = welch(rec, fs=fs, nperseg=nperseg, scaling="density")
    f_o, P_o = welch(obs, fs=fs, nperseg=nperseg, scaling="density")

    # Skip DC; period in years
    per_r = 1.0 / f_r[1:];  P_r = P_r[1:]
    per_o = 1.0 / f_o[1:];  P_o = P_o[1:]

    # Nyquist = fs/2 → period = 2/fs years
    nyq_yr = 2.0 / fs   # = 2 months in years

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.loglog(per_r, P_r, color=C_REC, lw=2.0, label=REC_LABEL)
    ax.loglog(per_o, P_o, color=C_OBS, lw=2.0, ls="--", label=OBS_LABEL)

    ax.set_xlabel("Period")
    ax.set_ylabel("PSD  [(Pg C yr⁻¹)² yr]")
    ax.set_title(
        "Power spectral density\n"
        "Global net ocean CO₂ uptake — reconstruction vs reference"
    )
    ax.legend(framealpha=0.85)
    ax.grid(True, which="major", color="#999999", ls="--", lw=0.8, alpha=0.8)
    ax.grid(True, which="minor", color="#CCCCCC", ls=":", lw=0.4, alpha=0.5)
    ax.invert_xaxis()
    ax.set_xlim(left=ax.get_xlim()[0], right=nyq_yr * 0.85)

    _set_period_ticks(ax, nyq_yr)
    _clean(ax, minor=False)

    out = cfg.FIG_DIR / "fig_spectra.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"[fig] {out}")


# ===========================================================================
# FIGURE 6 — CROSS-POWER SPECTRUM
# ===========================================================================

def plot_cps(J_rec, J_obs):
    from scipy.signal import coherence as sig_coherence

    fs  = 12.0
    rec = J_rec.values; obs = J_obs.values
    ok  = np.isfinite(rec) & np.isfinite(obs)
    rec = rec[ok] - np.mean(rec[ok])
    obs = obs[ok] - np.mean(obs[ok])

    n        = len(rec)
    nperseg  = max(n // 4, 24)
    noverlap = nperseg // 2
    nfft     = 2 ** int(np.ceil(np.log2(n)) + 3)

    f, Pxx = welch(obs, fs=fs, nperseg=nperseg, noverlap=noverlap,
                   nfft=nfft, scaling="density")
    _, Pyy = welch(rec, fs=fs, nperseg=nperseg, noverlap=noverlap,
                   nfft=nfft, scaling="density")
    _, Pxy = csd(obs, rec, fs=fs, nperseg=nperseg, noverlap=noverlap,
                 nfft=nfft, scaling="density")
    f_coh, Cxy = sig_coherence(obs, rec, fs=fs, nperseg=nperseg,
                                noverlap=noverlap, nfft=nfft)

    # Skip DC
    f     = f[1:];     Pxx  = Pxx[1:];  Pyy  = Pyy[1:];  Pxy  = Pxy[1:]
    f_coh = f_coh[1:]; Cxy  = Cxy[1:]

    period     = 1.0 / f
    period_coh = 1.0 / f_coh
    gain       = np.abs(Pxy) / Pxx
    phase_days = np.angle(Pxy) / (2 * np.pi * f) * 365.25

    nyq_yr = 2.0 / fs   # 2 months

    fig, axes = plt.subplots(3, 1, figsize=(11, 11))

    # ── Gain ─────────────────────────────────────────────────────────────────
    axes[0].semilogx(period, gain, color=C_REC, lw=1.8)
    axes[0].axhline(1, color="red", lw=2.0, ls="--")
    axes[0].set_ylabel("Gain")
    axes[0].set_title(
        f"{REC_LABEL}\n"
        f"vs {OBS_LABEL} — cross-power spectrum"
    )
    _clean(axes[0], minor=False)

    # ── Phase ─────────────────────────────────────────────────────────────────
    axes[1].semilogx(period, phase_days, color="#E69F00", lw=1.8)
    axes[1].axhline(0, color="red", lw=2.0, ls="--")
    axes[1].set_ylabel("Phase  [days]")
    _clean(axes[1], minor=False)

    # ── Coherence ─────────────────────────────────────────────────────────────
    axes[2].semilogx(period_coh, Cxy, color=C_CORR, lw=1.8)
    axes[2].axhline(0.5, color="red", lw=2.0, ls="--", label="0.5 threshold")
    axes[2].set_ylim(0, 1.05)
    axes[2].set_ylabel("Magnitude-squared coherence")
    axes[2].set_xlabel("Period")
    axes[2].legend(fontsize=10, framealpha=0.85)
    _clean(axes[2], minor=False)

    # Period ticks + darker major grid for all panels
    for ax, per_arr in zip(axes, [period, period, period_coh]):
        ax.grid(True, which="major", color="#999999", ls="--", lw=0.8, alpha=0.8)
        xlim_right = nyq_yr * 0.85
        ax.set_xlim(left=max(per_arr) * 1.05, right=xlim_right)
        _set_period_ticks(ax, nyq_yr)

    fig.tight_layout(h_pad=1.5)
    out = cfg.FIG_DIR / "fig_cps.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"[fig] {out}")


# ===========================================================================
# MAIN
# ===========================================================================

def main():
    print("[load] Data (cutoff: 2024-01-01) ...")
    fgco2_rec, fgco2_obs, J_rec, J_obs, ds_surf = load_data()
    print(f"  Record: {str(J_rec.time.values[0])[:10]} → "
          f"{str(J_rec.time.values[-1])[:10]}  ({len(J_rec.time)} months)")

    print("[compute] Global scalar metrics ...")
    metrics = scalar_metrics(fgco2_rec, fgco2_obs)
    print("\n  Global validation metrics:")
    for k, v in metrics.items():
        print(f"    {k:30s}: {v:.4f}")

    print("\n[compute] Pixel RMSD and bias ...")
    rmsd = pixel_rmsd(fgco2_rec, fgco2_obs)
    bias = pixel_bias(fgco2_rec, fgco2_obs)

    print("\n[plot] Time series + rolling skill ...")
    plot_timeseries(J_rec, J_obs)

    print("[plot] LOESS + breakpoint trend ...")
    plot_loess(J_rec, J_obs, frac=0.15)

    print("[plot] Validation map ...")
    plot_validation_map(rmsd, bias)

    print("[plot] Fay domain — time series ...")
    domain_metrics = plot_fay_ts(fgco2_rec, fgco2_obs)

    print("[plot] Fay domain — climatologies ...")
    plot_fay_clim(fgco2_rec, fgco2_obs)

    print("[plot] Power spectra ...")
    plot_spectra(J_rec, J_obs)

    print("[plot] Cross-power spectrum ...")
    plot_cps(J_rec, J_obs)

    all_metrics = {"global": metrics}
    all_metrics.update({f"domain_{k}": v for k, v in domain_metrics.items()})
    pd.DataFrame(all_metrics).T.to_csv(cfg.OUT_DIR / "validation_metrics.csv")
    print(f"\n[save] validation_metrics.csv")
    print(f"[done] Validation complete. Figures → {cfg.FIG_DIR}\n")


if __name__ == "__main__":
    main()
