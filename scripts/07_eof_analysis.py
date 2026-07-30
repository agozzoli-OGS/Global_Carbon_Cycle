"""
07_eof_analysis.py
==================
EOF (Empirical Orthogonal Function) analysis of the surface air-sea CO2
flux field for Stage 1.

Active analysis (v1.4.3):

    Tag           Temporal res.   Field
    annual_full   Annual mean     Full field (no mean removal)

The monthly_full, monthly_anom and annual_anom runs were deprecated in
v1.4.3 after evaluation showed none produced physically relevant results
beyond the annual full-field decomposition.

Figures produced:

    fig_eof_annual_full_scree.png
        Scree plot: variance explained per EOF (bars, first 10) +
        North et al. (1982) sampling error bars + cumulative variance line.
        Degenerate mode pairs highlighted in red.

    fig_eof_annual_full_modes.png
        3 rows x 2 columns composite:
            Left  : EOF spatial pattern (Robinson projection, white-centred
                    diverging colourmap, physical units [mol C m-2 yr-1])
            Right : PC timeseries (unit-std, dimensionless) + ONI and SAM
                    z-scores overlaid if index files present, Spearman rho.
        Rows 1-3 = EOF/PC 1, 2, 3.

Scaling convention
------------------
Unit-norm EOFs (standard oceanographic convention):
    1. SVD on area-weighted matrix: X_w = U Sigma V^T
    2. Unweight: EOF_unit = V^T / w_flat  (unit L2-norm, physical space)
    3. Project:  PC_raw = X_unweighted @ EOF_unit^T
    4. Normalise: PC_norm = PC_raw / std(PC_raw)  (dimensionless)
    5. Scale EOF: EOF_scaled = EOF_unit * std(PC_raw)
    -> EOF_scaled x PC_norm reconstructs the field in physical units.

Area weighting: sqrt(cos(lat)) per pixel before SVD, undone after.
North et al. (1982): delta_lambda = lambda * sqrt(2/N); degenerate pairs
shown as red bars on scree plot.

Usage
-----
    python scripts/07_eof_analysis.py

Changelog
---------
    v1.4.0 -- Initial implementation (annual anomaly only).
    v1.4.1 -- Four analyses; fixed PC scaling; 3x2 composite figure.
    v1.4.2 -- Fixed vmax NaN crash; suppressed empty-slice RuntimeWarning.
    v1.4.3 -- Deprecated monthly_full, monthly_anom, annual_anom.
              Only annual_full retained.
"""

import sys
from pathlib import Path
from dataclasses import dataclass

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mticker
import seaborn as sns
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent.parent))
import config as cfg

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    HAS_CARTOPY = True
except ImportError:
    HAS_CARTOPY = False
    print("[warn] cartopy not found -- maps will use plain axes.")

try:
    import cmocean.cm as cmo
    HAS_CMOCEAN = True
except ImportError:
    HAS_CMOCEAN = False

# ── Style ──────────────────────────────────────────────────────────────────────
sns.set_theme(
    style="ticks", context="paper", font_scale=1.1,
    rc={
        "axes.spines.top":    True,  "axes.spines.right":  True,
        "axes.linewidth":     1.5,
        "xtick.direction":    "inout", "ytick.direction":    "inout",
        "xtick.major.width":  1.2,    "ytick.major.width":  1.2,
        "xtick.minor.width":  0.84,   "ytick.minor.width":  0.84,
        "xtick.major.size":   5,      "ytick.major.size":   5,
        "xtick.minor.size":   3,      "ytick.minor.size":   3,
        "xtick.top":          True,   "ytick.right":        True,
        "font.size":          11,     "axes.labelsize":     11,
        "xtick.labelsize":    10,     "ytick.labelsize":    10,
        "legend.fontsize":    9,      "axes.titlesize":     11,
        "lines.linewidth":    1.8,
        "axes.grid":          True,   "grid.color":         "#CCCCCC",
        "grid.linestyle":     "--",   "grid.linewidth":     0.5,
        "grid.alpha":         0.6,    "axes.axisbelow":     True,
        "savefig.dpi":        300,    "savefig.bbox":       "tight",
        "savefig.pad_inches": 0.05,   "figure.facecolor":   "white",
        "axes.facecolor":     "white","legend.frameon":     True,
        "legend.framealpha":  0.85,   "legend.edgecolor":   "#AAAAAA",
    }
)

T_START  = "1993-01-01"
T_END    = "2023-12-31"
N_EOFS   = 3
N_SCREE  = 10
CMAP_EOF = cmo.balance if HAS_CMOCEAN else "RdBu_r"
C_PC     = "#0558AD"
C_ONI    = "#CC2936"
C_SAM    = "#2E8B57"

# Module-level lat/lon set in main(), read in fig_modes()
_LAT = None
_LON = None


@dataclass
class EOFRun:
    tag:      str   # filename identifier
    label:    str   # title label
    temporal: str   # "monthly" | "annual"
    anomaly:  bool  # True = remove time mean


RUNS = [
    EOFRun("annual_full",  "Annual mean — full field", "annual",  False),
    # monthly_full, monthly_anom, annual_anom deprecated in v1.4.3:
    # none produced physically relevant results beyond the annual full-field run.
]


# ===========================================================================
# HELPERS
# ===========================================================================

def _clean(ax, minor=True):
    for spine in ax.spines.values():
        spine.set_linewidth(mpl.rcParams["axes.linewidth"])
    if minor:
        ax.xaxis.set_minor_locator(mticker.AutoMinorLocator())
        ax.yaxis.set_minor_locator(mticker.AutoMinorLocator())


def _draw_eof_map(ax, lon, lat, eof_2d, vmax, title):
    """Draw one EOF spatial pattern on ax (with or without cartopy)."""
    levels = np.linspace(-vmax, vmax, 21)
    lon2d, lat2d = np.meshgrid(lon, lat)

    if HAS_CARTOPY:
        cf = ax.contourf(lon2d, lat2d, eof_2d, levels=levels,
                         cmap=CMAP_EOF, extend="both",
                         transform=ccrs.PlateCarree())
        ax.add_feature(cfeature.LAND,      facecolor="#d3d3d3", zorder=3)
        ax.add_feature(cfeature.COASTLINE, lw=0.6,              zorder=4)
        ax.set_global()
    else:
        cf = ax.contourf(lon2d, lat2d, eof_2d, levels=levels,
                         cmap=CMAP_EOF, extend="both")

    cb = plt.colorbar(cf, ax=ax, orientation="horizontal",
                      shrink=0.75, pad=0.04, aspect=30)
    cb.set_label("mol C m\u207b\u00b2 yr\u207b\u00b9", fontsize=9)
    cb.ax.tick_params(labelsize=8, width=0.8, length=2)
    ax.set_title(title, fontsize=10, pad=4)


def _draw_pc(ax, t_dates, pc, var_frac_k, idx_clim, k):
    """Draw one PC timeseries + optional ONI/SAM overlay."""
    ax.fill_between(t_dates, pc, 0, alpha=0.15, color=C_PC)
    ax.plot(t_dates, pc, color=C_PC, lw=1.8,
            label=f"PC {k+1}  ({100*var_frac_k:.1f}%)")
    ax.axhline(0, color="black", lw=0.7, ls=":")

    if idx_clim is not None:
        oni = idx_clim["oni"].reindex(t_dates, method="nearest").values
        sam = idx_clim["sam"].reindex(t_dates, method="nearest").values
        ok  = np.isfinite(pc) & np.isfinite(oni) & np.isfinite(sam)
        if ok.sum() >= 5:
            rho_oni, p_oni = stats.spearmanr(pc[ok], oni[ok])
            rho_sam, p_sam = stats.spearmanr(pc[ok], sam[ok])
            oni_z = (oni - np.nanmean(oni)) / (np.nanstd(oni) + 1e-12)
            sam_z = (sam - np.nanmean(sam)) / (np.nanstd(sam) + 1e-12)
            ax.plot(t_dates, oni_z, color=C_ONI, lw=1.4, ls="--", alpha=0.8,
                    label=f"ONI  \u03c1={rho_oni:+.2f} p={p_oni:.3f}")
            ax.plot(t_dates, sam_z, color=C_SAM, lw=1.4, ls="-.", alpha=0.8,
                    label=f"SAM  \u03c1={rho_sam:+.2f} p={p_sam:.3f}")

    ax.set_ylabel("Amplitude  [\u03c3]", fontsize=9)
    ax.tick_params(axis="x", labelsize=9, rotation=30)
    ax.tick_params(axis="y", labelsize=9)
    ax.legend(loc="upper left", fontsize=8, edgecolor="#AAAAAA",
              handlelength=1.5, ncol=1)
    _clean(ax, minor=False)


# ===========================================================================
# DATA LOADING
# ===========================================================================

def load_flux_fields():
    """Load monthly and annual flux fields and time axes."""
    print("[load] flux_3d.nc ...")
    ds   = xr.open_dataset(cfg.DATA_DIR / "flux_3d.nc")
    flux = ds["fgco2"].sel(time=slice(T_START, T_END))

    lat = flux.latitude.values
    lon = flux.longitude.values

    t_monthly = pd.to_datetime(flux.time.values)
    monthly   = flux.values

    flux_ann = flux.resample(time="1YE").mean()
    years    = flux_ann.time.dt.year.values
    t_annual = pd.to_datetime([f"{y}-07-01" for y in years])
    annual   = flux_ann.values

    print(f"       Monthly : n={len(t_monthly)}  "
          f"({t_monthly[0].strftime('%Y-%m')} -> {t_monthly[-1].strftime('%Y-%m')})")
    print(f"       Annual  : n={len(t_annual)}  ({years[0]} -> {years[-1]})")
    print(f"       Grid    : {len(lat)} lat x {len(lon)} lon")

    return dict(monthly=monthly, annual=annual,
                t_monthly=t_monthly, t_annual=t_annual,
                lat=lat, lon=lon)


def load_climate_indices(t_monthly, t_annual):
    """Load ONI and SAM, reindexed to monthly and annual time axes."""
    oni_file = cfg.DATA_DIR / "oni_monthly.txt"
    sam_file = cfg.DATA_DIR / "sam_monthly.txt"

    if not oni_file.exists() or not sam_file.exists():
        print("[info] Climate index files not found -- PC plots will omit ONI/SAM.")
        return None

    def _parse(fpath, col):
        rows = []
        with open(fpath) as f:
            for line in f:
                parts = line.split()
                if len(parts) == 13:
                    try:
                        yr = int(parts[0])
                        for m, v in enumerate([float(x) for x in parts[1:]], 1):
                            if v != -99.9:
                                rows.append({"year": yr, "month": m, col: v})
                    except ValueError:
                        pass
        df = pd.DataFrame(rows)
        df["time"] = pd.to_datetime(
            df["year"].astype(str) + "-" + df["month"].astype(str).str.zfill(2)
        )
        return df.set_index("time")[col].sort_index()

    oni_m = _parse(oni_file, "oni")
    sam_m = _parse(sam_file, "sam")

    # Monthly (reindex to flux time axis)
    oni_mon = oni_m.reindex(t_monthly, method="nearest")
    sam_mon = sam_m.reindex(t_monthly, method="nearest")

    # Annual (resample then reindex to mid-year dates)
    t_ye     = pd.date_range(T_START[:4], T_END[:4], freq="YE")
    oni_ye   = oni_m.resample("YE").mean().reindex(t_ye, method="nearest")
    sam_ye   = sam_m.resample("YE").mean().reindex(t_ye, method="nearest")
    oni_ann  = oni_ye.reindex(t_annual, method="nearest")
    sam_ann  = sam_ye.reindex(t_annual, method="nearest")

    return {
        "monthly": {"oni": oni_mon, "sam": sam_mon},
        "annual":  {"oni": oni_ann, "sam": sam_ann},
    }


# ===========================================================================
# EOF COMPUTATION
# ===========================================================================

def compute_eofs(data, lat):
    """
    Area-weighted SVD with unit-norm EOF scaling.

    Scaling: EOFs normalised to unit L2-norm, then scaled by std(PC)
    so that EOF_scaled x PC_norm reconstructs the field in physical units.
    PCs are dimensionless (unit variance). EOF maps show amplitude per
    one standard deviation of the PC [mol C m-2 yr-1].
    """
    n_time, n_lat, n_lon = data.shape

    # Area weights
    cos_lat = np.maximum(np.cos(np.deg2rad(lat)), 0.0)
    w_1d    = np.sqrt(cos_lat)
    w_2d    = np.broadcast_to(w_1d[:, None], (n_lat, n_lon))

    # Ocean mask — suppress empty-slice warning from all-NaN columns (e.g. polar rows)
    with np.errstate(all="ignore"):
        ocean_mask = np.isfinite(np.nanmean(data, axis=0))
    n_ocean    = ocean_mask.sum()
    print(f"       Ocean pixels : {n_ocean} / {n_lat * n_lon}")

    # Weighted matrix
    w_flat     = w_2d[ocean_mask]
    X          = data[:, ocean_mask]                          # unweighted
    X_weighted = np.nan_to_num(X * w_flat[None, :], nan=0.0)

    # SVD
    U, s, Vt = np.linalg.svd(X_weighted, full_matrices=False)

    # Variance
    eigenvalues = s ** 2
    var_frac    = eigenvalues / eigenvalues.sum()
    north_error = var_frac * np.sqrt(2.0 / n_time)

    # Unweight EOFs (unit L2-norm in physical space)
    eofs_unit = Vt / w_flat[None, :]                         # (k, n_ocean)

    # PCs via projection onto unit-norm EOFs
    pcs_raw   = X @ eofs_unit.T                              # (n_time, k)
    pc_std    = pcs_raw.std(axis=0)
    pc_std[pc_std == 0] = 1.0
    pcs_norm  = pcs_raw / pc_std[None, :]                   # unit variance

    # Scale EOFs by PC std -> physical units, reconstruction exact
    eofs_scaled = eofs_unit * pc_std[:, None]               # (k, n_ocean)

    # Rebuild 2-D grids
    eofs_2d = np.full((eofs_scaled.shape[0], n_lat, n_lon), np.nan)
    eofs_2d[:, ocean_mask] = eofs_scaled

    print("       Variance : " +
          "  ".join([f"EOF{i+1}={100*v:.1f}%"
                     for i, v in enumerate(var_frac[:min(5, len(var_frac))])]))

    return dict(eofs=eofs_2d, pcs=pcs_norm,
                var_frac=var_frac, north_error=north_error)


# ===========================================================================
# FIGURES
# ===========================================================================

def fig_scree(result, run):
    """Scree plot with North error bars and cumulative variance."""
    n     = min(N_SCREE, len(result["var_frac"]))
    idx   = np.arange(1, n + 1)
    frac  = result["var_frac"][:n] * 100
    err   = result["north_error"][:n] * 100
    cumul = np.cumsum(result["var_frac"][:n]) * 100

    fig, ax1 = plt.subplots(figsize=(9, 5))
    bars = ax1.bar(idx, frac, color=C_PC, alpha=0.70, zorder=3,
                   label="Variance explained")
    ax1.errorbar(idx, frac, yerr=err, fmt="none",
                 ecolor="#333333", elinewidth=1.3, capsize=4, zorder=4)

    for i in range(n - 1):
        if frac[i] - err[i] <= frac[i + 1] + err[i + 1]:
            for bi in [i, i + 1]:
                bars[bi].set_edgecolor("#CC2936")
                bars[bi].set_linewidth(2.0)

    ax1.set_xlabel("EOF")
    ax1.set_ylabel("Variance explained  [%]", color=C_PC)
    ax1.tick_params(axis="y", labelcolor=C_PC)
    ax1.set_xticks(idx)

    ax2 = ax1.twinx()
    ax2.plot(idx, cumul, color="#E69F00", lw=2.0, marker="o", ms=5,
             label="Cumulative variance")
    ax2.axhline(90, color="#E69F00", lw=0.8, ls=":", alpha=0.6)
    ax2.set_ylabel("Cumulative variance  [%]", color="#E69F00")
    ax2.tick_params(axis="y", labelcolor="#E69F00")
    ax2.set_ylim(0, 105)

    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, loc="center right",
               fontsize=9, edgecolor="#AAAAAA")
    ax1.set_title(
        f"Scree plot -- {run.label}\n"
        "Error bars = North et al. (1982)  |  Red outline = degenerate",
        fontsize=11,
    )
    _clean(ax1, minor=False)
    fig.tight_layout()

    out = cfg.FIG_DIR / f"fig_eof_{run.tag}_scree.png"
    fig.savefig(out, dpi=cfg.FIGURE_DPI)
    plt.close(fig)
    print(f"[fig] {out}")


def fig_modes(result, run, t_dates, idx_clim):
    """
    3-row x 2-col composite: left = EOF map, right = PC timeseries.
    """
    global _LAT, _LON
    eofs     = result["eofs"]
    pcs      = result["pcs"]
    var_frac = result["var_frac"]
    n_rows   = min(N_EOFS, eofs.shape[0])

    fig = plt.figure(figsize=(18, 4.8 * n_rows))
    gs  = gridspec.GridSpec(
        n_rows, 2, figure=fig,
        wspace=0.08, hspace=0.50,
        left=0.02, right=0.98,
        top=0.93,  bottom=0.05,
    )

    for k in range(n_rows):
        eof_map = eofs[k]
        pc      = pcs[:, k]
        finite_vals = eof_map[np.isfinite(eof_map)]
        vmax = float(np.percentile(np.abs(finite_vals), 97)) if len(finite_vals) > 0 else 1.0
        if not np.isfinite(vmax) or vmax == 0.0:
            vmax = 1.0

        # Left: map
        if HAS_CARTOPY:
            ax_map = fig.add_subplot(gs[k, 0], projection=ccrs.Robinson())
        else:
            ax_map = fig.add_subplot(gs[k, 0])

        _draw_eof_map(
            ax_map, _LON, _LAT, eof_map, vmax,
            title=f"EOF {k+1}  --  {100*var_frac[k]:.1f}% variance",
        )

        # Right: PC
        ax_pc = fig.add_subplot(gs[k, 1])
        _draw_pc(ax_pc, t_dates, pc, var_frac[k], idx_clim, k)
        ax_pc.set_title(f"PC {k+1}", fontsize=10, pad=4)

    fig.suptitle(
        f"EOF modes 1-{n_rows}  --  {run.label}\n"
        "Left: spatial pattern [mol C m\u207b\u00b2 yr\u207b\u00b9]  |  "
        "Right: PC (unit \u03c3) + ONI/SAM z-scores",
        fontsize=12, y=0.975,
    )

    out = cfg.FIG_DIR / f"fig_eof_{run.tag}_modes.png"
    fig.savefig(out, dpi=cfg.FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {out}")


# ===========================================================================
# MAIN
# ===========================================================================

def main():
    global _LAT, _LON
    cfg.FIG_DIR.mkdir(parents=True, exist_ok=True)

    fields = load_flux_fields()
    _LAT   = fields["lat"]
    _LON   = fields["lon"]

    print("[load] Climate indices ...")
    clim_all = load_climate_indices(fields["t_monthly"], fields["t_annual"])

    for run in RUNS:
        print(f"\n{'='*60}")
        print(f"[eof] {run.tag.upper()}  --  {run.label}")
        print(f"{'='*60}")

        if run.temporal == "monthly":
            raw      = fields["monthly"]
            t_dates  = fields["t_monthly"]
            idx_clim = clim_all["monthly"] if clim_all is not None else None
        else:
            raw      = fields["annual"]
            t_dates  = fields["t_annual"]
            idx_clim = clim_all["annual"] if clim_all is not None else None

        if run.anomaly:
            data = raw - np.nanmean(raw, axis=0, keepdims=True)
            print("       Anomaly: time-mean removed per pixel")
        else:
            data = raw.copy()
            print("       Full field: no mean removal")

        result = compute_eofs(data, _LAT)

        print("[fig] Scree ...")
        fig_scree(result, run)

        print("[fig] Modes composite ...")
        fig_modes(result, run, t_dates, idx_clim)

    print(f"\n[done] {len(RUNS)*2} figures saved to {cfg.FIG_DIR}\n")

if __name__ == "__main__":
    main()
