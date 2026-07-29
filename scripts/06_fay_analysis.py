"""
06_fay_analysis.py
==================
Fay (2014) biogeochemical domain analysis for Stage 1.

Separated from 05_plot_results.py because the 17-domain integration loop
is computationally heavy (each domain requires an area-weighted sum over the
full 3D flux field for every time step).

Run this script independently after 03_compute_flux.py has completed.
Requires: data/Time_Varying_Biomes.cmems.nc

Figures produced
----------------
    fig_fay_ts.png
        17-panel monthly timeseries of area-integrated flux per Fay domain
        [Pg C yr⁻¹] + global biome overview map.

    fig_fay_clim.png
        17-panel climatological seasonal cycle ±1σ per domain [Pg C yr⁻¹].

    fig_fay_trends.png
        17-panel Sen's slope per domain [Pg C yr⁻¹ per decade] with
        Mann-Kendall significance bars (filled = p ≤ 0.01, hatched = p > 0.01).

Usage
-----
    python scripts/06_fay_analysis.py

    fig_fay_climate_regression.png
        17-panel grid: per-domain Spearman ρ for J_net vs ONI and J_net vs
        SAM (separate bars), plus the multilinear combination ρ. Significance
        stars at α = 0.05. Final panel: biome overview map.

Changelog
---------
    v1.3.0 — Split from 05_plot_results.py for performance reasons.
              All Fay domain code moved here verbatim with minor cleanups.
    v1.3.1 — Added load_climate_indices() (ONI + SAM, shared logic with 05).
              Added fig_fay_climate_regression(): 17-panel per-domain Spearman
              correlation bars following Bellacicco et al. (2025).
"""

import sys
import requests
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
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

# ── Style ─────────────────────────────────────────────────────────────────────
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

REC_LABEL = "GLOBAL_MULTIYEAR_BGC_001_029"
T_START   = "1993-01-01"
T_END     = "2023-12-31"


# ===========================================================================
# HELPERS
# ===========================================================================

def _clean(ax, minor=True):
    for spine in ax.spines.values():
        spine.set_linewidth(mpl.rcParams["axes.linewidth"])
    if minor:
        ax.xaxis.set_minor_locator(mticker.AutoMinorLocator())
        ax.yaxis.set_minor_locator(mticker.AutoMinorLocator())


def compute_grid_cell_area(lat, lon):
    R       = cfg.EARTH_RADIUS_M
    d_lat   = float(np.abs(lat.diff("latitude").mean()))
    d_lon   = float(np.abs(lon.diff("longitude").mean()))
    area_1d = (R**2) * np.deg2rad(d_lon) * np.deg2rad(d_lat) * np.cos(np.deg2rad(lat))
    area_2d = area_1d.expand_dims({"longitude": lon}).transpose("latitude", "longitude")
    return area_2d.rename("cell_area")


def mann_kendall_pvalue(y):
    x  = np.arange(len(y))
    ok = np.isfinite(y)
    if ok.sum() < 5:
        return np.nan
    _, p = stats.kendalltau(x[ok], y[ok])
    return p


# ===========================================================================
# DATA LOADING
# ===========================================================================

def load_flux() -> xr.DataArray:
    ds = xr.open_dataset(cfg.DATA_DIR / "flux_3d.nc")
    return ds["fgco2"].sel(time=slice(T_START, T_END))


def load_fay_biomes():
    fay_file = cfg.DATA_DIR / "Time_Varying_Biomes.cmems.nc"
    if not fay_file.exists():
        raise FileNotFoundError(
            f"Fay biome file not found: {fay_file}\n"
            "Download and place in the data/ directory."
        )
    ds   = xr.open_dataset(fay_file)
    mask = (ds["MeanBiomes"].isel(year=0).values
            if "year" in ds["MeanBiomes"].dims
            else ds["MeanBiomes"].values)
    lon  = ds["lon"].values[0, :] if ds["lon"].ndim == 2 else ds["lon"].values
    lat  = ds["lat"].values[:, 0] if ds["lat"].ndim == 2 else ds["lat"].values
    return mask.astype(float), lon, lat


def domain_flux_timeseries(flux_da, biome_mask, domain_id, cell_area):
    """Area-integrated flux [Pg C yr⁻¹] for a single Fay domain."""
    mask_2d = xr.DataArray(
        (biome_mask == domain_id).astype(float),
        dims=["latitude", "longitude"],
        coords={"latitude": flux_da["latitude"], "longitude": flux_da["longitude"]},
    )
    total = (flux_da * cell_area * mask_2d).sum(
        dim=["latitude", "longitude"], skipna=True
    )
    return (total * cfg.MOL_C_TO_PG).values


# ===========================================================================
# SHARED MAP HELPER
# ===========================================================================

def _fay_global_biome_map(fig, subplot_spec, biome_mask, lon_fay, lat_fay, n=17):
    """Biome overview map for the final panel of Fay figures."""
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
    ax.set_title("Fay (2014) Mean Biomes", fontsize=10)
    plt.colorbar(mm, ax=ax, shrink=0.8, ticks=np.arange(1, n + 1),
                 orientation="horizontal", pad=0.06)


# ===========================================================================
# FIG A — FAY DOMAIN TIMESERIES
# ===========================================================================

def fig_fay_ts(flux_da, biome_mask, lon_fay, lat_fay, cell_area) -> None:
    """
    17-panel monthly timeseries of area-integrated flux per Fay domain
    [Pg C yr⁻¹], coloured by domain (jet colourmap matching the biome map).
    """
    n_domains  = 17
    cmap_jet   = plt.cm.jet
    dom_colors = {i: cmap_jet((i - 1) / (n_domains - 1))
                  for i in range(1, n_domains + 1)}

    flux_r = flux_da.sel(latitude=lat_fay, longitude=lon_fay, method="nearest")
    area_r = cell_area.sel(latitude=lat_fay, longitude=lon_fay, method="nearest")
    t_pd   = pd.to_datetime(flux_r.time.values)

    fig = plt.figure(figsize=(20, 24))
    gs  = gridspec.GridSpec(6, 3, figure=fig,
                            hspace=0.32, wspace=0.25,
                            top=0.93, bottom=0.04,
                            left=0.06, right=0.97)

    for i in range(1, n_domains + 1):
        print(f"  [fay-ts] domain {i:2d}/{n_domains} ...")
        row, col = divmod(i - 1, 3)
        ax = fig.add_subplot(gs[row, col])
        dc = dom_colors[i]

        ts = domain_flux_timeseries(flux_r, biome_mask, i, area_r)
        ax.plot(t_pd, ts, color=dc, lw=1.4, alpha=0.95)
        ax.axhline(0, color="gray", lw=0.6, ls=":")
        ax.set_ylabel("Pg C yr⁻¹", fontsize=11)
        ax.tick_params(axis="x", labelsize=10, rotation=30)
        ax.tick_params(axis="y", labelsize=10)
        ax.set_title(f"Domain {i}", fontsize=12, pad=5, fontweight="bold")
        _clean(ax, minor=False)

    _fay_global_biome_map(fig, gs[5, 2], biome_mask, lon_fay, lat_fay)
    fig.suptitle(f"Fay (2014) domain flux — monthly timeseries\n{REC_LABEL}",
                 fontsize=14, y=0.995)

    out = cfg.FIG_DIR / "fig_fay_ts.png"
    fig.savefig(out, dpi=cfg.FIGURE_DPI)
    plt.close(fig)
    print(f"[fig] {out}")


# ===========================================================================
# FIG B — FAY DOMAIN SEASONAL CLIMATOLOGIES
# ===========================================================================

def fig_fay_clim(flux_da, biome_mask, lon_fay, lat_fay, cell_area) -> None:
    """
    17-panel climatological seasonal cycle ±1σ per Fay domain [Pg C yr⁻¹].
    """
    n_domains   = 17
    cmap_jet    = plt.cm.jet
    dom_colors  = {i: cmap_jet((i - 1) / (n_domains - 1))
                   for i in range(1, n_domains + 1)}
    months_abbr = ["Jan","Feb","Mar","Apr","May","Jun",
                   "Jul","Aug","Sep","Oct","Nov","Dec"]
    m_idx = np.arange(1, 13)

    flux_r = flux_da.sel(latitude=lat_fay, longitude=lon_fay, method="nearest")
    area_r = cell_area.sel(latitude=lat_fay, longitude=lon_fay, method="nearest")
    t_pd   = pd.to_datetime(flux_r.time.values)

    fig = plt.figure(figsize=(20, 24))
    gs  = gridspec.GridSpec(6, 3, figure=fig,
                            hspace=0.32, wspace=0.25,
                            top=0.93, bottom=0.04,
                            left=0.06, right=0.97)

    for i in range(1, n_domains + 1):
        print(f"  [fay-clim] domain {i:2d}/{n_domains} ...")
        row, col = divmod(i - 1, 3)
        ax = fig.add_subplot(gs[row, col])
        dc = dom_colors[i]

        ts   = domain_flux_timeseries(flux_r, biome_mask, i, area_r)
        s    = pd.Series(ts, index=t_pd)
        clim = s.groupby(s.index.month).mean().values
        std  = s.groupby(s.index.month).std().values

        ax.plot(m_idx, clim, color=dc, lw=2.0, marker="o", ms=5)
        ax.fill_between(m_idx, clim - std, clim + std, color=dc, alpha=0.18)
        ax.axhline(0, color="gray", lw=0.6, ls=":")
        ax.set_xticks(m_idx)
        ax.set_xticklabels(months_abbr, fontsize=9, rotation=45)
        ax.set_ylabel("Pg C yr⁻¹", fontsize=11)
        ax.tick_params(axis="y", labelsize=10)
        ax.set_title(f"Domain {i}", fontsize=12, pad=5, fontweight="bold")
        _clean(ax, minor=False)

    _fay_global_biome_map(fig, gs[5, 2], biome_mask, lon_fay, lat_fay)

    handles = [
        plt.Line2D([0],[0], color="gray", lw=2.0, marker="o", ms=5,
                   label=f"{REC_LABEL}  (domain colour)"),
        mpatches.Patch(color="gray", alpha=0.2, label="±1σ across years"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=2,
               fontsize=11, bbox_to_anchor=(0.5, 0.975),
               framealpha=0.9, edgecolor="#555555")
    fig.suptitle(
        f"Fay (2014) domain flux — climatological seasonal cycle ±1σ\n{REC_LABEL}",
        fontsize=14, y=0.995,
    )

    out = cfg.FIG_DIR / "fig_fay_clim.png"
    fig.savefig(out, dpi=cfg.FIGURE_DPI)
    plt.close(fig)
    print(f"[fig] {out}")


# ===========================================================================
# FIG C — FAY DOMAIN TREND ATTRIBUTION
# ===========================================================================

def fig_fay_trends(flux_da, biome_mask, lon_fay, lat_fay, cell_area) -> None:
    """
    17-panel Sen's slope per domain [Pg C yr⁻¹ per decade] with
    Mann-Kendall significance. Filled = p ≤ 0.01; hatched = p > 0.01.
    """
    flux_r = flux_da.sel(latitude=lat_fay, longitude=lon_fay, method="nearest")
    area_r = cell_area.sel(latitude=lat_fay, longitude=lon_fay, method="nearest")
    t_pd   = pd.to_datetime(flux_r.time.values)

    n_domains  = 17
    dom_slopes = {}
    dom_pvals  = {}

    for i in range(1, n_domains + 1):
        ts_monthly = domain_flux_timeseries(flux_r, biome_mask, i, area_r)
        s          = pd.Series(ts_monthly, index=t_pd)
        s_annual   = s.resample("YE").mean().dropna()
        if len(s_annual) < 5:
            dom_slopes[i] = np.nan
            dom_pvals[i]  = np.nan
            continue
        x   = np.arange(len(s_annual), dtype=float)
        y   = s_annual.values
        res = theilslopes(y, x)
        dom_slopes[i] = res.slope * 10
        dom_pvals[i]  = mann_kendall_pvalue(y)

    cmap_jet   = plt.cm.jet
    dom_colors = {i: cmap_jet((i - 1) / (n_domains - 1))
                  for i in range(1, n_domains + 1)}

    fig = plt.figure(figsize=(20, 24))
    gs  = gridspec.GridSpec(6, 3, figure=fig,
                            hspace=0.32, wspace=0.25,
                            top=0.93, bottom=0.04,
                            left=0.06, right=0.97)

    for i in range(1, n_domains + 1):
        print(f"  [fay-trend] domain {i:2d}/{n_domains} ...")
        row, col = divmod(i - 1, 3)
        ax  = fig.add_subplot(gs[row, col])
        dc  = dom_colors[i]
        sl  = dom_slopes[i]
        pv  = dom_pvals[i]
        sig = np.isfinite(pv) and pv <= 0.01

        if np.isfinite(sl):
            bar_kw = dict(color=dc, alpha=0.85, edgecolor="k", linewidth=0.8)
            if not sig:
                bar_kw.update(hatch="///", alpha=0.5)
            ax.bar(0, sl, **bar_kw)

        ax.axhline(0, color="gray", lw=0.8, ls=":")
        ax.set_xlim(-0.5, 0.5)
        ax.set_xticks([])
        ax.set_ylabel("Pg C yr⁻¹ per decade", fontsize=11)
        ax.tick_params(axis="y", labelsize=10)
        sig_str = f"p = {pv:.3f}" if np.isfinite(pv) else "p = N/A"
        star    = " ✱" if sig else ""
        ax.set_title(f"Domain {i}  [{sig_str}]{star}",
                     fontsize=12, pad=5, fontweight="bold")
        _clean(ax, minor=False)

    _fay_global_biome_map(fig, gs[5, 2], biome_mask, lon_fay, lat_fay)

    handles = [
        mpatches.Patch(color="gray", alpha=0.85, edgecolor="k",
                       label="Significant trend (p ≤ 0.01)  ✱"),
        mpatches.Patch(color="gray", alpha=0.5, hatch="///", edgecolor="k",
                       label="Not significant (p > 0.01)"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=2,
               fontsize=11, bbox_to_anchor=(0.5, 0.975),
               framealpha=0.9, edgecolor="#555555")
    fig.suptitle(
        "Fay (2014) domain flux — Sen's slope per decade (Mann-Kendall, α = 0.01)\n"
        f"{REC_LABEL}",
        fontsize=14, y=0.995,
    )

    out = cfg.FIG_DIR / "fig_fay_trends.png"
    fig.savefig(out, dpi=cfg.FIGURE_DPI)
    plt.close(fig)
    print(f"[fig] {out}")


# ===========================================================================
# CLIMATE INDEX LOADER  (mirrors 05_plot_results.load_climate_indices)
# ===========================================================================

def load_climate_indices() -> dict | None:
    """
    Load ONI and SAM indices. Returns None if files unavailable (e.g. HPC
    compute node with no outbound internet). See 05_plot_results.py for
    full documentation and manual download instructions.

    Manual download (login node with internet):
        wget -O <DATA_DIR>/oni_monthly.txt \\
            https://psl.noaa.gov/data/correlation/oni.data
        wget -O <DATA_DIR>/sam_monthly.txt \\
            https://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/aao/monthly.aao.index.b79.current.ascii.table
    """
    T_START_YR = int(T_START[:4])
    T_END_YR   = int(T_END[:4])

    oni_file = cfg.DATA_DIR / "oni_monthly.txt"
    sam_file = cfg.DATA_DIR / "sam_monthly.txt"

    _DOWNLOAD_NOTE = (
        "\n[climate] *** File not found and download failed (no internet). ***\n"
        "[climate] Download both files on a login node and place in data/:\n"
        f"[climate]   wget -O {oni_file} \\\n"
        "[climate]       https://psl.noaa.gov/data/correlation/oni.data\n"
        f"[climate]   wget -O {sam_file} \\\n"
        "[climate]       https://www.cpc.ncep.noaa.gov/products/precip/CWlink/"
        "daily_ao_index/aao/monthly.aao.index.b79.current.ascii.table\n"
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
        print(f"[warn] ONI file {oni_file} parsed empty.")
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
        print(f"[warn] SAM file {sam_file} parsed empty.")
        return None
    sam_df = pd.DataFrame(sam_rows)
    sam_df["time"] = pd.to_datetime(
        sam_df["year"].astype(str) + "-" + sam_df["month"].astype(str).str.zfill(2)
    )
    sam_monthly = sam_df.set_index("time")["sam"].sort_index()

    t_range    = pd.date_range(f"{T_START_YR}", f"{T_END_YR}", freq="YE")
    oni_annual = oni_monthly.resample("YE").mean().reindex(t_range, method="nearest")
    sam_annual = sam_monthly.resample("YE").mean().reindex(t_range, method="nearest")

    return {
        "oni_monthly": oni_monthly,
        "sam_monthly": sam_monthly,
        "oni_annual":  oni_annual,
        "sam_annual":  sam_annual,
    }


# ===========================================================================
# FIG D — FAY DOMAIN CLIMATE INDEX REGRESSION
# ===========================================================================

def fig_fay_climate_regression(flux_da, biome_mask, lon_fay, lat_fay,
                                cell_area, idx: dict) -> None:
    """
    17-panel bar chart: per-domain Spearman ρ of annual J_net against ONI,
    SAM, and their OLS multilinear combination (α·ONI + β·SAM).

    Each panel shows three bars:
        Red   — ρ(J_net, ONI)
        Green — ρ(J_net, SAM)
        Blue  — ρ(J_net, α·ONI + β·SAM)   [multilinear]

    Filled bars = significant at α = 0.05 (Spearman p ≤ 0.05).
    Hatched bars = not significant.

    Scientific rationale
    --------------------
    Bellacicco et al. (2025) show that global PIP export is dominated by ONI
    with SAM providing secondary modulation, and that the signal is primarily
    driven by the Southern Ocean. Here we test which Fay domains at the
    surface show the same fingerprint, allowing direct comparison between
    surface and subsurface carbon cycle responses to climate variability.
    """
    flux_r = flux_da.sel(latitude=lat_fay, longitude=lon_fay, method="nearest")
    area_r = cell_area.sel(latitude=lat_fay, longitude=lon_fay, method="nearest")
    t_pd   = pd.to_datetime(flux_r.time.values)

    oni_ann = idx["oni_annual"]
    sam_ann = idx["sam_annual"]

    n_domains = 17
    C_ONI = "#CC2936"
    C_SAM = "#2E8B57"
    C_ML  = "#0558AD"

    cmap_jet   = plt.cm.jet
    dom_colors = {i: cmap_jet((i - 1) / (n_domains - 1))
                  for i in range(1, n_domains + 1)}

    fig = plt.figure(figsize=(20, 24))
    gs  = gridspec.GridSpec(6, 3, figure=fig,
                            hspace=0.40, wspace=0.30,
                            top=0.93, bottom=0.05,
                            left=0.06, right=0.97)

    for i in range(1, n_domains + 1):
        print(f"  [fay-climate] domain {i:2d}/{n_domains} ...")
        row, col = divmod(i - 1, 3)
        ax = fig.add_subplot(gs[row, col])

        # Build annual domain flux
        ts_monthly = domain_flux_timeseries(flux_r, biome_mask, i, area_r)
        s          = pd.Series(ts_monthly, index=t_pd)
        s_annual   = s.resample("YE").mean()

        # Align all three to common index
        combined = pd.DataFrame({
            "j":   s_annual,
            "oni": oni_ann,
            "sam": sam_ann,
        }).dropna()

        if len(combined) < 5:
            ax.set_title(f"Domain {i}  [n<5]", fontsize=11,
                         pad=4, fontweight="bold")
            _clean(ax, minor=False)
            continue

        j_v   = combined["j"].values
        oni_v = combined["oni"].values
        sam_v = combined["sam"].values

        # Normalise for multilinear
        oni_z = (oni_v - oni_v.mean()) / (oni_v.std() + 1e-9)
        sam_z = (sam_v - sam_v.mean()) / (sam_v.std() + 1e-9)
        X     = np.column_stack([oni_z, sam_z, np.ones(len(j_v))])
        coef, *_ = np.linalg.lstsq(X, j_v, rcond=None)
        j_pred   = X @ coef

        rho_oni, p_oni = stats.spearmanr(j_v, oni_v)
        rho_sam, p_sam = stats.spearmanr(j_v, sam_v)
        rho_ml,  p_ml  = stats.spearmanr(j_v, j_pred)

        # Draw bars
        positions = [-0.25, 0.0, 0.25]
        rhos      = [rho_oni, rho_sam, rho_ml]
        pvals     = [p_oni,   p_sam,   p_ml]
        colors    = [C_ONI,   C_SAM,   C_ML]
        labels    = ["ONI", "SAM", "ONI+SAM"]

        for xpos, rho, pv, col, lbl in zip(positions, rhos, pvals, colors, labels):
            sig    = np.isfinite(pv) and pv <= 0.05
            bar_kw = dict(color=col, alpha=0.85, edgecolor="k",
                          linewidth=0.7, width=0.22)
            if not sig:
                bar_kw.update(hatch="///", alpha=0.45)
            ax.bar(xpos, rho, **bar_kw)
            star = "✱" if sig else ""
            ax.text(xpos, rho + (0.04 if rho >= 0 else -0.07),
                    star, ha="center", va="bottom" if rho >= 0 else "top",
                    fontsize=10, color=col, fontweight="bold")

        ax.axhline(0, color="gray", lw=0.8, ls=":")
        ax.axhline( 0.4, color="gray", lw=0.4, ls=":", alpha=0.5)
        ax.axhline(-0.4, color="gray", lw=0.4, ls=":", alpha=0.5)
        ax.set_xlim(-0.5, 0.5)
        ax.set_ylim(-1.0, 1.0)
        ax.set_xticks(positions)
        ax.set_xticklabels(labels, fontsize=9)
        ax.set_ylabel("Spearman ρ", fontsize=10)
        ax.tick_params(axis="y", labelsize=9)
        ax.set_title(f"Domain {i}  (n={len(combined)})",
                     fontsize=11, pad=4, fontweight="bold")
        _clean(ax, minor=False)

    # Last panel: biome map
    _fay_global_biome_map(fig, gs[5, 2], biome_mask, lon_fay, lat_fay)

    # Legend
    handles = [
        mpatches.Patch(color=C_ONI, alpha=0.85, edgecolor="k",
                       label="ρ(J_net, ONI)"),
        mpatches.Patch(color=C_SAM, alpha=0.85, edgecolor="k",
                       label="ρ(J_net, SAM)"),
        mpatches.Patch(color=C_ML,  alpha=0.85, edgecolor="k",
                       label="ρ(J_net, α·ONI+β·SAM)"),
        mpatches.Patch(color="gray", alpha=0.45, hatch="///", edgecolor="k",
                       label="Not significant (p > 0.05)"),
        mpatches.Patch(color="white", edgecolor="white",
                       label="✱ = p ≤ 0.05"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=5,
               fontsize=10, bbox_to_anchor=(0.5, 0.975),
               framealpha=0.9, edgecolor="#555555")
    fig.suptitle(
        "Fay (2014) domain — Spearman ρ with ONI, SAM and ONI+SAM\n"
        f"{REC_LABEL}  ·  following Bellacicco et al. (2025, Nat. Commun.)",
        fontsize=13, y=0.997,
    )

    out = cfg.FIG_DIR / "fig_fay_climate_regression.png"
    fig.savefig(out, dpi=cfg.FIGURE_DPI)
    plt.close(fig)
    print(f"[fig] {out}")


# ===========================================================================
# MAIN
# ===========================================================================

def main():
    cfg.FIG_DIR.mkdir(parents=True, exist_ok=True)

    print("[load] flux_3d.nc ...")
    flux_da = load_flux()
    print(f"       time: {flux_da.time.values[0]} → {flux_da.time.values[-1]}")

    print("[load] Fay biome mask ...")
    biome_mask, lon_fay, lat_fay = load_fay_biomes()

    print("[compute] Grid cell areas ...")
    cell_area = compute_grid_cell_area(flux_da["latitude"], flux_da["longitude"])

    print("[load] Climate indices (ONI, SAM) ...")
    idx = load_climate_indices()

    print("\n[fig A] Fay domain timeseries ...")
    fig_fay_ts(flux_da, biome_mask, lon_fay, lat_fay, cell_area)

    print("\n[fig B] Fay domain seasonal climatologies ...")
    fig_fay_clim(flux_da, biome_mask, lon_fay, lat_fay, cell_area)

    print("\n[fig C] Fay domain trend attribution ...")
    fig_fay_trends(flux_da, biome_mask, lon_fay, lat_fay, cell_area)

    if idx is not None:
        print("\n[fig D] Fay domain climate index regression ...")
        fig_fay_climate_regression(flux_da, biome_mask, lon_fay, lat_fay,
                                   cell_area, idx)
    else:
        print("\n[skip] fig D — climate index files not available "
              "(see [climate] messages above).")

    print(f"\n[done] Fay figures saved to {cfg.FIG_DIR}\n")


if __name__ == "__main__":
    main()
