"""
04_validate.py
==============
Validation of the GLOBAL_MULTIYEAR_BGC_001_029 surface CO2 flux reconstruction
against the CMEMS MULTIOBS SOCAT-NN observation-based product.

Figures produced:

    fig_validation_ts.png
        TWO-curve global time series [Pg C yr⁻¹] — reconstruction vs MULTIOBS.

    fig_validation_loess.png
        Same two curves with LOESS smoothing overlaid to reveal multi-year
        variability without spectral analysis.

    fig_validation_map.png
        1×2 contourf map: RMSD (top) and Bias (bottom) of the reconstruction
        vs MULTIOBS. Fixed colourbar extents. Robinson projection.

    fig_fay_domains.png
        17-panel figure (one per Fay biome) + global biome map in panel 18.
        Each panel: full monthly time series (rec vs MULTIOBS, left y-axis)
        + climatological seasonal cycle (right y-axis) for that domain.

    fig_spectra.png
        Log-log power spectral density: reconstruction vs MULTIOBS overlaid.
        X-axis in years (period).

    fig_cps.png
        Cross-power spectrum decomposition: gain, phase [days], coherence.
        High nfft for smooth output.

    validation_metrics.csv
        Scalar statistics: bias, RMSD, Pearson r, σ_rec, σ_obs.
        Per-domain statistics for each Fay biome.

All figures at 300 dpi.

Changelog:
    v1.0.0 — Initial.
    v1.1.0 — Added improved reconstruction; 2×2 maps; Taylor diagram; 300 dpi.
    v1.2.0 — Dropped MULTIOBS pCO2 reconstruction; dropped Taylor diagram;
              single reconstruction (GLOBAL_MULTIYEAR_BGC_001_029) only;
              added LOESS multi-year variability plot;
              added Fay (2014) domain time series + climatology validation;
              added power spectra and cross-power spectrum (gain/phase/coherence);
              fixed MULTIOBS sign (no flip);
              product label updated to GLOBAL_MULTIYEAR_BGC_001_029.

Usage:
    python scripts/04_validate.py
    Requires: data/Time_Varying_Biomes.cmems.nc  (Fay biome mask on CMEMS grid)
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mticker
from matplotlib.colors import ListedColormap
from scipy import stats
from scipy.signal import welch, csd

sys.path.insert(0, str(Path(__file__).parent.parent))
import config as cfg

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    HAS_CARTOPY = True
except ImportError:
    HAS_CARTOPY = False
    print("[warn] cartopy not found — maps will use plain imshow.")

# Plot style
plt.rcParams.update({
    "font.family":  "sans-serif",
    "font.size":    10,
    "axes.labelsize": 11,
    "axes.titlesize": 11,
    "legend.fontsize": 9,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
})

# Label used in all figures for the reconstruction
REC_LABEL = "GLOBAL_MULTIYEAR_BGC_001_029"
OBS_LABEL = "MULTIOBS (SOCAT-NN)"


# ===========================================================================
# LOAD DATA
# ===========================================================================

def load_data():
    """
    Load reconstruction flux, MULTIOBS validation flux, and global integrals.

    Returns
    -------
    fgco2_rec : xr.DataArray   reconstruction [mol C m⁻² yr⁻¹]
    fgco2_obs : xr.DataArray   MULTIOBS flux  [mol C m⁻² yr⁻¹]
    J_rec     : xr.DataArray   reconstruction global integral [Pg C yr⁻¹]
    J_obs     : xr.DataArray   MULTIOBS global integral      [Pg C yr⁻¹]
    ds_surf   : xr.Dataset     processed_surface (for mask, coords)
    """
    ds_flux = xr.open_dataset(cfg.DATA_DIR / "flux_3d.nc",           chunks="auto")
    ds_surf = xr.open_dataset(cfg.DATA_DIR / "processed_surface.nc",  chunks="auto")
    ds_glob = xr.open_dataset(cfg.OUT_DIR  / "global_flux.nc")

    if "fgco2_obs" not in ds_surf:
        raise FileNotFoundError(
            "fgco2_obs not found in processed_surface.nc. "
            "Check MULTIOBS download (01_download_data.py)."
        )

    fgco2_rec = ds_flux["fgco2"]
    # v1.2.0: NO sign flip — use fgco2_obs as stored
    fgco2_obs = ds_surf["fgco2_obs"]

    J_rec = ds_glob["J_net_PgC"]

    # Align to common time window
    common_time = fgco2_rec.time[
        np.isin(fgco2_rec.time.values, fgco2_obs.time.values)
    ]
    fgco2_rec = fgco2_rec.sel(time=common_time)
    fgco2_obs = fgco2_obs.sel(time=common_time)
    J_rec     = J_rec.sel(time=common_time)

    # MULTIOBS global integral
    cell_area = compute_grid_cell_area(fgco2_obs["latitude"], fgco2_obs["longitude"])
    F_obs_mol = (fgco2_obs.where(ds_surf["ocean_mask"] == 1) * cell_area).sum(
        dim=["latitude", "longitude"], skipna=True
    )
    J_obs = (F_obs_mol * cfg.MOL_C_TO_PG).rename("J_obs_PgC")
    J_obs.attrs = {"units": "Pg C yr-1", "long_name": "MULTIOBS global CO2 uptake"}

    return fgco2_rec, fgco2_obs, J_rec, J_obs, ds_surf


# ===========================================================================
# SKILL METRICS
# ===========================================================================

def compute_grid_cell_area(lat: xr.DataArray, lon: xr.DataArray) -> xr.DataArray:
    R       = cfg.EARTH_RADIUS_M
    d_lat   = float(np.abs(lat.diff("latitude").mean()))
    d_lon   = float(np.abs(lon.diff("longitude").mean()))
    area_1d = (R**2) * np.deg2rad(d_lon) * np.deg2rad(d_lat) * np.cos(np.deg2rad(lat))
    area_2d = area_1d.expand_dims({"longitude": lon}).transpose("latitude", "longitude")
    area_2d.attrs = {"long_name": "Grid cell area", "units": "m2"}
    return area_2d.rename("cell_area")


def pixel_rmsd(rec: xr.DataArray, obs: xr.DataArray) -> xr.DataArray:
    rmsd = np.sqrt(((rec - obs)**2).mean(dim="time"))
    rmsd.attrs = {"long_name": "RMSD", "units": "mol C m-2 yr-1"}
    return rmsd


def pixel_bias(rec: xr.DataArray, obs: xr.DataArray) -> xr.DataArray:
    bias = (rec - obs).mean(dim="time")
    bias.attrs = {"long_name": "Bias (rec − obs)", "units": "mol C m-2 yr-1"}
    return bias


def scalar_metrics(rec: xr.DataArray, obs: xr.DataArray) -> dict:
    r = rec.values.ravel()
    o = obs.values.ravel()
    ok = np.isfinite(r) & np.isfinite(o)
    r, o = r[ok], o[ok]
    corr, pval = stats.pearsonr(r, o)
    return {
        "bias_mol_C_m2_yr": float(np.mean(r - o)),
        "rmsd_mol_C_m2_yr": float(np.sqrt(np.mean((r - o)**2))),
        "pearson_r":        float(corr),
        "pearson_pval":     float(pval),
        "std_rec":          float(np.std(r)),
        "std_obs":          float(np.std(o)),
    }


# ===========================================================================
# LOESS SMOOTHER
# ===========================================================================

def loess_smooth(x: np.ndarray, y: np.ndarray, frac: float = 0.15) -> np.ndarray:
    """
    LOWESS smoother (locally weighted scatterplot smoothing).
    frac : fraction of data used in each local regression window.
           0.15 ≈ ~18-month window on a monthly series — captures multi-year
           variability while removing the seasonal cycle.
    Returns smoothed y values at the same x positions.
    """
    from statsmodels.nonparametric.smoothers_lowess import lowess
    smoothed = lowess(y, x, frac=frac, it=3, return_sorted=False)
    return smoothed


# ===========================================================================
# FIGURE 1 — GLOBAL TIME SERIES (two curves, raw monthly)
# ===========================================================================

def plot_timeseries(J_rec: xr.DataArray, J_obs: xr.DataArray) -> None:
    """
    Global net ocean CO2 uptake [Pg C yr⁻¹]:
        blue  solid  — GLOBAL_MULTIYEAR_BGC_001_029 reconstruction
        red   dashed — MULTIOBS (SOCAT-NN) reference
    """
    fig, ax = plt.subplots(figsize=(14, 5))

    t_rec = pd.to_datetime(J_rec.time.values)
    t_obs = pd.to_datetime(J_obs.time.values)

    ax.plot(t_rec, J_rec.values,
            color="steelblue", lw=1.3, alpha=0.85,
            label=REC_LABEL)
    ax.plot(t_obs, J_obs.values,
            color="firebrick", lw=1.3, ls="--", alpha=0.85,
            label=f"{OBS_LABEL}  (reference)")

    ax.axhline(0, color="black", lw=0.7, ls=":")
    ax.set_xlabel("Year")
    ax.set_ylabel("Net ocean CO₂ uptake  [Pg C yr⁻¹]")
    ax.set_title("Global net ocean CO₂ uptake — reconstruction vs observation-based reference")
    ax.legend(loc="upper left")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()

    out = cfg.FIG_DIR / "fig_validation_ts.png"
    fig.savefig(out, dpi=cfg.FIGURE_DPI)
    plt.close(fig)
    print(f"[fig] {out}")


# ===========================================================================
# FIGURE 2 — LOESS MULTI-YEAR VARIABILITY
# ===========================================================================

def plot_loess(J_rec: xr.DataArray, J_obs: xr.DataArray, frac: float = 0.15) -> None:
    """
    Raw monthly time series with LOESS smooth overlaid for both curves.
    The LOESS smooth (frac≈0.15 → ~18-month window) reveals multi-year
    variability by suppressing the seasonal cycle without spectral analysis.
    """
    fig, ax = plt.subplots(figsize=(14, 5))

    t_rec = pd.to_datetime(J_rec.time.values)
    t_obs = pd.to_datetime(J_obs.time.values)
    x_rec = np.arange(len(t_rec), dtype=float)
    x_obs = np.arange(len(t_obs), dtype=float)

    # Raw monthly — faint
    ax.plot(t_rec, J_rec.values,
            color="steelblue", lw=0.7, alpha=0.3, label="_nolegend_")
    ax.plot(t_obs, J_obs.values,
            color="firebrick", lw=0.7, alpha=0.3, ls="--", label="_nolegend_")

    # LOESS smooth — bold
    rec_vals = J_rec.values.copy()
    obs_vals = J_obs.values.copy()
    ok_rec   = np.isfinite(rec_vals)
    ok_obs   = np.isfinite(obs_vals)

    smooth_rec = loess_smooth(x_rec[ok_rec], rec_vals[ok_rec], frac=frac)
    smooth_obs = loess_smooth(x_obs[ok_obs], obs_vals[ok_obs], frac=frac)

    ax.plot(t_rec[ok_rec], smooth_rec,
            color="steelblue", lw=2.2, label=f"{REC_LABEL}  (LOESS)")
    ax.plot(t_obs[ok_obs], smooth_obs,
            color="firebrick", lw=2.2, ls="--", label=f"{OBS_LABEL}  (LOESS)")

    ax.axhline(0, color="black", lw=0.7, ls=":")
    ax.set_xlabel("Year")
    ax.set_ylabel("Net ocean CO₂ uptake  [Pg C yr⁻¹]")
    ax.set_title(
        f"Multi-year variability — LOESS smoothing (frac = {frac:.2f})\n"
        "Pale lines: raw monthly  |  Bold lines: LOESS smooth"
    )
    ax.legend(loc="upper left")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()

    out = cfg.FIG_DIR / "fig_validation_loess.png"
    fig.savefig(out, dpi=cfg.FIGURE_DPI)
    plt.close(fig)
    print(f"[fig] {out}")


# ===========================================================================
# FIGURE 3 — VALIDATION MAP (1×2: RMSD top, Bias bottom)
# ===========================================================================

def plot_validation_map(
    rmsd: xr.DataArray,
    bias: xr.DataArray,
) -> None:
    """
    Two-panel contourf map:
        Top    — RMSD  [0 – VAL_RMSD_MAX mol C m⁻² yr⁻¹]  — Reds
        Bottom — Bias  [VAL_BIAS_MIN – VAL_BIAS_MAX]        — RdBu_r
    Fixed colourbar extents for reproducibility.
    """
    rmsd_levels = np.linspace(0,               cfg.VAL_RMSD_MAX, 21)
    bias_levels = np.linspace(cfg.VAL_BIAS_MIN, cfg.VAL_BIAS_MAX, 21)

    fig = plt.figure(figsize=(12, 10))
    gs  = gridspec.GridSpec(2, 1, hspace=0.35,
                            left=0.05, right=0.88, top=0.93, bottom=0.05)
    cax_rmsd = fig.add_axes([0.90, 0.55, 0.018, 0.35])
    cax_bias = fig.add_axes([0.90, 0.08, 0.018, 0.35])

    panels = [
        (rmsd, rmsd_levels, "Reds",        cax_rmsd,
         f"RMSD — {REC_LABEL} vs {OBS_LABEL}  [mol C m⁻² yr⁻¹]",
         "RMSD  [mol C m⁻² yr⁻¹]"),
        (bias, bias_levels, cfg.CMAP_FLUX,  cax_bias,
         f"Bias (rec − obs) — {REC_LABEL} vs {OBS_LABEL}  [mol C m⁻² yr⁻¹]",
         "Bias  [mol C m⁻² yr⁻¹]"),
    ]

    for row, (da, levels, cmap, cax, title, cbar_label) in enumerate(panels):
        lat = da["latitude"].values
        lon = da["longitude"].values
        Z   = da.values

        if HAS_CARTOPY:
            ax = fig.add_subplot(gs[row], projection=ccrs.Robinson())
            cf = ax.contourf(
                lon, lat, Z,
                levels=levels, cmap=cmap, extend="both",
                transform=ccrs.PlateCarree(),
            )
            ax.add_feature(cfeature.LAND,      facecolor="#d3d3d3", zorder=3)
            ax.add_feature(cfeature.COASTLINE, lw=0.3,              zorder=4)
            ax.set_global()
        else:
            ax = fig.add_subplot(gs[row])
            cf = ax.contourf(lon, lat, Z,
                             levels=levels, cmap=cmap, extend="both")
            ax.set_xlabel("Longitude")
            ax.set_ylabel("Latitude")

        ax.set_title(title, pad=6)
        plt.colorbar(cf, cax=cax, label=cbar_label)

    out = cfg.FIG_DIR / "fig_validation_map.png"
    fig.savefig(out, dpi=cfg.FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {out}")


# ===========================================================================
# FAY DOMAIN ANALYSIS
# ===========================================================================

def load_fay_biomes() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Load the Fay (2014) biome mask interpolated onto the CMEMS 0.25° grid.
    Expected file: data/Time_Varying_Biomes.cmems.nc
    Returns (biome_mask 2D, lon 1D, lat 1D).
    Biome IDs are integers 1–17.
    """
    fay_file = cfg.DATA_DIR / "Time_Varying_Biomes.cmems.nc"
    if not fay_file.exists():
        raise FileNotFoundError(
            f"Fay biome file not found: {fay_file}\n"
            "Run interpolator.py with input_type=cmems to generate it."
        )
    ds  = xr.open_dataset(fay_file)
    # Use MeanBiomes (static mask); isel year=0 if it has a year dimension
    if "year" in ds["MeanBiomes"].dims:
        mask = ds["MeanBiomes"].isel(year=0).values
    else:
        mask = ds["MeanBiomes"].values

    # Coordinates — handle both (y,x) with lon/lat 2D or 1D
    if ds["lon"].ndim == 2:
        lon = ds["lon"].values[0, :]
        lat = ds["lat"].values[:, 0]
    else:
        lon = ds["lon"].values
        lat = ds["lat"].values

    return mask.astype(float), lon, lat


def domain_mean_timeseries(
    flux: xr.DataArray,
    biome_mask: np.ndarray,
    domain_id: int,
) -> np.ndarray:
    """
    Area-weighted mean flux time series for a single Fay domain.
    Returns 1-D array [time] in mol C m⁻² yr⁻¹.
    """
    mask_2d = (biome_mask == domain_id).astype(float)   # (lat, lon)
    mask_da = xr.DataArray(
        mask_2d,
        dims=["latitude", "longitude"],
        coords={
            "latitude":  flux["latitude"],
            "longitude": flux["longitude"],
        },
    )
    # Area-weighted mean: Σ(F·A·mask) / Σ(A·mask)
    cell_area = compute_grid_cell_area(flux["latitude"], flux["longitude"])
    weight    = cell_area * mask_da
    total_w   = float(weight.sum())
    if total_w == 0:
        return np.full(flux.sizes["time"], np.nan)
    ts = (flux * weight).sum(dim=["latitude", "longitude"], skipna=True) / total_w
    return ts.values


def compute_climatology(ts: np.ndarray, time: pd.DatetimeIndex) -> np.ndarray:
    """
    Monthly climatology: mean of each calendar month across all years.
    Returns array of length 12 (Jan–Dec).
    """
    da = pd.Series(ts, index=time)
    return da.groupby(da.index.month).mean().values


def plot_fay_domains(
    fgco2_rec: xr.DataArray,
    fgco2_obs: xr.DataArray,
) -> None:
    """
    17-panel figure + global biome map (panel 18).
    Each of the 17 panels contains:
        Left  y-axis : full monthly time series — rec (blue) vs obs (red dashed)
        Right y-axis : climatological seasonal cycle — same colours, bold
    Panel 18: global MeanBiomes map.
    """
    print("[load] Fay biome mask ...")
    biome_mask, lon_fay, lat_fay = load_fay_biomes()
    n_domains = 17

    # Regrid flux arrays to the Fay grid if needed (should already be 0.25°)
    # Use nearest-neighbour reindex to snap to common grid
    fgco2_rec_r = fgco2_rec.sel(
        latitude=lat_fay,  longitude=lon_fay, method="nearest"
    )
    fgco2_obs_r = fgco2_obs.sel(
        latitude=lat_fay, longitude=lon_fay, method="nearest"
    )

    time_pd = pd.to_datetime(fgco2_rec_r.time.values)
    months  = ["J","F","M","A","M","J","J","A","S","O","N","D"]

    fig = plt.figure(figsize=(22, 26))
    nrows, ncols = 6, 3

    # Store per-domain metrics
    domain_metrics = {}

    for i in range(1, n_domains + 1):
        print(f"  [fay] domain {i:2d}/{n_domains} ...")

        ax = fig.add_subplot(nrows, ncols, i)
        ax2 = ax.twinx()

        ts_rec = domain_mean_timeseries(fgco2_rec_r, biome_mask, i)
        ts_obs = domain_mean_timeseries(fgco2_obs_r, biome_mask, i)

        # Full monthly time series
        ax.plot(time_pd, ts_rec,
                color="steelblue", lw=0.9, alpha=0.7, label=REC_LABEL)
        ax.plot(time_pd, ts_obs,
                color="firebrick", lw=0.9, alpha=0.7, ls="--", label=OBS_LABEL)
        ax.axhline(0, color="black", lw=0.5, ls=":")
        ax.set_ylabel("mol C m⁻² yr⁻¹", fontsize=7)
        ax.tick_params(axis="x", labelsize=7, rotation=30)
        ax.tick_params(axis="y", labelsize=7)

        # Climatology on twin axis
        clim_rec = compute_climatology(ts_rec, time_pd)
        clim_obs = compute_climatology(ts_obs, time_pd)
        m_idx    = np.arange(1, 13)

        ax2.plot(m_idx, clim_rec,
                 color="steelblue", lw=2.0, marker="o", ms=3,
                 label=f"{REC_LABEL} clim")
        ax2.plot(m_idx, clim_obs,
                 color="firebrick", lw=2.0, marker="o", ms=3, ls="--",
                 label=f"{OBS_LABEL} clim")
        ax2.set_xticks(m_idx)
        ax2.set_xticklabels(months, fontsize=7)
        ax2.tick_params(axis="y", labelsize=7, labelcolor="dimgray")
        ax2.set_ylabel("Climatology\nmol C m⁻² yr⁻¹", fontsize=6, color="dimgray")

        # Inset domain map
        pos         = ax.get_position()
        inset_w     = 0.22 * pos.width
        inset_h     = 0.28 * pos.height
        inset_x     = pos.x0 + 0.77 * pos.width
        inset_y     = pos.y0 + 0.68 * pos.height

        lon2d, lat2d = np.meshgrid(lon_fay, lat_fay)

        if HAS_CARTOPY:
            axins = fig.add_axes(
                [inset_x, inset_y, inset_w, inset_h],
                projection=ccrs.PlateCarree(),
            )
            masked = np.where(biome_mask == float(i), 1, np.nan)
            axins.pcolormesh(lon2d, lat2d, masked,
                             cmap=ListedColormap(["darkred"]),
                             transform=ccrs.PlateCarree())
            axins.add_feature(cfeature.LAND, color="lightgray", zorder=2)
            axins.coastlines(linewidth=0.3)
        else:
            axins = fig.add_axes([inset_x, inset_y, inset_w, inset_h])
            masked = np.where(biome_mask == float(i), 1, np.nan)
            axins.contourf(lon2d, lat2d, masked, colors=["darkred"])
            axins.axis("off")

        ax.set_title(f"Domain {i}", fontsize=9, pad=3)

        # Per-domain scalar metrics
        ok  = np.isfinite(ts_rec) & np.isfinite(ts_obs)
        if ok.sum() > 5:
            r, _ = stats.pearsonr(ts_rec[ok], ts_obs[ok])
            domain_metrics[i] = {
                "bias": float(np.mean(ts_rec[ok] - ts_obs[ok])),
                "rmsd": float(np.sqrt(np.mean((ts_rec[ok] - ts_obs[ok])**2))),
                "pearson_r": float(r),
            }

    # Panel 18: global biome map
    ax_map = fig.add_subplot(nrows, ncols, 18,
                             projection=ccrs.PlateCarree() if HAS_CARTOPY else None)
    import matplotlib as mpl
    cmap_jet = plt.cm.jet
    bounds   = np.arange(0.5, n_domains + 1.5, 1)
    norm     = mpl.colors.BoundaryNorm(bounds, cmap_jet.N)
    lon2d, lat2d = np.meshgrid(lon_fay, lat_fay)

    if HAS_CARTOPY:
        mm = ax_map.pcolormesh(lon2d, lat2d, biome_mask,
                               cmap=cmap_jet, norm=norm,
                               transform=ccrs.PlateCarree())
        ax_map.add_feature(cfeature.LAND,      color="gray")
        ax_map.add_feature(cfeature.COASTLINE, linewidth=0.4)
        ax_map.gridlines(draw_labels=False, linestyle="--", alpha=0.5)
    else:
        mm = ax_map.contourf(lon_fay, lat_fay, biome_mask,
                             levels=bounds, cmap=cmap_jet, norm=norm)
    ax_map.set_title("Fay (2014) Mean Biomes", fontsize=9)
    plt.colorbar(mm, ax=ax_map, shrink=0.7,
                 ticks=np.arange(1, n_domains + 1))

    # Figure-level legend
    handles = [
        plt.Line2D([0], [0], color="steelblue", lw=1.5, label=REC_LABEL),
        plt.Line2D([0], [0], color="firebrick", lw=1.5, ls="--", label=OBS_LABEL),
        plt.Line2D([0], [0], color="steelblue", lw=2.5, label=f"{REC_LABEL} climatology"),
        plt.Line2D([0], [0], color="firebrick", lw=2.5, ls="--",
                   label=f"{OBS_LABEL} climatology"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=2,
               fontsize=9, bbox_to_anchor=(0.5, 0.995))

    fig.suptitle(
        f"Fay (2014) domain validation — {REC_LABEL} vs {OBS_LABEL}\n"
        "Thin lines: monthly time series  |  Bold lines: climatological seasonal cycle",
        fontsize=10, y=1.00,
    )

    out = cfg.FIG_DIR / "fig_fay_domains.png"
    fig.savefig(out, dpi=cfg.FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {out}")

    return domain_metrics


# ===========================================================================
# FIGURE 5 — POWER SPECTRAL DENSITY (log-log, overlaid)
# ===========================================================================

def plot_spectra(J_rec: xr.DataArray, J_obs: xr.DataArray) -> None:
    """
    Log-log power spectral density of reconstruction vs MULTIOBS global integral.
    X-axis: period in years. Uses Welch's method (nperseg = 1/3 of record length).
    """
    fs   = 12.0    # monthly → 12 samples per year
    rec  = J_rec.values
    obs  = J_obs.values

    # Remove NaN
    ok   = np.isfinite(rec) & np.isfinite(obs)
    rec  = rec[ok] - np.mean(rec[ok])
    obs  = obs[ok] - np.mean(obs[ok])

    nperseg = max(len(rec) // 3, 24)   # at least 2 years

    f_rec, Pxx_rec = welch(rec, fs=fs, nperseg=nperseg, scaling="density")
    f_obs, Pxx_obs = welch(obs, fs=fs, nperseg=nperseg, scaling="density")

    # Convert frequency (cycles/yr) to period (yr); skip DC (f=0)
    period_rec = 1.0 / f_rec[1:]
    period_obs = 1.0 / f_obs[1:]
    Pxx_rec    = Pxx_rec[1:]
    Pxx_obs    = Pxx_obs[1:]

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.loglog(period_rec, Pxx_rec,
              color="steelblue", lw=1.8, label=REC_LABEL)
    ax.loglog(period_obs, Pxx_obs,
              color="firebrick", lw=1.8, ls="--", label=OBS_LABEL)

    # Reference slope lines
    p_ref = np.array([0.1, 20.0])
    ax.loglog(p_ref, 1e-3 * p_ref**(-5/3),
              color="gray", lw=0.8, ls=":", label="−5/3 slope")

    ax.set_xlabel("Period  [years]")
    ax.set_ylabel("PSD  [( Pg C yr⁻¹ )² yr]")
    ax.set_title("Power spectral density — global net ocean CO₂ uptake")
    ax.legend()
    ax.grid(True, which="both", alpha=0.2)
    ax.invert_xaxis()   # long periods on the left
    fig.tight_layout()

    out = cfg.FIG_DIR / "fig_spectra.png"
    fig.savefig(out, dpi=cfg.FIGURE_DPI)
    plt.close(fig)
    print(f"[fig] {out}")


# ===========================================================================
# FIGURE 6 — CROSS-POWER SPECTRUM (gain, phase, coherence)
# ===========================================================================

def plot_cps(J_rec: xr.DataArray, J_obs: xr.DataArray) -> None:
    """
    Cross-power spectrum: MULTIOBS as reference, reconstruction as target.
    Three panels: gain, phase [days], squared coherence.
    High nfft for smooth curves.
    """
    fs      = 12.0    # samples per year
    rec     = J_rec.values
    obs     = J_obs.values
    ok      = np.isfinite(rec) & np.isfinite(obs)
    rec     = rec[ok] - np.mean(rec[ok])
    obs     = obs[ok] - np.mean(obs[ok])

    n       = len(rec)
    nfft    = 2 ** int(np.ceil(np.log2(n)) + 2)   # next power of 2 × 4 — very smooth
    nperseg = min(n, nfft // 2)
    noverlap= nperseg // 2

    # Auto-spectra and cross-spectrum
    f, Pxx = welch(obs, fs=fs, nperseg=nperseg, noverlap=noverlap,
                   nfft=nfft, scaling="density")
    _, Pyy = welch(rec, fs=fs, nperseg=nperseg, noverlap=noverlap,
                   nfft=nfft, scaling="density")
    _, Pxy = csd(obs, rec, fs=fs, nperseg=nperseg, noverlap=noverlap,
                 nfft=nfft, scaling="density")

    # Skip DC
    f    = f[1:]
    Pxx  = Pxx[1:]
    Pyy  = Pyy[1:]
    Pxy  = Pxy[1:]

    period     = 1.0 / f               # years
    period_days= period * 365.25       # days

    gain       = np.abs(Pxy) / Pxx                        # dimensionless
    phase_days = np.angle(Pxy) / (2 * np.pi * f) * 365.25 # days
    coherence  = (np.abs(Pxy)**2) / (Pxx * Pyy)           # 0–1

    fig, axes = plt.subplots(3, 1, figsize=(11, 10), sharex=True)

    for ax in axes:
        ax.set_xscale("log")
        ax.grid(True, which="both", alpha=0.2)
        ax.invert_xaxis()

    axes[0].semilogx(period, gain, color="steelblue", lw=1.5)
    axes[0].axhline(1, color="gray", lw=0.8, ls="--")
    axes[0].set_ylabel("Gain")
    axes[0].set_title(
        f"Cross-power spectrum — {REC_LABEL}  vs  {OBS_LABEL}\n"
        f"(reference: {OBS_LABEL})"
    )

    axes[1].semilogx(period, phase_days, color="darkorange", lw=1.5)
    axes[1].axhline(0, color="gray", lw=0.8, ls="--")
    axes[1].set_ylabel("Phase  [days]")

    axes[2].semilogx(period, coherence, color="seagreen", lw=1.5)
    axes[2].axhline(0.5, color="gray", lw=0.8, ls="--", label="0.5 threshold")
    axes[2].set_ylim(0, 1.05)
    axes[2].set_ylabel("Squared coherence")
    axes[2].set_xlabel("Period  [years]")
    axes[2].legend(fontsize=8)

    # Mark annual and semi-annual periods
    for ax in axes:
        for p, lbl in [(1.0, "Annual"), (0.5, "Semi-annual")]:
            ax.axvline(p, color="red", lw=0.7, ls=":", alpha=0.6)
            ax.text(p, ax.get_ylim()[1] * 0.95, lbl,
                    fontsize=7, color="red", ha="center", va="top")

    fig.tight_layout()
    out = cfg.FIG_DIR / "fig_cps.png"
    fig.savefig(out, dpi=cfg.FIGURE_DPI)
    plt.close(fig)
    print(f"[fig] {out}")


# ===========================================================================
# MAIN
# ===========================================================================

def main():
    print("[load] Data ...")
    fgco2_rec, fgco2_obs, J_rec, J_obs, ds_surf = load_data()

    # --- Global scalar metrics ---
    print("[compute] Global scalar metrics ...")
    metrics = scalar_metrics(fgco2_rec, fgco2_obs)
    print("\n  Global validation metrics:")
    for k, v in metrics.items():
        print(f"    {k:30s}: {v:.4f}")

    # --- Spatial metrics ---
    print("\n[compute] Pixel RMSD and bias ...")
    rmsd = pixel_rmsd(fgco2_rec, fgco2_obs)
    bias = pixel_bias(fgco2_rec, fgco2_obs)

    # --- Figures ---
    print("\n[plot] Time series ...")
    plot_timeseries(J_rec, J_obs)

    print("[plot] LOESS multi-year variability ...")
    plot_loess(J_rec, J_obs, frac=0.15)

    print("[plot] Validation map ...")
    plot_validation_map(rmsd, bias)

    print("[plot] Fay domain validation ...")
    domain_metrics = plot_fay_domains(fgco2_rec, fgco2_obs)

    print("[plot] Power spectra ...")
    plot_spectra(J_rec, J_obs)

    print("[plot] Cross-power spectrum ...")
    plot_cps(J_rec, J_obs)

    # --- Save metrics ---
    all_metrics = {"global": metrics}
    all_metrics.update({f"domain_{k}": v for k, v in domain_metrics.items()})
    pd.DataFrame(all_metrics).T.to_csv(cfg.OUT_DIR / "validation_metrics.csv")
    print(f"\n[save] validation_metrics.csv")

    print(f"\n[done] Validation complete. Figures → {cfg.FIG_DIR}\n")


if __name__ == "__main__":
    main()
