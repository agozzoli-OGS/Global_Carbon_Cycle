"""
04_validate.py
==============
Cross-validation of the reconstructed air-sea CO2 flux (from 03_compute_flux.py)
against the observation-based CMEMS MULTIOBS surface carbon L4 product
(MULTIOBS_GLO_BGC_CARBON_SURFACE_MYNRT_015_008).

This step is critical: it tells us whether the numerical trend we will
report in the analysis is a real signal or a model artifact introduced by
the PISCES/FREEGLORYS2V4 physics chain.

Metrics computed:
    - RMSD (root mean squared difference) of the flux field, globally and by basin
    - Bias (mean difference: our reconstruction minus MULTIOBS)
    - Pearson r (temporal correlation at each pixel, and globally)
    - Taylor diagram statistics
    - Time series comparison of the global integral

Outputs:
    output/validation_metrics.csv      — scalar statistics
    output/figures/fig_validation_ts.png   — time series comparison
    output/figures/fig_validation_map.png  — spatial RMSD map
    output/figures/fig_taylor.png          — Taylor diagram

Usage:
    python scripts/04_validate.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent.parent))
import config as cfg


# ===========================================================================
# LOAD DATA
# ===========================================================================

def load_data() -> tuple[xr.DataArray, xr.DataArray, xr.DataArray]:
    """
    Load the reconstructed flux, the MULTIOBS validation flux, and the
    global integral time series. Return them aligned on the same time axis.

    Returns
    -------
    fgco2_rec : xr.DataArray
        Our reconstructed flux [mol C m⁻² yr⁻¹].
    fgco2_obs : xr.DataArray
        MULTIOBS observation-based flux [mol C m⁻² yr⁻¹].
    J_net_rec : xr.DataArray
        Global integral of reconstructed flux [Pg C yr⁻¹].
    """
    ds_flux  = xr.open_dataset(cfg.DATA_DIR / "flux_3d.nc", chunks="auto")
    ds_surf  = xr.open_dataset(cfg.DATA_DIR / "processed_surface.nc", chunks="auto")
    ds_glob  = xr.open_dataset(cfg.OUT_DIR / "global_flux.nc")

    fgco2_rec = ds_flux["fgco2"]
    J_net_rec = ds_glob["J_net_PgC"]

    if "fgco2_obs" not in ds_surf:
        raise FileNotFoundError(
            "fgco2_obs not found in processed_surface.nc.\n"
            "MULTIOBS download may have been skipped (see 01_download_data.py)."
        )
    fgco2_obs = ds_surf["fgco2_obs"]

    # Align on common time
    common_time = fgco2_rec.time[
        np.isin(fgco2_rec.time.values, fgco2_obs.time.values)
    ]
    fgco2_rec = fgco2_rec.sel(time=common_time)
    fgco2_obs = fgco2_obs.sel(time=common_time)

    return fgco2_rec, fgco2_obs, J_net_rec


# ===========================================================================
# SKILL METRICS
# ===========================================================================

def pixel_rmsd(rec: xr.DataArray, obs: xr.DataArray) -> xr.DataArray:
    """
    Root mean squared difference at every (lat, lon) pixel over time.

        RMSD(lat, lon) = sqrt( mean_t[ (rec - obs)² ] )
    """
    diff = rec - obs
    rmsd = np.sqrt((diff**2).mean(dim="time"))
    rmsd.attrs = {"long_name": "RMSD reconstructed vs MULTIOBS", "units": "mol C m-2 yr-1"}
    return rmsd


def pixel_bias(rec: xr.DataArray, obs: xr.DataArray) -> xr.DataArray:
    """
    Mean bias at every pixel: mean_t(rec - obs). Positive = our flux is too high."""
    bias = (rec - obs).mean(dim="time")
    bias.attrs = {"long_name": "Bias reconstructed - MULTIOBS", "units": "mol C m-2 yr-1"}
    return bias


def global_integral_obs(obs: xr.DataArray, mask: xr.DataArray) -> xr.DataArray:
    """
    Compute global integral of the MULTIOBS flux [Pg C yr⁻¹] for comparison.
    Uses the same grid cell area computation as 03_compute_flux.py.
    """
    from scripts_03 import compute_grid_cell_area

    cell_area = compute_grid_cell_area(obs["latitude"], obs["longitude"])
    F_global = (obs.where(mask == 1) * cell_area).sum(
        dim=["latitude", "longitude"], skipna=True
    )
    J_net_obs = -F_global * cfg.MOL_C_TO_PG
    J_net_obs.attrs = {"long_name": "MULTIOBS global CO2 uptake", "units": "Pg C yr-1"}
    return J_net_obs


def scalar_metrics(rec: xr.DataArray, obs: xr.DataArray) -> dict:
    """
    Compute global scalar skill metrics over the common time–space domain.

    Returns a dict with: bias, rmsd, r (Pearson), std_rec, std_obs.
    """
    r   = rec.values.ravel()
    o   = obs.values.ravel()
    ok  = np.isfinite(r) & np.isfinite(o)
    r, o = r[ok], o[ok]

    bias = float(np.mean(r - o))
    rmsd = float(np.sqrt(np.mean((r - o)**2)))
    corr, _ = stats.pearsonr(r, o)
    std_rec = float(np.std(r))
    std_obs = float(np.std(o))

    return {
        "bias_mol_C_m2_yr": bias,
        "rmsd_mol_C_m2_yr": rmsd,
        "pearson_r":        corr,
        "std_rec":          std_rec,
        "std_obs":          std_obs,
    }


# ===========================================================================
# PLOTS
# ===========================================================================

def plot_timeseries(
    J_net_rec: xr.DataArray,
    J_net_obs: xr.DataArray,
) -> None:
    """
    Figure: time series of global net ocean CO2 uptake — reconstructed vs MULTIOBS.

    Both curves in Pg C yr⁻¹ (positive = uptake).
    """
    fig, ax = plt.subplots(figsize=(12, 5))

    t_rec = pd.to_datetime(J_net_rec.time.values)
    t_obs = pd.to_datetime(J_net_obs.time.values)

    ax.plot(t_rec, J_net_rec.values, color="steelblue",  lw=1.5,
            label="This study (Wanninkhof 2014 + PISCES pCO₂)")
    ax.plot(t_obs, J_net_obs.values, color="firebrick",  lw=1.5, ls="--",
            label="CMEMS MULTIOBS (SOCAT NN)")

    ax.axhline(0, color="black", lw=0.8, ls=":")
    ax.set_xlabel("Year", fontsize=12)
    ax.set_ylabel("Net ocean CO₂ uptake  [Pg C yr⁻¹]", fontsize=12)
    ax.set_title("Global ocean CO₂ uptake: reconstructed vs observation-based validation product",
                 fontsize=12)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    out = cfg.FIG_DIR / "fig_validation_ts.png"
    fig.savefig(out, dpi=cfg.FIGURE_DPI)
    plt.close(fig)
    print(f"[fig] {out}")


def plot_rmsd_map(rmsd: xr.DataArray, bias: xr.DataArray) -> None:
    """
    Figure: 2-panel map of pixel-wise RMSD and bias between reconstructed
    and MULTIOBS flux. Highlights where the reconstruction diverges most.
    """
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    fig = plt.figure(figsize=(14, 8))
    gs  = gridspec.GridSpec(2, 1, hspace=0.3)

    proj = ccrs.Robinson()

    for i, (da, title, cmap, label) in enumerate([
        (rmsd, "RMSD  [mol C m⁻² yr⁻¹]", "Reds",   "RMSD"),
        (bias, "Bias (rec − obs)  [mol C m⁻² yr⁻¹]", cfg.CMAP_FLUX, "Bias"),
    ]):
        ax = fig.add_subplot(gs[i], projection=proj)
        im = da.plot(
            ax=ax, transform=ccrs.PlateCarree(),
            cmap=cmap, add_colorbar=True,
            cbar_kwargs={"label": label, "shrink": 0.6},
        )
        ax.add_feature(cfeature.LAND,  facecolor="lightgray", zorder=2)
        ax.add_feature(cfeature.COASTLINE, lw=0.4, zorder=3)
        ax.set_title(title, fontsize=11)
        ax.set_global()

    out = cfg.FIG_DIR / "fig_validation_map.png"
    fig.savefig(out, dpi=cfg.FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {out}")


def plot_taylor_diagram(metrics: dict) -> None:
    """
    Taylor diagram showing correlation, standard deviation ratio, and RMSD
    between reconstructed and MULTIOBS flux (global pixel-wise statistics).

    The target point is the MULTIOBS product (normalised std = 1, r = 1).
    Our reconstruction is plotted at (r, std_rec/std_obs).
    """
    r       = metrics["pearson_r"]
    std_rec = metrics["std_rec"]
    std_obs = metrics["std_obs"]

    fig, ax = plt.subplots(figsize=(6, 6), subplot_kw={"projection": "polar"})

    # Reference (MULTIOBS = perfect point at r=1, std_ratio=1)
    ax.plot(0, 1, marker="*", color="black", ms=14, label="MULTIOBS (reference)")

    # Our reconstruction
    theta = np.arccos(r)   # angle = arccos(correlation)
    rho   = std_rec / std_obs
    ax.plot(theta, rho, marker="o", color="steelblue", ms=10, label="This study")

    ax.set_thetamax(90)
    ax.set_rlim(0, max(1.5, rho * 1.1))
    ax.set_xlabel("Normalised standard deviation", labelpad=20)
    ax.set_title(f"Taylor diagram  (r = {r:.3f})", pad=20)
    ax.legend(loc="upper right", bbox_to_anchor=(1.3, 1.1))

    out = cfg.FIG_DIR / "fig_taylor.png"
    fig.savefig(out, dpi=cfg.FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {out}")


# ===========================================================================
# MAIN
# ===========================================================================

def main():
    fgco2_rec, fgco2_obs, J_net_rec = load_data()

    # --- Spatial metrics ---
    print("[compute] Pixel RMSD and bias ...")
    rmsd = pixel_rmsd(fgco2_rec, fgco2_obs)
    bias = pixel_bias(fgco2_rec, fgco2_obs)

    # --- Scalar metrics ---
    print("[compute] Global scalar metrics ...")
    metrics = scalar_metrics(fgco2_rec, fgco2_obs)
    print("\n  Validation metrics (global, all ocean pixels):")
    for k, v in metrics.items():
        print(f"    {k:30s}: {v:.4f}")

    pd.Series(metrics).to_csv(cfg.OUT_DIR / "validation_metrics.csv")
    print(f"\n[save] validation_metrics.csv")

    # --- Global integral of MULTIOBS ---
    ds_surf  = xr.open_dataset(cfg.DATA_DIR / "processed_surface.nc")
    J_net_obs = global_integral_obs(fgco2_obs, ds_surf["ocean_mask"])

    # --- Plots ---
    print("\n[plot] Time series comparison ...")
    plot_timeseries(J_net_rec, J_net_obs)

    print("[plot] RMSD and bias maps ...")
    plot_rmsd_map(rmsd, bias)

    print("[plot] Taylor diagram ...")
    plot_taylor_diagram(metrics)

    print("\n[done] Validation complete.\n")


if __name__ == "__main__":
    main()
