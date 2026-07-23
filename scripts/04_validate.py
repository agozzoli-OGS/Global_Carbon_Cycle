"""
04_validate.py
==============
Cross-validation of BOTH flux reconstructions against the CMEMS MULTIOBS
observation-based flux (SOCAT neural network).

Produces three comparison figures:

    fig_validation_ts.png
        THREE time series of global net ocean CO2 uptake [Pg C yr⁻¹]:
        (1) Original reconstruction  (v1.0.x — PISCES pCO2, monthly wind)
        (2) Improved reconstruction  (v1.1.0 — MULTIOBS pCO2 + wind variance)
        (3) CMEMS MULTIOBS           (observation-based reference)

    fig_validation_map.png
        2×2 panel map (contourf, fixed colourbars):
        Left column  — original reconstruction:  RMSD (top), Bias (bottom)
        Right column — improved reconstruction:  RMSD (top), Bias (bottom)
        RMSD colourbar: 0 – 5 mol C m⁻² yr⁻¹
        Bias colourbar: −5 – +5 mol C m⁻² yr⁻¹

    fig_taylor.png
        Taylor diagram with TWO markers:
        ★  MULTIOBS (reference)
        ●  Original reconstruction
        ▲  Improved reconstruction

    validation_metrics.csv
        Scalar statistics (RMSD, bias, r, σ) for both reconstructions.

All figures saved at 300 dpi.

Changelog:
    v1.1.0 — Added improved reconstruction to all three figures;
              maps changed to 2×2 contourf with fixed colourbars;
              time series extended to three curves;
              Taylor diagram extended to two model markers;
              DPI raised to 300.

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

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    HAS_CARTOPY = True
except ImportError:
    HAS_CARTOPY = False
    print("[warn] cartopy not found — maps will use plain imshow.")


# ===========================================================================
# LOAD DATA
# ===========================================================================

def load_data():
    """
    Load both reconstructed fluxes, the MULTIOBS validation flux, and
    the global integral time series.

    Returns
    -------
    fgco2_rec : xr.DataArray      original reconstruction [mol C m⁻² yr⁻¹]
    fgco2_imp : xr.DataArray|None improved reconstruction [mol C m⁻² yr⁻¹]
    fgco2_obs : xr.DataArray      MULTIOBS validation flux [mol C m⁻² yr⁻¹]
    J_orig    : xr.DataArray      original global integral [Pg C yr⁻¹]
    J_imp     : xr.DataArray|None improved global integral [Pg C yr⁻¹]
    """
    ds_flux = xr.open_dataset(cfg.DATA_DIR / "flux_3d.nc",        chunks="auto")
    ds_surf = xr.open_dataset(cfg.DATA_DIR / "processed_surface.nc", chunks="auto")
    ds_glob = xr.open_dataset(cfg.OUT_DIR  / "global_flux.nc")

    if "fgco2_obs" not in ds_surf:
        raise FileNotFoundError(
            "fgco2_obs not found in processed_surface.nc. "
            "MULTIOBS download may have been skipped (see 01_download_data.py)."
        )

    fgco2_rec = ds_flux["fgco2"]
    fgco2_imp = ds_flux["fgco2_improved"] if "fgco2_improved" in ds_flux else None
    # MULTIOBS sign convention is opposite — flip to match our reconstruction
    fgco2_obs = -ds_surf["fgco2_obs"]

    J_orig = ds_glob["J_net_PgC"]
    J_imp  = ds_glob["J_net_improved_PgC"] if "J_net_improved_PgC" in ds_glob else None

    # Align all to the common time window where MULTIOBS is available
    common_time = fgco2_rec.time[
        np.isin(fgco2_rec.time.values, fgco2_obs.time.values)
    ]
    fgco2_rec = fgco2_rec.sel(time=common_time)
    fgco2_obs = fgco2_obs.sel(time=common_time)
    if fgco2_imp is not None:
        fgco2_imp = fgco2_imp.sel(time=common_time)

    return fgco2_rec, fgco2_imp, fgco2_obs, J_orig, J_imp


# ===========================================================================
# SKILL METRICS
# ===========================================================================

def pixel_rmsd(rec: xr.DataArray, obs: xr.DataArray) -> xr.DataArray:
    """Per-pixel RMSD over time: sqrt(mean[(rec − obs)²])."""
    rmsd = np.sqrt(((rec - obs)**2).mean(dim="time"))
    rmsd.attrs = {"long_name": "RMSD", "units": "mol C m-2 yr-1"}
    return rmsd


def pixel_bias(rec: xr.DataArray, obs: xr.DataArray) -> xr.DataArray:
    """Per-pixel mean bias: mean(rec − obs)."""
    bias = (rec - obs).mean(dim="time")
    bias.attrs = {"long_name": "Bias (rec − obs)", "units": "mol C m-2 yr-1"}
    return bias


def scalar_metrics(rec: xr.DataArray, obs: xr.DataArray, label: str) -> dict:
    """Global scalar skill metrics over all ocean pixels and all times."""
    r = rec.values.ravel()
    o = obs.values.ravel()
    ok = np.isfinite(r) & np.isfinite(o)
    r, o = r[ok], o[ok]
    bias = float(np.mean(r - o))
    rmsd = float(np.sqrt(np.mean((r - o)**2)))
    corr, _ = stats.pearsonr(r, o)
    return {
        f"{label}_bias_mol_C_m2_yr": bias,
        f"{label}_rmsd_mol_C_m2_yr": rmsd,
        f"{label}_pearson_r":        corr,
        f"{label}_std_rec":          float(np.std(r)),
        f"{label}_std_obs":          float(np.std(o)),
    }


# ===========================================================================
# GRID-AREA HELPER (duplicated here for independence from 03_compute_flux.py)
# ===========================================================================

def compute_grid_cell_area(lat: xr.DataArray, lon: xr.DataArray) -> xr.DataArray:
    R         = cfg.EARTH_RADIUS_M
    d_lat     = float(np.abs(lat.diff("latitude").mean()))
    d_lon     = float(np.abs(lon.diff("longitude").mean()))
    area_1d   = (R**2) * np.deg2rad(d_lon) * np.deg2rad(d_lat) * np.cos(np.deg2rad(lat))
    area_2d   = area_1d.expand_dims({"longitude": lon}).transpose("latitude", "longitude")
    area_2d.attrs = {"long_name": "Grid cell area", "units": "m2"}
    return area_2d.rename("cell_area")


def global_integral_obs(obs: xr.DataArray, mask: xr.DataArray) -> xr.DataArray:
    """Global integral of MULTIOBS flux [Pg C yr⁻¹], positive = uptake."""
    cell_area = compute_grid_cell_area(obs["latitude"], obs["longitude"])
    F_global  = (obs.where(mask == 1) * cell_area).sum(
        dim=["latitude", "longitude"], skipna=True
    )
    J = F_global * cfg.MOL_C_TO_PG
    J.attrs = {"long_name": "MULTIOBS global CO2 uptake", "units": "Pg C yr-1"}
    return J


# ===========================================================================
# FIGURE 1 — THREE-CURVE TIME SERIES
# ===========================================================================

def plot_timeseries(
    J_orig: xr.DataArray,
    J_imp:  xr.DataArray | None,
    J_obs:  xr.DataArray,
) -> None:
    """
    Global net ocean CO2 uptake [Pg C yr⁻¹] — three time series:
        blue  solid   — original reconstruction  (v1.0.x)
        green solid   — improved reconstruction  (v1.1.0)
        red   dashed  — CMEMS MULTIOBS reference
    """
    fig, ax = plt.subplots(figsize=(14, 6))

    t_orig = pd.to_datetime(J_orig.time.values)
    t_obs  = pd.to_datetime(J_obs.time.values)

    ax.plot(t_orig, J_orig.values,
            color="steelblue", lw=1.4,
            label="Original reconstruction  (PISCES pCO₂, monthly wind)")

    if J_imp is not None:
        t_imp = pd.to_datetime(J_imp.time.values)
        ax.plot(t_imp, J_imp.values,
                color="seagreen", lw=1.4,
                label="Improved reconstruction  (MULTIOBS pCO₂ + wind variance corr.)")

    ax.plot(t_obs, J_obs.values,
            color="firebrick", lw=1.4, ls="--",
            label="CMEMS MULTIOBS  (SOCAT-NN, reference)")

    ax.axhline(0, color="black", lw=0.8, ls=":")
    ax.set_xlabel("Year", fontsize=12)
    ax.set_ylabel("Net ocean CO₂ uptake  [Pg C yr⁻¹]", fontsize=12)
    ax.set_title(
        "Global ocean CO₂ uptake: original vs improved reconstruction vs observation-based reference",
        fontsize=12,
    )
    ax.legend(fontsize=10, loc="upper left")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    out = cfg.FIG_DIR / "fig_validation_ts.png"
    fig.savefig(out, dpi=cfg.FIGURE_DPI)
    plt.close(fig)
    print(f"[fig] {out}")


# ===========================================================================
# FIGURE 2 — 2×2 CONTOURF COMPARISON MAP
# ===========================================================================

def plot_comparison_maps(
    rmsd_orig: xr.DataArray,
    bias_orig: xr.DataArray,
    rmsd_imp:  xr.DataArray | None,
    bias_imp:  xr.DataArray | None,
) -> None:
    """
    2×2 panel map comparing original (left) and improved (right) reconstructions.
    Top row: RMSD  [0 – VAL_RMSD_MAX mol C m⁻² yr⁻¹]   — colourmap Reds
    Bottom row: Bias [VAL_BIAS_MIN – VAL_BIAS_MAX]       — colourmap RdBu_r
    Both rows share the same fixed colourbar extent so results are directly
    comparable between the two reconstructions.
    All panels use contourf for a clean publication-quality look.
    """
    ncols = 2 if (rmsd_imp is not None and bias_imp is not None) else 1
    titles_top    = ["RMSD — Original  [mol C m⁻² yr⁻¹]",
                     "RMSD — Improved  [mol C m⁻² yr⁻¹]"]
    titles_bottom = ["Bias (rec − obs) — Original  [mol C m⁻² yr⁻¹]",
                     "Bias (rec − obs) — Improved  [mol C m⁻² yr⁻¹]"]

    rmsd_levels = np.linspace(0,                 cfg.VAL_RMSD_MAX, 21)
    bias_levels = np.linspace(cfg.VAL_BIAS_MIN,  cfg.VAL_BIAS_MAX, 21)

    fig = plt.figure(figsize=(7 * ncols + 1, 9))
    gs  = gridspec.GridSpec(
        2, ncols,
        hspace=0.35, wspace=0.15,
        left=0.05, right=0.88,
    )

    # Colour-bar axes: one per row, placed to the right of the grid
    cax_rmsd = fig.add_axes([0.90, 0.55, 0.02, 0.35])
    cax_bias = fig.add_axes([0.90, 0.10, 0.02, 0.35])

    def _make_panel(ax, data, levels, cmap, title):
        """Draw one contourf panel on ax (with or without cartopy)."""
        lat = data["latitude"].values
        lon = data["longitude"].values
        Z   = data.values

        if HAS_CARTOPY:
            ax2 = fig.add_subplot(ax.get_subplotspec(),
                                  projection=ccrs.Robinson())
            ax.remove()
            cf = ax2.contourf(
                lon, lat, Z,
                levels=levels, cmap=cmap, extend="both",
                transform=ccrs.PlateCarree(),
            )
            ax2.add_feature(cfeature.LAND,      facecolor="lightgray", zorder=3)
            ax2.add_feature(cfeature.COASTLINE, lw=0.3,                zorder=4)
            ax2.set_global()
            ax2.set_title(title, fontsize=10, pad=6)
            return cf, ax2
        else:
            cf = ax.contourf(lon, lat, Z, levels=levels, cmap=cmap, extend="both")
            ax.set_title(title, fontsize=10)
            ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude")
            return cf, ax

    # ---- Row 0: RMSD ----
    datasets_rmsd = [rmsd_orig] + ([rmsd_imp] if rmsd_imp is not None else [])
    cf_rmsd = None
    for col, (da, title) in enumerate(zip(datasets_rmsd, titles_top)):
        ax = fig.add_subplot(gs[0, col])
        cf_rmsd, _ = _make_panel(ax, da, rmsd_levels, "Reds", title)
    plt.colorbar(cf_rmsd, cax=cax_rmsd, label="RMSD  [mol C m⁻² yr⁻¹]")

    # ---- Row 1: Bias ----
    datasets_bias = [bias_orig] + ([bias_imp] if bias_imp is not None else [])
    cf_bias = None
    for col, (da, title) in enumerate(zip(datasets_bias, titles_bottom)):
        ax = fig.add_subplot(gs[1, col])
        cf_bias, _ = _make_panel(ax, da, bias_levels, cfg.CMAP_FLUX, title)
    plt.colorbar(cf_bias, cax=cax_bias, label="Bias (rec − obs)  [mol C m⁻² yr⁻¹]")

    out = cfg.FIG_DIR / "fig_validation_map.png"
    fig.savefig(out, dpi=cfg.FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {out}")


# ===========================================================================
# FIGURE 3 — TAYLOR DIAGRAM (two model markers)
# ===========================================================================

def plot_taylor_diagram(
    metrics_orig: dict,
    metrics_imp:  dict | None,
) -> None:
    """
    Taylor diagram with up to two model markers.
        ★  MULTIOBS (reference, r=1, std_ratio=1)
        ●  Original reconstruction
        ▲  Improved reconstruction (if available)

    Angle = arccos(r)  |  Radius = σ_rec / σ_obs (normalised std)
    """
    fig, ax = plt.subplots(figsize=(7, 7), subplot_kw={"projection": "polar"})

    def _plot_point(r, std_rec, std_obs, marker, color, label):
        theta = np.arccos(np.clip(r, -1, 1))
        rho   = std_rec / std_obs
        ax.plot(theta, rho, marker=marker, color=color, ms=12,
                linestyle="none", label=label)
        return rho

    # Reference
    ax.plot(0, 1, marker="*", color="black", ms=16,
            linestyle="none", label="MULTIOBS (reference)")

    r_orig = metrics_orig["orig_pearson_r"]
    rho_orig = _plot_point(
        r_orig,
        metrics_orig["orig_std_rec"],
        metrics_orig["orig_std_obs"],
        marker="o", color="steelblue",
        label=f"Original  (r = {r_orig:.3f})",
    )

    rho_max = rho_orig
    if metrics_imp is not None:
        r_imp = metrics_imp["imp_pearson_r"]
        rho_imp = _plot_point(
            r_imp,
            metrics_imp["imp_std_rec"],
            metrics_imp["imp_std_obs"],
            marker="^", color="seagreen",
            label=f"Improved  (r = {r_imp:.3f})",
        )
        rho_max = max(rho_max, rho_imp)

    ax.set_thetamax(90)
    ax.set_rlim(0, max(1.5, rho_max * 1.15))
    ax.set_xlabel("Normalised standard deviation", labelpad=20, fontsize=11)
    ax.set_title("Taylor diagram", pad=20, fontsize=12)
    ax.legend(loc="upper right", bbox_to_anchor=(1.45, 1.15), fontsize=10)

    out = cfg.FIG_DIR / "fig_taylor.png"
    fig.savefig(out, dpi=cfg.FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {out}")


# ===========================================================================
# MAIN
# ===========================================================================

def main():
    fgco2_rec, fgco2_imp, fgco2_obs, J_orig, J_imp = load_data()

    ds_surf = xr.open_dataset(cfg.DATA_DIR / "processed_surface.nc")

    # --- MULTIOBS global integral ---
    J_obs = global_integral_obs(fgco2_obs, ds_surf["ocean_mask"])

    # --- Spatial metrics ---
    print("[compute] Pixel RMSD and bias — original ...")
    rmsd_orig = pixel_rmsd(fgco2_rec, fgco2_obs)
    bias_orig = pixel_bias(fgco2_rec, fgco2_obs)

    rmsd_imp = bias_imp = None
    if fgco2_imp is not None:
        print("[compute] Pixel RMSD and bias — improved ...")
        rmsd_imp = pixel_rmsd(fgco2_imp, fgco2_obs)
        bias_imp = pixel_bias(fgco2_imp, fgco2_obs)

    # --- Scalar metrics ---
    print("[compute] Scalar metrics ...")
    m_orig = scalar_metrics(fgco2_rec, fgco2_obs, "orig")
    all_metrics = {**m_orig}
    if fgco2_imp is not None:
        m_imp = scalar_metrics(fgco2_imp, fgco2_obs, "imp")
        all_metrics.update(m_imp)

    print("\n  Validation metrics:")
    for k, v in all_metrics.items():
        print(f"    {k:40s}: {v:.4f}")

    pd.Series(all_metrics).to_csv(cfg.OUT_DIR / "validation_metrics.csv")
    print(f"\n[save] validation_metrics.csv")

    # --- Figures ---
    print("\n[plot] Three-curve time series ...")
    plot_timeseries(J_orig, J_imp, J_obs)

    print("[plot] 2×2 comparison maps ...")
    plot_comparison_maps(rmsd_orig, bias_orig, rmsd_imp, bias_imp)

    print("[plot] Taylor diagram ...")
    # Pass separate dicts with short-keyed names for the Taylor function
    taylor_orig = {
        "orig_pearson_r": m_orig["orig_pearson_r"],
        "orig_std_rec":   m_orig["orig_std_rec"],
        "orig_std_obs":   m_orig["orig_std_obs"],
    }
    taylor_imp = None
    if fgco2_imp is not None:
        taylor_imp = {
            "imp_pearson_r": m_imp["imp_pearson_r"],
            "imp_std_rec":   m_imp["imp_std_rec"],
            "imp_std_obs":   m_imp["imp_std_obs"],
        }
    plot_taylor_diagram(taylor_orig, taylor_imp)

    print(f"\n[done] Validation complete. Figures in {cfg.FIG_DIR}\n")


if __name__ == "__main__":
    main()
