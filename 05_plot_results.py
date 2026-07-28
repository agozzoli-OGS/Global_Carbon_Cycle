"""
05_plot_results.py
==================
Produces all main result figures for Stage 1. Reads from the computed
outputs in data/ and output/ — no computation happens here.

Figures produced:
    fig01_flux_vs_co2.png
        The primary figure from the meeting: global net ocean CO2 uptake
        (J_net, Pg C yr⁻¹) vs. time alongside atmospheric pCO2 (ppm).
        Shows whether the ocean sink is growing, saturating, or flat.

    fig02_annual_flux_map.png
        Mean spatial map of air-sea flux, averaged over the full record.
        Red = outgassing (source), blue = uptake (sink).

    fig03_flux_trend_map.png
        Pixel-wise linear trend (Sen's slope) of annual mean flux.
        Answers: where is the sink intensifying / weakening over time?

    fig04_delta_pco2_map.png
        Mean ΔpCO2 = pCO2_ocean − pCO2_atm map. The thermodynamic driver
        of the flux — positive regions are always potential sources.

    fig05_monthly_seasonal_cycle.png
        Global mean seasonal cycle (climatological monthly anomaly) of flux
        and atmospheric pCO2. Shows the competing seasonal signals.

Usage:
    python scripts/05_plot_results.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

sys.path.insert(0, str(Path(__file__).parent.parent))
import config as cfg

# Optional: cartopy for map projections
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

def load_all() -> dict:
    """Load all required output files and return as a dict of DataArrays."""
    ds_flux = xr.open_dataset(cfg.DATA_DIR / "flux_3d.nc")
    ds_glob = xr.open_dataset(cfg.OUT_DIR / "global_flux.nc")
    ds_surf = xr.open_dataset(cfg.DATA_DIR / "processed_surface.nc")

    # Atmospheric CO2 in ppm for the right-hand axis of fig01
    noaa = pd.read_csv(
        cfg.NOAA_CO2_FILE,
        comment="#",
        names=[
            "year",
            "month",
            "decimal",
            "average",
            "average_unc",
            "trend",
            "trend_unc",
        ],
    )

    # Fix NOAA column type
    noaa["average"] = pd.to_numeric(noaa["average"], errors="coerce")

    # Remove missing values
    noaa = noaa[noaa["average"] > 0]

    noaa["time"] = pd.to_datetime(
        noaa["year"].astype(str)
        + "-"
        + noaa["month"].astype(str).str.zfill(2)
    )
    noaa = noaa.set_index("time").sort_index()

    return {
        "flux_3d":      ds_flux["fgco2"],
        "J_net":        ds_glob["J_net_PgC"],
        "spco2_ocean":  ds_surf["spco2_ocean"],
        "spco2_atm":    ds_surf["spco2_atm"],
        "ocean_mask":   ds_surf["ocean_mask"],
        "noaa_co2_ppm": noaa["average"],
    }


# ===========================================================================
# FIG 01 — GLOBAL FLUX VS ATMOSPHERIC CO2 (the "key figure" from the meeting)
# ===========================================================================

def fig01_flux_vs_co2(data: dict) -> None:
    """
    Twin-axis time series:
        Left  axis:  J_net [Pg C yr⁻¹] — annual mean global ocean uptake
        Right axis:  Atmospheric CO2 [ppm] — NOAA GML global mean

    This is the central figure requested in the meeting:
    'I want a graph where I put atmospheric CO2 on one axis and the net
     flux on the other — is the ocean uptake growing, flat, or saturating?'

    Monthly values are shown in pale colours; annual means as solid lines
    to reduce visual noise while preserving sub-annual variability context.
    """
    fig, ax1 = plt.subplots(figsize=(14, 6))

    J_net   = data["J_net"]
    noaa    = data["noaa_co2_ppm"]

    time_monthly = pd.to_datetime(J_net.time.values)

    # Annual mean of J_net
    J_annual = (
        xr.DataArray(J_net.values, coords={"time": time_monthly}, dims="time")
        .resample(time="1YE").mean()
    )
    t_annual = pd.to_datetime(J_annual.time.values)

    # Align NOAA to the same annual axis
    noaa_annual = noaa.resample("YE").mean()

    # --- Left axis: flux ---
    ax1.fill_between(
        time_monthly, J_net.values, 0,
        alpha=0.2, color="steelblue", label="_nolegend_",
    )
    ax1.plot(time_monthly, J_net.values, color="steelblue",
             lw=0.8, alpha=0.5, label="_nolegend_")
    ax1.plot(t_annual, J_annual.values, color="steelblue",
             lw=2.2, label="Net ocean uptake J_net (annual mean)")

    ax1.axhline(0, color="black", lw=0.8, ls=":")
    ax1.set_xlabel("Year", fontsize=12)
    ax1.set_ylabel("Net ocean CO₂ uptake  [Pg C yr⁻¹]\n(positive = ocean sink)",
                   color="steelblue", fontsize=12)
    ax1.tick_params(axis="y", labelcolor="steelblue")
    ax1.set_ylim(-5, 5)

    # --- Right axis: atmospheric CO2 ---
    ax2 = ax1.twinx()
    ax2.plot(noaa_annual.index, noaa_annual.values, color="firebrick",
             lw=2.2, ls="--", label="Atmospheric CO₂ (NOAA GML)")
    ax2.set_ylabel("Atmospheric CO₂  [ppm]", color="firebrick", fontsize=12)
    ax2.tick_params(axis="y", labelcolor="firebrick")

    # --- Legend ---
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper left", fontsize=10)

    ax1.set_title(
        "Global ocean CO₂ uptake vs. atmospheric CO₂ (1993–2026)\n"
        "Wanninkhof (2014) k  ·  Weiss (1974) K0  ·  PISCES surface pCO₂",
        fontsize=12,
    )
    ax1.grid(True, alpha=0.3)
    fig.tight_layout()

    out = cfg.FIG_DIR / "fig01_flux_vs_co2.png"
    fig.savefig(out, dpi=cfg.FIGURE_DPI)
    plt.close(fig)
    print(f"[fig] {out}")


# ===========================================================================
# FIG 02 — MEAN SPATIAL FLUX MAP
# ===========================================================================

def fig02_mean_flux_map(data: dict) -> None:
    """
    Map of the time-mean air-sea CO2 flux [mol C m⁻² yr⁻¹].
    Red = outgassing (source), blue = ocean uptake (sink).

    Classic result: equatorial Pacific is red (upwelling brings CO2-rich
    deep water to surface); Southern Ocean and North Atlantic are intense blue.
    """
    flux_mean = data["flux_3d"].mean(dim="time").where(data["ocean_mask"] == 1)

    vmax = float(np.nanpercentile(np.abs(flux_mean.values), 97))

    if HAS_CARTOPY:
        fig, ax = plt.subplots(
            figsize=(14, 7),
            subplot_kw={"projection": ccrs.Robinson()},
        )
        im = flux_mean.plot(
            ax=ax, transform=ccrs.PlateCarree(),
            cmap=cfg.CMAP_FLUX, vmin=-vmax, vmax=vmax,
            add_colorbar=True,
            cbar_kwargs={
                "label": "Air-sea CO₂ flux  [mol C m⁻² yr⁻¹]\n(+ = uptake, − = outgassing)",
                "shrink": 0.6,
                "orientation": "horizontal",
                "pad": 0.04,
            },
        )
        ax.add_feature(cfeature.LAND, facecolor="lightgray", zorder=2)
        ax.add_feature(cfeature.COASTLINE, lw=0.4, zorder=3)
        ax.set_global()
    else:
        fig, ax = plt.subplots(figsize=(14, 6))
        im = ax.imshow(
            flux_mean.values[::-1],   # flip lat axis
            cmap=cfg.CMAP_FLUX, vmin=-vmax, vmax=vmax,
            extent=[-180, 180, -90, 90], aspect="auto",
        )
        plt.colorbar(im, ax=ax, label="Air-sea CO₂ flux [mol C m⁻² yr⁻¹]")

    ax.set_title(
        "Time-mean air-sea CO₂ flux (full record)\n"
        "Blue = uptake (sink), Red = outgassing (source)",
        fontsize=12,
    )

    out = cfg.FIG_DIR / "fig02_annual_flux_map.png"
    fig.savefig(out, dpi=cfg.FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {out}")


# ===========================================================================
# FIG 03 — TREND MAP (Sen's slope per pixel)
# ===========================================================================

def fig03_trend_map(data: dict) -> None:
    """
    Map of the per-pixel linear trend in annual mean air-sea flux.

    Uses scipy.stats.theilslopes (= Sen's slope / Theil-Sen estimator) —
    the robust, median-based trend estimator discussed in the action plan
    (resistant to outliers and non-normal error distributions).

    Units: mol C m⁻² yr⁻¹  per  decade  (trend × 10)

    Red pixels: outgassing is increasing (or uptake is weakening) — BAD
    Blue pixels: uptake is intensifying — GOOD
    """
    from scipy.stats import theilslopes

    print("[compute] Per-pixel Sen's slope (this may take a few minutes) ...")

    # Work with annual means to reduce noise
    flux_annual = data["flux_3d"].resample(time="1YE").mean()
    times       = np.arange(len(flux_annual.time))  # integer year index

    lat = flux_annual.latitude.values
    lon = flux_annual.longitude.values
    flux_np = flux_annual.values  # (years, lat, lon)

    slope_arr = np.full((len(lat), len(lon)), np.nan)

    for i in range(len(lat)):
        for j in range(len(lon)):
            y = flux_np[:, i, j]
            if np.sum(np.isfinite(y)) < 5:   # need at least 5 years
                continue
            ok = np.isfinite(y)
            res = theilslopes(y[ok], times[ok])
            slope_arr[i, j] = res.slope   # units: change per year

    # Convert to change per decade
    slope_da = xr.DataArray(
        slope_arr * 10,
        coords={"latitude": lat, "longitude": lon},
        dims=["latitude", "longitude"],
        attrs={
            "long_name": "Sen's slope of annual air-sea CO2 flux",
            "units":     "mol C m-2 yr-1 per decade",
        },
    )
    slope_da = slope_da.where(data["ocean_mask"] == 1)

    vmax = float(np.nanpercentile(np.abs(slope_da.values), 97))

    if HAS_CARTOPY:
        fig, ax = plt.subplots(
            figsize=(14, 7),
            subplot_kw={"projection": ccrs.Robinson()},
        )
        slope_da.plot(
            ax=ax, transform=ccrs.PlateCarree(),
            cmap=cfg.CMAP_TREND, vmin=-vmax, vmax=vmax,
            add_colorbar=True,
            cbar_kwargs={
                "label": "Trend in CO₂ flux  [mol C m⁻² yr⁻¹ per decade]",
                "shrink": 0.6, "orientation": "horizontal", "pad": 0.04,
            },
        )
        ax.add_feature(cfeature.LAND, facecolor="lightgray", zorder=2)
        ax.add_feature(cfeature.COASTLINE, lw=0.4, zorder=3)
        ax.set_global()
    else:
        fig, ax = plt.subplots(figsize=(14, 6))
        im = ax.imshow(
            slope_da.values[::-1],
            cmap=cfg.CMAP_TREND, vmin=-vmax, vmax=vmax,
            extent=[-180, 180, -90, 90], aspect="auto",
        )
        plt.colorbar(im, ax=ax, label="Trend [mol C m⁻² yr⁻¹ per decade]")

    ax.set_title(
        "Per-pixel trend in air-sea CO₂ flux (Sen's slope, per decade)\n"
        "Positive trend = increasing uptake (sink strengthening)\n"
        "Negative trend = increasing outgassing (sink weakening)",
        fontsize=12,
    )

    out = cfg.FIG_DIR / "fig03_flux_trend_map.png"
    fig.savefig(out, dpi=cfg.FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {out}")


# ===========================================================================
# FIG 04 — ΔpCO2 MAP
# ===========================================================================

def fig04_delta_pco2_map(data: dict) -> None:
    """
    Map of the time-mean ΔpCO2 = pCO2_ocean − pCO2_atm [µatm].

    ΔpCO2 is the thermodynamic driver of the flux:
        ΔpCO2 > 0 → ocean is a potential CO2 source
        ΔpCO2 < 0 → ocean is a potential CO2 sink

    Converted from atm to µatm (× 1e6) for conventional units.
    """
    dpco2 = (data["spco2_ocean"] - data["spco2_atm"]).mean(dim="time")
    dpco2 = dpco2 * 1e6   # atm → µatm
    dpco2 = dpco2.where(data["ocean_mask"] == 1)
    dpco2.attrs["units"] = "µatm"

    vmax = float(np.nanpercentile(np.abs(dpco2.values), 97))

    if HAS_CARTOPY:
        fig, ax = plt.subplots(
            figsize=(14, 7),
            subplot_kw={"projection": ccrs.Robinson()},
        )
        dpco2.plot(
            ax=ax, transform=ccrs.PlateCarree(),
            cmap=cfg.CMAP_FLUX, vmin=-vmax, vmax=vmax,
            add_colorbar=True,
            cbar_kwargs={
                "label": "ΔpCO₂ = pCO₂(ocean) − pCO₂(atm)  [µatm]\n(+ = source driver, − = sink driver)",
                "shrink": 0.6, "orientation": "horizontal", "pad": 0.04,
            },
        )
        ax.add_feature(cfeature.LAND, facecolor="lightgray", zorder=2)
        ax.add_feature(cfeature.COASTLINE, lw=0.4, zorder=3)
        ax.set_global()
    else:
        fig, ax = plt.subplots(figsize=(14, 6))
        im = ax.imshow(
            dpco2.values[::-1],
            cmap=cfg.CMAP_FLUX, vmin=-vmax, vmax=vmax,
            extent=[-180, 180, -90, 90], aspect="auto",
        )
        plt.colorbar(im, ax=ax, label="ΔpCO₂ [µatm]")

    ax.set_title(
        "Time-mean ΔpCO₂ = pCO₂(ocean) − pCO₂(atm)\n"
        "Thermodynamic driver of the air-sea flux",
        fontsize=12,
    )

    out = cfg.FIG_DIR / "fig04_delta_pco2_map.png"
    fig.savefig(out, dpi=cfg.FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {out}")


# ===========================================================================
# FIG 05 — SEASONAL CYCLE
# ===========================================================================

def fig05_seasonal_cycle(data: dict) -> None:
    """
    Climatological seasonal cycle of:
        - Global mean J_net [Pg C yr⁻¹] (monthly climatology)
        - Global mean atmospheric CO2 [ppm]

    Shows that the biological pump and solubility pump have opposite-signed
    seasonal cycles in different hemispheres — the global integral is their
    net residual.
    """
    J_net  = data["J_net"]
    noaa   = data["noaa_co2_ppm"]

    # Monthly climatology of J_net
    J_da = xr.DataArray(
        J_net.values,
        coords={"time": pd.to_datetime(J_net.time.values)},
        dims="time",
    )
    J_clim = J_da.groupby("time.month").mean()

    # Monthly climatology of NOAA CO2
    noaa.index = pd.to_datetime(noaa.index)
    noaa_clim = noaa.groupby(noaa.index.month).mean()

    fig, ax1 = plt.subplots(figsize=(10, 5))
    months = np.arange(1, 13)
    month_labels = ["Jan","Feb","Mar","Apr","May","Jun",
                    "Jul","Aug","Sep","Oct","Nov","Dec"]

    ax1.bar(months, J_clim.values, color="steelblue", alpha=0.7,
            label="J_net seasonal cycle")
    ax1.axhline(0, color="black", lw=0.8)
    ax1.set_xlabel("Month")
    ax1.set_ylabel("Global net uptake  [Pg C yr⁻¹]", color="steelblue")
    ax1.tick_params(axis="y", labelcolor="steelblue")
    ax1.set_xticks(months)
    ax1.set_xticklabels(month_labels)

    ax2 = ax1.twinx()
    ax2.plot(months, noaa_clim.values, color="firebrick", lw=2.2, marker="o",
             label="Atm CO₂ (ppm)")
    ax2.set_ylabel("Atmospheric CO₂  [ppm]", color="firebrick")
    ax2.tick_params(axis="y", labelcolor="firebrick")

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="lower right", fontsize=10)
    ax1.set_title("Climatological seasonal cycle: global ocean uptake vs. atmospheric CO₂",
                  fontsize=12)
    ax1.grid(True, alpha=0.3)
    fig.tight_layout()

    out = cfg.FIG_DIR / "fig05_monthly_seasonal_cycle.png"
    fig.savefig(out, dpi=cfg.FIGURE_DPI)
    plt.close(fig)
    print(f"[fig] {out}")


# ===========================================================================
# MAIN
# ===========================================================================

def main():
    print("[load] Output files ...")
    data = load_all()

    print("\n[figure 1] Global flux vs atmospheric CO2 ...")
    fig01_flux_vs_co2(data)

    print("[figure 2] Mean spatial flux map ...")
    fig02_mean_flux_map(data)

    print("[figure 3] Trend map (Sen's slope) ...")
    fig03_trend_map(data)

    print("[figure 4] ΔpCO2 map ...")
    fig04_delta_pco2_map(data)

    print("[figure 5] Seasonal cycle ...")
    fig05_seasonal_cycle(data)

    print(f"\n[done] All figures saved to {cfg.FIG_DIR}\n")


if __name__ == "__main__":
    main()
