"""
03_compute_flux.py
==================
Core physics module for Stage 1. Computes TWO air-sea CO2 flux reconstructions:

    ORIGINAL  (v1.0.x):
        F = k(u_monthly) · K0(SST,SSS) · (pCO2_atm − pCO2_PISCES)

    IMPROVED  (v1.1.0) — two changes stacked:
        (1) pCO2 driver replaced: PISCES → MULTIOBS SOCAT-NN surface pCO2
        (2) Wind variance correction applied:
            k_corr = a · (<u>² + σ²_u) · (Sc/660)^(-0.5)
            where σ²_u = sub-monthly wind speed variance from ERA5 daily winds

Sign convention (both reconstructions):
    F > 0  →  ocean UPTAKE   (flux from atmosphere into ocean, sink)
    F < 0  →  ocean OUTGASSING (flux to atmosphere, source)

All physical functions are pure (no I/O). main() does all loading and saving.

Outputs:
    data/flux_3d.nc         — (time, lat, lon) fields:
        fgco2               — original reconstruction   [mol C m⁻² yr⁻¹]
        fgco2_improved      — improved reconstruction   [mol C m⁻² yr⁻¹]
        k, K0, Sc           — gas exchange intermediates (from original)
    output/global_flux.nc   — time series:
        J_net_PgC           — original global integral  [Pg C yr⁻¹]
        J_net_improved_PgC  — improved global integral  [Pg C yr⁻¹]

References:
    Wanninkhof (2014) DOI:10.4319/lom.2014.12.351
    Weiss (1974)      DOI:10.1016/0304-4203(74)90015-2

Changelog:
    v1.1.0 — Added compute_flux_improved() using MULTIOBS pCO2 + variance
              correction; both reconstructions saved to the same output files.

Usage:
    python scripts/03_compute_flux.py
"""

import sys
from pathlib import Path

import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).parent.parent))
import config as cfg


# ===========================================================================
# PHYSICAL FUNCTIONS  (pure — no I/O)
# ===========================================================================

def schmidt_number_co2(sst_degC: xr.DataArray) -> xr.DataArray:
    """
    Schmidt number of CO2 in seawater (Wanninkhof 2014, Table 1).

    Sc = A - B·T + C·T² - D·T³ + E·T⁴      (T in °C)
    """
    T  = sst_degC
    Sc = (
        cfg.SC_A
        - cfg.SC_B * T
        + cfg.SC_C * T**2
        - cfg.SC_D * T**3
        + cfg.SC_E * T**4
    )
    Sc.attrs = {
        "long_name": "Schmidt number of CO2 in seawater",
        "units":     "dimensionless",
        "reference": "Wanninkhof (2014) Table 1",
    }
    return Sc.rename("schmidt_number")


def gas_transfer_velocity(
    wind_speed: xr.DataArray,
    Sc: xr.DataArray,
    wind_variance: xr.DataArray | None = None,
) -> xr.DataArray:
    """
    Gas transfer velocity k (Wanninkhof 2014).

    Without variance correction (original, v1.0.x):
        k = a · u² · (Sc/660)^(-0.5)

    With variance correction (improved, v1.1.0):
        k = a · (u² + σ²_u) · (Sc/660)^(-0.5)

    where σ²_u = sub-monthly wind speed variance from ERA5 daily data.
    This accounts for the fact that Wanninkhof's coefficient was calibrated
    against the full wind speed distribution, not monthly means. Using
    monthly-mean u alone underestimates k, particularly in the Southern
    Ocean and storm tracks where sub-monthly wind bursts dominate.

    Parameters
    ----------
    wind_speed : xr.DataArray
        Monthly-mean 10 m scalar wind speed [m/s].
    Sc : xr.DataArray
        Schmidt number of CO2 (from schmidt_number_co2).
    wind_variance : xr.DataArray, optional
        Sub-monthly variance of wind speed σ²_u [m²/s²].
        If None, the uncorrected formula is used.

    Returns
    -------
    xr.DataArray
        Gas transfer velocity k [m/s].
    """
    u2 = wind_speed**2
    if wind_variance is not None:
        u2 = u2 + wind_variance    # <u²> = <u>² + σ²_u

    k_cmhr = cfg.WANNINKHOF_A * u2 * (Sc / cfg.SC_REF) ** (-0.5)
    k_ms   = k_cmhr * cfg.CMHR_TO_MS

    note = "" if wind_variance is None else " + variance correction (σ²_u)"
    k_ms.attrs = {
        "long_name": "Gas transfer velocity for CO2",
        "units":     "m s-1",
        "reference": f"Wanninkhof (2014){note}",
    }
    return k_ms.rename("k")


def co2_solubility_K0(
    sst_degC: xr.DataArray,
    sss: xr.DataArray,
) -> xr.DataArray:
    """
    CO2 solubility K0 (Weiss 1974).

    ln(K0) = A1 + A2·(100/T) + A3·ln(T/100)
             + S·[B1 + B2·(T/100) + B3·(T/100)²]    (T in K, S in PSU)

    Returns K0 in mol m⁻³ atm⁻¹ (converted from mol L⁻¹ atm⁻¹ × 1000).
    """
    T     = sst_degC + 273.15
    ln_K0 = (
        cfg.WEISS_A1
        + cfg.WEISS_A2 * (100.0 / T)
        + cfg.WEISS_A3 * np.log(T / 100.0)
        + sss * (
            cfg.WEISS_B1
            + cfg.WEISS_B2 * (T / 100.0)
            + cfg.WEISS_B3 * (T / 100.0) ** 2
        )
    )
    K0 = np.exp(ln_K0) * 1000.0   # mol/(L·atm) → mol/(m³·atm)
    K0.attrs = {
        "long_name": "CO2 solubility K0",
        "units":     "mol m-3 atm-1",
        "reference": "Weiss (1974) Eq. 12",
    }
    return K0.rename("K0")


def compute_flux(
    k: xr.DataArray,
    K0: xr.DataArray,
    spco2_atm: xr.DataArray,
    spco2_ocean: xr.DataArray,
    ocean_mask: xr.DataArray,
    label: str = "fgco2",
) -> xr.DataArray:
    """
    Air-sea CO2 flux via the bulk formula:

        F = k · K0 · (pCO2_atm − pCO2_ocean)    [mol m⁻² s⁻¹]

    Converted to mol m⁻² yr⁻¹. Land pixels set to NaN.

    Sign convention: F > 0 = ocean uptake; F < 0 = outgassing.

    Parameters
    ----------
    label : str
        Name for the output DataArray ('fgco2' or 'fgco2_improved').
    """
    delta_pco2 = spco2_atm - spco2_ocean   # positive = uptake
    F_per_yr   = k * K0 * delta_pco2 * cfg.S_TO_YR
    F_per_yr   = F_per_yr.where(ocean_mask == 1)
    F_per_yr.attrs = {
        "long_name":        f"Air-sea CO2 flux ({label})",
        "units":            "mol C m-2 yr-1",
        "sign_convention":  "positive = ocean uptake; negative = outgassing",
        "references":       "Wanninkhof (2014); Weiss (1974)",
    }
    return F_per_yr.rename(label)


# ===========================================================================
# GRID-AREA HELPER
# ===========================================================================

def compute_grid_cell_area(
    lat: xr.DataArray,
    lon: xr.DataArray,
) -> xr.DataArray:
    """
    Area of each 0.25° grid cell in m² via spherical geometry:
        A(lat) = R² · Δlon_rad · Δlat_rad · cos(lat_rad)
    """
    R         = cfg.EARTH_RADIUS_M
    d_lat     = float(np.abs(lat.diff("latitude").mean()))
    d_lon     = float(np.abs(lon.diff("longitude").mean()))
    d_lat_rad = np.deg2rad(d_lat)
    d_lon_rad = np.deg2rad(d_lon)
    lat_rad   = np.deg2rad(lat)

    area_1d = (R**2) * d_lon_rad * d_lat_rad * np.cos(lat_rad)
    area_2d = area_1d.expand_dims({"longitude": lon}).transpose("latitude", "longitude")
    area_2d.attrs = {"long_name": "Grid cell area", "units": "m2"}
    return area_2d.rename("cell_area")


# ===========================================================================
# GLOBAL INTEGRATION
# ===========================================================================

def global_integral(
    flux: xr.DataArray,
    cell_area: xr.DataArray,
    ocean_mask: xr.DataArray,
    out_name: str = "J_net_PgC",
) -> xr.DataArray:
    """
    Integrate flux over the global ocean surface at each time step.

        J(t) = Σ F(t,lat,lon) · A(lat,lon)

    Returns a 1-D time series in Pg C yr⁻¹ (positive = uptake).

    Parameters
    ----------
    out_name : str
        Name for the output DataArray.
    """
    masked   = flux.where(ocean_mask == 1)
    F_mol    = (masked * cell_area).sum(dim=["latitude", "longitude"], skipna=True)
    J_PgC    = F_mol * cfg.MOL_C_TO_PG
    J_PgC.attrs = {
        "long_name":       f"Global net ocean CO2 uptake ({out_name})",
        "units":           "Pg C yr-1",
        "sign_convention": "positive = ocean uptake (sink)",
    }
    return J_PgC.rename(out_name)


# ===========================================================================
# MAIN
# ===========================================================================

def main():
    """
    Compute both the original and improved flux reconstructions and save.

    Outputs
    -------
    data/flux_3d.nc
        fgco2              — original reconstruction
        fgco2_improved     — improved reconstruction (if inputs available)
        k, K0, Sc          — gas-exchange intermediates
    output/global_flux.nc
        J_net_PgC          — original global integral
        J_net_improved_PgC — improved global integral (if available)
        *_annual variants  — annual resamples of each
    """
    flux_3d_file    = cfg.DATA_DIR / "flux_3d.nc"
    global_out_file = cfg.OUT_DIR  / "global_flux.nc"

    print("[load] processed_surface.nc ...")
    ds = xr.open_dataset(cfg.DATA_DIR / "processed_surface.nc", chunks="auto")

    if "wind_speed" not in ds:
        raise FileNotFoundError(
            "wind_speed not in processed_surface.nc — "
            "download ERA5 monthly winds first (01_download_data.py)."
        )

    # -----------------------------------------------------------------------
    # Shared intermediates (same for both reconstructions)
    # -----------------------------------------------------------------------
    print("[compute] Schmidt number ...")
    Sc = schmidt_number_co2(ds["sst"])

    print("[compute] CO2 solubility K0 ...")
    K0 = co2_solubility_K0(ds["sst"], ds["sss"])

    # -----------------------------------------------------------------------
    # ORIGINAL reconstruction (v1.0.x)
    #   pCO2 driver: PISCES BGC hindcast
    #   k: no variance correction (monthly-mean wind only)
    # -----------------------------------------------------------------------
    print("[compute] k — original (no variance correction) ...")
    k_orig = gas_transfer_velocity(ds["wind_speed"], Sc, wind_variance=None)

    print("[compute] Flux — original ...")
    flux_orig = compute_flux(
        k_orig, K0,
        ds["spco2_atm"], ds["spco2_ocean"],
        ds["ocean_mask"],
        label="fgco2",
    )

    # -----------------------------------------------------------------------
    # IMPROVED reconstruction (v1.1.0)
    #   pCO2 driver: MULTIOBS SOCAT-NN  (replaces PISCES)
    #   k: wind variance correction applied (if ERA5 daily available)
    # -----------------------------------------------------------------------
    has_obs_pco2  = "spco2_ocean_obs" in ds
    has_wind_var  = "wind_variance"   in ds
    has_improved  = has_obs_pco2   # minimum requirement for improved flux

    flux_improved = None
    k_improved    = None

    if has_improved:
        print("[compute] k — improved (variance correction: "
              + ("YES" if has_wind_var else "NO — daily wind missing") + ") ...")
        wind_var  = ds["wind_variance"] if has_wind_var else None
        k_improved = gas_transfer_velocity(ds["wind_speed"], Sc, wind_variance=wind_var)

        print("[compute] Flux — improved ...")
        flux_improved = compute_flux(
            k_improved, K0,
            ds["spco2_atm"], ds["spco2_ocean_obs"],
            ds["ocean_mask"],
            label="fgco2_improved",
        )
        if not has_wind_var:
            flux_improved.attrs["note"] = (
                "pCO2 driver: MULTIOBS SOCAT-NN; "
                "wind variance correction NOT applied (ERA5 daily wind not available)"
            )
        else:
            flux_improved.attrs["note"] = (
                "pCO2 driver: MULTIOBS SOCAT-NN; "
                "wind variance correction applied (ERA5 daily σ²_u)"
            )
    else:
        print("[warn] spco2_ocean_obs not found — improved reconstruction skipped.")
        print("       Download MULTIOBS and re-run 02_preprocess.py.")

    # -----------------------------------------------------------------------
    # Save 3D flux fields
    # -----------------------------------------------------------------------
    flux_vars = {"fgco2": flux_orig, "k": k_orig, "K0": K0, "Sc": Sc}
    if flux_improved is not None:
        flux_vars["fgco2_improved"] = flux_improved
        if k_improved is not None:
            flux_vars["k_improved"] = k_improved

    ds_flux = xr.Dataset(flux_vars, attrs={
        "title":   "Stage 1 — Air-sea CO2 flux (3D fields)",
        "version": "1.1.0",
        "project": "Ocean Carbon Cycle — Net Flux Project",
    })
    print(f"[save] {flux_3d_file} ...")
    ds_flux.to_netcdf(flux_3d_file)

    # -----------------------------------------------------------------------
    # Global integrals
    # -----------------------------------------------------------------------
    cell_area = compute_grid_cell_area(ds["latitude"], ds["longitude"])

    print("[compute] Global integral — original ...")
    J_orig = global_integral(flux_orig, cell_area, ds["ocean_mask"], "J_net_PgC")

    global_vars = {"J_net_PgC": J_orig}

    if flux_improved is not None:
        print("[compute] Global integral — improved ...")
        J_improved = global_integral(
            flux_improved, cell_area, ds["ocean_mask"], "J_net_improved_PgC"
        )
        global_vars["J_net_improved_PgC"] = J_improved

    # Add annual resamples
    ds_global = xr.Dataset(global_vars)
    for var in list(ds_global.data_vars):
        annual      = ds_global[var].resample(time="1YE").mean()
        annual_name = var + "_annual"
        ds_global[annual_name] = annual.reindex(
            time=ds_global.time, method="nearest"
        )

    ds_global.attrs = {
        "title":   "Stage 1 — Global surface CO2 flux time series",
        "version": "1.1.0",
        "project": "Ocean Carbon Cycle — Net Flux Project",
    }
    print(f"[save] {global_out_file} ...")
    ds_global.to_netcdf(global_out_file)

    print("\n[done] Flux computation complete.")
    print(f"  Original  J_net range: "
          f"{float(J_orig.min()):.3f} – {float(J_orig.max()):.3f} Pg C yr⁻¹")
    if flux_improved is not None:
        print(f"  Improved  J_net range: "
              f"{float(J_improved.min()):.3f} – {float(J_improved.max()):.3f} Pg C yr⁻¹")


if __name__ == "__main__":
    main()
