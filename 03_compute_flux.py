"""
03_compute_flux.py
==================
Core physics module for Stage 1. Computes the air-sea CO2 flux reconstruction
using the standard bulk parameterisation:

    F = k · K0 · (pCO2_atm − pCO2_ocean)

pCO2 source : GLOBAL_MULTIYEAR_BGC_001_029  (CMEMS BGC hindcast, pure model)
k            : Wanninkhof (2014) quadratic wind parameterisation
               — with wind variance correction (σ²_u) from ERA5 daily winds
                 when available, monthly-mean only otherwise
K0           : Weiss (1974) solubility, driven by GLORYS12 SST/SSS

Sign convention:
    F > 0  →  ocean UPTAKE   (flux from atmosphere into ocean, sink)
    F < 0  →  ocean OUTGASSING (flux to atmosphere, source)

Outputs:
    data/flux_3d.nc         — (time, lat, lon):
        fgco2               — reconstruction [mol C m⁻² yr⁻¹]
        k, K0, Sc           — gas-exchange intermediates
    output/global_flux.nc   — time series:
        J_net_PgC           — global integral [Pg C yr⁻¹]
        J_net_PgC_annual    — annual resample

References:
    Wanninkhof (2014) DOI:10.4319/lom.2014.12.351
    Weiss (1974)      DOI:10.1016/0304-4203(74)90015-2

Changelog:
    v1.0.0 — Initial implementation, monthly-mean wind only.
    v1.1.0 — Added wind variance correction; added improved MULTIOBS reconstruction.
    v1.2.0 — Dropped MULTIOBS reconstruction; single reconstruction only
              (GLOBAL_MULTIYEAR_BGC_001_029 pCO2 + wind variance correction).
              Both original and variance-corrected k applied to same pCO2 source.

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

    Without variance correction:
        k = a · u² · (Sc/660)^(-0.5)

    With variance correction (preferred when ERA5 daily available):
        k = a · (u² + σ²_u) · (Sc/660)^(-0.5)

    where σ²_u = sub-monthly wind speed variance from ERA5 daily data.
    The correction accounts for the fact that Wanninkhof's coefficient was
    calibrated against the full wind speed distribution, not monthly means:
        ⟨u²⟩ = ⟨u⟩² + σ²_u

    Parameters
    ----------
    wind_speed    : monthly-mean 10 m scalar wind speed [m/s]
    Sc            : Schmidt number of CO2
    wind_variance : sub-monthly variance σ²_u [m²/s²], or None
    """
    u2 = wind_speed**2
    if wind_variance is not None:
        u2 = u2 + wind_variance

    k_cmhr = cfg.WANNINKHOF_A * u2 * (Sc / cfg.SC_REF) ** (-0.5)
    k_ms   = k_cmhr * cfg.CMHR_TO_MS

    note = " + variance correction σ²_u" if wind_variance is not None else ""
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
    ln(K0) = A1 + A2·(100/T) + A3·ln(T/100) + S·[B1 + B2·(T/100) + B3·(T/100)²]
    T in Kelvin, S in PSU. Returns K0 in mol m⁻³ atm⁻¹.
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
) -> xr.DataArray:
    """
    Air-sea CO2 flux: F = k · K0 · (pCO2_atm − pCO2_ocean) [mol m⁻² s⁻¹]
    Converted to mol m⁻² yr⁻¹. Land pixels set to NaN.
    Sign: F > 0 = ocean uptake; F < 0 = outgassing.
    """
    delta_pco2 = spco2_atm - spco2_ocean
    F_per_yr   = k * K0 * delta_pco2 * cfg.S_TO_YR
    F_per_yr   = F_per_yr.where(ocean_mask == 1)
    F_per_yr.attrs = {
        "long_name":       "Air-sea CO2 flux — GLOBAL_MULTIYEAR_BGC_001_029",
        "units":           "mol C m-2 yr-1",
        "sign_convention": "positive = ocean uptake; negative = outgassing",
        "pco2_source":     "GLOBAL_MULTIYEAR_BGC_001_029 (CMEMS BGC hindcast)",
        "k_source":        "Wanninkhof (2014); ERA5 winds",
        "K0_source":       "Weiss (1974); GLORYS12 SST/SSS",
    }
    return F_per_yr.rename("fgco2")


# ===========================================================================
# GRID-AREA HELPER
# ===========================================================================

def compute_grid_cell_area(
    lat: xr.DataArray,
    lon: xr.DataArray,
) -> xr.DataArray:
    """
    Area of each 0.25° grid cell in m²:
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
    Integrate flux over the global ocean surface:
        J(t) = Σ F(t,lat,lon) · A(lat,lon)    [mol C yr⁻¹] → Pg C yr⁻¹
    Positive = ocean uptake.
    """
    masked = flux.where(ocean_mask == 1)
    F_mol  = (masked * cell_area).sum(dim=["latitude", "longitude"], skipna=True)
    J_PgC  = F_mol * cfg.MOL_C_TO_PG
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
    Single reconstruction: GLOBAL_MULTIYEAR_BGC_001_029 pCO2
    + Wanninkhof (2014) k (with wind variance correction if available)
    + Weiss (1974) K0 from GLORYS12 SST/SSS.

    Outputs
    -------
    data/flux_3d.nc       : fgco2, k, K0, Sc
    output/global_flux.nc : J_net_PgC, J_net_PgC_annual
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

    # --- Shared intermediates ---
    print("[compute] Schmidt number ...")
    Sc = schmidt_number_co2(ds["sst"])

    print("[compute] CO2 solubility K0 ...")
    K0 = co2_solubility_K0(ds["sst"], ds["sss"])

    # --- Gas transfer velocity (with variance correction if available) ---
    has_wind_var = "wind_variance" in ds
    wind_var     = ds["wind_variance"] if has_wind_var else None
    print(f"[compute] k — wind variance correction: "
          f"{'YES (ERA5 daily σ²_u)' if has_wind_var else 'NO (monthly mean only)'}")
    k = gas_transfer_velocity(ds["wind_speed"], Sc, wind_variance=wind_var)

    # --- Flux ---
    print("[compute] Air-sea CO2 flux F ...")
    flux = compute_flux(
        k, K0,
        ds["spco2_atm"], ds["spco2_ocean"],
        ds["ocean_mask"],
    )
    if has_wind_var:
        flux.attrs["wind_note"] = "Wind variance correction applied (ERA5 daily σ²_u)"
    else:
        flux.attrs["wind_note"] = "Monthly-mean wind only (ERA5 daily not available)"

    # --- Save 3D flux ---
    ds_flux = xr.Dataset(
        {"fgco2": flux, "k": k, "K0": K0, "Sc": Sc},
        attrs={
            "title":   "Stage 1 — Air-sea CO2 flux (GLOBAL_MULTIYEAR_BGC_001_029)",
            "version": "1.2.0",
            "project": "Ocean Carbon Cycle — Net Flux Project",
        },
    )
    print(f"[save] {flux_3d_file} ...")
    ds_flux.to_netcdf(flux_3d_file)

    # --- Grid cell areas ---
    cell_area = compute_grid_cell_area(ds["latitude"], ds["longitude"])

    # --- Global integral ---
    print("[compute] Global surface flux integral ...")
    J = global_integral(flux, cell_area, ds["ocean_mask"], "J_net_PgC")

    # Annual resample
    J_annual = J.resample(time="1YE").mean()

    ds_global = xr.Dataset(
        {
            "J_net_PgC":        J,
            "J_net_PgC_annual": J_annual.reindex(time=J.time, method="nearest"),
        },
        attrs={
            "title":   "Stage 1 — Global surface CO2 flux time series",
            "version": "1.2.0",
            "project": "Ocean Carbon Cycle — Net Flux Project",
        },
    )
    print(f"[save] {global_out_file} ...")
    ds_global.to_netcdf(global_out_file)

    print("\n[done] Flux computation complete.")
    print(f"  J_net range: {float(J.min()):.3f} – {float(J.max()):.3f} Pg C yr⁻¹")


if __name__ == "__main__":
    main()
