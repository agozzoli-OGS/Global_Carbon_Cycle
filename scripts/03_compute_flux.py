"""
03_compute_flux.py
==================
Core physics module for Stage 1. Computes the air-sea CO2 flux F at every
grid point and time step using the standard bulk parameterization:

    F = k · K0 · (pCO2_atm − pCO2_ocean)

where:
    k    = gas transfer velocity [m/s]       (Wanninkhof 2014)
    K0   = CO2 solubility       [mol/m³/atm] (Weiss 1974)
    pCO2 = partial pressure     [atm]

Sign convention (OCMIP / CMEMS standard):
    F > 0  →  ocean OUTGASSING (flux to atmosphere, source)
    F < 0  →  ocean UPTAKE     (flux from atmosphere, sink)

All functions are pure — they take xr.DataArray inputs and return
xr.DataArray outputs. No I/O happens here; that is done by main() at the
bottom.

References:
    Wanninkhof, R. (2014). Relationship between wind speed and gas exchange
        over the ocean revisited. Limnology and Oceanography: Methods, 12(6),
        351–362. https://doi.org/10.4319/lom.2014.12.351
    Weiss, R.F. (1974). Carbon dioxide in water and seawater: the solubility
        of a non-ideal gas. Marine Chemistry, 2(3), 203–215.
        https://doi.org/10.1016/0304-4203(74)90015-2

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
# PHYSICAL FUNCTIONS
# ===========================================================================

def schmidt_number_co2(sst_degC: xr.DataArray) -> xr.DataArray:
    """
    Compute the Schmidt number of CO2 in seawater as a function of SST.

    The Schmidt number Sc is the ratio of the kinematic viscosity of seawater
    to the diffusion coefficient of CO2. It is used to normalise the gas
    transfer velocity k to a reference Schmidt number of 660 (CO2 at 20°C).

    Polynomial fit from Wanninkhof (2014), Table 1:
        Sc = A - B·T + C·T² - D·T³ + E·T⁴
    where T is in °C.

    Parameters
    ----------
    sst_degC : xr.DataArray
        Sea surface temperature in degrees Celsius.

    Returns
    -------
    xr.DataArray
        Schmidt number (dimensionless), same shape as input.
    """
    T = sst_degC
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
) -> xr.DataArray:
    """
    Compute the gas transfer velocity k using the Wanninkhof (2014)
    quadratic wind-speed parameterisation.

    Equation:
        k = a · u² · (Sc / 660)^(-0.5)

    where:
        a    = 0.251 cm hr⁻¹ (m s⁻¹)⁻²   [Wanninkhof 2014 coefficient]
        u    = 10 m wind speed [m/s]
        Sc   = Schmidt number of CO2       [dimensionless]
        660  = Schmidt number of CO2 at 20°C in seawater [reference]

    The factor (Sc/660)^(-0.5) normalises k to the reference condition.

    k is first computed in cm/hr, then converted to m/s.

    Parameters
    ----------
    wind_speed : xr.DataArray
        10 m scalar wind speed [m/s].
    Sc : xr.DataArray
        Schmidt number of CO2 (from schmidt_number_co2).

    Returns
    -------
    xr.DataArray
        Gas transfer velocity k [m/s].
    """
    # k in cm/hr
    k_cmhr = cfg.WANNINKHOF_A * (wind_speed**2) * (Sc / cfg.SC_REF) ** (-0.5)

    # Convert to m/s
    k_ms = k_cmhr * cfg.CMHR_TO_MS

    k_ms.attrs = {
        "long_name": "Gas transfer velocity for CO2",
        "units":     "m s-1",
        "reference": "Wanninkhof (2014), Eq. 3; a=0.251 cm/hr/(m/s)²",
    }
    return k_ms.rename("k")


def co2_solubility_K0(
    sst_degC: xr.DataArray,
    sss: xr.DataArray,
) -> xr.DataArray:
    """
    Compute the solubility of CO2 in seawater (K0) using Weiss (1974).

    K0 is defined as the ratio of dissolved CO2 concentration to its
    partial pressure:
        [CO2(aq)] = K0 · pCO2

    Equation (Weiss 1974, Eq. 12 in natural log form):
        ln(K0) = A1 + A2·(100/T) + A3·ln(T/100)
                 + S·[B1 + B2·(T/100) + B3·(T/100)²]

    where T is in Kelvin and S is salinity in PSU.

    K0 from Weiss (1974) is in mol / (L · atm).
    We convert to mol / (m³ · atm) by multiplying by 1000.

    Parameters
    ----------
    sst_degC : xr.DataArray
        Sea surface temperature [°C].
    sss : xr.DataArray
        Sea surface salinity [PSU].

    Returns
    -------
    xr.DataArray
        K0, CO2 solubility [mol m⁻³ atm⁻¹].
    """
    # Temperature in Kelvin
    T = sst_degC + 273.15

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

    # K0 in mol/(L·atm) → mol/(m³·atm)  [×1000]
    K0 = np.exp(ln_K0) * 1000.0

    K0.attrs = {
        "long_name": "CO2 solubility K0",
        "units":     "mol m-3 atm-1",
        "reference": "Weiss (1974) Eq. 12; converted from mol/(L·atm)×1000",
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
    Compute the air-sea CO2 flux F at every grid point and time step.

    Equation:
        F = k · K0 · (pCO2_atm − pCO2_ocean)    [mol m⁻² s⁻¹]

    Sign convention:
        F > 0  →  pCO2_ocean > pCO2_atm  →  outgassing (source)
        F < 0  →  pCO2_ocean < pCO2_atm  →  ocean uptake (sink)

    The flux is then converted to mol m⁻² yr⁻¹ for reporting.
    Land pixels (ocean_mask == 0) are set to NaN.

    Parameters
    ----------
    k : xr.DataArray
        Gas transfer velocity [m/s].
    K0 : xr.DataArray
        CO2 solubility [mol m⁻³ atm⁻¹].
    spco2_atm : xr.DataArray
        Atmospheric pCO2 [atm], broadcast to (time, lat, lon).
    spco2_ocean : xr.DataArray
        Surface ocean pCO2 [atm].
    ocean_mask : xr.DataArray
        Ocean mask (1=ocean, 0=land).

    Returns
    -------
    xr.DataArray
        Air-sea CO2 flux F [mol C m⁻² yr⁻¹].
        Positive = outgassing, negative = uptake.
    """
    # Core bulk formula in mol m⁻² s⁻¹
    delta_pco2 = spco2_ocean - spco2_atm   # positive = outgassing
    F_per_s = k * K0 * delta_pco2

    # Convert s⁻¹ → yr⁻¹
    F_per_yr = F_per_s * cfg.S_TO_YR

    # Mask land
    F_per_yr = F_per_yr.where(ocean_mask == 1)

    F_per_yr.attrs = {
        "long_name":   "Air-sea CO2 flux",
        "units":       "mol C m-2 yr-1",
        "sign_convention": (
            "positive = ocean outgassing (source); "
            "negative = ocean uptake (sink)"
        ),
        "references": (
            "Wanninkhof (2014) doi:10.4319/lom.2014.12.351; "
            "Weiss (1974) doi:10.1016/0304-4203(74)90015-2"
        ),
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
    Compute the area of each 0.25° grid cell in m².

    For a regular lat/lon grid with spacing Δlat × Δlon (both in degrees):
        A(lat) = R² · |Δlat_rad| · |Δlon_rad| · cos(lat_rad)

    where R = 6.371×10⁶ m (Earth mean radius).

    Parameters
    ----------
    lat : xr.DataArray
        1-D latitude coordinate [degrees].
    lon : xr.DataArray
        1-D longitude coordinate [degrees].

    Returns
    -------
    xr.DataArray
        Grid cell area [m²], shape (lat, lon).
    """
    R = cfg.EARTH_RADIUS_M

    # Grid spacing (assume uniform)
    d_lat = float(np.abs(lat.diff("latitude").mean()))
    d_lon = float(np.abs(lon.diff("longitude").mean()))

    d_lat_rad = np.deg2rad(d_lat)
    d_lon_rad = np.deg2rad(d_lon)
    lat_rad   = np.deg2rad(lat)

    # Area in m²: R² * dlon * dlat * cos(lat)
    area_1d = (R**2) * d_lon_rad * d_lat_rad * np.cos(lat_rad)

    # Broadcast to (lat, lon) — same value for all longitudes at a given lat
    area_2d = area_1d.expand_dims({"longitude": lon})
    area_2d = area_2d.transpose("latitude", "longitude")

    area_2d.attrs = {
        "long_name": "Grid cell area",
        "units":     "m2",
    }
    return area_2d.rename("cell_area")


# ===========================================================================
# GLOBAL INTEGRATION
# ===========================================================================

def global_integral(
    flux: xr.DataArray,
    cell_area: xr.DataArray,
    ocean_mask: xr.DataArray,
) -> xr.DataArray:
    """
    Integrate the air-sea CO2 flux over the global ocean surface.

    Computes, at each time step:
        J(t) = Σ_{lat,lon} F(t, lat, lon) · A(lat, lon)    [mol C yr⁻¹]

    Then converts to Pg C yr⁻¹ for human-readable reporting.

    The 'J_in − J_out' framing in the meeting corresponds to:
        J_net = -J(t)  [sign flip: positive J = uptake INTO the ocean]

    Returns both the raw integral (mol C yr⁻¹) and the sign-flipped
    'J_net' in Pg C yr⁻¹ (positive = ocean uptake).

    Parameters
    ----------
    flux : xr.DataArray
        Flux field [mol C m⁻² yr⁻¹], shape (time, lat, lon).
    cell_area : xr.DataArray
        Grid cell area [m²], shape (lat, lon).
    ocean_mask : xr.DataArray
        Ocean mask (1=ocean, 0=land), shape (lat, lon).

    Returns
    -------
    xr.Dataset
        Dataset with variables:
            F_global_mol     [mol C yr⁻¹]  — raw integral (+ = outgassing)
            J_net_PgC        [Pg C yr⁻¹]   — ocean uptake (+ = uptake)
    """
    # Weighted area sum over lat/lon (NaN pixels don't contribute)
    masked_flux = flux.where(ocean_mask == 1)
    F_global_mol = (masked_flux * cell_area).sum(
        dim=["latitude", "longitude"], skipna=True
    )

    # Convert mol C → Pg C
    F_global_PgC = F_global_mol * cfg.MOL_C_TO_PG

    # J_net (positive = ocean uptake): flip sign
    J_net_PgC = -F_global_PgC

    F_global_mol.attrs = {
        "long_name": "Global ocean–atmosphere CO2 flux integral",
        "units": "mol C yr-1",
        "sign_convention": "positive = outgassing",
    }
    J_net_PgC.attrs = {
        "long_name": "Net global ocean CO2 uptake (J_net)",
        "units": "Pg C yr-1",
        "sign_convention": "positive = ocean uptake (sink)",
    }

    return xr.Dataset(
        {"F_global_mol": F_global_mol, "J_net_PgC": J_net_PgC}
    )


# ===========================================================================
# MAIN
# ===========================================================================

def main():
    """
    Load preprocessed data, compute flux, integrate globally, save output.

    Output files:
        data/flux_3d.nc        — full (time, lat, lon) flux field
        output/global_flux.nc  — time series of global integral + J_net
    """
    flux_3d_file   = cfg.DATA_DIR / "flux_3d.nc"
    global_out_file = cfg.OUT_DIR / "global_flux.nc"

    # Load preprocessed surface data
    print("[load] preprocessed_surface.nc ...")
    ds = xr.open_dataset(cfg.DATA_DIR / "processed_surface.nc", chunks="auto")

    # Check that wind speed is available
    if "wind_speed" not in ds:
        raise FileNotFoundError(
            "wind_speed not found in processed_surface.nc.\n"
            "Download ERA5 wind data first — see 01_download_data.py."
        )

    # --- Compute physical quantities ---
    print("[compute] Schmidt number ...")
    Sc = schmidt_number_co2(ds["sst"])

    print("[compute] Gas transfer velocity k ...")
    k = gas_transfer_velocity(ds["wind_speed"], Sc)

    print("[compute] CO2 solubility K0 ...")
    K0 = co2_solubility_K0(ds["sst"], ds["sss"])

    print("[compute] Air-sea CO2 flux F ...")
    flux = compute_flux(
        k            = k,
        K0           = K0,
        spco2_atm    = ds["spco2_atm"],
        spco2_ocean  = ds["spco2_ocean"],
        ocean_mask   = ds["ocean_mask"],
    )

    # --- Save full 3D flux field ---
    ds_flux = xr.Dataset(
        {"fgco2": flux, "k": k, "K0": K0, "Sc": Sc},
        attrs={
            "title":   "Stage 1 — Air-sea CO2 flux (full 3D field)",
            "project": "Ocean Carbon Cycle — Net Flux Project",
        },
    )
    print(f"[save] {flux_3d_file} ...")
    ds_flux.to_netcdf(flux_3d_file)

    # --- Compute grid cell areas ---
    print("[compute] Grid cell areas ...")
    cell_area = compute_grid_cell_area(ds["latitude"], ds["longitude"])

    # --- Global integral ---
    print("[compute] Global surface flux integral ...")
    ds_global = global_integral(flux, cell_area, ds["ocean_mask"])

    # Add annual means alongside the monthly time series
    ds_global_annual = ds_global.resample(time="1YE").mean()

    ds_global["F_global_mol_annual"] = ds_global_annual["F_global_mol"].reindex(
        time=ds_global.time, method="nearest"
    )
    ds_global["J_net_PgC_annual"] = ds_global_annual["J_net_PgC"].reindex(
        time=ds_global.time, method="nearest"
    )

    print(f"[save] {global_out_file} ...")
    ds_global.to_netcdf(global_out_file)

    print("\n[done] Flux computation complete.")
    print(f"       Monthly global J_net range: "
          f"{float(ds_global['J_net_PgC'].min()):.3f} to "
          f"{float(ds_global['J_net_PgC'].max()):.3f} Pg C yr⁻¹")


if __name__ == "__main__":
    main()
