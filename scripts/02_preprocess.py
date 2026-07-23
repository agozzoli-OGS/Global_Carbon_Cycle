"""
02_preprocess.py
================
Load, regrid, unit-convert, and harmonize all downloaded datasets onto a
common 0.25° monthly grid ready for flux computation.

Steps:
    1. Load BGC hindcast surface pCO2 (already 0.25°).
    2. Load physical reanalysis SST and SSS (0.083°) and regrid to 0.25°.
    3. Load ERA5 wind (0.25° or 0.5°), compute scalar speed, regrid to 0.25°.
    4. Load NOAA atmospheric CO2 CSV, parse to a DataArray on the same
       time axis as the CMEMS products.
    5. Load MULTIOBS validation flux (already 0.25°).
    6. Harmonize time axes (cftime → np.datetime64), align all on a single
       monthly time coordinate.
    7. Apply land mask (ocean-only pixels from BGC hindcast's spco2 fill values).
    8. Apply unit conversions (Pa → atm, °C → K, etc.).
    9. Save one merged NetCDF per variable family to data/.

Outputs:
    data/processed_surface.nc   — merged xarray Dataset with:
        spco2_ocean [atm]       — surface ocean pCO2
        spco2_atm   [atm]       — atmospheric pCO2 (broadcast to lat/lon)
        sst         [°C]        — sea surface temperature
        sss         [PSU]       — sea surface salinity
        wind_speed  [m/s]       — 10 m scalar wind speed
        fgco2_obs   [mol/m²/yr] — validation flux from MULTIOBS
        ocean_mask              — 1 = ocean, 0 = land/ice

Usage:
    python scripts/02_preprocess.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, str(Path(__file__).parent.parent))
import config as cfg


# ===========================================================================
# LOADERS
# ===========================================================================

def load_spco2_ocean() -> xr.DataArray:
    """
    Load surface ocean pCO2 from the BGC hindcast and convert Pa → atm.

    The CMEMS variable 'spco2' is in Pa.
    1 atm = 101325 Pa.
    """
    print("[load] BGC hindcast pCO2 ...")
    ds = xr.open_dataset(cfg.DATA_DIR / "bgc_hindcast_spco2.nc", chunks="auto")

    da = ds[cfg.BGC_PCO2_VAR]

    # Squeeze out the depth dimension (we downloaded only the surface level)
    if "depth" in da.dims:
        da = da.isel(depth=0, drop=True)

    # Unit conversion: Pa → atm
    da = da * cfg.PA_TO_ATM
    da.attrs["units"] = "atm"
    da.attrs["long_name"] = "Surface ocean pCO2"
    return da.rename("spco2_ocean")


def load_sst_sss() -> xr.Dataset:
    """
    Load SST and SSS from GLORYS12 physical reanalysis (0.083°) and
    regrid to the BGC hindcast's 0.25° grid using bilinear interpolation.

    SST is in °C (thetao), SSS is in PSU (so).
    Both are squeezed to surface (first depth level).
    """
    print("[load] Physical reanalysis SST/SSS ...")
    ds_phy = xr.open_dataset(cfg.DATA_DIR / "phy_reanalysis_sst_sss.nc", chunks="auto")

    sst = ds_phy[cfg.PHY_TEMP_VAR]
    sss = ds_phy[cfg.PHY_SAL_VAR]

    if "depth" in sst.dims:
        sst = sst.isel(depth=0, drop=True)
        sss = sss.isel(depth=0, drop=True)

    # Load the BGC target grid from the already-downloaded pCO2 file
    ds_bgc = xr.open_dataset(cfg.DATA_DIR / "bgc_hindcast_spco2.nc", chunks="auto")
    target_lat = ds_bgc["latitude"]
    target_lon = ds_bgc["longitude"]

    # Regrid: interp to the coarser 0.25° BGC grid
    # (bilinear interpolation is appropriate for physical fields on a regular grid)
    print("[regrid] SST/SSS 0.083° → 0.25° ...")
    sst = sst.interp(latitude=target_lat, longitude=target_lon, method="linear")
    sss = sss.interp(latitude=target_lat, longitude=target_lon, method="linear")

    sst.attrs["units"] = "degC"
    sss.attrs["units"] = "PSU"

    out = xr.Dataset({"sst": sst.rename("sst"), "sss": sss.rename("sss")})
    return out


def load_wind_speed() -> xr.DataArray:
    """
    Load ERA5 10 m wind components and compute scalar wind speed.

    Returns the magnitude |u| = sqrt(u10² + v10²) in m/s, regridded to
    the 0.25° BGC grid if the ERA5 native resolution differs.

    If the ERA5 file does not exist yet (CDS download not done), returns
    a placeholder DataArray of NaN — the rest of the pipeline will raise a
    clear error in 03_compute_flux.py.
    """
    wind_file = cfg.ERA5_WIND_FILE
    if not wind_file.exists():
        print("[warn] ERA5 wind file not found — returning NaN placeholder.")
        print("       Run download step with CDS API first (see 01_download_data.py).")
        return None

    print("[load] ERA5 10 m wind ...")
    ds_era5 = xr.open_dataset(wind_file, chunks="auto")

    u10 = ds_era5[cfg.ERA5_U10_VAR]
    v10 = ds_era5[cfg.ERA5_V10_VAR]

    # Scalar wind speed [m/s]
    wspd = np.sqrt(u10**2 + v10**2)
    wspd.attrs["units"] = "m/s"
    wspd.attrs["long_name"] = "10 m scalar wind speed"

    # Regrid to 0.25° BGC grid if needed
    ds_bgc = xr.open_dataset(cfg.DATA_DIR / "bgc_hindcast_spco2.nc", chunks="auto")
    wspd = wspd.interp(
        latitude=ds_bgc["latitude"],
        longitude=ds_bgc["longitude"],
        method="linear",
    )

    return wspd.rename("wind_speed")


def load_atmospheric_co2() -> xr.DataArray:
    """
    Load NOAA GML global average monthly CO2 (ppm) and convert to atm.

    The CSV file has comment lines starting with '#'. The data columns are:
        year, month, decimal, average, average_unc, trend, trend_unc

    We parse year + month → a monthly datetime index, then convert:
        ppm → atm  (1 ppm = 1e-6 atm)

    Returns a 1-D DataArray indexed by 'time' (numpy datetime64[ns]).
    """
    print("[load] NOAA GML atmospheric CO2 ...")
    df = pd.read_csv(
        cfg.NOAA_CO2_FILE,
        comment="#",
        names=["year", "month", "decimal", "average", "average_unc", "trend", "trend_unc"],
    )

    # Drop rows where average is flagged as missing (-9.99 or -99.99)
    df = df[df["average"] > 0]

    # Build monthly datetime index
    df["time"] = pd.to_datetime(
        df["year"].astype(str) + "-" + df["month"].astype(str).str.zfill(2)
    )
    df = df.set_index("time").sort_index()

    # Convert ppm → atm
    co2_atm = df["average"] * cfg.UATM_TO_ATM   # 1 ppm ≈ 1 µatm

    da = xr.DataArray(
        co2_atm.values,
        coords={"time": co2_atm.index.values.astype("datetime64[ns]")},
        dims=["time"],
        name="spco2_atm",
        attrs={"units": "atm", "long_name": "Atmospheric CO2 mole fraction (global mean)"},
    )
    return da


def load_multiobs_flux() -> xr.DataArray:
    """
    Load the observation-based air-sea CO2 flux from the CMEMS MULTIOBS
    surface carbon L4 product (SOCAT-trained neural network).

    Used only for validation — never mixed into the primary computation.

    Variable fgco2 units may be kg C / m² / s or mol C / m² / yr:
    check the NetCDF attributes and convert to mol/m²/yr.
        1 kg C / m² / s = (1000/12.011) mol / m² / s * 3.156e7 s/yr
                        ≈ 2.628e9 mol / m² / yr
    The conversion is applied if units attribute contains 'kg'.
    """
    f = cfg.DATA_DIR / "multiobs_surface_carbon.nc"
    if not f.exists():
        print("[warn] MULTIOBS validation file not found — skipping.")
        return None

    print("[load] MULTIOBS surface carbon flux ...")
    ds = xr.open_dataset(f, chunks="auto")
    da = ds[cfg.MULTIOBS_FLUX_VAR]

    # Auto-detect units and convert to mol / m² / yr
    units = da.attrs.get("units", "")
    if "kg" in units.lower():
        print(f"[unit] Converting fgco2 from {units} to mol/m²/yr ...")
        # mol/m²/yr = (kg C/m²/s) * (1000 g/kg) / (12.011 g/mol) * (s/yr)
        da = da * (1000.0 / cfg.MOL_C_TO_G) * cfg.S_TO_YR
        da.attrs["units"] = "mol C m-2 yr-1"
    elif "mol" in units.lower() and "s" in units.lower():
        print(f"[unit] Converting fgco2 from {units} to mol/m²/yr ...")
        da = da * cfg.S_TO_YR
        da.attrs["units"] = "mol C m-2 yr-1"

    return da.rename("fgco2_obs")


# ===========================================================================
# HARMONIZE & MERGE
# ===========================================================================

def harmonize_time(da: xr.DataArray | xr.Dataset, name: str) -> xr.DataArray | xr.Dataset:
    """
    Ensure the time coordinate is numpy datetime64[ns] (not cftime objects).
    CMEMS files sometimes use cftime, which causes issues with xr.merge.
    """
    if hasattr(da, "time"):
        try:
            da["time"] = da.indexes["time"].to_datetimeindex()
        except Exception:
            pass   # already datetime64 or not a standard calendar
    return da


def build_ocean_mask(spco2_ocean: xr.DataArray) -> xr.DataArray:
    """
    Derive an ocean mask from the BGC hindcast pCO2 field.
    Land/ice pixels are NaN in the PISCES output; ocean pixels have valid values.
    We take the mask as the 'any valid time step' union, then it is static.

    Returns a 2-D DataArray (lat, lon) with 1 = ocean, 0 = land/ice.
    """
    # True wherever at least one month is non-NaN
    valid = spco2_ocean.notnull().any(dim="time")
    mask = valid.astype(np.float32)
    mask.attrs["long_name"] = "Ocean mask (1=ocean, 0=land/ice)"
    return mask.rename("ocean_mask")


def broadcast_atm_co2(
    co2_atm: xr.DataArray,
    reference: xr.DataArray,
) -> xr.DataArray:
    """
    Broadcast the 1-D atmospheric CO2 time series to the (time, lat, lon) grid
    of the CMEMS product, so it can participate in pixel-wise arithmetic.

    The NOAA time series and CMEMS time axes may not align exactly
    (different day-of-month anchors). We interpolate the 1-D series to
    the CMEMS monthly time coordinate before broadcasting.
    """
    # Interpolate NOAA monthly series to CMEMS time axis (linear OK for CO2)
    co2_interp = co2_atm.interp(time=reference.time, method="linear")

    # Expand to (time, lat, lon) — no copy of data, just coordinates
    co2_broad = co2_interp.expand_dims(
        {"latitude": reference.latitude, "longitude": reference.longitude}
    )
    co2_broad = co2_broad.transpose("time", "latitude", "longitude")
    co2_broad.attrs.update(co2_atm.attrs)
    return co2_broad.rename("spco2_atm")


# ===========================================================================
# MAIN
# ===========================================================================

def main():
    out_file = cfg.DATA_DIR / "processed_surface.nc"
    if out_file.exists():
        print(f"[skip] {out_file.name} already exists. Delete to reprocess.")
        return

    # --- Load each component ---
    spco2_ocean = load_spco2_ocean()
    phy         = load_sst_sss()
    wind_speed  = load_wind_speed()
    co2_atm     = load_atmospheric_co2()
    fgco2_obs   = load_multiobs_flux()

    # --- Harmonize time axes ---
    spco2_ocean = harmonize_time(spco2_ocean, "spco2_ocean")
    phy         = harmonize_time(phy, "phy")
    if wind_speed is not None:
        wind_speed = harmonize_time(wind_speed, "wind_speed")

    # --- Ocean mask ---
    ocean_mask = build_ocean_mask(spco2_ocean)

    # --- Broadcast atmospheric CO2 to the CMEMS lat/lon grid ---
    spco2_atm = broadcast_atm_co2(co2_atm, spco2_ocean)

    # --- Align time axes: use the BGC hindcast as reference ---
    ref_time = spco2_ocean.time

    sst = phy["sst"].interp(time=ref_time, method="nearest")
    sss = phy["sss"].interp(time=ref_time, method="nearest")

    if wind_speed is not None:
        wind_speed = wind_speed.interp(time=ref_time, method="nearest")

    if fgco2_obs is not None:
        fgco2_obs = fgco2_obs.interp(time=ref_time, method="nearest")

    # --- Build merged Dataset ---
    ds_vars = {
        "spco2_ocean": spco2_ocean,
        "spco2_atm":   spco2_atm,
        "sst":         sst,
        "sss":         sss,
        "ocean_mask":  ocean_mask,
    }
    if wind_speed is not None:
        ds_vars["wind_speed"] = wind_speed
    if fgco2_obs is not None:
        ds_vars["fgco2_obs"] = fgco2_obs

    ds_out = xr.Dataset(ds_vars)

    # --- Global attributes ---
    ds_out.attrs = {
        "title":       "Stage 1 — Surface CO2 flux input data",
        "project":     "Ocean Carbon Cycle — Net Flux Project",
        "created_by":  "02_preprocess.py",
        "conventions": "CF-1.8",
        "references":  (
            "BGC hindcast: GLOBAL_MULTIYEAR_BGC_001_029 | "
            "Physics: GLOBAL_MULTIYEAR_PHY_001_030 | "
            "Atm CO2: NOAA GML https://gml.noaa.gov/ccgg/trends/"
        ),
    }

    # --- Save ---
    print(f"\n[save] Writing → {out_file} ...")
    ds_out.to_netcdf(out_file)
    print("[done] Preprocessing complete.\n")
    print(ds_out)


if __name__ == "__main__":
    main()
