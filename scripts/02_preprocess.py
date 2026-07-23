"""
02_preprocess.py
================
Load, regrid, unit-convert, and harmonize all downloaded datasets onto a
common 0.25° monthly grid ready for flux computation.

Steps:
    1. Load BGC hindcast surface pCO2 (already 0.25°) — driver for the
       original (v1.0.x) reconstruction.
    2. Load MULTIOBS surface pCO2 (SOCAT-NN) and regrid to 0.25° — driver
       for the improved (v1.1.0) reconstruction.
    3. Load physical reanalysis SST and SSS (0.083°) and regrid to 0.25°.
    4. Load ERA5 monthly wind, compute scalar speed and regrid to 0.25°.
    5. [NEW v1.1.0] Load ERA5 daily wind, compute sub-monthly variance σ²_u
       per month per pixel, save as wind_variance for the variance correction.
    6. Load NOAA atmospheric CO2 CSV, parse to a DataArray.
    7. Load MULTIOBS validation flux (fgco2) for cross-validation.
    8. Harmonize time axes, align all on the BGC hindcast time coordinate.
    9. Build ocean mask from PISCES NaN pattern.
    10. Save one merged NetCDF to data/processed_surface.nc.

Outputs:
    data/processed_surface.nc   — merged xarray Dataset with:
        spco2_ocean     [atm]           — PISCES surface pCO2 (original driver)
        spco2_ocean_obs [atm]           — MULTIOBS/SOCAT surface pCO2 (improved driver)
        spco2_atm       [atm]           — atmospheric pCO2 broadcast to lat/lon
        sst             [°C]            — sea surface temperature
        sss             [PSU]           — sea surface salinity
        wind_speed      [m/s]           — monthly-mean 10 m scalar wind speed
        wind_variance   [m²/s²]         — sub-monthly wind speed variance σ²_u
        fgco2_obs       [mol/m²/yr]     — MULTIOBS validation flux
        ocean_mask                      — 1=ocean, 0=land/ice

Changelog:
    v1.1.0 — Added spco2_ocean_obs (MULTIOBS pCO2 driver for improved reconstruction);
              added wind_variance (sub-monthly σ²_u from ERA5 daily winds);
              added MULTIOBS spatial regrid step.

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
    Load surface ocean pCO2 from the PISCES BGC hindcast and convert Pa → atm.
    This is the pCO2 driver for the ORIGINAL (v1.0.x) reconstruction.

    The CMEMS variable 'spco2' is in Pa.  1 atm = 101325 Pa.
    """
    print("[load] BGC hindcast pCO2 (PISCES) ...")
    ds = xr.open_dataset(cfg.DATA_DIR / "bgc_hindcast_spco2.nc", chunks="auto")
    da = ds[cfg.BGC_PCO2_VAR]

    if "depth" in da.dims:
        da = da.isel(depth=0, drop=True)

    da = da * cfg.PA_TO_ATM
    da.attrs["units"]     = "atm"
    da.attrs["long_name"] = "Surface ocean pCO2 (PISCES BGC hindcast)"
    da.attrs["source"]    = "GLOBAL_MULTIYEAR_BGC_001_029"
    return da.rename("spco2_ocean")


def load_spco2_ocean_obs() -> xr.DataArray | None:
    """
    Load surface ocean pCO2 from the CMEMS MULTIOBS product (SOCAT neural network).
    This is the pCO2 driver for the IMPROVED (v1.1.0) reconstruction.

    The MULTIOBS spco2 units may be Pa or µatm — auto-detected from attributes.
    Returns None if the MULTIOBS file has not been downloaded yet.
    """
    f = cfg.DATA_DIR / "multiobs_surface_carbon.nc"
    if not f.exists():
        print("[warn] MULTIOBS file not found — improved pCO2 driver unavailable.")
        return None

    print("[load] MULTIOBS surface pCO2 (SOCAT-NN) ...")
    ds = xr.open_dataset(f, chunks="auto")

    if cfg.MULTIOBS_PCO2_VAR not in ds:
        print(f"[warn] Variable '{cfg.MULTIOBS_PCO2_VAR}' not in MULTIOBS file — skipping.")
        return None

    da = ds[cfg.MULTIOBS_PCO2_VAR]

    # Auto-detect units and convert to atm
    units = da.attrs.get("units", "").lower()
    if "pa" in units and "µ" not in units and "u" not in units:
        # Pascals → atm
        print(f"[unit] MULTIOBS spco2: Pa → atm")
        da = da * cfg.PA_TO_ATM
    elif "uatm" in units or "µatm" in units or "ppm" in units:
        # µatm ≈ ppm → atm
        print(f"[unit] MULTIOBS spco2: µatm → atm")
        da = da * cfg.UATM_TO_ATM
    else:
        print(f"[warn] MULTIOBS spco2 units='{units}' — assuming µatm, converting ×1e-6")
        da = da * cfg.UATM_TO_ATM

    da.attrs["units"]     = "atm"
    da.attrs["long_name"] = "Surface ocean pCO2 (MULTIOBS SOCAT-NN)"
    da.attrs["source"]    = "MULTIOBS_GLO_BGC_CARBON_SURFACE_MYNRT_015_008"
    return da.rename("spco2_ocean_obs")


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

    ds_bgc     = xr.open_dataset(cfg.DATA_DIR / "bgc_hindcast_spco2.nc", chunks="auto")
    target_lat = ds_bgc["latitude"]
    target_lon = ds_bgc["longitude"]

    print("[regrid] SST/SSS 0.083° → 0.25° ...")
    sst = sst.interp(latitude=target_lat, longitude=target_lon, method="linear")
    sss = sss.interp(latitude=target_lat, longitude=target_lon, method="linear")

    sst.attrs["units"] = "degC"
    sss.attrs["units"] = "PSU"

    return xr.Dataset({"sst": sst.rename("sst"), "sss": sss.rename("sss")})


def load_wind_speed() -> xr.DataArray | None:
    """
    Load ERA5 monthly 10 m wind components and compute scalar wind speed.
    Returns the magnitude |u| = sqrt(u10² + v10²) in m/s, regridded to 0.25°.
    Returns None if the ERA5 monthly file has not been downloaded yet.
    """
    wind_file = cfg.ERA5_WIND_FILE
    if not wind_file.exists():
        print("[warn] ERA5 monthly wind file not found — returning None.")
        print("       Run 01_download_data.py with CDS credentials first.")
        return None

    print("[load] ERA5 monthly 10 m wind ...")
    ds_era5 = xr.open_dataset(wind_file, chunks="auto")

    u10 = ds_era5[cfg.ERA5_U10_VAR]
    v10 = ds_era5[cfg.ERA5_V10_VAR]

    # ERA5 CDS sometimes uses "valid_time" instead of "time"
    if "valid_time" in u10.dims:
        u10 = u10.rename({"valid_time": "time"})
        v10 = v10.rename({"valid_time": "time"})

    wspd = np.sqrt(u10**2 + v10**2)
    wspd.attrs["units"]     = "m/s"
    wspd.attrs["long_name"] = "10 m scalar wind speed (monthly mean)"

    ds_bgc = xr.open_dataset(cfg.DATA_DIR / "bgc_hindcast_spco2.nc", chunks="auto")
    wspd = wspd.interp(
        latitude=ds_bgc["latitude"],
        longitude=ds_bgc["longitude"],
        method="linear",
    )
    return wspd.rename("wind_speed")


def load_wind_variance() -> xr.DataArray | None:
    """
    Compute the sub-monthly wind speed variance σ²_u from ERA5 daily winds.

    The Wanninkhof (2014) gas transfer velocity k = a·u²·(Sc/660)^(-0.5) was
    calibrated against the *full wind speed distribution*, not monthly means.
    Applying it to monthly-mean u systematically underestimates k, because:
        <u²> = <u>² + σ²_u
    where σ²_u is the sub-monthly variance of the scalar wind speed.

    The corrected piston velocity uses:
        k_corr = a · (<u>² + σ²_u) · (Sc/660)^(-0.5)

    This function:
        1. Loads ERA5 daily 10 m wind components (u10, v10).
        2. Computes daily scalar speed: wspd_daily = sqrt(u10² + v10²).
        3. Groups by year-month, computes: σ²_u = Var(wspd_daily) per pixel.
        4. Regrids to the 0.25° BGC grid.

    Returns None if the ERA5 daily file has not been downloaded.

    Reference:
        Wanninkhof (2014) Eq. 3 and discussion; Boutin & Etcheto (1991) for
        the variance correction derivation.
    """
    daily_file = cfg.ERA5_WIND_DAILY_FILE
    if not daily_file.exists():
        print("[warn] ERA5 daily wind file not found — wind variance correction unavailable.")
        print("       Download ERA5 daily winds via CDS to enable this improvement.")
        return None

    print("[load] ERA5 daily winds for sub-monthly variance correction ...")
    ds_daily = xr.open_dataset(daily_file, chunks="auto")

    u10 = ds_daily[cfg.ERA5_U10_VAR]
    v10 = ds_daily[cfg.ERA5_V10_VAR]

    if "valid_time" in u10.dims:
        u10 = u10.rename({"valid_time": "time"})
        v10 = v10.rename({"valid_time": "time"})

    # Daily scalar wind speed [m/s]
    wspd_daily = np.sqrt(u10**2 + v10**2)

    # Monthly variance of daily wind speed: σ²_u [m²/s²]
    print("[compute] Sub-monthly wind variance σ²_u ...")
    wind_var = wspd_daily.resample(time="1ME").var(dim="time")
    wind_var.attrs["units"]     = "m2 s-2"
    wind_var.attrs["long_name"] = "Sub-monthly wind speed variance (σ²_u)"
    wind_var.attrs["note"]      = (
        "Used for Wanninkhof (2014) variance correction: "
        "k ∝ (<u>² + σ²_u) instead of <u>²"
    )

    # Regrid to 0.25° BGC grid
    ds_bgc = xr.open_dataset(cfg.DATA_DIR / "bgc_hindcast_spco2.nc", chunks="auto")
    wind_var = wind_var.interp(
        latitude=ds_bgc["latitude"],
        longitude=ds_bgc["longitude"],
        method="linear",
    )
    return wind_var.rename("wind_variance")


def load_atmospheric_co2() -> xr.DataArray:
    """
    Load NOAA GML global average monthly CO2 (ppm) and convert to atm.
    """
    print("[load] NOAA GML atmospheric CO2 ...")
    df = pd.read_csv(cfg.NOAA_CO2_FILE, comment="#")
    df.columns = df.columns.str.strip()

    for col in ["year", "month", "decimal", "average", "average_unc", "trend", "trend_unc"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=["average"])
    df = df[df["average"] > 0]

    df["time"] = pd.to_datetime(
        df["year"].astype(int).astype(str)
        + "-"
        + df["month"].astype(int).astype(str).str.zfill(2)
    )
    df = df.set_index("time").sort_index()

    co2_atm = df["average"] * cfg.UATM_TO_ATM
    return xr.DataArray(
        co2_atm.values,
        coords={"time": co2_atm.index.values.astype("datetime64[ns]")},
        dims=["time"],
        name="spco2_atm",
        attrs={"units": "atm", "long_name": "Atmospheric CO2 mole fraction (global mean)"},
    )


def load_multiobs_flux() -> xr.DataArray | None:
    """
    Load the MULTIOBS pre-computed air-sea CO2 flux for validation.
    Converts to mol C m⁻² yr⁻¹ if needed (auto-detects kg or mol s⁻¹ input).
    Returns None if the file has not been downloaded.
    """
    f = cfg.DATA_DIR / "multiobs_surface_carbon.nc"
    if not f.exists():
        print("[warn] MULTIOBS validation file not found — skipping.")
        return None

    print("[load] MULTIOBS air-sea CO2 flux (validation) ...")
    ds    = xr.open_dataset(f, chunks="auto")
    da    = ds[cfg.MULTIOBS_FLUX_VAR]
    units = da.attrs.get("units", "")

    if "kg" in units.lower():
        print(f"[unit] fgco2: {units} → mol C m⁻² yr⁻¹")
        da = da * (1000.0 / cfg.MOL_C_TO_G) * cfg.S_TO_YR
        da.attrs["units"] = "mol C m-2 yr-1"
    elif "mol" in units.lower() and "s" in units.lower():
        print(f"[unit] fgco2: {units} → mol C m⁻² yr⁻¹")
        da = da * cfg.S_TO_YR
        da.attrs["units"] = "mol C m-2 yr-1"

    return da.rename("fgco2_obs")


# ===========================================================================
# HARMONIZE & MERGE HELPERS
# ===========================================================================

def harmonize_time(da, name: str):
    """Convert cftime → numpy datetime64[ns] if needed."""
    if hasattr(da, "time"):
        try:
            da["time"] = da.indexes["time"].to_datetimeindex()
        except Exception:
            pass
    return da


def build_ocean_mask(spco2_ocean: xr.DataArray) -> xr.DataArray:
    """
    2-D ocean mask (lat, lon) derived from PISCES NaN pattern.
    1 = ocean pixel with at least one valid month; 0 = land/ice.
    """
    valid = spco2_ocean.notnull().any(dim="time")
    mask  = valid.astype(np.float32)
    mask.attrs["long_name"] = "Ocean mask (1=ocean, 0=land/ice)"
    return mask.rename("ocean_mask")


def broadcast_atm_co2(co2_atm: xr.DataArray, reference: xr.DataArray) -> xr.DataArray:
    """
    Interpolate the 1-D NOAA CO2 series to the CMEMS time axis, then
    broadcast to (time, lat, lon) so it can be used in pixel-wise arithmetic.
    """
    co2_interp = co2_atm.interp(time=reference.time, method="linear")
    co2_broad  = co2_interp.expand_dims(
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

    # --- Load all components ---
    spco2_ocean     = load_spco2_ocean()
    spco2_ocean_obs = load_spco2_ocean_obs()   # NEW v1.1.0 — may be None
    phy             = load_sst_sss()
    wind_speed      = load_wind_speed()
    wind_variance   = load_wind_variance()     # NEW v1.1.0 — may be None
    co2_atm         = load_atmospheric_co2()
    fgco2_obs       = load_multiobs_flux()

    # --- Harmonize time axes ---
    spco2_ocean = harmonize_time(spco2_ocean, "spco2_ocean")
    phy         = harmonize_time(phy, "phy")
    if wind_speed is not None:
        wind_speed = harmonize_time(wind_speed, "wind_speed")
    if wind_variance is not None:
        wind_variance = harmonize_time(wind_variance, "wind_variance")
    if spco2_ocean_obs is not None:
        spco2_ocean_obs = harmonize_time(spco2_ocean_obs, "spco2_ocean_obs")

    # --- Ocean mask (from PISCES NaN pattern) ---
    ocean_mask = build_ocean_mask(spco2_ocean)

    # --- Broadcast atmospheric CO2 ---
    spco2_atm = broadcast_atm_co2(co2_atm, spco2_ocean)

    # --- Align all to BGC hindcast time axis ---
    ref_time = spco2_ocean.time

    sst = phy["sst"].interp(time=ref_time, method="nearest")
    sss = phy["sss"].interp(time=ref_time, method="nearest")

    if wind_speed is not None:
        wind_speed = wind_speed.interp(time=ref_time, method="nearest")

    if wind_variance is not None:
        wind_variance = wind_variance.interp(time=ref_time, method="nearest")

    if fgco2_obs is not None:
        print("[regrid] MULTIOBS flux → BGC 0.25° grid ...")
        fgco2_obs = fgco2_obs.interp(
            latitude=spco2_ocean.latitude,
            longitude=spco2_ocean.longitude,
            method="linear",
        ).interp(time=ref_time, method="nearest")

    if spco2_ocean_obs is not None:
        print("[regrid] MULTIOBS pCO2 → BGC 0.25° grid ...")
        spco2_ocean_obs = spco2_ocean_obs.interp(
            latitude=spco2_ocean.latitude,
            longitude=spco2_ocean.longitude,
            method="linear",
        ).interp(time=ref_time, method="nearest")

    # --- Assemble merged Dataset ---
    ds_vars = {
        "spco2_ocean": spco2_ocean,
        "spco2_atm":   spco2_atm,
        "sst":         sst,
        "sss":         sss,
        "ocean_mask":  ocean_mask,
    }
    if spco2_ocean_obs is not None:
        ds_vars["spco2_ocean_obs"] = spco2_ocean_obs
    if wind_speed is not None:
        ds_vars["wind_speed"] = wind_speed
    if wind_variance is not None:
        ds_vars["wind_variance"] = wind_variance
    if fgco2_obs is not None:
        ds_vars["fgco2_obs"] = fgco2_obs

    # Print dimension summary for quick sanity check
    for name, var in ds_vars.items():
        print(f"  {name:20s}  {dict(var.sizes)}")

    ds_out = xr.Dataset(ds_vars)
    ds_out.attrs = {
        "title":       "Stage 1 — Surface CO2 flux input data",
        "version":     "1.1.0",
        "project":     "Ocean Carbon Cycle — Net Flux Project",
        "created_by":  "02_preprocess.py",
        "conventions": "CF-1.8",
        "references":  (
            "BGC hindcast: GLOBAL_MULTIYEAR_BGC_001_029 | "
            "Physics: GLOBAL_MULTIYEAR_PHY_001_030 | "
            "MULTIOBS: MULTIOBS_GLO_BGC_CARBON_SURFACE_MYNRT_015_008 | "
            "Atm CO2: NOAA GML https://gml.noaa.gov/ccgg/trends/"
        ),
    }

    print(f"\n[save] Writing → {out_file} ...")
    ds_out.to_netcdf(out_file)
    print("[done] Preprocessing complete.\n")
    print(ds_out)


if __name__ == "__main__":
    main()
