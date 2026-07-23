"""
01_download_data.py
===================
Downloads all data required for Stage 1 (surface net air-sea CO2 flux).

What this script does:
    1. Downloads monthly surface pCO2 from the Copernicus BGC hindcast
       (GLOBAL_MULTIYEAR_BGC_001_029) over the full record 1993–2026.
    2. Downloads monthly SST and SSS from the Copernicus physical reanalysis
       (GLOBAL_MULTIYEAR_PHY_001_030) — needed for K0 and Schmidt number.
    3. Downloads the observation-based flux product
       (MULTIOBS_GLO_BGC_CARBON_SURFACE_MYNRT_015_008) for cross-validation.
    4. Downloads the NOAA GML global average atmospheric CO2 monthly series.

Data sources:
    - CMEMS BGC hindcast:
        https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_BGC_001_029
    - CMEMS physical reanalysis (GLORYS12V1):
        https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030
    - CMEMS surface carbon L4:
        https://data.marine.copernicus.eu/product/MULTIOBS_GLO_BGC_CARBON_SURFACE_MYNRT_015_008
    - NOAA GML global CO2:
        https://gml.noaa.gov/ccgg/trends/gl_data.html

Authentication:
    Set your CMEMS credentials either:
        (a) via `copernicusmarine login` in terminal (stores ~/.copernicusmarine), OR
        (b) via environment variables COPERNICUSMARINE_SERVICE_USERNAME and
            COPERNICUSMARINE_SERVICE_PASSWORD

Usage:
    python scripts/01_download_data.py [--test]

    --test    : Download only 2010 (one year) as a quick sanity check.

Requirements:
    copernicusmarine >= 1.3
    requests
"""

import argparse
import sys
from pathlib import Path

import copernicusmarine
import requests

# ── bring project config in regardless of where the script is called from ──
sys.path.insert(0, str(Path(__file__).parent.parent))
import config as cfg


# ===========================================================================
# HELPERS
# ===========================================================================

def _time_range(test: bool) -> tuple[str, str]:
    """Return (start, end) date strings depending on test mode."""
    if test:
        return cfg.TEST_TIME_START, cfg.TEST_TIME_END
    return cfg.TIME_START, cfg.TIME_END


def download_bgc_hindcast(start: str, end: str) -> None:
    """
    Download monthly surface pCO2 from the CMEMS global BGC hindcast.

    Product : GLOBAL_MULTIYEAR_BGC_001_029
    Variable: spco2  (surface partial pressure of CO2, Pa)
    Depth   : surface level only (first layer)
    Grid    : 0.25° global
    """
    out_file = cfg.DATA_DIR / "bgc_hindcast_spco2.nc"
    if out_file.exists():
        print(f"[skip] {out_file.name} already exists.")
        return

    print(f"[download] BGC hindcast pCO2  {start} → {end} ...")
    copernicusmarine.subset(
        dataset_id     = cfg.BGC_HINDCAST_DATASET,
        variables      = [cfg.BGC_PCO2_VAR],
        minimum_longitude = cfg.BBOX[0],
        maximum_longitude = cfg.BBOX[1],
        minimum_latitude  = cfg.BBOX[2],
        maximum_latitude  = cfg.BBOX[3],
        minimum_depth  = cfg.DEPTH_MIN,
        maximum_depth  = cfg.DEPTH_MAX,
        start_datetime = start,
        end_datetime   = end,
        output_filename = str(out_file),
        force_download  = True,
    )
    print(f"[ok] Saved → {out_file}")


def download_physical_reanalysis(start: str, end: str) -> None:
    """
    Download monthly SST and SSS from CMEMS GLORYS12V1 physical reanalysis.

    Product : GLOBAL_MULTIYEAR_PHY_001_030
    Variables: thetao (potential temperature, °C), so (practical salinity, PSU)
    Depth   : surface level only
    Grid    : 0.083° → will be used as-is; regrid to 0.25° in 02_preprocess.py
    """
    out_file = cfg.DATA_DIR / "phy_reanalysis_sst_sss.nc"
    if out_file.exists():
        print(f"[skip] {out_file.name} already exists.")
        return

    print(f"[download] Physical reanalysis SST/SSS  {start} → {end} ...")
    copernicusmarine.subset(
        dataset_id     = cfg.PHY_REANALYSIS_DATASET,
        variables      = [cfg.PHY_TEMP_VAR, cfg.PHY_SAL_VAR],
        minimum_longitude = cfg.BBOX[0],
        maximum_longitude = cfg.BBOX[1],
        minimum_latitude  = cfg.BBOX[2],
        maximum_latitude  = cfg.BBOX[3],
        minimum_depth  = cfg.DEPTH_MIN,
        maximum_depth  = cfg.DEPTH_MAX,
        start_datetime = start,
        end_datetime   = end,
        output_filename = str(out_file),
        force_download  = True,
    )
    print(f"[ok] Saved → {out_file}")


def download_multiobs_surface(start: str, end: str) -> None:
    """
    Download the observation-based surface carbon L4 product for validation.

    Product : MULTIOBS_GLO_BGC_CARBON_SURFACE_MYNRT_015_008
    Variables: fgco2 (air-sea CO2 flux), spco2 (surface pCO2)
    Note    : This is a SOCAT-trained neural-network product — used only
              for cross-validation of our reconstructed flux, NOT as primary data.

    ⚠️  Confirm the exact dataset name on the CMEMS portal before running —
        MULTIOBS product dataset names sometimes change between catalogue updates.
        Product page: https://data.marine.copernicus.eu/product/
                      MULTIOBS_GLO_BGC_CARBON_SURFACE_MYNRT_015_008/description
    """
    out_file = cfg.DATA_DIR / "multiobs_surface_carbon.nc"
    if out_file.exists():
        print(f"[skip] {out_file.name} already exists.")
        return

    print(f"[download] MULTIOBS surface carbon  {start} → {end} ...")
    copernicusmarine.subset(
        dataset_id     = cfg.MULTIOBS_SURFACE_DATASET,
        variables      = [cfg.MULTIOBS_FLUX_VAR, cfg.MULTIOBS_PCO2_VAR],
        minimum_longitude = cfg.BBOX[0],
        maximum_longitude = cfg.BBOX[1],
        minimum_latitude  = cfg.BBOX[2],
        maximum_latitude  = cfg.BBOX[3],
        start_datetime = start,
        end_datetime   = end,
        output_filename = str(out_file),
        force_download  = True,
    )
    print(f"[ok] Saved → {out_file}")


def download_noaa_co2() -> None:
    """
    Download NOAA GML global average monthly atmospheric CO2 (ppm).

    Source: https://gml.noaa.gov/ccgg/trends/gl_data.html
    File  : co2_mm_gl.csv  (global mean marine surface CO2, monthly)

    Column description (from NOAA header):
        year, month, decimal_date, average, trend, ...
        'average' = monthly mean CO2 in ppm

    The file has a multi-line comment header starting with '#' — pandas
    will skip these in 02_preprocess.py using comment='#'.

    Note: 1 ppm ≈ 1 µatm at standard atmosphere pressure (1 atm total).
    """
    out_file = cfg.NOAA_CO2_FILE
    if out_file.exists():
        print(f"[skip] {out_file.name} already exists.")
        return

    url = "https://gml.noaa.gov/webdata/ccgg/trends/co2/co2_mm_gl.csv"
    print(f"[download] NOAA GML CO2 from {url} ...")
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    out_file.write_bytes(r.content)
    print(f"[ok] Saved → {out_file}")


def download_era5_wind() -> None:
    """
    ERA5 10-metre wind components (u10, v10) — monthly means.

    ERA5 is downloaded via the CDS API (separate from CMEMS):
        https://cds.climate.copernicus.eu/

    Prerequisites:
        pip install cdsapi
        Set up ~/.cdsapirc with your CDS UID and API key:
            https://cds.climate.copernicus.eu/api-how-to

    The function below shows the correct CDS API call. Uncomment and run
    after setting up your CDS credentials.

    Product : ERA5 monthly averaged reanalysis on single levels
    Dataset : reanalysis-era5-single-levels-monthly-means
    Variables: 10m_u_component_of_wind, 10m_v_component_of_wind
    Years   : 1993–2025 (download in blocks if the request is too large)

    Output file: data/era5_wind10m_monthly.nc
    """
    out_file = cfg.ERA5_WIND_FILE
    if out_file.exists():
        print(f"[skip] {out_file.name} already exists.")
        return

    print("[info] ERA5 wind requires CDS API — see docstring for setup.")
    print("[info] Uncomment the block below once ~/.cdsapirc is configured.\n")

    # ── uncomment after CDS setup ──────────────────────────────────────────
    # import cdsapi
    # c = cdsapi.Client()
    # c.retrieve(
    #     "reanalysis-era5-single-levels-monthly-means",
    #     {
    #         "product_type": "monthly_averaged_reanalysis",
    #         "variable": [
    #             "10m_u_component_of_wind",
    #             "10m_v_component_of_wind",
    #         ],
    #         "year": [str(y) for y in range(1993, 2026)],
    #         "month": [f"{m:02d}" for m in range(1, 13)],
    #         "time": "00:00",
    #         "format": "netcdf",
    #         "area": [90, -180, -90, 180],   # global
    #     },
    #     str(out_file),
    # )
    # print(f"[ok] Saved → {out_file}")
    # ──────────────────────────────────────────────────────────────────────


# ===========================================================================
# MAIN
# ===========================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Download all Stage 1 data from CMEMS, NOAA, and ERA5."
    )
    parser.add_argument(
        "--test",
        action="store_true",
        help="Download one year only (2010) for a quick sanity check.",
    )
    args = parser.parse_args()

    start, end = _time_range(args.test)
    if args.test:
        print(f"\n[TEST MODE] Downloading {start} → {end} only.\n")

    download_bgc_hindcast(start, end)
    download_physical_reanalysis(start, end)
    download_multiobs_surface(start, end)
    download_noaa_co2()
    download_era5_wind()

    print("\n[done] All Stage 1 data downloaded.\n")


if __name__ == "__main__":
    main()
