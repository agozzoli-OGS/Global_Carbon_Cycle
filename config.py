"""
config.py
=========
Central configuration file for Stage 1 of the Ocean Carbon Cycle project.
All paths, constants, dataset parameters, and tunable settings live here.
No computation happens in this file — import it from every other script.

Project: Ocean Carbon Cycle — Net Flux & Layered Accumulation
Stage:   1 — Surface net air-sea CO2 flux, J(0 m)
Authors: [your names]
Date:    2026-07
"""

# ===========================================================================
# PATHS
# ===========================================================================

from pathlib import Path

# Root of the stage 1 working directory
ROOT = Path(__file__).parent

# Where raw downloaded NetCDF files are stored
DATA_DIR = ROOT / "data"

# Where figures and output CSVs are saved
OUT_DIR  = ROOT / "output"
FIG_DIR  = OUT_DIR / "figures"

# Make directories if they don't exist yet
for d in [DATA_DIR, OUT_DIR, FIG_DIR]:
    d.mkdir(parents=True, exist_ok=True)


# ===========================================================================
# TIME RANGE
# ===========================================================================

# Full record of GLOBAL_MULTIYEAR_BGC_001_029
TIME_START = "1993-01-01"
TIME_END   = "2026-04-30"

# For quick test runs, narrow to a single year
TEST_TIME_START = "2010-01-01"
TEST_TIME_END   = "2010-12-31"


# ===========================================================================
# COPERNICUS MARINE SERVICE — PRODUCT & VARIABLE IDs
# ===========================================================================

# --- Primary BGC hindcast (pCO2 at the ocean surface) ---
BGC_HINDCAST_ID      = "GLOBAL_MULTIYEAR_BGC_001_029"
BGC_HINDCAST_DATASET = "cmems_mod_glo_bgc_my_0.25deg_P1M-m"   # monthly means
BGC_PCO2_VAR         = "spco2"      # surface partial pressure of CO2 (Pa or uatm — check units on download)

# --- Physical reanalysis (SST, SSS — needed for K0 and k) ---
PHY_REANALYSIS_ID      = "GLOBAL_MULTIYEAR_PHY_001_030"
PHY_REANALYSIS_DATASET = "cmems_mod_glo_phy_my_0.083deg_P1M-m"  # monthly means
PHY_TEMP_VAR           = "thetao"   # potential temperature (deg C)
PHY_SAL_VAR            = "so"       # practical salinity (PSU)

# --- Cell thickness (from BGC hindcast, needed for volume integrals later) ---
THICKNESS_VAR = "thkcello"

# --- Validation: observation-based surface carbon L4 product ---
MULTIOBS_SURFACE_ID      = "MULTIOBS_GLO_BGC_CARBON_SURFACE_MYNRT_015_008"
MULTIOBS_SURFACE_DATASET = "dataset-carbon-rep-monthly"   # confirm exact dataset name on CMEMS portal
MULTIOBS_FLUX_VAR        = "fgco2"    # air-sea CO2 flux (mol C / m2 / yr or kg C / m2 / s — check units)
MULTIOBS_PCO2_VAR        = "spco2"    # surface pCO2 from SOCAT-NN (for cross-check)

# --- Spatial extent (global ocean) ---
# Copernicus toolbox uses (lon_min, lon_max, lat_min, lat_max)
BBOX = (-180.0, 180.0, -90.0, 90.0)

# --- Depth: for Stage 1 we only need the surface level ---
DEPTH_MIN = 0.0   # m
DEPTH_MAX = 1.0   # m — take just the first model layer


# ===========================================================================
# ATMOSPHERIC CO2 — NOAA GML MARINE BOUNDARY LAYER REFERENCE
# ===========================================================================

# Download from: https://gml.noaa.gov/ccgg/trends/gl_data.html
# File: co2_mm_mlo.csv  (Mauna Loa monthly means, ppm)
# OR global marine boundary layer product:
# https://gml.noaa.gov/ccgg/mbl/index.html
# File: co2_GHGreference.GLOBALAV.monthly.txt

NOAA_CO2_FILE = DATA_DIR / "co2_mm_gl.csv"   # global average monthly file
NOAA_CO2_COLNAME = "average"                  # column holding the ppm value


# ===========================================================================
# ERA5 WIND SPEED (for gas-transfer velocity k)
# ===========================================================================

# Source: Copernicus Climate Data Store (CDS) — separate from CMEMS
# Product: ERA5 monthly averaged reanalysis on single levels
# Variable: 10 metre U/V wind components  (u10, v10)
# Units: m/s
# Download page: https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels-monthly-means
# After download, compute scalar speed: wspd = sqrt(u10**2 + v10**2)

ERA5_WIND_FILE = DATA_DIR / "era5_wind10m_monthly.nc"
ERA5_U10_VAR   = "u10"
ERA5_V10_VAR   = "v10"


# ===========================================================================
# GAS TRANSFER VELOCITY — WANNINKHOF (2014)
# ===========================================================================

# Wanninkhof, R. (2014). Relationship between wind speed and gas exchange
# over the ocean revisited. Limnology and Oceanography: Methods, 12(6), 351–362.
# DOI: 10.4319/lom.2014.12.351
#
# k = a · <u^2> · (Sc/660)^(-0.5)
# where:
#   a  = 0.251  [cm/hr / (m/s)^2]  -- Wanninkhof 2014 coefficient
#   u  = 10 m wind speed [m/s]
#   Sc = Schmidt number for CO2 in seawater (function of SST)
#   k  is returned in cm/hr; convert to m/s for flux calculation

WANNINKHOF_A = 0.251    # cm hr^-1 (m s^-1)^-2
SC_REF       = 660.0    # Schmidt number of CO2 in seawater at 20°C (reference)


# ===========================================================================
# CO2 SOLUBILITY — WEISS (1974)
# ===========================================================================

# Weiss, R.F. (1974). Carbon dioxide in water and seawater: the solubility
# of a non-ideal gas. Marine Chemistry, 2(3), 203–215.
# DOI: 10.1016/0304-4203(74)90015-2
#
# ln(K0) = A1 + A2*(100/T) + A3*ln(T/100) + S*(B1 + B2*(T/100) + B3*(T/100)^2)
# where T is temperature in Kelvin, S is salinity in PSU
# K0 is in mol / (L · atm); convert to mol / (m3 · Pa) for SI

WEISS_A1 = -58.0931
WEISS_A2 =  90.5069
WEISS_A3 =  22.2940
WEISS_B1 =   0.027766
WEISS_B2 =  -0.025888
WEISS_B3 =   0.0050578


# ===========================================================================
# UNIT CONVERSIONS
# ===========================================================================

# pCO2: the BGC hindcast spco2 is in Pa — convert to atm for Weiss K0 formula
PA_TO_ATM   = 1.0 / 101325.0

# pCO2: NOAA CO2 is in ppm = uatm (at 1 atm total pressure, 1 ppm ≈ 1 uatm)
# Convert uatm to atm:
UATM_TO_ATM = 1e-6

# k: Wanninkhof gives k in cm/hr — convert to m/s:
CMHR_TO_MS  = 1.0 / (100.0 * 3600.0)

# Flux F = k * K0 * delta_pCO2 will be in mol / m2 / s
# Convert to mol / m2 / yr:
S_TO_YR     = 365.25 * 24.0 * 3600.0

# Convert mol C to g C (for reporting):
MOL_C_TO_G  = 12.011    # g/mol

# Convert mol C to Pg C (for global integral reporting):
# 1 Pg = 1e15 g
MOL_C_TO_PG = MOL_C_TO_G / 1e15


# ===========================================================================
# SCHMIDT NUMBER POLYNOMIAL — WANNINKHOF (1992/2014) for CO2 in seawater
# ===========================================================================

# Sc(CO2) = A - B*T + C*T^2 - D*T^3 + E*T^4
# T in degrees Celsius
# Coefficients from Wanninkhof (2014), Table 1

SC_A = 2116.8
SC_B =  136.25
SC_C =    4.7353
SC_D =    0.092307
SC_E =    0.0007555


# ===========================================================================
# OCEAN AREA (for global integration)
# ===========================================================================

# Earth mean radius (m)
EARTH_RADIUS_M = 6.371e6

# Approximate total ocean area (m2) — used as a sanity check
TOTAL_OCEAN_AREA_M2 = 3.619e14


# ===========================================================================
# PLOTTING DEFAULTS
# ===========================================================================

FIGURE_DPI = 150
CMAP_FLUX   = "RdBu_r"   # diverging: red=outgassing, blue=uptake
CMAP_TREND  = "RdBu_r"
CMAP_PCO2   = "viridis"
