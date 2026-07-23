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

Changelog:
    v1.1.0 — Added ERA5 daily wind file path and variance correction constant;
              added MULTIOBS_PCO2_VAR entry for improved reconstruction pCO2 driver;
              updated FIGURE_DPI to 300.
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
BGC_PCO2_VAR         = "spco2"      # surface partial pressure of CO2 (Pa)

# --- Physical reanalysis (SST, SSS — needed for K0 and k) ---
PHY_REANALYSIS_ID      = "GLOBAL_MULTIYEAR_PHY_001_030"
PHY_REANALYSIS_DATASET = "cmems_mod_glo_phy_my_0.083deg_P1M-m"  # monthly means
PHY_TEMP_VAR           = "thetao"   # potential temperature (deg C)
PHY_SAL_VAR            = "so"       # practical salinity (PSU)

# --- Cell thickness (from BGC hindcast, needed for volume integrals later) ---
THICKNESS_VAR = "thkcello"

# --- Validation & improved reconstruction: observation-based surface carbon L4 ---
# Used as:
#   (a) validation product: fgco2 (pre-computed flux, sign-flipped in 04_validate.py)
#   (b) improved pCO2 driver: spco2 (replaces PISCES spco2 in the improved reconstruction)
MULTIOBS_SURFACE_ID      = "MULTIOBS_GLO_BGC_CARBON_SURFACE_MYNRT_015_008"
MULTIOBS_SURFACE_DATASET = "dataset-carbon-rep-monthly"   # confirm exact name on CMEMS portal
MULTIOBS_FLUX_VAR        = "fgco2"    # air-sea CO2 flux (mol C m⁻² s⁻¹ or kg C m⁻² s⁻¹)
MULTIOBS_PCO2_VAR        = "spco2"    # surface pCO2 from SOCAT-NN [Pa or uatm — check attrs]

# --- Spatial extent (global ocean) ---
# Copernicus toolbox uses (lon_min, lon_max, lat_min, lat_max)
BBOX = (-180.0, 180.0, -90.0, 90.0)

# --- Depth: for Stage 1 we only need the surface level ---
DEPTH_MIN = 0.0   # m
DEPTH_MAX = 1.0   # m — take just the first model layer


# ===========================================================================
# ATMOSPHERIC CO2 — NOAA GML MARINE BOUNDARY LAYER REFERENCE
# ===========================================================================

NOAA_CO2_FILE    = DATA_DIR / "co2_mm_gl.csv"
NOAA_CO2_COLNAME = "average"    # column holding the ppm value


# ===========================================================================
# ERA5 WIND SPEED (for gas-transfer velocity k)
# ===========================================================================

# Monthly-mean wind file (used by v1.0.x reconstruction)
# Source: CDS — reanalysis-era5-single-levels-monthly-means
# Variables: u10, v10 [m/s]
ERA5_WIND_FILE = DATA_DIR / "era5_wind10m_monthly.nc"
ERA5_U10_VAR   = "u10"
ERA5_V10_VAR   = "v10"

# Daily wind file (used for sub-monthly variance correction in v1.1.0)
# Source: CDS — reanalysis-era5-single-levels  (daily means or 6-hourly)
# Same variables u10, v10 — aggregated to monthly mean and variance in 02_preprocess.py
# Download note: request per year to keep file sizes manageable.
ERA5_WIND_DAILY_FILE = DATA_DIR / "era5_wind10m_daily.nc"


# ===========================================================================
# GAS TRANSFER VELOCITY — WANNINKHOF (2014)
# ===========================================================================

# k = a · (u² + σ²_u) · (Sc/660)^(-0.5)
# where σ²_u is the sub-monthly wind speed variance (variance correction, v1.1.0)
#
# Wanninkhof, R. (2014). Relationship between wind speed and gas exchange
# over the ocean revisited. Limnology and Oceanography: Methods, 12(6), 351–362.
# DOI: 10.4319/lom.2014.12.351

WANNINKHOF_A = 0.251    # cm hr^-1 (m s^-1)^-2
SC_REF       = 660.0    # reference Schmidt number (CO2 at 20°C in seawater)


# ===========================================================================
# CO2 SOLUBILITY — WEISS (1974)
# ===========================================================================

# ln(K0) = A1 + A2*(100/T) + A3*ln(T/100) + S*(B1 + B2*(T/100) + B3*(T/100)^2)
# T in Kelvin, S in PSU; K0 in mol/(L·atm) → ×1000 → mol/(m³·atm)
#
# Weiss, R.F. (1974). Marine Chemistry, 2(3), 203–215.
# DOI: 10.1016/0304-4203(74)90015-2

WEISS_A1 = -58.0931
WEISS_A2 =  90.5069
WEISS_A3 =  22.2940
WEISS_B1 =   0.027766
WEISS_B2 =  -0.025888
WEISS_B3 =   0.0050578


# ===========================================================================
# UNIT CONVERSIONS
# ===========================================================================

PA_TO_ATM   = 1.0 / 101325.0       # Pa → atm
UATM_TO_ATM = 1e-6                  # ppm (≈µatm) → atm
CMHR_TO_MS  = 1.0 / (100.0 * 3600.0)   # cm/hr → m/s
S_TO_YR     = 365.25 * 24.0 * 3600.0   # s → yr
MOL_C_TO_G  = 12.011                # g/mol C
MOL_C_TO_PG = MOL_C_TO_G / 1e15    # mol C → Pg C


# ===========================================================================
# SCHMIDT NUMBER POLYNOMIAL — WANNINKHOF (2014), Table 1
# ===========================================================================

# Sc(CO2) = A - B*T + C*T^2 - D*T^3 + E*T^4,  T in °C

SC_A = 2116.8
SC_B =  136.25
SC_C =    4.7353
SC_D =    0.092307
SC_E =    0.0007555


# ===========================================================================
# OCEAN AREA (for global integration)
# ===========================================================================

EARTH_RADIUS_M      = 6.371e6    # m
TOTAL_OCEAN_AREA_M2 = 3.619e14   # m² (sanity check only)


# ===========================================================================
# PLOTTING DEFAULTS
# ===========================================================================

FIGURE_DPI  = 300               # increased from 150 in v1.1.0
CMAP_FLUX   = "RdBu_r"         # diverging: red=outgassing, blue=uptake
CMAP_TREND  = "RdBu_r"
CMAP_PCO2   = "viridis"

# Fixed colorbar extents for validation comparison maps (v1.1.0)
# Ensures both old and improved reconstruction maps share the same scale
VAL_RMSD_MAX  = 5.0    # mol C m⁻² yr⁻¹ — upper bound for RMSD colourbar
VAL_BIAS_MIN  = -5.0   # mol C m⁻² yr⁻¹ — lower bound for bias colourbar
VAL_BIAS_MAX  =  5.0   # mol C m⁻² yr⁻¹ — upper bound for bias colourbar
