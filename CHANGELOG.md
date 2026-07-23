# Changelog

All notable changes to this project are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
Versioning follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html):
`MAJOR.MINOR.PATCH` — major = breaking change, minor = new feature, patch = bug fix.

---

## [1.0.0] — 2026-07-23

### Summary

First complete implementation of **Stage 1: Surface Net Air-Sea CO₂ Flux**.

This version establishes the full pipeline from raw Copernicus Marine Service
data downloads through to the key result figures, including cross-validation
against an independent observation-based product. It constitutes the primary
deliverable discussed in the project meeting of 22 July 2026 and implements
exactly the analysis described as "Figure 1" in that meeting: a time series
of the global net ocean CO₂ uptake plotted alongside atmospheric CO₂, covering
the full available reanalysis record (1993–2026).

---

### Added

#### `config.py`
- Central configuration module containing all project-wide constants, file
  paths, CMEMS product/dataset/variable identifiers, physical constants
  (Wanninkhof 2014 coefficient `a`, Weiss 1974 solubility coefficients,
  Schmidt number polynomial coefficients), unit conversion factors, and
  plotting defaults.
- Single source of truth — all other scripts import from here.

#### `scripts/01_download_data.py`
- Automated download of all raw data required for Stage 1:
  - **Surface ocean pCO₂** (`spco2`) from the CMEMS Global Ocean Biogeochemistry
    Hindcast (`GLOBAL_MULTIYEAR_BGC_001_029`) — 1993–2026, 0.25°, monthly.
  - **Sea surface temperature and salinity** (`thetao`, `so`) from the CMEMS
    GLORYS12V1 Physical Reanalysis (`GLOBAL_MULTIYEAR_PHY_001_030`) — 1993–present,
    1/12°, monthly.
  - **Observation-based air-sea CO₂ flux** (`fgco2`, `spco2`) from the CMEMS
    MULTIOBS surface carbon L4 product (`MULTIOBS_GLO_BGC_CARBON_SURFACE_MYNRT_015_008`)
    — SOCAT neural-network product, used exclusively for cross-validation.
  - **Global mean atmospheric CO₂** (ppm) from NOAA GML
    (`co2_mm_gl.csv`) — downloaded directly over HTTP.
  - **ERA5 10 m wind components** (`u10`, `v10`) from the Copernicus CDS
    (`reanalysis-era5-single-levels-monthly-means`) — provided as a commented
    template requiring separate CDS API credentials.
- `--test` flag for single-year (2010) download to verify credentials and
  connectivity before committing to the full 30-year download.
- Idempotent: skips any file that already exists on disk.

#### `scripts/02_preprocess.py`
- Loads all raw downloaded files into `xarray` Datasets with Dask-backed
  lazy evaluation (`chunks="auto"`).
- Regrids GLORYS12 SST/SSS (1/12°) to the BGC hindcast's 0.25° grid via
  bilinear interpolation (`xr.DataArray.interp`).
- Regrids ERA5 wind speed to 0.25° on the same grid.
- Harmonizes time axes: converts CMEMS cftime objects to numpy datetime64
  to allow `xr.merge` across products.
- Unit conversions:
  - `spco2`: Pa → atm (`× 1/101325`)
  - `spco2_atm` (NOAA): ppm → atm (`× 1e-6`)
  - `fgco2_obs` (MULTIOBS): auto-detects kg C m⁻² s⁻¹ or mol C m⁻² s⁻¹
    and converts to mol C m⁻² yr⁻¹
- Derives ocean/land/ice mask from NaN pattern of the PISCES pCO₂ field.
- Broadcasts 1-D atmospheric CO₂ time series to full (time, lat, lon) grid
  for pixel-wise arithmetic.
- Saves merged output to `data/processed_surface.nc`.

#### `scripts/03_compute_flux.py`
- Implements the full bulk air-sea CO₂ flux parameterisation
  $F = k \cdot K_0 \cdot \Delta pCO_2$:

  - **`schmidt_number_co2(sst)`** — Wanninkhof (2014) Table 1 polynomial,
    4th-order in SST (°C).
  - **`gas_transfer_velocity(wind_speed, Sc)`** — Wanninkhof (2014) quadratic
    wind parameterisation, $k = a \cdot u_{10}^2 \cdot (Sc/660)^{-0.5}$,
    converted cm hr⁻¹ → m s⁻¹.
  - **`co2_solubility_K0(sst, sss)`** — Weiss (1974) natural-log formulation
    in T (K) and S (PSU), converted mol L⁻¹ atm⁻¹ → mol m⁻³ atm⁻¹ (×1000).
  - **`compute_flux(k, K0, spco2_atm, spco2_ocean, ocean_mask)`** — assembles
    the bulk formula; converts mol m⁻² s⁻¹ → mol m⁻² yr⁻¹; applies ocean mask.
  - **`compute_grid_cell_area(lat, lon)`** — spherical geometry, area element
    $R^2 \cos\varphi \, d\varphi \, d\lambda$, returns (lat, lon) DataArray in m².
  - **`global_integral(flux, cell_area, ocean_mask)`** — spatially integrates
    the flux over the global ocean at each time step; returns both raw
    mol C yr⁻¹ and sign-flipped Pg C yr⁻¹ (`J_net`, positive = uptake).
- Saves full (time, lat, lon) flux field to `data/flux_3d.nc`.
- Saves global integral time series to `output/global_flux.nc`.
- All physical functions are pure (no I/O, no side effects) and
  accept/return `xr.DataArray` — independently testable and replaceable.

#### `scripts/04_validate.py`
- Cross-validates the reconstructed bulk flux against the CMEMS MULTIOBS
  observation-based flux (SOCAT neural network) over their period of overlap.
- Metrics computed:
  - Pixel-wise RMSD and bias (mean difference, reconstructed − MULTIOBS)
  - Global scalar Pearson correlation coefficient
  - Standard deviation ratio
- Computes MULTIOBS global integral for direct time series comparison.
- Produces:
  - `output/validation_metrics.csv` — scalar statistics
  - `fig_validation_ts.png` — time series of J_net: reconstructed vs MULTIOBS
  - `fig_validation_map.png` — two-panel map of RMSD and bias
  - `fig_taylor.png` — Taylor diagram

#### `scripts/05_plot_results.py`
- Produces all main result figures (no computation — reads pre-computed files):
  - **`fig01_flux_vs_co2.png`** — the primary figure from the project meeting:
    global net ocean CO₂ uptake J_net [Pg C yr⁻¹] on the left axis vs.
    atmospheric CO₂ [ppm] on the right axis, both plotted over 1993–2026.
    Monthly values shown in pale fill; annual means as solid lines.
  - **`fig02_annual_flux_map.png`** — time-mean air-sea CO₂ flux map (Robinson
    projection). Red = outgassing, blue = uptake.
  - **`fig03_flux_trend_map.png`** — per-pixel Sen's slope (Theil-Sen estimator)
    of annual mean flux, in mol C m⁻² yr⁻¹ per decade. Identifies where the
    sink is intensifying or weakening.
  - **`fig04_delta_pco2_map.png`** — time-mean ΔpCO₂ = pCO₂(ocean) − pCO₂(atm)
    in µatm. The thermodynamic driver map.
  - **`fig05_monthly_seasonal_cycle.png`** — climatological monthly mean of
    J_net and atmospheric CO₂, showing the seasonal signal in the global flux.
- Cartopy used for map projections when available; falls back to plain imshow.

#### `README.md`
- Full scientific documentation including:
  - Mathematical formulation of the bulk flux equation with all sub-equations
    (gas transfer velocity, Schmidt number, solubility K0, global integral,
    Sen's slope)
  - Complete data source table with product IDs, URLs, temporal coverage,
    resolution, and model details
  - Software requirements and installation instructions
  - Credentials setup for CMEMS and CDS
  - Step-by-step run instructions
  - Output file inventory
  - Full bibliography (10 references)
  - Known limitations and caveats (5 items)
  - Code structure and design principles
  - Git setup recommendations including `.gitignore` template
  - License and citation guidance

#### `requirements.txt`
- Pinned minimum versions for all Python dependencies.

---

### Physical constants and parameterisations (v1.0.0)

| Quantity | Formula / value | Reference |
|---|---|---|
| Gas transfer velocity coefficient | $a = 0.251$ cm hr⁻¹ (m s⁻¹)⁻² | Wanninkhof (2014) |
| Schmidt number reference | $Sc_{ref} = 660$ (CO₂ at 20°C in seawater) | Wanninkhof (2014) |
| Schmidt number polynomial | 4th-order in SST (°C), 5 coefficients | Wanninkhof (2014) Table 1 |
| Solubility K₀ | Natural-log formula in T (K) and S (PSU) | Weiss (1974) |
| K₀ unit conversion | mol L⁻¹ atm⁻¹ × 1000 → mol m⁻³ atm⁻¹ | — |
| Earth radius | 6.371 × 10⁶ m | — |
| C molar mass | 12.011 g mol⁻¹ | — |

---

### Known issues and planned improvements (tracked for v1.1.0)

- [ ] **Sub-monthly wind variance correction** — the Wanninkhof (2014) coefficient
  was calibrated against global bomb-¹⁴C inventories and implicitly accounts
  for the variance of the wind speed distribution. Applying it to monthly-mean
  winds underestimates $k$. A correction term $\sigma_u^2$ (sub-monthly wind
  variance from ERA5 6-hourly data) should be added.
- [ ] **Sea-ice open-water fraction** — the current implementation uses the PISCES
  NaN mask as a proxy for ice cover. An explicit sea-ice concentration field
  (e.g., from GLORYS12 or NSIDC) should be used to scale the flux by the
  open-water fraction in partially ice-covered pixels.
- [ ] **Physical consistency** — SST/SSS come from GLORYS12 (ERA5-forced, assimilated)
  while surface pCO₂ comes from the BGC hindcast (ERA-Interim-forced, free-running).
  These are not dynamically consistent. Resolution pending clarification from
  Mercator Ocean on availability of the PISCES-internal T/S fields.
- [ ] **Sen's slope loop performance** — the per-pixel Theil-Sen estimator in
  `05_plot_results.py` currently runs in a Python loop over the lat/lon grid.
  Should be replaced with `xr.apply_ufunc` + `dask` parallelism for acceptable
  runtime on the full global domain.
- [ ] **MULTIOBS dataset name** — the exact dataset string for
  `MULTIOBS_GLO_BGC_CARBON_SURFACE_MYNRT_015_008` in `config.py` is a
  placeholder (`"dataset-carbon-rep-monthly"`) and must be confirmed against
  the current CMEMS catalogue before running.

---

[1.0.0]: https://github.com/[your-repo]/releases/tag/v1.0.0

---

## [1.1.0] — 2026-07-23

### Summary

Introduced a second, **improved flux reconstruction** alongside the original,
addressing the two dominant error sources identified during validation against
the CMEMS MULTIOBS product (r = 0.445, normalised σ = 0.22 in v1.0.0):

1. **pCO₂ driver substitution** — replaces the free-running PISCES model pCO₂
   with the observation-constrained MULTIOBS SOCAT-NN surface pCO₂ field.
   The PISCES pCO₂ lacked sufficient regional variance (dominant cause of the
   low normalised standard deviation), whereas SOCAT covers tens of millions of
   in-situ fCO₂ measurements and correctly resolves regional source/sink
   contrasts (equatorial Pacific outgassing, Southern Ocean and North Atlantic
   uptake).

2. **Wind variance correction** — adds the sub-monthly wind speed variance σ²_u
   (computed from ERA5 daily winds) to the monthly-mean u² in the Wanninkhof
   (2014) gas transfer velocity: k ∝ (⟨u⟩² + σ²_u) instead of ⟨u⟩² alone.
   Monthly averaging suppresses wind speed variance, causing k — and therefore
   the flux — to be underestimated, particularly in the Southern Ocean and
   storm-track regions where sub-monthly wind bursts dominate.

Both reconstructions are preserved in all output files and all validation
figures. The original (v1.0.x) reconstruction is not removed or replaced.

---

### Added

#### `config.py`
- `ERA5_WIND_DAILY_FILE` — path for ERA5 daily wind download (new CDS request).
- `VAL_RMSD_MAX`, `VAL_BIAS_MIN`, `VAL_BIAS_MAX` — fixed colorbar extents for
  the 2×2 validation maps (0–5 and −5–+5 mol C m⁻² yr⁻¹ respectively).
- `FIGURE_DPI` raised from 150 to **300**.

#### `scripts/02_preprocess.py`
- `load_spco2_ocean_obs()` — loads MULTIOBS `spco2` (SOCAT-NN surface pCO₂),
  auto-detects units (Pa or µatm), converts to atm, regrids to the BGC 0.25°
  grid, and saves as `spco2_ocean_obs` in `processed_surface.nc`.
- `load_wind_variance()` — loads ERA5 daily 10 m winds, computes the
  sub-monthly scalar wind speed variance σ²_u per month per pixel via
  `resample(time="1ME").var()`, regrids to 0.25°, saves as `wind_variance`
  in `processed_surface.nc`. Returns None gracefully if the daily wind file
  has not been downloaded.
- Both new variables are optional — if their source files are absent, the
  preprocessing completes with the original variable set and downstream scripts
  fall back to the original reconstruction only.

#### `scripts/03_compute_flux.py`
- `gas_transfer_velocity()` extended — accepts an optional `wind_variance`
  argument; when provided, computes k ∝ (⟨u⟩² + σ²_u) instead of ⟨u⟩².
- `compute_flux()` made generic via a `label` parameter so the same function
  produces both `fgco2` and `fgco2_improved` without duplication.
- `main()` now computes and saves both reconstructions:
  - `data/flux_3d.nc`: gains `fgco2_improved` and `k_improved` fields.
  - `output/global_flux.nc`: gains `J_net_improved_PgC` and its annual resample.
  - Gracefully skips the improved reconstruction if `spco2_ocean_obs` is absent.

#### `scripts/04_validate.py`  (full rework of all three figures)
- **Time series** (`fig_validation_ts.png`) — extended from 2 to **3 curves**:
  original (blue), improved (green), MULTIOBS reference (red dashed).
- **Validation maps** (`fig_validation_map.png`) — reworked from a 2×1 layout
  (RMSD + bias for original only) to a **2×2 contourf layout**:
  - Left column: original reconstruction RMSD (top) and bias (bottom).
  - Right column: improved reconstruction RMSD (top) and bias (bottom).
  - Both RMSD panels share a fixed colourbar 0–5 mol C m⁻² yr⁻¹ (Reds).
  - Both bias panels share a fixed colourbar −5–+5 mol C m⁻² yr⁻¹ (RdBu_r).
  - All panels rendered with `contourf` (21 levels) for publication quality.
- **Taylor diagram** (`fig_taylor.png`) — extended from 1 to **2 model markers**:
  ● original (blue circle) and ▲ improved (green triangle).
- `validation_metrics.csv` — now contains metrics for both reconstructions
  (prefixed `orig_` and `imp_`).
- All figures saved at **300 dpi**.

### Changed

- `config.py`: `FIGURE_DPI` 150 → 300.
- `05_plot_results.py`: no code changes — DPI increase is inherited from config.

### Known issues carried forward

- Sub-monthly wind variance correction is inactive until the ERA5 daily wind
  file (`data/era5_wind10m_daily.nc`) is downloaded via the CDS API. The
  improved reconstruction then applies only the pCO₂ driver substitution.
- MULTIOBS `spco2` NRT record ends before the BGC hindcast (Apr 2026) — the
  improved reconstruction time series will be shorter than the original.
- Physical consistency caveat (GLORYS12 vs. FREEGLORYS2V4) unchanged.

---

[1.1.0]: https://github.com/[your-repo]/releases/tag/v1.1.0
[1.0.0]: https://github.com/[your-repo]/releases/tag/v1.0.0
