# Stage 1 — Surface Net Air-Sea CO₂ Flux

**Ocean Carbon Cycle · Net Flux & Layered Accumulation Project**

---

## Overview

Stage 1 estimates the **net air-sea CO₂ flux at the ocean surface** at every 0.25° grid point and every month across the full Copernicus BGC hindcast record (1993–2024), then:

1. Integrates it over the **global ocean surface** to produce a time series of net uptake in Pg C yr⁻¹.
2. Plots it against **atmospheric CO₂** (NOAA GML) to answer the central question: *is the ocean's carbon sink growing, flattening, or saturating under rising CO₂?*
3. Maps it spatially to identify **regional sources and sinks**.
4. Estimates **per-pixel trends** (Sen's slope) over the full record.
5. **Cross-validates** the reconstructed flux against the CMEMS MULTIOBS SOCAT-NN observation-based product, globally and per Fay (2014) biome domain.
6. Characterises **multi-year variability** via LOESS smoothing with piecewise trend detection.
7. Characterises **spectral properties** and cross-spectral coherence between reconstruction and reference.

The reconstruction uses exclusively **`GLOBAL_MULTIYEAR_BGC_001_029`** as the pCO₂ source — a free-running NEMO-PISCES biogeochemical hindcast with no BGC data assimilation. GLORYS12V1 supplies SST and SSS; ERA5 supplies wind speed. The MULTIOBS product (`MULTIOBS_GLO_BGC_CARBON_SURFACE_MYNRT_015_008`) is used strictly as validation — never as a primary input.

---

## Mathematical formulation

### Bulk air-sea CO₂ flux

The air-sea CO₂ flux per unit area is given by:

$$\boxed{F = k \cdot K_0 \cdot \left( pCO_2^{\,\mathrm{atm}} - pCO_2^{\,\mathrm{ocn}} \right)}$$

| Symbol | Definition | Units |
|--------|-----------|-------|
| $F$ | Net CO₂ flux | mol C m⁻² yr⁻¹ |
| $k$ | Gas transfer velocity | m s⁻¹ |
| $K_0$ | CO₂ solubility | mol m⁻³ atm⁻¹ |
| $pCO_2^{\mathrm{atm}}$ | Atmospheric pCO₂ | atm |
| $pCO_2^{\mathrm{ocn}}$ | Surface ocean pCO₂ | atm |

**Sign convention:** $F > 0$ = ocean uptake (sink); $F < 0$ = outgassing (source).

---

### Gas transfer velocity — Wanninkhof (2014)

The gas transfer velocity $k$ is parameterised as a quadratic function of the 10 m wind speed:

$$k = a \cdot \left(\langle u_{10} \rangle^2 + \sigma_u^2\right) \cdot \left(\frac{Sc}{660}\right)^{-1/2}$$

where:

| Symbol | Value / definition |
|--------|--------------------|
| $a$ | 0.251 cm hr⁻¹ (m s⁻¹)⁻² — Wanninkhof (2014) coefficient calibrated against global bomb-¹⁴C inventories |
| $\langle u_{10} \rangle$ | Monthly-mean 10 m scalar wind speed [m s⁻¹] from ERA5 |
| $\sigma_u^2$ | Sub-monthly wind speed variance [m² s⁻²] from ERA5 daily winds (optional — see §Wind variance correction) |
| $Sc$ | Schmidt number of CO₂ in seawater (see below) |
| 660 | Reference Schmidt number of CO₂ at 20°C in seawater |

The variance correction term $\sigma_u^2$ accounts for the fact that Wanninkhof's coefficient was calibrated against the full wind speed distribution, not monthly means. Applying the formula to $\langle u \rangle^2$ alone underestimates $k$ because $\langle u^2 \rangle = \langle u \rangle^2 + \sigma_u^2$. When the ERA5 daily wind file is not available, the correction term drops out and only $\langle u_{10} \rangle^2$ is used.

$k$ is computed in cm hr⁻¹ and converted to m s⁻¹:
$$k \,[\text{m s}^{-1}] = k \,[\text{cm hr}^{-1}] \times \frac{1}{100 \times 3600}$$

---

### Schmidt number — Wanninkhof (2014), Table 1

The Schmidt number of CO₂ in seawater as a function of sea surface temperature $T$ (in °C):

$$Sc = A - B \cdot T + C \cdot T^2 - D \cdot T^3 + E \cdot T^4$$

| $A$ | $B$ | $C$ | $D$ | $E$ |
|-----|-----|-----|-----|-----|
| 2116.8 | 136.25 | 4.7353 | 0.092307 | 0.0007555 |

The factor $(Sc/660)^{-1/2}$ normalises $k$ to the reference condition (CO₂ at 20°C in seawater).

---

### CO₂ solubility — Weiss (1974)

The solubility of CO₂ in seawater $K_0$ is computed using the empirical relation of Weiss (1974):

$$\ln K_0 = A_1 + A_2 \cdot \frac{100}{T} + A_3 \cdot \ln\!\left(\frac{T}{100}\right) + S \cdot \left[ B_1 + B_2 \cdot \frac{T}{100} + B_3 \cdot \left(\frac{T}{100}\right)^2 \right]$$

where $T$ is temperature in **Kelvin** and $S$ is salinity in **PSU**.

| $A_1$ | $A_2$ | $A_3$ | $B_1$ | $B_2$ | $B_3$ |
|-------|-------|-------|-------|-------|-------|
| −58.0931 | 90.5069 | 22.2940 | 0.027766 | −0.025888 | 0.0050578 |

Weiss (1974) gives $K_0$ in mol L⁻¹ atm⁻¹; we multiply by 1000 to convert to mol m⁻³ atm⁻¹.

---

### Global surface integral

The global net ocean CO₂ uptake at time $t$ is:

$$J_{\mathrm{net}}(t) = \sum_{i,j} F(t, i, j) \cdot A(i, j)$$

where $A(i,j) = R^2 \cos\varphi \,\Delta\varphi\,\Delta\lambda$ is the area of each 0.25° grid cell on the sphere ($R = 6.371 \times 10^6$ m, $\Delta\varphi = \Delta\lambda = 0.25°$). Positive $J_{\mathrm{net}}$ = ocean uptake.

The result is converted from mol C yr⁻¹ to Pg C yr⁻¹:

$$J_{\mathrm{net}} \,[\text{Pg C yr}^{-1}] = J_{\mathrm{net}} \,[\text{mol C yr}^{-1}] \times \frac{12.011 \,\text{g mol}^{-1}}{10^{15} \,\text{g Pg}^{-1}}$$

---

### Per-pixel trend estimation — Sen's slope

The trend at each pixel is estimated using the **Theil-Sen estimator** (Sen's slope):

$$\hat{\beta} = \mathrm{median}\!\left(\frac{y_j - y_i}{t_j - t_i}\right), \quad \forall \; i < j$$

This is the median of all pairwise slopes. It is preferred over OLS because it is robust to outliers (e.g., ENSO-driven anomalous years), makes no assumption of normality on the residuals, and is consistent under serial correlation. Units: mol C m⁻² yr⁻¹ per decade.

---

### Validation metrics

**Rolling RMSD and Pearson r** are computed on raw monthly global integrals within a 12-month rolling window:

$$\text{RMSD}(t) = \sqrt{\frac{1}{w}\sum_{i=t-w+1}^{t}(J_{\mathrm{rec},i} - J_{\mathrm{obs},i})^2}, \qquad r(t) = \text{Pearson}(J_{\mathrm{rec}}, J_{\mathrm{obs}})_{[t-w+1,t]}$$

**LOESS multi-year variability** uses locally weighted scatterplot smoothing (`frac = 0.15`, ≈18-month window) to suppress the seasonal cycle and reveal interannual variability. PELT breakpoint detection (`ruptures`, RBF cost, `pen = 3`) identifies regime changes in the reconstruction LOESS smooth; piecewise linear trends are fitted to each segment and reported in Pg C yr⁻² with arrows anchored to the trend line.

**Spectral analysis** uses Welch's method (`nperseg = n/3`) for PSD and `scipy.signal.coherence` with `nperseg = n/4` for the cross-power spectrum (ensuring multiple averaging segments and physically meaningful coherence). `nfft = 2^(ceil(log2(n))+3)` for smooth interpolation between Fourier bins.

---

## Data sources

### Surface ocean pCO₂ (primary driver)

| | |
|---|---|
| **Product** | CMEMS Global Ocean Biogeochemistry Hindcast |
| **Product ID** | `GLOBAL_MULTIYEAR_BGC_001_029` |
| **URL** | https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_BGC_001_029/description |
| **Variable** | `spco2` — surface partial pressure of CO₂ [Pa → converted to atm] |
| **Model** | NEMO-PISCES (Mercator Ocean International) |
| **Assimilation** | None — free-running hindcast (no BGC data assimilation) |
| **Record** | 1993-01 to 2026-04 |
| **Resolution** | 0.25°, 75 vertical levels, monthly means |
| **Physics forcing** | FREEGLORYS2V4 / ERA-Interim |

> ⚠️ This is a **hindcast**, not a reanalysis in the strict sense. The physical ocean state is observation-constrained (via FREEGLORYS2V4), but all biogeochemical fields — pCO₂, nutrients, pH — are produced by a free-running PISCES simulation with no assimilation of BGC observations.

### Physical fields — SST, SSS

| | |
|---|---|
| **Product** | CMEMS Global Ocean Physics Reanalysis (GLORYS12V1) |
| **Product ID** | `GLOBAL_MULTIYEAR_PHY_001_030` |
| **URL** | https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/description |
| **Variables** | `thetao` (potential temperature, °C), `so` (practical salinity, PSU) |
| **Assimilation** | Yes — along-track altimetry, satellite SST, sea-ice concentration, in-situ T/S profiles (reduced-order Kalman filter) |
| **Record** | 1993–present |
| **Resolution** | 1/12° → bilinearly regridded to 0.25° |
| **Physics forcing** | ERA5 |

> ⚠️ GLORYS12V1 and `GLOBAL_MULTIYEAR_BGC_001_029` are **not dynamically consistent** — they use different physics (ERA5 vs ERA-Interim) and different ocean simulations. The PISCES-internal T/S fields are not publicly distributed. This is a known limitation (see §Known limitations).

### Wind speed — monthly mean

| | |
|---|---|
| **Product** | ERA5 monthly averaged reanalysis on single levels |
| **Dataset ID** | `reanalysis-era5-single-levels-monthly-means` |
| **URL** | https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels-monthly-means |
| **Variables** | `u10`, `v10` → scalar speed $\sqrt{u_{10}^2 + v_{10}^2}$ [m s⁻¹] |
| **Record** | 1940–present |
| **Resolution** | 0.25° |
| **Access** | Copernicus CDS API (separate credentials from CMEMS — see §Credentials) |

### Wind speed — daily (optional, for variance correction)

| | |
|---|---|
| **Product** | ERA5 reanalysis on single levels |
| **Dataset ID** | `reanalysis-era5-single-levels` |
| **Variables** | `u10`, `v10` at 12:00 UTC (daily representative snapshot) |
| **Purpose** | Sub-monthly wind variance $\sigma_u^2$ for the Wanninkhof (2014) correction: $\langle u^2 \rangle = \langle u \rangle^2 + \sigma_u^2$ |
| **Size warning** | ~40–50 GB for the full 1993–2024 record. Downloaded year-by-year to avoid CDS timeout |

### Atmospheric CO₂

| | |
|---|---|
| **Source** | NOAA Global Monitoring Laboratory — Marine Boundary Layer Reference |
| **File** | `co2_mm_gl.csv` |
| **Direct download** | https://gml.noaa.gov/webdata/ccgg/trends/co2/co2_mm_gl.csv |
| **Product page** | https://gml.noaa.gov/ccgg/trends/gl_data.html |
| **Units** | ppm → converted to atm (1 ppm ≈ 1 µatm at 1 atm total pressure) |
| **Record** | 1979–present, monthly |

### Validation — observation-based flux (L4)

| | |
|---|---|
| **Product** | CMEMS MULTIOBS Global Ocean Surface Carbon L4 |
| **Product ID** | `MULTIOBS_GLO_BGC_CARBON_SURFACE_MYNRT_015_008` |
| **URL** | https://data.marine.copernicus.eu/product/MULTIOBS_GLO_BGC_CARBON_SURFACE_MYNRT_015_008/description |
| **Variables used** | `fgco2` (air-sea CO₂ flux) |
| **Method** | Ensemble neural network trained on SOCAT in-situ fCO₂ observations; gap-filled globally using satellite covariates (SST, Chl, MLD) |
| **Type** | L4 — observation-constrained, gridded |
| **Record** | 1985–present |
| **Resolution** | 0.25°, monthly |

> This product is used **only for validation** — never as a primary input. The key distinction vs the hindcast: MULTIOBS is constrained by tens of millions of SOCAT in-situ fCO₂ measurements; the hindcast pCO₂ is purely model-generated.

### Biome domains

| | |
|---|---|
| **Source** | Fay & McKinley (2014) time-varying biome classification |
| **File** | `Time_Varying_Biomes.nc` → interpolated to CMEMS 0.25° grid → `Time_Varying_Biomes.cmems.nc` |
| **Domains** | 17 biomes covering the global ocean |
| **Used in** | `04_validate.py` — per-domain time series and climatology validation |

---

## Software requirements

```bash
pip install copernicusmarine cdsapi xarray[io] dask numpy scipy pandas \
            matplotlib cartopy seaborn statsmodels ruptures cmocean requests netCDF4
```

| Package | Role | Optional? |
|---------|------|-----------|
| `copernicusmarine` | CMEMS data download | No |
| `cdsapi` | ERA5 data download via CDS | No |
| `xarray` + `dask` | Out-of-core NetCDF handling | No |
| `numpy`, `scipy`, `pandas` | Numerical computation | No |
| `matplotlib`, `cartopy` | Plotting and map projections | No |
| `seaborn` | Plot styling | No |
| `statsmodels` | LOWESS smoother (`lowess` in `04_validate.py`) | No |
| `ruptures` | PELT breakpoint detection in LOESS figure | Yes — degrades gracefully |
| `cmocean` | Ocean-specific colormaps for maps | Yes — falls back to matplotlib |
| `requests` | NOAA CO₂ CSV download | No |
| `netCDF4` | NetCDF backend for xarray | No |

Recommended Python: ≥ 3.11.

---

## Credentials

**CMEMS** (scripts `01`, `02`, `03`):
```bash
copernicusmarine login
# credentials stored in ~/.copernicusmarine/
```
Register at: https://data.marine.copernicus.eu/register

**CDS / ERA5** (script `01`, wind downloads):

Create `~/.cdsapirc`:
```
url: https://cds.climate.copernicus.eu/api
key: <YOUR_UID>:<YOUR_API_KEY>
```
Register at: https://cds.climate.copernicus.eu/user/register

---

## Running the pipeline

```bash
# 1. Download all data (--test for one year only as a sanity check)
python scripts/01_download_data.py [--test]

# 2. Preprocess: regrid, unit-convert, harmonise, merge → processed_surface.nc
python scripts/02_preprocess.py

# 3. Compute flux and global integral → flux_3d.nc, global_flux.nc
python scripts/03_compute_flux.py

# 4. Full validation suite → figures + validation_metrics.csv
python scripts/04_validate.py

# 5. Main result figures (trend maps, seasonal cycle, etc.)
python scripts/05_plot_results.py
```

**Before running `04_validate.py`**, generate the Fay biome mask on the CMEMS 0.25° grid. On a login node this will be killed by the memory manager — use the provided SLURM script:

```bash
# From the biomes directory, on a compute node:
sbatch run_interpolator_cmems.slurm
# Then move the output:
mv Time_Varying_Biomes.cmems.nc /path/to/stage1/data/
```

Or interactively if memory permits:
```bash
python interpolator.py cmems
mv Time_Varying_Biomes.cmems.nc /path/to/stage1/data/
```

**Expected runtimes** (full 1993–2024 record, global 0.25°, modern laptop):

| Step | Runtime |
|------|---------|
| Download (CMEMS + NOAA) | 30–90 min |
| Download (ERA5 monthly) | 20–40 min |
| Download (ERA5 daily, optional) | Several hours — 40–50 GB |
| Preprocess (without daily wind) | 10–20 min |
| Preprocess (with daily wind variance) | 1–3 hours — memory intensive |
| Flux computation | 15–30 min |
| Validation | 30–60 min |

> ⚠️ The ERA5 daily wind variance computation (`load_wind_variance()` in `02_preprocess.py`) loads a ~50 GB file and computes monthly variance per pixel. On systems with < 16 GB RAM this may swap heavily. Run overnight or on a compute node. A standalone pre-computation script is recommended for memory-constrained systems.

---

## Repository structure

```
stage1/
├── README.md
├── CHANGELOG.md
├── config.py                        ← central configuration
├── requirements.txt
├── data/                            ← not tracked by git
│   ├── bgc_hindcast_spco2.nc
│   ├── phy_reanalysis_sst_sss.nc
│   ├── multiobs_surface_carbon.nc
│   ├── era5_wind10m_monthly.nc
│   ├── era5_wind10m_daily.nc        ← optional (wind variance correction)
│   ├── co2_mm_gl.csv
│   ├── Time_Varying_Biomes.cmems.nc ← generated by interpolator.py cmems
│   ├── processed_surface.nc         ← output of 02_preprocess.py
│   └── flux_3d.nc                   ← output of 03_compute_flux.py
├── output/
│   ├── global_flux.nc
│   ├── validation_metrics.csv       ← global + per-domain Fay metrics
│   └── figures/
│       ├── fig_validation_ts.png
│       ├── fig_validation_loess.png
│       ├── fig_validation_map.png
│       ├── fig_fay_ts.png
│       ├── fig_fay_clim.png
│       ├── fig_spectra.png
│       └── fig_cps.png
└── scripts/
    ├── 01_download_data.py
    ├── 02_preprocess.py
    ├── 03_compute_flux.py
    ├── 04_validate.py
    └── 05_plot_results.py
```

### Script dependency chain

```
01_download_data.py
        │
        ▼
02_preprocess.py ────────────────────────────────────────────┐
        │                                                     │
        ▼                                                     │
03_compute_flux.py                                           │
        │                                                     │
        ├──────────────────────────┐                         │
        ▼                          ▼                         │
04_validate.py          05_plot_results.py ◄─────────────────┘
```

### Design principles

- **`config.py` is the single source of truth.** All product IDs, variable names, physical constants, file paths, and unit conversion factors live there. If a CMEMS product is updated and a variable name changes, fix it in one place.
- **Scripts are stateless and idempotent.** Each script reads files, does its work, and writes output. Re-running produces the same result; delete an output file to force recomputation.
- **No computation in `05_plot_results.py`.** It reads pre-computed NetCDF files only — figure aesthetics can be iterated without rerunning the physics.
- **Pure functions in `03_compute_flux.py`.** All physical parameterisations (`schmidt_number_co2`, `gas_transfer_velocity`, `co2_solubility_K0`, `compute_flux`) take `xr.DataArray` inputs and return `xr.DataArray` outputs — independently testable and swappable for alternative formulations.
- **Dask-backed lazy loading.** All `xr.open_dataset` calls use `chunks="auto"` so the full global 3D archive is never loaded into RAM at once.

---

## Output figures

| File | Description |
|------|-------------|
| `fig_validation_ts.png` | 2-row: global time series J_net [Pg C yr⁻¹] (top) + 12-month rolling RMSD and Pearson r (bottom) |
| `fig_validation_loess.png` | Raw monthly + LOESS smooth (frac=0.15) + PELT breakpoints + piecewise trend slopes [Pg C yr⁻²] |
| `fig_validation_map.png` | Spatial RMSD (top) and bias (bottom) maps — contourf, fixed colourbars (RMSD: 0–5, bias: ±5 mol C m⁻² yr⁻¹), Robinson projection |
| `fig_fay_ts.png` | 17-panel monthly time series per Fay (2014) domain — reconstruction (domain colour) vs MULTIOBS (rose red) |
| `fig_fay_clim.png` | 17-panel climatological seasonal cycle ±1σ per domain |
| `fig_spectra.png` | Log-log PSD (Welch) — reconstruction vs MULTIOBS; period ticks at 30d, 90d, 180d, 1yr, 5yr, 10yr |
| `fig_cps.png` | Cross-power spectrum: gain, phase [days], magnitude-squared coherence |

---

## Known limitations and caveats

1. **Physical consistency:** `GLOBAL_MULTIYEAR_BGC_001_029` was forced by FREEGLORYS2V4/ERA-Interim, while SST/SSS here come from GLORYS12V1 (ERA5). These are not dynamically consistent — they use different atmospheric forcing and different ocean simulations. The PISCES-internal T/S fields are not publicly distributed; obtain them from Mercator Ocean International if exact consistency is required.

2. **Wind variance correction:** Active only when `era5_wind10m_daily.nc` is present in `data/`. Without it, `load_wind_variance()` returns `None` gracefully and $k$ is computed from monthly-mean $\langle u \rangle^2$ only, which underestimates $k$ in high-variability regions (Southern Ocean, storm tracks) where sub-monthly wind bursts dominate.

3. **Sea-ice masking:** The current implementation uses the PISCES pCO₂ NaN mask as a proxy for ice-covered pixels, rather than an explicit sea-ice concentration field. Grid cells with partial ice cover may have their flux underestimated. An explicit open-water fraction from NSIDC or GLORYS12 sea-ice output would improve this.

4. **Hindcast, not reanalysis:** Because `GLOBAL_MULTIYEAR_BGC_001_029` assimilates no BGC observations, its surface pCO₂ reflects whatever PISCES produces — with known regional biases in the equatorial Pacific (too-weak upwelling signal) and Southern Ocean. This directly limits flux accuracy independently of the $k$ and $K_0$ formulations.

5. **Record length and pseudo-steady-state:** A 30-year window may not fully capture carbon cycle timescales of 100–1000 years. Results should be interpreted as a snapshot of the current flux regime, not a long-term equilibrium. Rolling sub-window trends (see `fig_validation_loess.png`) are provided as a sensitivity check.

6. **MULTIOBS sign convention:** The `fgco2` variable in `MULTIOBS_GLO_BGC_CARBON_SURFACE_MYNRT_015_008` sign convention has varied between product versions. The current code uses `fgco2_obs = ds_surf["fgco2_obs"]` with no sign flip (v1.2.0+). Verify against the product QUID before updating the MULTIOBS dataset.

---

## Extending the project

This codebase is structured to grow stage by stage:

1. Create a new `stage2/` or `stage3/` directory at the same level as `stage1/`.
2. Each stage has its own `config.py`, `scripts/`, `data/`, and `output/` — keeping stages self-contained.
3. Shared utility functions (e.g., `compute_grid_cell_area`, ocean basin masks) should eventually be factored into a common `utils/` package at the repo root.
4. Any new physical parameterisation should be added as a pure function in the relevant `03_compute_*.py` with full docstring including the reference equation and citation.

Planned subsequent stages:
- **Stage 2** — export flux at 100 m and 500 m (gravitational pump)
- **Stage 3** — layered total carbon content (0–100 m / 100–500 m / >500 m) and mass-balance consistency check

---

## Bibliography

**Gas transfer velocity:**
> Wanninkhof, R. (2014). Relationship between wind speed and gas exchange over the ocean revisited. *Limnology and Oceanography: Methods*, 12(6), 351–362. https://doi.org/10.4319/lom.2014.12.351

**CO₂ solubility:**
> Weiss, R. F. (1974). Carbon dioxide in water and seawater: the solubility of a non-ideal gas. *Marine Chemistry*, 2(3), 203–215. https://doi.org/10.1016/0304-4203(74)90015-2

**BGC hindcast model:**
> Aumont, O., Ethé, C., Tagliabue, A., Bopp, L., & Gehlen, M. (2015). PISCES-v2: an ocean biogeochemical model for carbon and ecosystem studies. *Geoscientific Model Development*, 8, 2465–2513. https://doi.org/10.5194/gmd-8-2465-2015

**Physical reanalysis (GLORYS12):**
> Lellouche, J.-M. et al. (2021). The Copernicus Global 1/12° Oceanic and Sea Ice GLORYS12 Reanalysis. *Frontiers in Earth Science*, 9, 698876. https://doi.org/10.3389/feart.2021.698876

**Validation product — SOCAT:**
> Bakker, D. C. E. et al. (2016). A multi-decade record of high-quality fCO₂ data in version 3 of the Surface Ocean CO₂ Atlas (SOCAT). *Earth System Science Data*, 8, 383–413. https://doi.org/10.5194/essd-8-383-2016

**Key motivating reference — physical injection pump:**
> Bellacicco, M., Marullo, S., Dall'Olmo, G., Iudicone, D., & Buongiorno Nardelli, B. (2025). The oceanic physical injection pump of organic carbon. *Nature Communications*, 16, 7100. https://doi.org/10.1038/s41467-025-62363-z

**Carbon pump conceptual framework:**
> Boyd, P. W., Claustre, H., Levy, M., Siegel, D. A., & Weber, T. (2019). Multi-faceted particle pumps drive carbon sequestration in the ocean. *Nature*, 568, 327–335. https://doi.org/10.1038/s41586-019-1098-2

**ERA5 reanalysis:**
> Hersbach, H. et al. (2020). The ERA5 global reanalysis. *Quarterly Journal of the Royal Meteorological Society*, 146(730), 1999–2049. https://doi.org/10.1002/qj.3803

**Biome domains:**
> Fay, A. R. & McKinley, G. A. (2014). Global open-ocean biomes: mean and temporal variability. *Earth System Science Data*, 6, 273–284. https://doi.org/10.5194/essd-6-273-2014

**Robust trend estimation:**
> Sen, P. K. (1968). Estimates of the regression coefficient based on Kendall's tau. *Journal of the American Statistical Association*, 63(324), 1379–1389. https://doi.org/10.1080/01621459.1968.10480934

---

## License and citation

**Data licenses:** All CMEMS products are distributed under the [Copernicus Marine Service licence](https://marine.copernicus.eu/user-corner/service-commitments-and-licence). ERA5 data is subject to the [Copernicus Climate Change Service licence](https://cds.climate.copernicus.eu/api/v2/terms/static/licence-to-use-copernicus-products.pdf). NOAA GML data is public domain.

If you use this code in a publication, please cite the underlying data products and parameterisations listed above, and acknowledge the Copernicus Marine Service:

> *"This study has been conducted using E.U. Copernicus Marine Service Information; https://doi.org/10.48670/moi-00019"*
