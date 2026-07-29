# Stage 1 — Surface Net Air-Sea CO₂ Flux

**Ocean Carbon Cycle · Net Flux & Layered Accumulation Project**

---

## Overview

Stage 1 estimates the **net air-sea CO₂ flux at the ocean surface** at every 0.25° grid point and every month across the full Copernicus BGC hindcast record (1993–2024), then:

1. Integrates it over the **global ocean surface** to produce a time series of net uptake in Pg C yr⁻¹.
2. Plots it against **atmospheric CO₂** (NOAA GML) to answer the central question: *is the ocean's carbon sink growing, flattening, or saturating under rising CO₂?*
3. Maps it spatially to identify **regional sources and sinks**.
4. Estimates **per-pixel trends** (Sen's slope) with **Mann-Kendall significance testing** (α = 0.05) over the full record.
5. Characterises **multi-year variability** via LOESS smoothing with PELT piecewise trend detection.
6. Analyses **ocean sink saturation** by regressing annual J_net against annual atmospheric CO₂.
7. Disaggregates flux by **Fay (2014) biogeochemical domain** — timeseries, seasonal climatologies, and trend attribution per province.
8. **Cross-validates** the reconstructed flux against the CMEMS MULTIOBS SOCAT-NN observation-based product, globally and per Fay domain.
9. Characterises **spectral properties** and cross-spectral coherence between reconstruction and reference.

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

$$k = a \cdot \left(\langle u_{10} \rangle^2 + \sigma_u^2\right) \cdot \left(\frac{Sc}{660}\right)^{-1/2}$$

| Symbol | Value / definition |
|--------|--------------------|
| $a$ | 0.251 cm hr⁻¹ (m s⁻¹)⁻² |
| $\langle u_{10} \rangle$ | Monthly-mean 10 m scalar wind speed [m s⁻¹] from ERA5 |
| $\sigma_u^2$ | Sub-monthly wind speed variance [m² s⁻²] from ERA5 daily winds (optional) |
| $Sc$ | Schmidt number of CO₂ in seawater |
| 660 | Reference Schmidt number of CO₂ at 20°C |

---

### Schmidt number — Wanninkhof (2014), Table 1

$$Sc = A - B \cdot T + C \cdot T^2 - D \cdot T^3 + E \cdot T^4 \qquad (T \text{ in °C})$$

| $A$ | $B$ | $C$ | $D$ | $E$ |
|-----|-----|-----|-----|-----|
| 2116.8 | 136.25 | 4.7353 | 0.092307 | 0.0007555 |

---

### CO₂ solubility — Weiss (1974)

$$\ln K_0 = A_1 + A_2 \cdot \frac{100}{T} + A_3 \cdot \ln\!\left(\frac{T}{100}\right) + S \cdot \left[ B_1 + B_2 \cdot \frac{T}{100} + B_3 \cdot \left(\frac{T}{100}\right)^2 \right]$$

$T$ in Kelvin, $S$ in PSU. Converted mol L⁻¹ atm⁻¹ → mol m⁻³ atm⁻¹ (×1000).

| $A_1$ | $A_2$ | $A_3$ | $B_1$ | $B_2$ | $B_3$ |
|-------|-------|-------|-------|-------|-------|
| −58.0931 | 90.5069 | 22.2940 | 0.027766 | −0.025888 | 0.0050578 |

---

### Global surface integral

$$J_{\mathrm{net}}(t) = \sum_{i,j} F(t, i, j) \cdot A(i, j)$$

where $A(i,j) = R^2 \cos\varphi \,\Delta\varphi\,\Delta\lambda$ ($R = 6.371 \times 10^6$ m, $\Delta\varphi = \Delta\lambda = 0.25°$).

$$J_{\mathrm{net}} \,[\text{Pg C yr}^{-1}] = J_{\mathrm{net}} \,[\text{mol C yr}^{-1}] \times \frac{12.011}{10^{15}}$$

---

### Per-pixel trend estimation — Sen's slope + Mann-Kendall test

The trend at each pixel is estimated using the **Theil-Sen estimator**:

$$\hat{\beta} = \mathrm{median}\!\left(\frac{y_j - y_i}{t_j - t_i}\right), \quad \forall \; i < j$$

Statistical significance is assessed with the **Mann-Kendall test** (two-sided, α = 0.05). The MK test statistic $S$ is computed via `scipy.stats.kendalltau` against a monotone integer index — mathematically equivalent to the standard MK test. Pixels with p > 0.05 are marked with hatching in `fig03`.

Working on annual means (rather than monthly) reduces serial autocorrelation and avoids inflating significance through short-lag correlation.

---

### LOESS multi-year variability

LOESS smoothing (`frac = 0.15`, ≈ 18-month window, 3 robustness iterations via `statsmodels.lowess`) suppresses the dominant seasonal cycle and reveals interannual variability without harmonic assumptions.

PELT breakpoint detection (`ruptures.Pelt`, RBF cost, `pen = 3`) identifies regime changes in the LOESS smooth. Piecewise linear segments are fitted by OLS; slopes annotated in Pg C yr⁻².

---

### Ocean sink saturation analysis

Annual J_net is regressed against annual atmospheric CO₂ (OLS):

$$J_{\mathrm{net}}^{\mathrm{annual}} = \beta_0 + \beta_1 \cdot [\mathrm{CO}_2]_{\mathrm{atm}}^{\mathrm{annual}} + \varepsilon$$

The slope $\beta_1$ [Pg C yr⁻¹ ppm⁻¹] quantifies ocean uptake sensitivity to atmospheric CO₂:
- $\beta_1 > 0$ significant → sink growing with rising CO₂ (no saturation)
- $\beta_1 \approx 0$ or not significant → sink stagnant (saturation or decoupling)
- Points coloured by year expose temporal drift in the relationship

---

### Climate index regression — Bellacicco et al. (2025)

Annual J_net is regressed against the Oceanic Niño Index (ONI) and Southern
Annular Mode (SAM) index via Spearman rank correlation (separate and multilinear):

**Separate correlations:**

$$\rho_\mathrm{ONI} = \mathrm{Spearman}(J_\mathrm{net}^\mathrm{annual},\, \mathrm{ONI}^\mathrm{annual})$$
$$\rho_\mathrm{SAM} = \mathrm{Spearman}(J_\mathrm{net}^\mathrm{annual},\, \mathrm{SAM}^\mathrm{annual})$$

**Multilinear combination (OLS):**

$$\hat{J} = \alpha \cdot \tilde{\mathrm{ONI}} + \beta \cdot \tilde{\mathrm{SAM}} + \gamma$$

where $\tilde{\cdot}$ denotes z-score normalisation. Spearman $\rho$ is then
computed between $J_\mathrm{net}$ and $\hat{J}$.

Bellacicco et al. (2025) found for PIP export: $\rho_\mathrm{ONI} = 0.57$,
$\rho_\mathrm{ONI+SAM} = 0.617$ ($p = 0.001$). The analysis is repeated here
for surface J_net globally and per Fay (2014) biogeochemical province to test
whether the ENSO/SAM teleconnection manifests already at the air-sea interface.

### Fay domain integration

Flux is area-integrated over each Fay (2014) biogeochemical province:

$$J_{\mathrm{domain}}^{(d)}(t) = \sum_{i,j \in d} F(t,i,j) \cdot A(i,j) \times \frac{12.011}{10^{15}} \quad [\text{Pg C yr}^{-1}]$$

Sen's slope and Mann-Kendall p-value computed on annual domain integrals allow attribution of the global trend to specific provinces.

---

### Validation metrics

**Rolling RMSD and Pearson r** on raw monthly global integrals, 12-month window.

**LOESS multi-year variability** + PELT breakpoints as above.

**Spectral analysis** via Welch's method (`nperseg = n/3`) for PSD; `scipy.signal.coherence` with `nperseg = n/4` for cross-power spectrum.

---

## Data sources

### Surface ocean pCO₂ (primary driver)

| | |
|---|---|
| **Product** | CMEMS Global Ocean Biogeochemistry Hindcast |
| **ID** | `GLOBAL_MULTIYEAR_BGC_001_029` |
| **Variable** | `spco2` [Pa] → converted to [atm] |
| **Depth** | Surface level only |
| **Resolution** | 0.25°, monthly, 1993–2024 |

### Physical driver fields (SST, SSS)

| | |
|---|---|
| **Product** | CMEMS Global Physical Reanalysis (GLORYS12V1) |
| **ID** | `GLOBAL_MULTIYEAR_PHY_001_030` |
| **Variables** | `thetao` [°C], `so` [PSU] |
| **Resolution** | 1/12°, monthly → regridded to 0.25° by bilinear interpolation |

### Wind speed

| | |
|---|---|
| **Product** | ERA5 monthly + daily reanalysis |
| **Variables** | `u10`, `v10` → $\|u\| = \sqrt{u_{10}^2 + v_{10}^2}$ [m s⁻¹] |
| **Variance correction** | σ²_u from daily data (optional; see `01_download_data.py`) |

### Atmospheric CO₂

| | |
|---|---|
| **Source** | NOAA Global Monitoring Laboratory |
| **File** | `co2_mm_gl.csv` |
| **Variable** | `average` [ppm] → converted to [atm] (×10⁻⁶) |

### Validation product

| | |
|---|---|
| **Product** | CMEMS MULTIOBS SOCAT Neural Network L4 |
| **ID** | `MULTIOBS_GLO_BGC_CARBON_SURFACE_MYNRT_015_008` |
| **Variables** | `fgco2` (flux), `spco2` (surface pCO₂) |
| **Use** | Cross-validation only — never as primary input |

### Biogeochemical provinces

| | |
|---|---|
| **Source** | Fay & McKinley (2014) time-varying biome mask |
| **File** | `data/Time_Varying_Biomes.cmems.nc` |
| **Variable** | `MeanBiomes` (17 open-ocean provinces) |

### Climate indices

| | |
|---|---|
| **ONI** | Oceanic Niño Index — NOAA PSL (ERSSTv5 Niño-3.4 3-month running mean) |
| **URL** | https://psl.noaa.gov/data/correlation/oni.data |
| **File** | `data/oni_monthly.txt` (auto-downloaded on first run) |

| | |
|---|---|
| **SAM** | Southern Annular Mode — Marshall (2003) station-based, NOAA CPC update |
| **URL** | https://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/aao/monthly.aao.index.b79.current.ascii.table |
| **File** | `data/sam_monthly.txt` (auto-downloaded on first run) |

---

## File structure

```
stage1/
├── config.py
├── data/
│   ├── bgc_hindcast_spco2.nc
│   ├── phy_reanalysis_sst_sss.nc
│   ├── multiobs_surface_carbon.nc
│   ├── era5_wind10m_monthly.nc
│   ├── era5_wind10m_daily.nc        ← optional (wind variance correction)
│   ├── co2_mm_gl.csv
│   ├── Time_Varying_Biomes.cmems.nc ← required for Fay figures
│   ├── oni_monthly.txt              ← auto-downloaded on first run
│   ├── sam_monthly.txt              ← auto-downloaded on first run
│   ├── processed_surface.nc         ← output of 02_preprocess.py
│   └── flux_3d.nc                   ← output of 03_compute_flux.py
├── output/
│   ├── global_flux.nc
│   ├── validation_metrics.csv
│   └── figures/
│       ├── fig01_flux_vs_co2.png
│       ├── fig02_mean_flux_map.png
│       ├── fig03_trend_significance_map.png
│       ├── fig04_delta_pco2_map.png
│       ├── fig05_seasonal_cycle.png
│       ├── fig06_sink_saturation.png
│       ├── fig07_climate_regression_global.png
│       ├── fig08_fay_ts.png
│       ├── fig09_fay_clim.png
│       ├── fig10_fay_trends.png
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
    ├── 05_plot_results.py
    └── plot_style.py
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

- **`config.py` is the single source of truth.** All product IDs, variable names, physical constants, file paths, and unit conversion factors live there.
- **Scripts are stateless and idempotent.** Re-running produces the same result; delete an output file to force recomputation.
- **No physics computation in `05_plot_results.py`.** It reads pre-computed NetCDF files only. The exception is lightweight statistical operations (Sen's slope, Mann-Kendall, OLS regression) applied directly on the already-computed flux arrays for plotting purposes.
- **Pure functions in `03_compute_flux.py`.** All physical parameterisations take `xr.DataArray` inputs and return `xr.DataArray` outputs.
- **Dask-backed lazy loading.** All `xr.open_dataset` calls use `chunks="auto"`.
- **Optional dependencies degrade gracefully.** `cartopy`, `cmocean`, `statsmodels`, `ruptures` each have fallback paths; the script never crashes due to a missing optional package.

---

## Output figures

### Primary result figures (`05_plot_results.py`)

| File | Description |
|------|-------------|
| `fig01_flux_vs_co2.png` | Twin-axis: J_net [Pg C yr⁻¹] (left) + atmospheric CO₂ [ppm] (right), monthly + annual mean |
| `fig02_mean_flux_map.png` | Time-mean air-sea flux map — blue = uptake, red = outgassing, Robinson projection |
| `fig03_trend_significance_map.png` | Per-pixel Sen's slope [mol C m⁻² yr⁻¹ per decade] + Mann-Kendall significance hatching (/// = p > 0.05) |
| `fig04_delta_pco2_map.png` | Time-mean ΔpCO₂ = pCO₂(ocean) − pCO₂(atm) [µatm] — thermodynamic driver |
| `fig05_seasonal_cycle.png` | Climatological seasonal cycle: J_net (bars) + atmospheric CO₂ (line) |
| `fig06_sink_saturation.png` | Annual J_net vs annual CO₂ scatter — ocean sink sensitivity / saturation analysis |
| `fig07_climate_regression_global.png` | 3-panel: (a) J_net + normalised ONI/SAM timeseries; (b) separate Spearman scatter vs ONI and SAM; (c) multilinear α·ONI+β·SAM scatter coloured by year — following Bellacicco et al. (2025). **Auto-skipped if index files absent (HPC).** |
| `fig08_fay_ts.png` | 17-panel monthly flux timeseries per Fay (2014) domain [Pg C yr⁻¹] + biome map |
| `fig09_fay_clim.png` | 17-panel climatological seasonal cycle ±1σ per Fay domain [Pg C yr⁻¹] |
| `fig10_fay_trends.png` | 17-panel Sen's slope per domain [Pg C yr⁻¹ per decade] + Mann-Kendall significance bars |

### Validation figures (`04_validate.py`)

| File | Description |
|------|-------------|
| `fig_validation_ts.png` | 2-row: global J_net timeseries (top) + 12-month rolling RMSD and Pearson r (bottom) |
| `fig_validation_loess.png` | LOESS + PELT breakpoints for reconstruction vs MULTIOBS |
| `fig_validation_map.png` | Spatial RMSD (top) and bias (bottom) maps, Robinson projection |
| `fig_fay_ts.png` | 17-panel monthly timeseries per Fay domain [Pg C yr⁻¹] + biome map |
| `fig_fay_clim.png` | 17-panel climatological seasonal cycle ±1σ per domain |
| `fig_fay_trends.png` | 17-panel Sen's slope per domain + Mann-Kendall significance (α=0.01) |
| `fig_fay_climate_regression.png` | 17-panel Spearman ρ bars: J_net vs ONI (red), SAM (green), α·ONI+β·SAM (blue); significance at α=0.05 — following Bellacicco et al. (2025). **Auto-skipped if index files absent (HPC).** |
| `fig_spectra.png` | Power spectral density (Welch) — reconstruction vs MULTIOBS |
| `fig_cps.png` | Cross-power spectrum: gain, phase [days], magnitude-squared coherence |

---

## Known limitations and caveats

1. **Physical consistency:** `GLOBAL_MULTIYEAR_BGC_001_029` was forced by FREEGLORYS2V4/ERA-Interim, while SST/SSS here come from GLORYS12V1 (ERA5). These are not dynamically consistent.

2. **Wind variance correction:** Active only when `era5_wind10m_daily.nc` is present. Without it, k is computed from monthly-mean ⟨u⟩² only, underestimating k in high-variability regions.

3. **Sea-ice masking:** Uses PISCES pCO₂ NaN mask as a proxy for ice-covered pixels rather than an explicit sea-ice concentration field.

4. **Hindcast, not reanalysis:** `GLOBAL_MULTIYEAR_BGC_001_029` assimilates no BGC observations. Known biases include a too-weak equatorial Pacific upwelling signal and Southern Ocean regional offsets.

5. **Autocorrelation in trends:** Mann-Kendall significance is computed on annual means to reduce serial autocorrelation. Monthly values would require pre-whitening or effective sample size correction.

6. **Saturation analysis (fig07):** The annual scatter regression is a diagnostic tool only. A significant positive slope does not prove linear forcing response — confounders include SST-driven solubility changes and ENSO-driven interannual variability.

7. **Record length:** A ~30-year window may not fully capture carbon cycle timescales of 100–1000 years. Results represent a snapshot of the current flux regime.

8. **HPC compute nodes (no internet):** `load_climate_indices()` auto-downloads ONI and SAM on first run. HPC compute nodes typically block outbound connections — the function degrades gracefully, prints `wget` commands for manual download from a login node, and skips fig07 / fig_fay_climate_regression without aborting the run. Manual download:
```bash
wget -O data/oni_monthly.txt https://psl.noaa.gov/data/correlation/oni.data
wget -O data/sam_monthly.txt https://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/aao/monthly.aao.index.b79.current.ascii.table
```

8. **MULTIOBS sign convention:** Verify against the product QUID before updating the MULTIOBS dataset version.

---

## Extending the project

This codebase is structured to grow stage by stage:

1. Create a new `stage2/` or `stage3/` directory at the same level as `stage1/`.
2. Each stage has its own `config.py`, `scripts/`, `data/`, and `output/` — keeping stages self-contained.
3. Shared utility functions (e.g., `compute_grid_cell_area`, ocean basin masks) should eventually be factored into a common `utils/` package.
4. Any new physical parameterisation should be added as a pure function with full docstring and citation.

Planned subsequent stages:
- **Stage 2** — export flux at 100 m and 500 m (gravitational pump)
- **Stage 3** — layered total carbon content (0–100 m / 100–500 m / >500 m) and mass-balance consistency check
- **Stage 4** — robust trend estimation across depth layers
- **Stage 5** — spatial mapping and hotspot clustering
- **Stage 6** — proposed extensions (thermal vs non-thermal flux decomposition)
- **Stage 7** — validation and synthesis

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

**Key motivating reference — physical injection pump & climate index regression:**
> Bellacicco, M., Marullo, S., Dall'Olmo, G., Iudicone, D., & Buongiorno Nardelli, B. (2025). The oceanic physical injection pump of organic carbon. *Nature Communications*, 16, 7100. https://doi.org/10.1038/s41467-025-62363-z

**SAM index:**
> Marshall, G. J. (2003). Trends in the Southern Annular Mode from observations and reanalyses. *Journal of Climate*, 16(24), 4134–4143. https://doi.org/10.1175/1520-0442(2003)016<4134:TITSAM>2.0.CO;2

**ONI index:**
> NOAA CPC (2024). Oceanic Niño Index (ONI). NOAA Physical Sciences Laboratory. https://psl.noaa.gov/data/correlation/oni.data

**Carbon pump conceptual framework:**
> Boyd, P. W., Claustre, H., Levy, M., Siegel, D. A., & Weber, T. (2019). Multi-faceted particle pumps drive carbon sequestration in the ocean. *Nature*, 568, 327–335. https://doi.org/10.1038/s41586-019-1098-2

**ERA5 reanalysis:**
> Hersbach, H. et al. (2020). The ERA5 global reanalysis. *Quarterly Journal of the Royal Meteorological Society*, 146(730), 1999–2049. https://doi.org/10.1002/qj.3803

**Biome domains:**
> Fay, A. R. & McKinley, G. A. (2014). Global open-ocean biomes: mean and temporal variability. *Earth System Science Data*, 6, 273–284. https://doi.org/10.5194/essd-6-273-2014

**Robust trend estimation:**
> Sen, P. K. (1968). Estimates of the regression coefficient based on Kendall's tau. *Journal of the American Statistical Association*, 63(324), 1379–1389. https://doi.org/10.1080/01621459.1968.10480934

**Mann-Kendall significance test:**
> Mann, H. B. (1945). Nonparametric tests against trend. *Econometrica*, 13(3), 245–259. https://doi.org/10.2307/1907187
> Kendall, M. G. (1975). *Rank Correlation Methods* (4th ed.). Griffin, London.

**LOESS smoothing:**
> Cleveland, W. S. (1979). Robust locally weighted regression and smoothing scatterplots. *Journal of the American Statistical Association*, 74(368), 829–836. https://doi.org/10.1080/01621459.1979.10481038

---

## License and citation

**Data licenses:** All CMEMS products are distributed under the [Copernicus Marine Service licence](https://marine.copernicus.eu/user-corner/service-commitments-and-licence). ERA5 data is subject to the [Copernicus Climate Change Service licence](https://cds.climate.copernicus.eu/api/v2/terms/static/licence-to-use-copernicus-products.pdf). NOAA GML data is public domain.

If you use this code in a publication, please cite the underlying data products and parameterisations listed above, and acknowledge the Copernicus Marine Service:

> *"This study has been conducted using E.U. Copernicus Marine Service Information; https://doi.org/10.48670/moi-00019"*
