# Stage 1 — Surface Net Air-Sea CO₂ Flux

**Ocean Carbon Cycle · Net Flux & Layered Accumulation Project**

---

## Overview

Stage 1 estimates the net air-sea CO₂ flux at every 0.25° grid point monthly across the Copernicus BGC hindcast (1993–2023), then integrates, maps, and analyses it. The reconstruction uses `GLOBAL_MULTIYEAR_BGC_001_029` (NEMO-PISCES, no BGC assimilation) for pCO₂; GLORYS12V1 for SST/SSS; ERA5 for wind. The MULTIOBS product is used for validation only.

---

## Physics

### Air-sea CO₂ flux

$$F = k \cdot K_0 \cdot \left( pCO_2^{\mathrm{atm}} - pCO_2^{\mathrm{ocn}} \right) \quad [\text{mol C m}^{-2}\text{ yr}^{-1}]$$

$F > 0$ = ocean uptake; $F < 0$ = outgassing.

**Gas transfer velocity** (Wanninkhof 2014):

$$k = a \cdot \left(\langle u_{10} \rangle^2 + \sigma_u^2\right) \cdot \left(\frac{Sc}{660}\right)^{-1/2}$$

$a = 0.251$ cm hr⁻¹ (m s⁻¹)⁻². $\sigma_u^2$ from ERA5 daily winds (optional).

**Schmidt number** (Wanninkhof 2014, Table 1): $Sc = A - BT + CT^2 - DT^3 + ET^4$, coefficients (2116.8, 136.25, 4.7353, 0.092307, 0.0007555).

**CO₂ solubility** (Weiss 1974): $\ln K_0 = f(T[\text{K}], S[\text{PSU}])$, coefficients from Eq. 12.

**Global integral**: $J_\mathrm{net}(t) = \sum_{i,j} F(t,i,j) \cdot R^2 \cos\varphi \,\Delta\varphi\,\Delta\lambda \times (12.011/10^{15})$ [Pg C yr⁻¹].

### Statistical methods

**Sen's slope** (Theil-Sen estimator) for per-pixel trends: robust to outliers, consistent with Mann-Kendall. **Mann-Kendall** (two-sided, α = 0.01) for significance: rank-based, distribution-free p-value. Applied to annual means to reduce autocorrelation.

**LOESS** (frac = 0.15, ≈ 18-month window) to reveal interannual variability. **PELT** (RBF cost, pen = 3) detects breakpoints on the LOESS smooth.

**Spearman ρ** for climate index correlations: rank-based, appropriate for small n (31 years) and potentially non-linear ENSO/flux relationships.

**EOF** (area-weighted SVD): decomposes the annual flux anomaly field into orthogonal spatial patterns (EOFs) and their temporal amplitudes (PCs). Pixels weighted by √cos(lat). North et al. (1982) sampling error used to identify degenerate modes.

---

## Script pipeline

```
01_download_data.py   →   02_preprocess.py   →   03_compute_flux.py
                                                        │
                          ┌─────────────────────────────┤
                          │                             │
                   04_validate.py          05_plot_results.py
                                           06_fay_analysis.py
                                           07_eof_analysis.py
```

| Script | Purpose |
|--------|---------|
| `01_download_data.py` | CMEMS, NOAA, ERA5 data acquisition |
| `02_preprocess.py` | Regrid, unit convert, harmonize, build `processed_surface.nc` |
| `03_compute_flux.py` | Flux computation, global integral, `flux_3d.nc`, `global_flux.nc` |
| `04_validate.py` | Validation vs MULTIOBS — timeseries, maps, spectra, Fay domains |
| `05_plot_results.py` | 7 primary result figures (global) |
| `06_fay_analysis.py` | 4 Fay domain figures (heavy — run separately) |
| `07_eof_analysis.py` | EOF decomposition of annual flux anomaly field |

---

## Data sources

| Dataset | Product / URL | Variables |
|---------|--------------|-----------|
| BGC hindcast | `GLOBAL_MULTIYEAR_BGC_001_029` | `spco2` |
| Physical reanalysis | `GLOBAL_MULTIYEAR_PHY_001_030` (GLORYS12V1) | `thetao`, `so` |
| Wind | ERA5 monthly + daily | `u10`, `v10` |
| Atmospheric CO₂ | NOAA GML `co2_mm_gl.csv` | `average` [ppm] |
| Validation | `MULTIOBS_GLO_BGC_CARBON_SURFACE_MYNRT_015_008` | `fgco2`, `spco2` |
| Fay domains | `Time_Varying_Biomes.cmems.nc` | `MeanBiomes` |
| ONI | https://psl.noaa.gov/data/correlation/oni.data | — |
| SAM | https://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/aao/monthly.aao.index.b79.current.ascii.table | — |

ONI and SAM are auto-downloaded on first run. On HPC nodes without internet, download manually from a login node:
```bash
wget -O data/oni_monthly.txt https://psl.noaa.gov/data/correlation/oni.data
wget -O data/sam_monthly.txt https://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/aao/monthly.aao.index.b79.current.ascii.table
```

---

## Output figures

### `05_plot_results.py`

| File | Description |
|------|-------------|
| `fig01_flux_vs_co2.png` | J_net + LOESS + B1/B4 breakpoint trends vs atmospheric CO₂ |
| `fig02_mean_flux_map.png` | Time-mean flux map, ±2.5 mol C m⁻² yr⁻¹ |
| `fig03_trend_significance_map.png` | Sen's slope ±0.5 per decade + MK significance (α=0.01) |
| `fig04_delta_pco2_map.png` | Time-mean ΔpCO₂ ±100 µatm |
| `fig05_seasonal_cycle.png` | Climatological seasonal cycle |
| `fig06_sink_saturation.png` | Annual J_net vs CO₂ — saturation diagnostic |
| `fig07_climate_regression_global.png` | Spearman ρ vs ONI/SAM + multilinear combination |

### `06_fay_analysis.py`

| File | Description |
|------|-------------|
| `fig_fay_ts.png` | 17-panel monthly flux timeseries per Fay domain |
| `fig_fay_clim.png` | 17-panel seasonal climatology ±1σ |
| `fig_fay_trends.png` | 17-panel Sen's slope + MK significance |
| `fig_fay_climate_regression.png` | 17-panel Spearman ρ vs ONI, SAM, ONI+SAM |

### `07_eof_analysis.py`

| File | Description |
|------|-------------|
| `fig_eof_annual_full_scree.png` | Scree + North error + cumulative variance — annual full field |
| `fig_eof_annual_full_modes.png` | 3×2 composite (EOF map + PC) modes 1–3 — annual full field |

### `04_validate.py`

| File | Description |
|------|-------------|
| `fig_validation_ts.png` | Global timeseries + rolling RMSD/r |
| `fig_validation_loess.png` | LOESS + PELT breakpoints, reconstruction vs MULTIOBS |
| `fig_validation_map.png` | Spatial RMSD and bias maps |
| `fig_fay_ts.png` | Fay domain timeseries, reconstruction vs MULTIOBS |
| `fig_fay_clim.png` | Fay domain seasonal cycle ±1σ |
| `fig_spectra.png` | Welch PSD comparison |
| `fig_cps.png` | Cross-power spectrum: gain, phase, coherence |

---

## Caveats

1. `GLOBAL_MULTIYEAR_BGC_001_029` forced by FREEGLORYS2V4/ERA-Interim; SST/SSS from GLORYS12V1 (ERA5). Not dynamically consistent.
2. Wind variance correction active only when `era5_wind10m_daily.nc` present.
3. Sea-ice mask derived from PISCES NaN pattern, not explicit ice concentration.
4. No BGC data assimilation — known biases in equatorial Pacific and Southern Ocean.
5. MK significance on annual means reduces but does not eliminate autocorrelation effects.
6. EOF degeneracy: modes where North et al. error bars overlap eigenvalues cannot be individually interpreted — check the scree plot. With n=31 annual means the North error is ~25% per eigenvalue; with monthly data (n≈372) it is ~7%, making more modes individually resolvable. The `vmax` colourbar falls back to 1.0 for any EOF mode whose map is entirely NaN.
7. 30-year record captures current flux regime but not centennial-scale carbon cycle dynamics.

---

## Planned stages

- **Stage 2** — Export flux at 100 m and 500 m
- **Stage 3** — Layered carbon content (0–100 m / 100–500 m / >500 m) + mass balance
- **Stage 4** — Robust trend estimation across depth layers
- **Stage 5** — Spatial mapping and hotspot clustering
- **Stage 6** — Thermal vs non-thermal flux decomposition
- **Stage 7** — Validation and synthesis

---

## Bibliography

> Wanninkhof (2014) https://doi.org/10.4319/lom.2014.12.351  
> Weiss (1974) https://doi.org/10.1016/0304-4203(74)90015-2  
> Aumont et al. (2015) https://doi.org/10.5194/gmd-8-2465-2015  
> Lellouche et al. (2021) https://doi.org/10.3389/feart.2021.698876  
> Bakker et al. (2016) https://doi.org/10.5194/essd-8-383-2016  
> Bellacicco et al. (2025) https://doi.org/10.1038/s41467-025-62363-z  
> Marshall (2003) https://doi.org/10.1175/1520-0442(2003)016<4134:TITSAM>2.0.CO;2  
> Boyd et al. (2019) https://doi.org/10.1038/s41586-019-1098-2  
> Hersbach et al. (2020) https://doi.org/10.1002/qj.3803  
> Fay & McKinley (2014) https://doi.org/10.5194/essd-6-273-2014  
> Sen (1968) https://doi.org/10.1080/01621459.1968.10480934  
> Mann (1945) https://doi.org/10.2307/1907187  
> Cleveland (1979) https://doi.org/10.1080/01621459.1979.10481038  
> North et al. (1982) https://doi.org/10.1175/1520-0493(1982)110<0699:SEITEO>2.0.CO;2  
> Björnsson & Venegas (1997) McGill CCGCR Report 97-1

Data: *"This study has been conducted using E.U. Copernicus Marine Service Information; https://doi.org/10.48670/moi-00019"*
