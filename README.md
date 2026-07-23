# Stage 1 — Surface Net Air-Sea CO₂ Flux

**Ocean Carbon Cycle · Net Flux & Layered Accumulation Project**

---

## 1. What this stage computes

Stage 1 estimates the **net air-sea CO₂ flux at the ocean surface** — the quantity the meeting labels *J(0 m)* — at every 0.25° grid point and every month across the full Copernicus reanalysis record (1993–2026), then:

1. Integrates it over the **global ocean surface** to get a single time series of net uptake in Pg C yr⁻¹.
2. Plots it against **atmospheric CO₂** (NOAA GML) to answer the central question: *is the ocean's carbon sink growing, flattening, or saturating under rising CO₂?*
3. Maps it spatially to identify **regional sources and sinks**.
4. Estimates **per-pixel trends** (Sen's slope) over the full record.
5. **Cross-validates** the reconstructed flux against an independent observation-based product (CMEMS MULTIOBS / SOCAT neural network).

---

## 2. Mathematical formulation

### 2.1 Bulk air-sea CO₂ flux parameterisation

The air-sea CO₂ flux per unit area is given by:

$$\boxed{F = k \cdot K_0 \cdot \left( pCO_2^{\mathrm{ocn}} - pCO_2^{\mathrm{atm}} \right)}$$

| Symbol | Name | Units |
|--------|------|-------|
| $F$ | Net CO₂ flux (positive = outgassing) | mol C m⁻² yr⁻¹ |
| $k$ | Gas transfer velocity | m s⁻¹ |
| $K_0$ | CO₂ solubility | mol m⁻³ atm⁻¹ |
| $pCO_2^{\mathrm{ocn}}$ | Surface ocean pCO₂ | atm |
| $pCO_2^{\mathrm{atm}}$ | Atmospheric pCO₂ | atm |

**Sign convention:** $F > 0$ = ocean outgassing (CO₂ source); $F < 0$ = ocean uptake (CO₂ sink). This matches the OCMIP and CMEMS convention.

---

### 2.2 Gas transfer velocity — Wanninkhof (2014)

The gas transfer velocity $k$ is parameterised as a quadratic function of the 10 m wind speed following Wanninkhof (2014):

$$k = a \cdot u_{10}^{2} \cdot \left(\frac{Sc}{660}\right)^{-1/2}$$

where:

| Symbol | Value / definition |
|--------|--------------------|
| $a$ | 0.251 cm hr⁻¹ (m s⁻¹)⁻² |
| $u_{10}$ | 10 m scalar wind speed [m s⁻¹] (from ERA5) |
| $Sc$ | Schmidt number of CO₂ in seawater (see §2.3) |
| 660 | Reference Schmidt number of CO₂ at 20°C in seawater |

The factor $(Sc/660)^{-1/2}$ normalises $k$ to the reference condition.

$k$ is computed in cm hr⁻¹ and converted to m s⁻¹:
$$k \,[\text{m s}^{-1}] = k \,[\text{cm hr}^{-1}] \times \frac{1}{100 \times 3600}$$

---

### 2.3 Schmidt number — Wanninkhof (2014)

The Schmidt number of CO₂ in seawater as a function of sea surface temperature $T$ (in °C):

$$Sc = A - B \cdot T + C \cdot T^2 - D \cdot T^3 + E \cdot T^4$$

Coefficients from Wanninkhof (2014), Table 1:

| $A$ | $B$ | $C$ | $D$ | $E$ |
|-----|-----|-----|-----|-----|
| 2116.8 | 136.25 | 4.7353 | 0.092307 | 0.0007555 |

---

### 2.4 CO₂ solubility — Weiss (1974)

The solubility of CO₂ in seawater $K_0$ is computed using the empirical relation of Weiss (1974):

$$\ln K_0 = A_1 + A_2 \cdot \frac{100}{T} + A_3 \cdot \ln\!\left(\frac{T}{100}\right) + S \cdot \left[ B_1 + B_2 \cdot \frac{T}{100} + B_3 \cdot \left(\frac{T}{100}\right)^2 \right]$$

where $T$ is temperature in **Kelvin** and $S$ is salinity in **PSU**.

Coefficients (Weiss 1974, Table I):

| $A_1$ | $A_2$ | $A_3$ | $B_1$ | $B_2$ | $B_3$ |
|-------|-------|-------|-------|-------|-------|
| −58.0931 | 90.5069 | 22.2940 | 0.027766 | −0.025888 | 0.0050578 |

Weiss (1974) gives $K_0$ in mol L⁻¹ atm⁻¹; we multiply by 1000 to convert to mol m⁻³ atm⁻¹.

---

### 2.5 Global surface integral

The global net ocean CO₂ uptake at time $t$ is:

$$J_{\mathrm{net}}(t) = -\iint_{\mathrm{ocean}} F(t, \varphi, \lambda) \; dA$$

where $dA = R^2 \cos\varphi \; d\varphi \; d\lambda$ is the area element on the sphere (R = 6.371 × 10⁶ m).

On the discrete 0.25° grid:

$$J_{\mathrm{net}}(t) = -\sum_{i,j} F(t, i, j) \cdot A(i, j)$$

with the sign flip so that **$J_{\mathrm{net}} > 0$** = ocean uptake (sink).

The result is converted from mol C yr⁻¹ to Pg C yr⁻¹:

$$J_{\mathrm{net}} \,[\text{Pg C yr}^{-1}] = J_{\mathrm{net}} \,[\text{mol C yr}^{-1}] \times \frac{12.011 \,\text{g mol}^{-1}}{10^{15} \,\text{g Pg}^{-1}}$$

---

### 2.6 Per-pixel trend estimation — Sen's slope

The trend at each pixel is estimated using the **Theil-Sen estimator** (Sen's slope):

$$\hat{\beta} = \mathrm{median}\!\left(\frac{y_j - y_i}{t_j - t_i}\right), \quad \forall \; i < j$$

This is the median of all pairwise slopes. It is preferred over OLS because:
- It is **robust to outliers** (e.g., ENSO-driven anomalous years).
- It makes **no assumption of normality** on the residuals.
- It is consistent under **serial correlation** (geophysical time series are typically autocorrelated).

Units: mol C m⁻² yr⁻¹ per decade (slope × 10).

---

## 3. Data sources

### 3.1 Primary — surface ocean pCO₂

| Field | Source | Product ID | URL |
|-------|--------|-----------|-----|
| Surface pCO₂ (`spco2`) | CMEMS Global Ocean Biogeochemistry Hindcast | `GLOBAL_MULTIYEAR_BGC_001_029` | https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_BGC_001_029/description |
| Record | 1993-01 to 2026-04 | Resolution | 0.25°, monthly |
| Model | NEMO-PISCES (Mercator Ocean) | Assimilation | None |

### 3.2 Physical fields — SST, SSS

| Field | Source | Product ID | URL |
|-------|--------|-----------|-----|
| SST (`thetao`), SSS (`so`) | CMEMS GLORYS12V1 Physical Reanalysis | `GLOBAL_MULTIYEAR_PHY_001_030` | https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/description |
| Record | 1993–present | Resolution | 1/12° → regridded to 0.25° |

### 3.3 Wind speed

| Field | Source | Product ID | URL |
|-------|--------|-----------|-----|
| u10, v10 → `wind_speed` | ERA5 monthly averaged reanalysis | `reanalysis-era5-single-levels-monthly-means` | https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels-monthly-means |
| Record | 1940–present | Resolution | 0.25° |
| Access | Copernicus CDS API (separate from CMEMS — requires CDS account) | | |

### 3.4 Atmospheric CO₂

| Field | Source | URL |
|-------|--------|-----|
| Global mean monthly CO₂ (ppm) | NOAA GML Marine Boundary Layer Reference | https://gml.noaa.gov/ccgg/trends/gl_data.html |
| File | `co2_mm_gl.csv` | Direct download | https://gml.noaa.gov/webdata/ccgg/trends/co2/co2_mm_gl.csv |
| Record | 1979–present | Frequency | Monthly |

> **Note:** 1 ppm ≈ 1 µatm at 1 atm total atmospheric pressure. For precise work, multiply by the mean atmospheric pressure at sea level (101325 Pa / atm) — the difference is < 0.1% and negligible here.

### 3.5 Validation — observation-based flux

| Field | Source | Product ID | URL |
|-------|--------|-----------|-----|
| Air-sea CO₂ flux (`fgco2`), surface pCO₂ | CMEMS MULTIOBS surface carbon L4 | `MULTIOBS_GLO_BGC_CARBON_SURFACE_MYNRT_015_008` | https://data.marine.copernicus.eu/product/MULTIOBS_GLO_BGC_CARBON_SURFACE_MYNRT_015_008/description |
| Method | Ensemble neural network trained on SOCAT in-situ fCO₂ data | | |
| Record | 1985–present | Resolution | 0.25°, monthly |

---

## 4. Software requirements

```
python          >= 3.11
copernicusmarine >= 1.3        # CMEMS download
cdsapi          >= 0.7         # ERA5 download (CDS)
xarray          >= 2024.1
dask            >= 2024.1      # out-of-core computation
numpy           >= 1.26
scipy           >= 1.12
pandas          >= 2.2
matplotlib      >= 3.8
cartopy         >= 0.23        # map projections (optional but recommended)
requests        >= 2.31        # NOAA CSV download
```

Install:
```bash
pip install copernicusmarine cdsapi xarray[io] dask numpy scipy pandas \
            matplotlib cartopy requests
```

---

## 5. Credentials setup

### CMEMS (for scripts 01, 02, 03)
```bash
copernicusmarine login
# enter your username and password when prompted
# credentials are stored in ~/.copernicusmarine/
```
Register at: https://data.marine.copernicus.eu/register

### CDS / ERA5 (for wind speed in script 01)
Create `~/.cdsapirc`:
```
url: https://cds.climate.copernicus.eu/api
key: <YOUR_UID>:<YOUR_API_KEY>
```
Get your key at: https://cds.climate.copernicus.eu/user/register

---

## 6. How to run

Run scripts in order. Each script checks for its expected input and skips or raises a clear error if something is missing.

```bash
# Step 0 — download everything (add --test for one year only)
python scripts/01_download_data.py [--test]

# Step 1 — preprocess: regrid, unit-convert, harmonize, merge
python scripts/02_preprocess.py

# Step 2 — compute flux F, integrate globally
python scripts/03_compute_flux.py

# Step 3 — cross-validate against MULTIOBS
python scripts/04_validate.py

# Step 4 — produce all result figures
python scripts/05_plot_results.py
```

Expected runtime on a modern laptop (full 1993–2026 record, global 0.25°):
- Download: 30–90 min (depends on connection)
- Preprocess: 10–20 min
- Flux compute: 15–30 min
- Validation: 20–40 min (pixel-wise loops)
- Plotting: 10–20 min

---

## 7. Output files

```
output/
├── global_flux.nc            # time series of global J_net [Pg C yr-1]
├── validation_metrics.csv    # scalar RMSD, bias, r vs MULTIOBS
└── figures/
    ├── fig01_flux_vs_co2.png         # KEY FIGURE — J_net vs atm CO2
    ├── fig02_annual_flux_map.png     # mean spatial flux
    ├── fig03_flux_trend_map.png      # Sen's slope per pixel
    ├── fig04_delta_pco2_map.png      # ΔpCO2 driver map
    ├── fig05_monthly_seasonal_cycle.png
    ├── fig_validation_ts.png         # validation time series
    ├── fig_validation_map.png        # RMSD / bias maps
    └── fig_taylor.png                # Taylor diagram
data/
├── bgc_hindcast_spco2.nc     # raw CMEMS download
├── phy_reanalysis_sst_sss.nc # raw CMEMS download
├── multiobs_surface_carbon.nc
├── era5_wind10m_monthly.nc
├── co2_mm_gl.csv             # NOAA GML
├── processed_surface.nc      # merged, regridded, unit-converted
└── flux_3d.nc                # full (time, lat, lon) flux field
```

---

## 8. Bibliography

**Flux parameterisation — gas transfer velocity:**
> Wanninkhof, R. (2014). Relationship between wind speed and gas exchange over the ocean revisited. *Limnology and Oceanography: Methods*, 12(6), 351–362. https://doi.org/10.4319/lom.2014.12.351

**CO₂ solubility:**
> Weiss, R. F. (1974). Carbon dioxide in water and seawater: the solubility of a non-ideal gas. *Marine Chemistry*, 2(3), 203–215. https://doi.org/10.1016/0304-4203(74)90015-2

**Schmidt number coefficients:**
> Wanninkhof, R. (1992). Relationship between wind speed and gas exchange over the ocean. *Journal of Geophysical Research: Oceans*, 97(C5), 7373–7382. https://doi.org/10.1029/92JC00188

**BGC reanalysis model:**
> Aumont, O., Ethé, C., Tagliabue, A., Bopp, L., & Gehlen, M. (2015). PISCES-v2: an ocean biogeochemical model for carbon and ecosystem studies. *Geoscientific Model Development*, 8, 2465–2513. https://doi.org/10.5194/gmd-8-2465-2015

**Physical reanalysis (GLORYS12):**
> Lellouche, J.-M. et al. (2021). The Copernicus Global 1/12° Oceanic and Sea Ice GLORYS12 Reanalysis. *Frontiers in Earth Science*, 9, 698876. https://doi.org/10.3389/feart.2021.698876

**Validation product — SOCAT neural network:**
> Bakker, D. C. E. et al. (2016). A multi-decade record of high-quality fCO₂ data in version 3 of the Surface Ocean CO₂ Atlas (SOCAT). *Earth System Science Data*, 8, 383–413. https://doi.org/10.5194/essd-8-383-2016

**Key motivating reference — physical injection pump:**
> Bellacicco, M., Marullo, S., Dall'Olmo, G., Iudicone, D., & Buongiorno Nardelli, B. (2025). The oceanic physical injection pump of organic carbon. *Nature Communications*, 16, 7100. https://doi.org/10.1038/s41467-025-62363-z

**Carbon pump conceptual framework:**
> Boyd, P. W., Claustre, H., Levy, M., Siegel, D. A., & Weber, T. (2019). Multi-faceted particle pumps drive carbon sequestration in the ocean. *Nature*, 568, 327–335. https://doi.org/10.1038/s41586-019-1098-2

**ERA5 reanalysis:**
> Hersbach, H. et al. (2020). The ERA5 global reanalysis. *Quarterly Journal of the Royal Meteorological Society*, 146(730), 1999–2049. https://doi.org/10.1002/qj.3803

**Robust trend estimation:**
> Sen, P. K. (1968). Estimates of the regression coefficient based on Kendall's tau. *Journal of the American Statistical Association*, 63(324), 1379–1389. https://doi.org/10.1080/01621459.1968.10480934

---

## 9. Known limitations and caveats

1. **pCO₂ consistency:** The BGC hindcast (`001_029`) was forced by FREEGLORYS2V4/ERA-Interim physics, while the SST/SSS used here comes from GLORYS12V1 (ERA5 forcing). These are not dynamically consistent. An alternative is to use only the T/S fields archived internally by the PISCES run itself — check with Mercator Ocean if available.

2. **Wind speed:** The Wanninkhof (2014) coefficient $a = 0.251$ was calibrated against global bomb-¹⁴C inventories and assumes a specific (quadratic) wind speed distribution. Applying it to monthly-mean wind speeds underestimates $k$ relative to applying it to high-frequency (6-hourly) data; a variance correction factor $c$ should be added: $k = a(u^2 + \sigma_u^2)(\ldots)$ where $\sigma_u^2$ is the sub-monthly wind variance. ERA5 daily or 6-hourly data can supply this.

3. **Sea-ice masking:** Grid cells with sea ice cover should have their flux zeroed or reduced proportionally to the open-water fraction. The current implementation uses the PISCES pCO₂ NaN mask (which PISCES sets to zero under ice) — verify this is handled correctly in the CMEMS product.

4. **DOC:** Dissolved organic carbon is not included in Stage 1 (only the air-sea flux). It becomes relevant in Stage 3 (layered carbon content). No 3D time-evolving DOC product currently exists in the global CMEMS catalogue.

5. **Record length and pseudo-steady-state assumption:** As discussed in the meeting, a 30-year window may not fully capture carbon cycle timescales of 100–1000 years. Results should be interpreted as a snapshot of the current flux regime, not a long-term equilibrium.

---

## 10. Code structure

```
stage1/
├── README.md               ← this file
├── CHANGELOG.md            ← version history
├── config.py               ← central configuration (paths, constants, product IDs)
├── requirements.txt        ← Python dependencies
├── data/                   ← raw and processed NetCDF files (not tracked by git — see §11)
│   ├── bgc_hindcast_spco2.nc
│   ├── phy_reanalysis_sst_sss.nc
│   ├── multiobs_surface_carbon.nc
│   ├── era5_wind10m_monthly.nc
│   ├── co2_mm_gl.csv
│   ├── processed_surface.nc    ← output of 02_preprocess.py
│   └── flux_3d.nc              ← output of 03_compute_flux.py
├── output/
│   ├── global_flux.nc          ← global integral time series
│   ├── validation_metrics.csv
│   └── figures/
│       ├── fig01_flux_vs_co2.png
│       ├── fig02_annual_flux_map.png
│       ├── fig03_flux_trend_map.png
│       ├── fig04_delta_pco2_map.png
│       ├── fig05_monthly_seasonal_cycle.png
│       ├── fig_validation_ts.png
│       ├── fig_validation_map.png
│       └── fig_taylor.png
└── scripts/
    ├── 01_download_data.py     ← download all raw data from CMEMS, NOAA, ERA5
    ├── 02_preprocess.py        ← regrid, unit-convert, harmonize, merge
    ├── 03_compute_flux.py      ← core physics: Sc, K0, k, F, global integral
    ├── 04_validate.py          ← cross-validation vs MULTIOBS (RMSD, bias, Taylor)
    └── 05_plot_results.py      ← all figures
```

### Script dependency chain

```
01_download_data.py
        │
        ▼
02_preprocess.py  ──────────────────────────────────────────────────────┐
        │                                                                │
        ▼                                                                │
03_compute_flux.py                                                       │
        │                                                                │
        ├──────────────────────────────────────────────────┐            │
        ▼                                                  ▼            │
04_validate.py                                    05_plot_results.py ◄──┘
```

### Design principles

- **`config.py` is the single source of truth.** All product IDs, variable names, physical constants, file paths, and unit conversion factors live there. If a CMEMS product is updated and a variable name changes, you fix it in one place.
- **Scripts are stateless and idempotent.** Each script reads files, does its work, and writes output files. Re-running a script produces the same result. If the output already exists, the script skips the computation (delete the file to force a rerun).
- **No computation in plot scripts.** `05_plot_results.py` reads pre-computed NetCDF files only. This separation means you can iterate on figure aesthetics without rerunning the heavy physics.
- **Pure functions in `03_compute_flux.py`.** All physical parameterisations (`schmidt_number_co2`, `gas_transfer_velocity`, `co2_solubility_K0`, `compute_flux`) are standalone functions that take `xr.DataArray` inputs and return `xr.DataArray` outputs — easy to test in isolation or swap for alternative formulations.
- **Dask-backed lazy loading.** All `xr.open_dataset` calls use `chunks="auto"` so that the full 30-year global 3D archive is never loaded into RAM at once.

---

## 11. Git setup recommendations

### `.gitignore`

Data files (NetCDF, CSV) and figure outputs are large and should not be tracked by git. Add this `.gitignore` at the repo root:

```gitignore
# Data files
data/
output/figures/
output/*.nc
output/*.csv
*.nc
*.h5
*.hdf5

# Python
__pycache__/
*.pyc
*.pyo
.eggs/
*.egg-info/
dist/
build/

# Environments
.venv/
env/
venv/

# Credentials — NEVER commit these
.cdsapirc
.copernicusmarine/

# OS
.DS_Store
Thumbs.db
```

### Recommended branch strategy

```
main          ← stable, tagged releases only
dev           ← integration branch for ongoing work
stage1/       ← feature branches per stage (stage1/flux-validation, etc.)
```

### Tagging releases

```bash
git tag -a v1.0.0 -m "Stage 1: surface air-sea CO2 flux (1993-2026)"
git push origin v1.0.0
```

---

## 12. Contributing and extending

This codebase is structured to grow stage by stage. When adding Stage 2 (export flux at 100 m / 500 m) or Stage 3 (layered carbon content):

1. Create a new `stage2/` or `stage3/` directory at the same level as `stage1/`.
2. Each stage has its own `config.py`, `scripts/`, `data/`, and `output/` — keeping stages self-contained.
3. Shared utility functions (e.g., `compute_grid_cell_area`, ocean basin masks) should eventually be factored into a common `utils/` package at the repo root.
4. Any new physical parameterisation should be added as a pure function in the relevant `03_compute_*.py` with full docstring including the reference equation and citation.

---

## 13. License and citation

> **Data licenses:** All CMEMS products are distributed under the [Copernicus Marine Service licence](https://marine.copernicus.eu/user-corner/service-commitments-and-licence). ERA5 data is subject to the [Copernicus Climate Change Service licence](https://cds.climate.copernicus.eu/api/v2/terms/static/licence-to-use-copernicus-products.pdf). NOAA GML data is public domain.

If you use this code in a publication, please cite the underlying data products and parameterisations as listed in §8, and acknowledge the Copernicus Marine Service:
> *"This study has been conducted using E.U. Copernicus Marine Service Information; https://doi.org/10.48670/moi-00019"*
