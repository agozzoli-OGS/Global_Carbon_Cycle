# Stage 1 — Surface Net Air-Sea CO₂ Flux: Technical Reference

**Ocean Carbon Cycle · Net Flux & Layered Accumulation Project**  
*Detailed mathematical, methodological, and implementation reference.*  
*For a concise operational summary see README.md.*

---

## Table of Contents

1. [Scientific context](#1-scientific-context)
2. [Data sources and products](#2-data-sources-and-products)
3. [Complete mathematical formulation](#3-complete-mathematical-formulation)
4. [Statistical methods](#4-statistical-methods)
5. [File structure and design principles](#5-file-structure-and-design-principles)
6. [Script pipeline — detailed description](#6-script-pipeline--detailed-description)
7. [Output figures — detailed description](#7-output-figures--detailed-description)
8. [Known limitations and caveats](#8-known-limitations-and-caveats)
9. [Extending the project](#9-extending-the-project)
10. [Full bibliography](#10-full-bibliography)
11. [Data licences and citation](#11-data-licences-and-citation)

---

## 1. Scientific context

The ocean is the largest active carbon sink on Earth, absorbing approximately 2–3 Pg C yr⁻¹ on a net annual basis — roughly 25–30% of anthropogenic CO₂ emissions. This uptake occurs through two coupled mechanisms: the **solubility pump** (CO₂ dissolves more readily in cold water, which sinks at high latitudes) and the **biological pump** (photosynthesis fixes dissolved CO₂ into organic matter that sinks before respiration returns it to the atmosphere). The balance between uptake and outgassing is set locally by the partial pressure difference ΔpCO₂ = pCO₂(ocean) − pCO₂(atm) and by how efficiently the gas can cross the air-sea interface, parameterised through the gas transfer velocity k.

Stage 1 focuses exclusively on the **surface interface flux** — what we call J(0 m) — using the standard bulk parameterisation. This is the entry point into the ocean carbon cycle: everything that happens at depth (export, accumulation, remineralisation) is downstream of this surface exchange. Understanding its spatial patterns, temporal trends, interannual variability, and dominant modes of variability is a prerequisite for interpreting the deeper stages.

The central scientific questions addressed here are:

- Is the global ocean sink growing proportionally with rising atmospheric CO₂, or is it beginning to saturate?
- Where is the sink intensifying or weakening, and is that trend statistically robust?
- What fraction of interannual variability is explained by large-scale climate modes (ENSO, SAM)?
- What are the dominant spatial patterns of flux variability, and how do they relate to known physical drivers?

The motivating paper is **Bellacicco et al. (2025, *Nature Communications*)**, which analyses the Physical Injection Pump (PIP) — the mechanical export of organic carbon below the mixed layer by physical processes. That paper finds strong interannual correlations between the PIP and ENSO/SAM indices. Stage 1 tests whether the same teleconnection already manifests at the air-sea interface, enabling a direct surface-to-depth comparison.

---

## 2. Data sources and products

### 2.1 Surface ocean pCO₂ — primary reconstruction input

| | |
|---|---|
| **Product** | CMEMS Global Ocean Biogeochemistry Hindcast |
| **ID** | `GLOBAL_MULTIYEAR_BGC_001_029` |
| **Model** | NEMO-PISCES (free-running, no BGC data assimilation) |
| **Variable** | `spco2` — surface partial pressure of CO₂ [Pa] → converted to [atm] |
| **Resolution** | 0.25°, monthly, 1993–2024 |
| **Depth** | Surface level only |
| **Nature** | Pure model output. All pCO₂ variability is model-generated. |

This is the sole pCO₂ source for the reconstruction. It is explicitly a **hindcast** (not a reanalysis) — the biogeochemical state is not constrained by any observational assimilation. Known model biases include an underestimate of equatorial Pacific outgassing (weak upwelling) and regional offsets in the Southern Ocean.

### 2.2 Sea surface temperature and salinity — K₀ and Sc drivers

| | |
|---|---|
| **Product** | CMEMS Global Physical Reanalysis (GLORYS12V1) |
| **ID** | `GLOBAL_MULTIYEAR_PHY_001_030` |
| **Variables** | `thetao` (potential temperature, °C), `so` (practical salinity, PSU) |
| **Resolution** | 1/12°, monthly → regridded to 0.25° in `02_preprocess.py` |
| **Nature** | Physical reanalysis with data assimilation (altimetry, SST, Argo) |

**Important**: `GLOBAL_MULTIYEAR_BGC_001_029` was forced by FREEGLORYS2V4 and ERA-Interim atmospheric fields. The SST/SSS used here come from GLORYS12V1, which is forced by ERA5. These are therefore **not dynamically consistent** — the pCO₂ field was generated under slightly different physical forcing than the SST/SSS used to compute K₀. This is a known limitation inherent to combining these products.

### 2.3 Wind speed — gas transfer velocity

| | |
|---|---|
| **Product** | ERA5 monthly averaged reanalysis |
| **Dataset** | `reanalysis-era5-single-levels-monthly-means` |
| **Variables** | `u10`, `v10` (10 m wind components) → scalar speed $\|u\| = \sqrt{u_{10}^2 + v_{10}^2}$ |
| **Resolution** | 0.25°, monthly |

| | |
|---|---|
| **Product** | ERA5 daily reanalysis (optional, for wind variance correction) |
| **Dataset** | `reanalysis-era5-single-levels` |
| **Variables** | `u10`, `v10` |
| **Resolution** | 0.25°, daily (noon snapshot as proxy for daily mean) |
| **Output** | `era5_wind10m_daily.nc` — ~30–50 GB for full 1993–2025 record |

The daily wind data is required only for the wind variance correction σ²_u. Without it, the reconstruction still runs using monthly-mean winds only.

### 2.4 Atmospheric CO₂

| | |
|---|---|
| **Source** | NOAA Global Monitoring Laboratory (GML) |
| **File** | `co2_mm_gl.csv` — global marine surface monthly mean |
| **Variable** | `average` [ppm] → converted to [atm] (×10⁻⁶) |
| **URL** | https://gml.noaa.gov/webdata/ccgg/trends/co2/co2_mm_gl.csv |
| **Note** | 1 ppm ≈ 1 µatm at standard total atmospheric pressure (1 atm) |

Used both for the atmospheric end of the ΔpCO₂ calculation (after spatial broadcast) and for plotting and saturation analysis.

### 2.5 Validation product — MULTIOBS

| | |
|---|---|
| **Product** | CMEMS MULTIOBS SOCAT Neural Network L4 |
| **ID** | `MULTIOBS_GLO_BGC_CARBON_SURFACE_MYNRT_015_008` |
| **Nature** | SOCAT-trained neural network — observation-constrained |
| **Variables** | `fgco2` (flux), `spco2` (surface pCO₂) |
| **Use** | Validation only — never enters the reconstruction |

This product is used strictly as a reference to quantify reconstruction skill (RMSD, bias, correlation, spectral coherence). It is NOT used as a primary input at any stage.

### 2.6 Biogeochemical provinces — Fay domains

| | |
|---|---|
| **Source** | Fay & McKinley (2014) |
| **File** | `data/Time_Varying_Biomes.cmems.nc` |
| **Variable** | `MeanBiomes` — integer domain ID (1–17) per pixel |
| **Basis** | SST, chlorophyll, MLD seasonality clustering |

17 open-ocean biomes with distinct carbon-cycle regimes. Used in `06_fay_analysis.py` to disaggregate the global flux into province-level timeseries, climatologies, and trend attribution.

### 2.7 Climate indices

| Index | Source | File | Format |
|-------|--------|------|--------|
| ONI | NOAA PSL — ERSSTv5 Niño-3.4 3-month running mean | `data/oni_monthly.txt` | year Jan…Dec, missing = −99.9 |
| SAM | Marshall (2003), NOAA CPC update | `data/sam_monthly.txt` | year Jan…Dec |

Both files are auto-downloaded on first run from internet-connected nodes. On HPC compute nodes without internet, download manually:
```bash
wget -O data/oni_monthly.txt https://psl.noaa.gov/data/correlation/oni.data
wget -O data/sam_monthly.txt https://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/aao/monthly.aao.index.b79.current.ascii.table
```

---

## 3. Complete mathematical formulation

### 3.1 Bulk air-sea CO₂ flux

$$\boxed{F = k \cdot K_0 \cdot \left( pCO_2^{\,\mathrm{atm}} - pCO_2^{\,\mathrm{ocn}} \right)}$$

| Symbol | Definition | Units |
|--------|-----------|-------|
| $F$ | Net CO₂ flux | mol C m⁻² s⁻¹ → converted to mol C m⁻² yr⁻¹ |
| $k$ | Gas transfer velocity | m s⁻¹ |
| $K_0$ | CO₂ solubility | mol m⁻³ atm⁻¹ |
| $pCO_2^{\mathrm{atm}}$ | Atmospheric pCO₂ (spatially uniform per time step) | atm |
| $pCO_2^{\mathrm{ocn}}$ | Surface ocean pCO₂ from BGC hindcast | atm |

**Sign convention:** $F > 0$ = ocean uptake (flux from atmosphere into ocean, sink); $F < 0$ = outgassing (flux to atmosphere, source).

This convention is consistent with CMEMS MULTIOBS and the GCB (Global Carbon Budget) reporting convention.

**Unit conversion:** F is computed in mol m⁻² s⁻¹ (natural SI units from k [m s⁻¹] × K₀ [mol m⁻³ atm⁻¹] × ΔpCO₂ [atm]) and converted to mol m⁻² yr⁻¹ by multiplying by `cfg.S_TO_YR` = 3.1557 × 10⁷ s yr⁻¹.

---

### 3.2 Gas transfer velocity — Wanninkhof (2014)

The Wanninkhof (2014) quadratic wind parameterisation:

$$k = a \cdot \langle u^2 \rangle \cdot \left(\frac{Sc}{660}\right)^{-1/2}$$

where $\langle u^2 \rangle$ is the mean square wind speed. When only monthly-mean winds are available:

$$\langle u^2 \rangle \approx \langle u \rangle^2$$

When ERA5 daily winds are available, the **wind variance correction** is applied:

$$\langle u^2 \rangle = \langle u \rangle^2 + \sigma_u^2$$

where $\sigma_u^2$ is the sub-monthly variance of daily wind speed computed per grid cell per month. This correction arises because Wanninkhof's coefficient $a$ was calibrated against the full wind speed distribution using bomb-¹⁴C and ³He/SF₆ tracer techniques, not against monthly means. Since $\langle u^2 \rangle \geq \langle u \rangle^2$ always, using monthly means alone systematically **underestimates** $k$ — particularly in high-wind, high-variability regions such as the Southern Ocean and North Atlantic.

| Symbol | Value / definition |
|--------|--------------------|
| $a$ | 0.251 cm hr⁻¹ (m s⁻¹)⁻² — Wanninkhof (2014) Table 1 |
| $\langle u_{10} \rangle$ | Monthly-mean 10 m scalar wind speed [m s⁻¹] from ERA5 |
| $\sigma_u^2$ | Sub-monthly variance of $\|u_{10}\|$ [m² s⁻²] from ERA5 daily (optional) |
| $Sc$ | Schmidt number of CO₂ in seawater (dimensionless) |
| 660 | Reference Schmidt number of CO₂ at 20 °C in seawater |

The exponent −1/2 arises from the assumption that gas transfer scales as $Sc^{-1/2}$ in the turbulent/rough-surface regime.

**Unit conversion of $k$:** computed in cm hr⁻¹, converted to m s⁻¹ by `cfg.CMHR_TO_MS` = 1/(100 × 3600).

---

### 3.3 Schmidt number — Wanninkhof (2014), Table 1

The Schmidt number Sc is the ratio of the kinematic viscosity of seawater to the diffusivity of CO₂. It depends strongly on temperature:

$$Sc = A - B \cdot T + C \cdot T^2 - D \cdot T^3 + E \cdot T^4 \qquad (T \text{ in °C})$$

Coefficients from Wanninkhof (2014), Table 1, for CO₂ in seawater:

| $A$ | $B$ | $C$ | $D$ | $E$ |
|-----|-----|-----|-----|-----|
| 2116.8 | 136.25 | 4.7353 | 0.092307 | 0.0007555 |

Sc decreases with temperature (CO₂ diffuses faster in warm water), which partially counteracts the decrease in K₀ with warming — the net effect on flux depends on the balance between reduced solubility and increased k.

---

### 3.4 CO₂ solubility — Weiss (1974)

$$\ln K_0 = A_1 + A_2 \cdot \frac{100}{T} + A_3 \cdot \ln\!\left(\frac{T}{100}\right) + S \cdot \left[ B_1 + B_2 \cdot \frac{T}{100} + B_3 \cdot \left(\frac{T}{100}\right)^2 \right]$$

$T$ in Kelvin, $S$ in PSU (practical salinity units). The formula returns K₀ in mol L⁻¹ atm⁻¹; multiplied by 1000 to convert to mol m⁻³ atm⁻¹.

Coefficients from Weiss (1974), Equation 12:

| $A_1$ | $A_2$ | $A_3$ | $B_1$ | $B_2$ | $B_3$ |
|-------|-------|-------|-------|-------|-------|
| −58.0931 | 90.5069 | 22.2940 | 0.027766 | −0.025888 | 0.0050578 |

K₀ decreases strongly with increasing temperature (CO₂ is less soluble in warm water). This is the primary physical driver of the positive ΔpCO₂ signal in the tropics and the reduced uptake expected under ocean warming. The salinity dependence is secondary but not negligible — higher salinity slightly reduces solubility (salting-out effect).

---

### 3.5 Global surface integral

The global net ocean CO₂ uptake at time $t$:

$$J_{\mathrm{net}}(t) = \sum_{i,j} F(t, i, j) \cdot A(i, j)$$

Grid cell area at latitude $\varphi$ and longitude $\lambda$:

$$A(i,j) = R^2 \cdot \Delta\varphi_\mathrm{rad} \cdot \Delta\lambda_\mathrm{rad} \cdot \cos\varphi$$

where $R = 6.371 \times 10^6$ m (Earth radius), $\Delta\varphi = \Delta\lambda = 0.25°$ for the BGC hindcast grid.

Unit conversion from mol C yr⁻¹ to Pg C yr⁻¹:

$$J_{\mathrm{net}} \,[\text{Pg C yr}^{-1}] = J_{\mathrm{net}} \,[\text{mol C yr}^{-1}] \times \frac{12.011 \,\text{g mol}^{-1}}{10^{15} \,\text{g Pg}^{-1}}$$

Land pixels (where `ocean_mask == 0`) are excluded via NaN masking before summation. `skipna=True` is used in the xarray sum to handle any remaining isolated NaN pixels (e.g. under sea ice) without biasing the integral.

---

### 3.6 Wind variance correction — physical basis

For any random variable $u$:

$$\langle u^2 \rangle = \langle u \rangle^2 + \mathrm{Var}(u) = \langle u \rangle^2 + \sigma_u^2$$

Since $k \propto \langle u^2 \rangle$, using the monthly mean $\langle u \rangle$ in place of $\langle u^2 \rangle$ underestimates k by exactly $a \cdot \sigma_u^2 \cdot (Sc/660)^{-1/2}$. In monthly ERA5 data, $\sigma_u^2$ is typically 1–5 m² s⁻² depending on region and season, corresponding to a k correction of 5–20% in high-wind regions. This matters most in the Southern Ocean (mean wind speeds 8–12 m s⁻¹, high variance) where the underestimate would be largest.

The variance σ²_u is computed from daily ERA5 winds as:

$$\sigma_u^2(t_\mathrm{month}, i, j) = \mathrm{Var}_{d \in \mathrm{month}}\!\left(\sqrt{u_{10,d}^2 + v_{10,d}^2}\right)$$

using `resample("1ME").var()` in xarray, which applies Bessel's correction (n−1 denominator) automatically.

---

## 4. Statistical methods

### 4.1 Sen's slope (Theil-Sen estimator)

Per-pixel trend estimation. The slope is the median of all pairwise slopes between data points:

$$\hat{\beta} = \mathrm{median}\!\left(\frac{y_j - y_i}{t_j - t_i}\right), \quad \forall \; i < j$$

For $n$ data points this involves $n(n-1)/2$ pairwise slopes. Implementation via `scipy.stats.theilslopes`.

**Why Sen's slope and not OLS regression?**
OLS minimises the sum of squared residuals, making it sensitive to outliers — a single anomalous year (e.g. a strong El Niño) can pull the OLS slope substantially. Sen's slope uses the median and is therefore robust to outliers and does not assume normality of residuals. It is also the natural companion to the Mann-Kendall test (both are rank-based), ensuring consistency between the trend estimate and its significance test.

Slopes reported in **mol C m⁻² yr⁻¹ per decade** (multiplied by 10 after computing the per-year slope from annual means).

---

### 4.2 Mann-Kendall significance test

For a timeseries $y_1, y_2, \ldots, y_n$, the MK statistic:

$$S = \sum_{i < j} \mathrm{sgn}(y_j - y_i)$$

where sgn = +1 if $y_j > y_i$, −1 if $y_j < y_i$, 0 if equal. Under the null hypothesis of no trend, $S$ has mean 0 and a known variance depending on $n$ and the number of ties.

The standardised test statistic:

$$Z = \begin{cases} (S-1)/\sqrt{\mathrm{Var}(S)} & S > 0 \\ 0 & S = 0 \\ (S+1)/\sqrt{\mathrm{Var}(S)} & S < 0 \end{cases}$$

The two-sided p-value is $p = 2 \cdot (1 - \Phi(|Z|))$ where $\Phi$ is the standard normal CDF.

**Implementation:** computed via `scipy.stats.kendalltau(x_indices, y)` against a monotone integer sequence $x = [0, 1, \ldots, n-1]$. This is mathematically equivalent to the standard MK test — Kendall's tau is a normalised version of $S$, and `kendalltau` returns the associated two-sided p-value directly.

**Why α = 0.01 (not the conventional 0.05)?**  
At α = 0.05 on a global 0.25° ocean grid (~500,000 ocean pixels), the expected number of false positives by chance is ~25,000 pixels — enough to produce misleading spatial patterns in the trend map. Setting α = 0.01 reduces this to ~5,000. The trade-off is reduced power (some genuinely weak trends, particularly in high-variability tropical regions, will be missed). Given the 30-year record with strong ENSO-driven interannual noise, the stricter threshold is the more defensible scientific choice.

**Why annual means?**  
Monthly timeseries have strong autocorrelation at lag 1 (seasonal cycle) and lag 12. Autocorrelation inflates the effective sample size, making p-values anti-conservative — trends appear significant when they are not. Working on annual means removes the seasonal autocorrelation and reduces lag-1 autocorrelation substantially, making the MK p-values more reliable without requiring pre-whitening.

---

### 4.3 LOESS smoothing

LOESS (Locally Weighted Scatterplot Smoothing, Cleveland 1979). At each point $x_0$, fits a polynomial using only the $k$ nearest neighbours (fraction `frac = 0.15` of the data), weighted by a tricubic kernel:

$$w_i = \left(1 - \left|\frac{x_i - x_0}{h}\right|^3\right)^3$$

where $h$ is the bandwidth (distance to the `frac`-th nearest neighbour). Implementation: `statsmodels.nonparametric.smoothers_lowess.lowess` with `frac=0.15`, `it=3` (3 robustness iterations using bisquare weights to down-weight outliers).

**Why LOESS and not a running mean or harmonic filter?**  
A 12-month running mean removes the seasonal cycle but introduces edge effects and cannot adapt to non-stationary variance. A harmonic (Fourier) filter assumes stationarity. LOESS makes no such assumptions — it adapts locally and handles irregular sampling and non-stationarity naturally. The `frac=0.15` corresponds to roughly an 18-month window on a 30-year monthly record, which suppresses the annual and semi-annual harmonics while preserving interannual signals.

---

### 4.4 PELT breakpoint detection

PELT (Pruned Exact Linear Time, Killick et al. 2012). Given a signal $y_1, \ldots, y_n$, finds the partition into $m$ segments that minimises:

$$\mathcal{C}(y) = \sum_{k=1}^{m} \left[ \mathrm{cost}(y_{\tau_{k-1}+1:\tau_k}) \right] + \mathrm{pen} \cdot m$$

where $\tau_0 = 0 < \tau_1 < \ldots < \tau_m = n$ are the breakpoint positions and pen = 3 is the penalty per breakpoint.

**Cost function (RBF):** for a segment, the RBF cost measures how much the signal deviates from the segment mean — effectively a sum of squared distances from the local mean. A flat, stationary segment has low cost; a segment that is rising, falling, or oscillating has high cost.

**The penalty** `pen = 3` is the price each breakpoint must pay. A breakpoint is placed only if it reduces the total segment cost by more than 3 units. Lower penalty → more breakpoints (risk of detecting noise); higher penalty → fewer breakpoints (risk of missing real regime shifts). `pen = 3` is a moderate setting empirically suited to 30-year monthly ocean timeseries.

**Why PELT and not simpler methods?** PELT is exact — it finds the globally optimal partition for the given penalty, not just a locally optimal one. Its pruning rule allows it to run in O(n) time in practice (linear in the length of the series), making it feasible on the LOESS smooth without approximation.

**B1 and B4:** PELT typically detects 4–7 breakpoints on a 30-year record. We retain only B1 (the first) and B4 (the fourth) to define three interpretable multi-decadal segments. B1 and B4 tend to correspond to the transitions entering and exiting the mid-record variability regime. If fewer than 4 breakpoints are detected, the first and last are used instead.

---

### 4.5 Ocean sink saturation analysis

Annual J_net regressed against annual atmospheric CO₂ (OLS):

$$J_{\mathrm{net}}^{\mathrm{annual}} = \beta_0 + \beta_1 \cdot [\mathrm{CO}_2]_{\mathrm{atm}}^{\mathrm{annual}} + \varepsilon$$

The slope $\beta_1$ [Pg C yr⁻¹ ppm⁻¹] is the **ocean uptake sensitivity** to atmospheric CO₂. Interpretation:

- $\beta_1 > 0$, significant → the ocean sink is growing in proportion to rising CO₂ (no saturation)
- $\beta_1 \approx 0$ or p > 0.05 → the sink is stagnating relative to CO₂ (partial saturation or decoupling)
- $\beta_1 < 0$ → the sink is weakening despite rising CO₂ (full saturation or reversal)

Points are coloured by year to reveal whether the relationship has been stable or has drifted over the record. A temporal drift (e.g., high-β₁ years clustering in the early record, low-β₁ years in the recent record) would suggest a changing sink efficiency independent of any trend in the annual mean.

**Caveats:** This is a diagnostic, not a mechanistic test. Confounders include SST-driven changes in K₀ (warming reduces solubility, partially opposing the uptake increase from rising ΔpCO₂), biological feedbacks, and ENSO-driven interannual variability. The regression should be interpreted alongside the LOESS/breakpoint analysis and the climate index regression.

---

### 4.6 Climate index regression — Bellacicco et al. (2025)

Spearman rank correlation of annual J_net against ONI and SAM, separately and as a multilinear combination.

**Spearman rank correlation:**

$$\rho = 1 - \frac{6 \sum d_i^2}{n(n^2 - 1)}$$

where $d_i$ is the difference in ranks of the $i$-th pair. Two-sided p-value from the t-distribution with $n-2$ degrees of freedom.

**Why Spearman and not Pearson?**
Three reasons: (1) With only 31 annual data points, Pearson is sensitive to normality assumptions and to individual extreme years. (2) The ENSO–flux relationship may be non-linear — a strong El Niño can have disproportionate effects. Spearman captures any monotone relationship. (3) Direct comparability with Bellacicco et al. (2025), who use Spearman — their reported values (ρ_ONI = 0.57, ρ_combined = 0.617) are directly comparable to ours.

**Multilinear combination (OLS):**

$$\hat{J} = \alpha \cdot \tilde{\mathrm{ONI}} + \beta \cdot \tilde{\mathrm{SAM}} + \gamma$$

where $\tilde{\cdot}$ denotes z-score normalisation (zero mean, unit variance). The Spearman ρ between $J_\mathrm{net}$ and $\hat{J}$ quantifies how well the combined predictor captures the observed flux variability. OLS coefficients α and β give the relative contribution of each index.

**Per-domain (Fay) analysis:** the same regression is performed for each of the 17 Fay domains using area-integrated domain flux [Pg C yr⁻¹]. This allows identification of which provinces drive the global ENSO/SAM signal, enabling direct comparison with Bellacicco et al.'s finding that Southern Ocean domains dominate the PIP interannual variability.

---

### 4.7 EOF analysis — area-weighted SVD

Empirical Orthogonal Function decomposition of the annual-mean surface flux
**full field** (no mean removal). Annual means are used to suppress the seasonal
cycle without harmonic assumptions; the full field (rather than anomalies) is
retained because evaluation showed the anomaly and monthly variants did not
produce physically distinct results for this dataset — the dominant spatial
pattern aligns closely with the climatological flux structure and is scientifically
meaningful in that form.

**Pre-processing:**
1. Annual means of `fgco2` from `flux_3d.nc`, sliced to 1993–2023
2. **No time-mean removal** — full field decomposition
3. Weight each pixel by $\sqrt{\cos\varphi}$ — area-weighting so high-latitude pixels (which cover less ocean area) do not dominate the variance decomposition
4. Build ocean mask from NaN pattern of the time-mean; flatten to $\mathbf{X}$ of shape (n_years × n_ocean_pixels) after masking land

**SVD decomposition:**

$$\mathbf{X}_w = \mathbf{U} \cdot \boldsymbol{\Sigma} \cdot \mathbf{V}^T$$

where $\mathbf{X}_w$ is the area-weighted matrix.

**Scaling convention (unit-norm EOFs — standard oceanographic):**

1. Unweight: $\mathrm{EOF}_k^\mathrm{unit} = \mathbf{V}^T_k / w$ — unit L2-norm in physical space
2. Project: $\mathrm{PC}_k^\mathrm{raw} = \mathbf{X}_\mathrm{unweighted} \cdot \mathrm{EOF}_k^\mathrm{unit}$
3. Normalise: $\mathrm{PC}_k^\mathrm{norm} = \mathrm{PC}_k^\mathrm{raw} / \sigma_k$ — unit variance, dimensionless
4. Scale: $\mathrm{EOF}_k^\mathrm{scaled} = \mathrm{EOF}_k^\mathrm{unit} \times \sigma_k$ — physical units [mol C m⁻² yr⁻¹]

Result: $\mathrm{EOF}_k^\mathrm{scaled} \times \mathrm{PC}_k^\mathrm{norm}$ exactly reconstructs the field. EOF maps show the typical flux pattern per one standard deviation of the PC. This matches the `eofs` package convention.

**Variance explained:**

$$f_k = \frac{\sigma_k^2}{\sum_j \sigma_j^2}$$

**North et al. (1982) sampling error:**

$$\delta\lambda_k = \lambda_k \cdot \sqrt{\frac{2}{N}}$$

where $N$ = number of annual time steps (31 for the 1993–2023 record). The North error is approximately 25% per eigenvalue at this sample size — degeneracy between neighbouring modes is common beyond EOF2–3 and should be checked on the scree plot before interpreting modes individually.

**PC–climate index correlation:** PC timeseries are correlated with z-score ONI and SAM (Spearman ρ) to attribute each mode to known climate drivers. A strong correlation between PC1 and ONI would indicate that ENSO modulates the leading mode of surface flux spatial variability.

**Deprecated analyses (v1.4.3):** monthly full-field (dominated by seasonal cycle), monthly anomaly (noise-dominated at 31-year record length), and annual anomaly (structurally redundant with annual full-field for this dataset) were evaluated and removed.

---

## 5. File structure and design principles

### 5.1 Directory layout

```
stage1/
├── config.py                         ← single source of truth for all constants
├── data/
│   ├── bgc_hindcast_spco2.nc
│   ├── phy_reanalysis_sst_sss.nc
│   ├── multiobs_surface_carbon.nc
│   ├── era5_wind10m_monthly.nc
│   ├── era5_wind10m_daily.nc         ← optional (wind variance correction)
│   ├── co2_mm_gl.csv
│   ├── Time_Varying_Biomes.cmems.nc  ← required for Fay analysis
│   ├── oni_monthly.txt               ← auto-downloaded or manual wget
│   ├── sam_monthly.txt               ← auto-downloaded or manual wget
│   ├── processed_surface.nc          ← output of 02_preprocess.py
│   └── flux_3d.nc                    ← output of 03_compute_flux.py
├── output/
│   ├── global_flux.nc
│   ├── validation_metrics.csv
│   └── figures/
│       └── [all figure files]
└── scripts/
    ├── 01_download_data.py
    ├── 02_preprocess.py
    ├── 03_compute_flux.py
    ├── 04_validate.py
    ├── 05_plot_results.py
    ├── 06_fay_analysis.py
    ├── 07_eof_analysis.py
    └── plot_style.py
```

### 5.2 Design principles

**`config.py` is the single source of truth.** All product IDs, variable names, physical constants, file paths, unit conversion factors, and colourmap choices live in `config.py`. Scripts import `config as cfg` and reference `cfg.WANNINKHOF_A`, `cfg.WEISS_A1`, `cfg.BBOX`, `cfg.DATA_DIR`, etc. Nothing is hardcoded inside script bodies except figure-specific aesthetic parameters.

**Scripts are stateless and idempotent.** Every script checks if its output files already exist before running. Delete an output file to force recomputation. Re-running a script that has already succeeded is a no-op for all file-creation operations.

**No physics computation in plotting scripts.** `05_plot_results.py`, `06_fay_analysis.py`, and `07_eof_analysis.py` read pre-computed NetCDF files only. The exception is lightweight statistical operations applied to the flux arrays for plotting purposes (Sen's slope, Mann-Kendall, OLS regression, SVD). No gas exchange parameterisation or physical unit conversion happens in these scripts.

**Pure functions in `03_compute_flux.py`.** Every physical parameterisation (`schmidt_number_co2`, `gas_transfer_velocity`, `co2_solubility_K0`, `compute_flux`) takes `xr.DataArray` inputs and returns `xr.DataArray` outputs with full metadata (long_name, units, reference). These functions have no side effects and can be called independently for testing.

**Dask-backed lazy loading.** All `xr.open_dataset` calls use `chunks="auto"`. Computation is deferred until an explicit `.compute()` or `.to_netcdf()` call. This allows the scripts to handle datasets larger than available RAM by streaming chunks.

**Optional dependencies degrade gracefully.** `cartopy` (map projections), `cmocean` (colormaps), `statsmodels` (LOESS), and `ruptures` (PELT breakpoints) are all optional. The scripts detect their absence at import time and fall back to plain axes, `RdBu_r`, skipping LOESS/breakpoints, respectively. A clear `[warn]` or `[skip]` message is printed; the script never hard-crashes due to a missing optional package.

**Climate index downloads degrade gracefully on HPC.** `load_climate_indices()` wraps all network calls in `try/except`. On failure it prints `wget` commands and returns `None`. All climate figures are skipped when `None` is returned; the rest of the pipeline continues.

---

## 6. Script pipeline — detailed description

### `01_download_data.py`

Downloads all Stage 1 data from CMEMS (via the `copernicusmarine` Python client), NOAA GML (via `requests`), and ERA5 (via `cdsapi`). All functions check for existing files before downloading. Supports `--test` flag to download 2010 only for a quick sanity check.

Functions: `download_bgc_hindcast`, `download_physical_reanalysis`, `download_multiobs_surface`, `download_noaa_co2`, `download_era5_wind` (monthly), `download_era5_wind_daily` (yearly chunks, merged).

### `02_preprocess.py`

Harmonises all downloaded data onto a common 0.25° monthly grid. Key operations:

- Regrid GLORYS12 SST/SSS from 1/12° to 0.25° by bilinear interpolation
- Convert BGC pCO₂ from Pa to atm
- Broadcast NOAA atmospheric CO₂ (scalar per month) to the full 2D grid
- Compute wind speed magnitude from ERA5 u10/v10 components
- Compute wind variance σ²_u from daily ERA5 (if available)
- Build ocean mask from PISCES pCO₂ NaN pattern (NaN = land or sea ice)
- Time-align all fields to a common monthly axis
- Output: `data/processed_surface.nc` containing: `sst`, `sss`, `wind_speed`, `wind_variance` (optional), `spco2_ocean`, `spco2_atm`, `ocean_mask`

### `03_compute_flux.py`

Core physics. Reads `processed_surface.nc`, applies the Wanninkhof/Weiss parameterisation, computes the global integral. Outputs:

- `data/flux_3d.nc`: `fgco2` [mol C m⁻² yr⁻¹], `k` [m s⁻¹], `K0` [mol m⁻³ atm⁻¹], `Sc` [dimensionless]
- `output/global_flux.nc`: `J_net_PgC` [Pg C yr⁻¹] (monthly), `J_net_PgC_annual` (annual resample)

### `04_validate.py`

Cross-validation of the reconstruction against MULTIOBS. Produces 7 validation figures covering global timeseries skill, rolling metrics, spatial RMSD/bias, Fay domain skill, and spectral properties.

### `05_plot_results.py`

Primary result figures (7 figures). Global-scale analysis: flux evolution, maps, trends, seasonal cycle, saturation, climate indices. Time window applied in `load_all()` so all figures use the same 1993–2023 slice automatically.

### `06_fay_analysis.py`

Fay domain analysis (4 figures). Separated from `05_plot_results.py` because the 17-domain area-integration loop (17 × n_time × n_lat × n_lon weighted sums) is computationally heavy. Run independently after `03_compute_flux.py`.

### `07_eof_analysis.py`

EOF decomposition (7 figures: 1 scree + 3 spatial + 3 PC). Self-contained — loads `flux_3d.nc` directly, performs area-weighted SVD, produces spatial patterns and PC timeseries. PC plots optionally overlay ONI/SAM if index files are present.

---

## 7. Output figures — detailed description

### `05_plot_results.py`

**`fig01_flux_vs_co2.png`**  
Twin-axis: left = J_net [Pg C yr⁻¹] (raw monthly faint + LOESS bold + B1/B4 breakpoints + per-segment OLS trend slopes annotated above trend lines); right = atmospheric CO₂ [ppm] (NOAA GML annual mean, dashed red, slope annotated below). B1/B4 labels on top twin x-axis. y-axis fixed at −0.2 → 3.0 Pg C yr⁻¹. The central figure of the analysis — directly addresses the saturation question visually.

**`fig02_mean_flux_map.png`**  
Time-mean air-sea CO₂ flux [mol C m⁻² yr⁻¹], Robinson projection. White-centred diverging colourmap (`cmocean.balance`). Fixed range ±2.5 mol C m⁻² yr⁻¹. Blue = uptake, red = outgassing.

**`fig03_trend_significance_map.png`**  
Per-pixel Sen's slope of annual mean flux [mol C m⁻² yr⁻¹ per decade]. Fixed range ±0.5 per decade. Mann-Kendall significance at α = 0.01: non-significant ocean pixels hatched (///); a gray contour outlines the boundary between significant and non-significant regions; land drawn on top (zorder > hatching) so hatching never appears on land.

**`fig04_delta_pco2_map.png`**  
Time-mean ΔpCO₂ = pCO₂(ocean) − pCO₂(atm) [µatm]. Range ±100 µatm. Shows the thermodynamic driver — positive (red) regions are net sources, negative (blue) regions are net sinks on thermodynamic grounds. Note that flux also depends on k, so a region with negative ΔpCO₂ in low-wind conditions may have small actual flux.

**`fig05_seasonal_cycle.png`**  
Climatological monthly mean (1993–2023): J_net as bars (left axis), atmospheric CO₂ as a line (right axis). Reveals the seasonal asymmetry driven by the biological pump (Northern Hemisphere spring bloom draws down pCO₂) and the hemispheric asymmetry (Southern Ocean summer uptake peaks December–February).

**`fig06_sink_saturation.png`**  
Annual J_net vs annual atmospheric CO₂ scatter, points coloured by year (plasma colourmap). OLS regression line with β₁ [Pg C yr⁻¹ ppm⁻¹], Pearson r, and p-value annotated. Temporal colour gradient shows whether the relationship has been stable or drifting.

**`fig07_climate_regression_global.png`**  
Three-panel climate index regression (auto-skipped on HPC if index files absent):
- (a) Annual J_net bars + z-score ONI (red dashed) and SAM (green dash-dot) on twin axis
- (b) Scatter J_net vs ONI (circles) and vs SAM (triangles) with separate OLS lines, Spearman ρ and p annotated
- (c) Scatter J_net vs OLS multilinear α·ONI + β·SAM, coloured by year, Spearman ρ annotated

### `06_fay_analysis.py`

**`fig_fay_ts.png`**  
17-panel 6×3 grid: monthly area-integrated flux [Pg C yr⁻¹] per Fay domain. Each domain coloured by its jet-colourmap index (matching the 18th biome overview map panel). Area-integrated (not area-mean) so values can be summed to recover J_net.

**`fig_fay_clim.png`**  
17-panel: climatological monthly mean ±1σ (across years) per domain. Reveals biological vs thermodynamic seasonal controls per province. Filled band is ±1σ, not error on the mean.

**`fig_fay_trends.png`**  
17-panel: Sen's slope bar per domain [Pg C yr⁻¹ per decade]. Mann-Kendall significance at α = 0.01: filled = significant, hatched (///) = not. ✱ on significant bars. Identifies which provinces drive the global trend.

**`fig_fay_climate_regression.png`**  
17-panel: three bars per domain — ρ(J_net, ONI) red, ρ(J_net, SAM) green, ρ(J_net, α·ONI+β·SAM) blue. Reference lines at ±0.4. Filled = p ≤ 0.05, hatched = p > 0.05. ✱ on significant bars. Auto-skipped if index files absent.

### `07_eof_analysis.py`

Active analysis: **annual full field** (v1.4.3). Two figures produced.

**`fig_eof_annual_full_scree.png`**  
Scree plot for the first 10 EOFs. Bars = variance explained [%] with North et al. (1982) sampling error bars. Line = cumulative variance on right axis with 90% reference. Red bar outlines mark degenerate mode pairs (error bars of neighbouring eigenvalues overlap — those modes cannot be individually interpreted as distinct physical patterns and should not be over-interpreted).

**`fig_eof_annual_full_modes.png`**  
3 rows × 2 columns composite figure for EOF/PC modes 1–3:
- *Left column*: spatial pattern of EOF N [mol C m⁻² yr⁻¹]. White-centred diverging colourmap (`cmocean.balance`). Colourbar scaled to 97th percentile of |EOF| per mode (falls back to ±1.0 if map is entirely NaN). Robinson projection, land drawn on top (zorder 3/4).
- *Right column*: PC N timeseries (unit variance, dimensionless) as filled line. If `data/oni_monthly.txt` and `data/sam_monthly.txt` are present: z-score ONI (red dashed) and SAM (green dash-dot) overlaid; Spearman ρ and p-value annotated in legend. Auto-skipped gracefully if index files absent.

*Deprecated figures (v1.4.3, no longer produced):* `fig_eof_monthly_full_*.png`, `fig_eof_monthly_anom_*.png`, `fig_eof_annual_anom_*.png`.

### `04_validate.py`

**`fig_validation_ts.png`**  
2-row: (a) monthly global J_net — reconstruction (blue) vs MULTIOBS (red); (b) 12-month rolling RMSD (pink) and rolling Pearson r (green) on twin axes.

**`fig_validation_loess.png`**  
LOESS smooth (frac=0.15) of both reconstruction and MULTIOBS. PELT breakpoints detected on reconstruction LOESS; same breakpoints applied to both curves. Per-segment OLS slopes annotated: reconstruction below trend line, MULTIOBS above. Segment colours consistent across both curves.

**`fig_validation_map.png`**  
2-panel Robinson projection: pixel-wise RMSD (cmocean.amp) and bias = reconstruction − MULTIOBS (cmocean.balance). Fixed colour limits from config.

**`fig_fay_ts.png`** / **`fig_fay_clim.png`**  
Same structure as `06_fay_analysis.py` equivalents but with MULTIOBS overlaid. Allows direct per-domain skill assessment of the reconstruction.

**`fig_spectra.png`**  
Welch power spectral density (nperseg = n/3) of global J_net — reconstruction vs MULTIOBS. Log-log axes, period on x-axis. Nyquist period set as rightmost tick.

**`fig_cps.png`**  
Cross-power spectrum: (a) gain, (b) phase [days], (c) magnitude-squared coherence. `scipy.signal.coherence` with nperseg = n/4. Coherence threshold at 0.5 (red). Gain threshold at 1 (red). Phase reference at 0 (red).

---

## 8. Known limitations and caveats

1. **Physical consistency:** `GLOBAL_MULTIYEAR_BGC_001_029` was forced by FREEGLORYS2V4/ERA-Interim (not ERA5). The SST and SSS used here come from GLORYS12V1 (ERA5-forced). These are not dynamically consistent — the pCO₂ field was generated under slightly different physical forcing than the physical fields used to compute K₀ and Sc. The magnitude of this inconsistency is unknown but likely small for SST-driven K₀ errors.

2. **Wind variance correction:** Active only when `era5_wind10m_daily.nc` is present (~30–50 GB). Without it, k is computed from monthly-mean ⟨u⟩² only, systematically underestimating k in high-wind, high-variability regions (Southern Ocean, North Atlantic storm track). The bias in the global integral is estimated at 5–10%.

3. **Sea-ice masking:** Uses the PISCES pCO₂ NaN pattern as a proxy for ice-covered pixels rather than an explicit sea-ice concentration field. Some marginal ice zone pixels may be incorrectly included or excluded, affecting high-latitude flux estimates.

4. **Hindcast, not reanalysis:** `GLOBAL_MULTIYEAR_BGC_001_029` assimilates no BGC observations. Known model biases include underestimation of equatorial Pacific outgassing (weak Ekman upwelling in the model) and regional offsets in Southern Ocean pCO₂ of order 10–20 µatm.

5. **Autocorrelation in trend significance:** Mann-Kendall significance is computed on annual means to reduce serial autocorrelation. Some residual autocorrelation may remain, particularly in regions with strong decadal variability (e.g. the Southern Ocean). For a fully rigorous analysis, effective sample size correction or block bootstrap would be preferred.

6. **Saturation analysis:** The β₁ regression is diagnostic only. A flat relationship could reflect saturation, or it could reflect compensation between opposing mechanisms (e.g., warming reducing K₀ while rising ΔpCO₂ increases flux). Disentangling these requires the thermal/non-thermal decomposition planned for Stage 6.

7. **EOF degeneracy:** Modes where North et al. (1982) error bars overlap cannot be individually interpreted. In a 30-year record (n=31 annual means), the North error is large (~25% for each eigenvalue), making degeneracy common beyond EOF2–3. Only EOF1 is reliably interpretable in isolation.

8. **Record length:** The 1993–2023 window (~30 years) captures the satellite era and is sufficient to characterise interannual variability and the recent trend. However, carbon cycle timescales span centuries to millennia. Results represent the current flux regime and should not be extrapolated to centennial projections.

9. **MULTIOBS sign convention:** Verify the sign of `fgco2` in the MULTIOBS product against the current product QUID before updating to a new MULTIOBS dataset version.

10. **HPC compute nodes:** `load_climate_indices()` fails gracefully when outbound internet is blocked. Climate figures are skipped without aborting the run. See Section 2.7 for manual download instructions.

---

## 9. Extending the project

This codebase is structured to grow stage by stage, with each stage as a self-contained directory.

**Adding a new stage:**
1. Create `stage2/` (or `stage3/`, etc.) at the same level as `stage1/`.
2. Each stage has its own `config.py`, `scripts/`, `data/`, and `output/`.
3. Shared utility functions (`compute_grid_cell_area`, basin masks, `mann_kendall_pvalue`) should eventually be factored into a common `utils/` package importable by all stages.
4. Any new physical parameterisation should be a pure function with full docstring, units table, and DOI-linked reference.

**Adding a new analysis to Stage 1:**
1. Create `scripts/08_new_analysis.py` following the existing pattern: load from `flux_3d.nc` or `processed_surface.nc`, compute, save to `output/figures/`.
2. Update `README.md` (output figures table) and `CHANGELOG.md`.
3. Optional dependencies should be wrapped in `try/except ImportError`.

**Planned subsequent stages:**
- **Stage 2** — Export flux at 100 m and 500 m (gravitational/biological pump)
- **Stage 3** — Layered total carbon content (0–100 m / 100–500 m / >500 m) + mass-balance consistency check (surface flux integral vs depth layer accumulation)
- **Stage 4** — Robust trend estimation across depth layers
- **Stage 5** — Spatial mapping and hotspot/refugia clustering (k-means on T, pH, Ω)
- **Stage 6** — Thermal vs non-thermal flux decomposition (solubility-driven vs biology-driven pCO₂ changes)
- **Stage 7** — Validation and synthesis across all stages

---

## 10. Full bibliography

**Gas transfer velocity:**
> Wanninkhof, R. (2014). Relationship between wind speed and gas exchange over the ocean revisited. *Limnology and Oceanography: Methods*, 12(6), 351–362. https://doi.org/10.4319/lom.2014.12.351

**CO₂ solubility:**
> Weiss, R. F. (1974). Carbon dioxide in water and seawater: the solubility of a non-ideal gas. *Marine Chemistry*, 2(3), 203–215. https://doi.org/10.1016/0304-4203(74)90015-2

**BGC hindcast model (PISCES-v2):**
> Aumont, O., Ethé, C., Tagliabue, A., Bopp, L., & Gehlen, M. (2015). PISCES-v2: an ocean biogeochemical model for carbon and ecosystem studies. *Geoscientific Model Development*, 8, 2465–2513. https://doi.org/10.5194/gmd-8-2465-2015

**Physical reanalysis (GLORYS12V1):**
> Lellouche, J.-M. et al. (2021). The Copernicus Global 1/12° Oceanic and Sea Ice GLORYS12 Reanalysis. *Frontiers in Earth Science*, 9, 698876. https://doi.org/10.3389/feart.2021.698876

**Validation product — SOCAT:**
> Bakker, D. C. E. et al. (2016). A multi-decade record of high-quality fCO₂ data in version 3 of the Surface Ocean CO₂ Atlas (SOCAT). *Earth System Science Data*, 8, 383–413. https://doi.org/10.5194/essd-8-383-2016

**Key motivating reference — physical injection pump and climate index regression:**
> Bellacicco, M., Marullo, S., Dall'Olmo, G., Iudicone, D., & Buongiorno Nardelli, B. (2025). The oceanic physical injection pump of organic carbon. *Nature Communications*, 16, 7100. https://doi.org/10.1038/s41467-025-62363-z

**Carbon pump conceptual framework:**
> Boyd, P. W., Claustre, H., Levy, M., Siegel, D. A., & Weber, T. (2019). Multi-faceted particle pumps drive carbon sequestration in the ocean. *Nature*, 568, 327–335. https://doi.org/10.1038/s41586-019-1098-2

**ERA5 reanalysis:**
> Hersbach, H. et al. (2020). The ERA5 global reanalysis. *Quarterly Journal of the Royal Meteorological Society*, 146(730), 1999–2049. https://doi.org/10.1002/qj.3803

**Biogeochemical provinces:**
> Fay, A. R. & McKinley, G. A. (2014). Global open-ocean biomes: mean and temporal variability. *Earth System Science Data*, 6, 273–284. https://doi.org/10.5194/essd-6-273-2014

**SAM index:**
> Marshall, G. J. (2003). Trends in the Southern Annular Mode from observations and reanalyses. *Journal of Climate*, 16(24), 4134–4143. https://doi.org/10.1175/1520-0442(2003)016<4134:TITSAM>2.0.CO;2

**ONI index:**
> NOAA CPC (2024). Oceanic Niño Index (ONI). NOAA Physical Sciences Laboratory. https://psl.noaa.gov/data/correlation/oni.data

**Robust trend estimation:**
> Sen, P. K. (1968). Estimates of the regression coefficient based on Kendall's tau. *Journal of the American Statistical Association*, 63(324), 1379–1389. https://doi.org/10.1080/01621459.1968.10480934

**Mann-Kendall significance test:**
> Mann, H. B. (1945). Nonparametric tests against trend. *Econometrica*, 13(3), 245–259. https://doi.org/10.2307/1907187  
> Kendall, M. G. (1975). *Rank Correlation Methods* (4th ed.). Griffin, London.

**LOESS smoothing:**
> Cleveland, W. S. (1979). Robust locally weighted regression and smoothing scatterplots. *Journal of the American Statistical Association*, 74(368), 829–836. https://doi.org/10.1080/01621459.1979.10481038

**PELT breakpoint detection:**
> Killick, R., Fearnhead, P., & Eckley, I. A. (2012). Optimal detection of changepoints with a linear computational cost. *Journal of the American Statistical Association*, 107(500), 1590–1598. https://doi.org/10.1080/01621459.2012.737745

**EOF sampling error:**
> North, G. R., Bell, T. L., Cahalan, R. F., & Moeng, F. J. (1982). Sampling errors in the estimation of empirical orthogonal functions. *Monthly Weather Review*, 110(7), 699–706. https://doi.org/10.1175/1520-0493(1982)110<0699:SEITEO>2.0.CO;2

**EOF methodology:**
> Björnsson, H. & Venegas, S. A. (1997). A manual for EOF and SVD analyses of climate data. McGill University, CCGCR Report No. 97-1.

---

## 11. Data licences and citation

**CMEMS products:** Distributed under the [Copernicus Marine Service licence](https://marine.copernicus.eu/user-corner/service-commitments-and-licence). Free for scientific and commercial use with attribution.

**ERA5:** Subject to the [Copernicus Climate Change Service licence](https://cds.climate.copernicus.eu/api/v2/terms/static/licence-to-use-copernicus-products.pdf).

**NOAA GML atmospheric CO₂:** Public domain.

**ONI and SAM indices:** Public domain (NOAA PSL and NOAA CPC).

**Required acknowledgement for CMEMS data:**

> *"This study has been conducted using E.U. Copernicus Marine Service Information; https://doi.org/10.48670/moi-00019"*

Please also cite the individual data product DOIs and the parameterisation references (Wanninkhof 2014, Weiss 1974) in any publication using this codebase.
