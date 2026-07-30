# Changelog

All notable changes to this project are documented here.  
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).  
Versioning follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [1.4.3] — 2026-07-29

### Changed — `07_eof_analysis.py` — deprecated three EOF analyses

After evaluation, the `monthly_full`, `monthly_anom`, and `annual_anom` runs
did not produce physically relevant results beyond the `annual_full` run:

- **Monthly full field** — EOF1 explains >55% of variance but is dominated
  by the seasonal cycle, not interannual variability. The pattern is not
  scientifically distinct from the time-mean flux map.
- **Monthly anomaly** — removes the seasonal cycle but the residual monthly
  noise dominates over interannual signal given the short (31-year) record.
- **Annual anomaly** — structurally equivalent to the annual full-field run
  for this dataset: the time-mean spatial pattern aligns closely with EOF1
  of the full field, so removing it simply reorders the modes without
  revealing additional information.

`RUNS` reduced to a single entry: `EOFRun("annual_full", ...)`.
Module docstring and `[done]` message updated accordingly.
Output reduced from 8 figures to 2: `fig_eof_annual_full_scree.png` and
`fig_eof_annual_full_modes.png`.

---

## [1.4.2] — 2026-07-29

### Fixed — `07_eof_analysis.py` — two runtime errors

**`ValueError: lower_level and upper_level cannot be NaN` in `contourf`**
- Root cause: `np.nanpercentile(np.abs(eof_map[np.isfinite(eof_map)]), 97)`
  returns NaN when the boolean index produces an empty array (can occur for
  any EOF mode whose map is entirely NaN after the ocean mask is applied,
  e.g. a degenerate or near-zero mode in the full-field run).
- Fix: extract finite values first, check the array is non-empty before
  calling `np.percentile`, and fall back to `vmax = 1.0` if the result is
  still not finite or is zero. `contourf` now always receives valid symmetric
  levels regardless of the EOF content.

**`RuntimeWarning: Mean of empty slice` in `compute_eofs`**
- Root cause: `np.nanmean(data, axis=0)` triggers a warning for grid columns
  that are all-NaN across every time step (e.g. polar rows entirely under
  sea ice or land). This is expected and harmless — `np.isfinite(NaN) = False`
  correctly excludes those pixels from the ocean mask — but produces noisy
  output.
- Fix: wrapped the `np.nanmean` call in `np.errstate(all="ignore")` to
  suppress the warning without masking genuine numerical issues elsewhere.

---

## [1.4.1] — 2026-07-29

### Changed — `07_eof_analysis.py` — full rework

Complete rewrite of the EOF analysis script. The v1.4.0 implementation ran
only one analysis (annual anomaly). v1.4.1 runs four analyses, fixes the
PC scaling convention, and replaces the separate spatial + PC figures with
a single composite per analysis.

#### Four analyses

An `EOFRun` dataclass (`tag`, `label`, `temporal`, `anomaly`) describes each
configuration. `RUNS` is a list of four instances iterated in `main()`:

| Tag | Temporal | Field |
|-----|----------|-------|
| `monthly_full` | Monthly | Full field (no mean removal) |
| `monthly_anom` | Monthly | Anomaly (time-mean removed per pixel) |
| `annual_full`  | Annual mean | Full field |
| `annual_anom`  | Annual mean | Anomaly |

Annual-mean runs suppress the seasonal cycle; full-field runs capture the
dominant climatological spatial structure; anomaly runs isolate interannual
variability. Together the four runs allow comparison of how the leading modes
change depending on temporal averaging and mean-field inclusion.

#### PC scaling convention fixed

v1.4.0 computed PCs by normalising to unit variance and scaling the EOF maps
to match — correct in principle but implemented ambiguously. v1.4.1 makes the
convention explicit and standard:

1. SVD on area-weighted matrix: `X_w = U Σ Vᵀ`
2. Unweight EOFs: `EOF_unit = Vt / w_flat` — unit L2-norm in physical space
3. Project original (unweighted) field: `PC_raw = X @ EOF_unit.T`
4. Normalise: `PC_norm = PC_raw / std(PC_raw)` — unit variance, dimensionless
5. Scale EOF: `EOF_scaled = EOF_unit × std(PC_raw)` — physical units

Result: `EOF_scaled × PC_norm` exactly reconstructs the field. EOF maps show
the typical anomaly pattern [mol C m⁻² yr⁻¹] per one standard deviation of
the PC. This matches the `eofs` package convention and is standard in
physical oceanography.

#### Climate indices split by temporal resolution

`load_climate_indices()` now accepts both `t_monthly` and `t_annual` time
axes and returns a dict with `"monthly"` and `"annual"` sub-dicts, each
reindexed to the appropriate axis. Monthly PC plots get monthly ONI/SAM;
annual PC plots get annual ONI/SAM.

#### Composite figure (3 rows × 2 columns)

Separate `fig_eof_spatial_N.png` and `fig_eof_pc_N.png` files replaced by a
single `fig_eof_{tag}_modes.png` per analysis:
- Rows 1–3: EOF 1, 2, 3 (top to bottom)
- Left column: spatial pattern (Robinson projection map with colourbar)
- Right column: corresponding PC timeseries + ONI/SAM overlay + Spearman ρ

Built with `GridSpec(n_rows=3, ncols=2)`. Map axes created with
`projection=ccrs.Robinson()` when cartopy is available; plain axes otherwise.
`_draw_eof_map()` and `_draw_pc()` are pure helpers that operate on a
pre-created axes object.

#### Scree plot per analysis

`fig_eof_variance.png` replaced by `fig_eof_{tag}_scree.png` (one per run).
Structure unchanged: variance bars + North error + cumulative line + red
outlines on degenerate pairs.

#### Output files: 8 total (4 analyses × 2 figures)

```
fig_eof_monthly_full_scree.png   fig_eof_monthly_full_modes.png
fig_eof_monthly_anom_scree.png   fig_eof_monthly_anom_modes.png
fig_eof_annual_full_scree.png    fig_eof_annual_full_modes.png
fig_eof_annual_anom_scree.png    fig_eof_annual_anom_modes.png
```

---

## [1.4.0] — 2026-07-29

### Added — `07_eof_analysis.py` — EOF decomposition of surface CO₂ flux

New standalone script performing Empirical Orthogonal Function analysis
on the annual-mean surface air-sea CO₂ flux anomaly field (time × lat × lon).

#### Method

1. Annual means of `fgco2` from `flux_3d.nc`, sliced to 1993–2023.
2. Time-mean removed at each pixel → anomaly field.
3. Area-weighting by √cos(lat) applied before SVD so high-latitude pixels
   (which cover less ocean area) do not dominate the variance decomposition.
4. Ocean mask built from NaN pattern of the time-mean; land excluded.
5. SVD on the weighted anomaly matrix (n_years × n_ocean_pixels).
6. EOFs unweighted after decomposition → physical units [mol C m⁻² yr⁻¹].
7. PCs normalised to unit variance (dimensionless); amplitude carried by EOF map.
8. North et al. (1982) sampling error computed for each eigenvalue:
   δλ = λ · √(2/N), used to identify degenerate modes on the scree plot.

#### Figures produced

**`fig_eof_variance.png`** — scree plot.
- Bars: variance explained [%] per EOF with North et al. error bars.
- Line: cumulative variance on right axis; 90% reference line.
- Red outlines on bars where North criterion indicates degeneracy
  (neighbouring eigenvalue error bars overlap — modes not individually
  interpretable as distinct physical patterns).

**`fig_eof_spatial_N.png`** (N = 1, 2, 3) — spatial EOF patterns.
- White-centred diverging colourmap (`cmocean.balance` / `RdBu_r`).
- Colourbar scale set to 97th percentile of |EOF| for each mode.
- Robinson projection; land on top of data (zorder=3/4).
- Variance explained annotated in title.

**`fig_eof_pc_N.png`** (N = 1, 2, 3) — PC timeseries.
- Normalised PC (unit variance) as filled line.
- If `data/oni_monthly.txt` and `data/sam_monthly.txt` are present:
  z-score ONI and SAM overlaid; Spearman ρ and p-value annotated
  per index in the legend. Allows direct attribution of flux modes
  to known climate variability (ENSO, SAM).
- If index files absent: PC plotted alone, no crash.

#### Design decisions

- Annual means chosen over monthly to avoid EOF patterns being dominated
  by the seasonal cycle rather than interannual variability.
- `scipy.linalg.svd` (via `numpy.linalg.svd`) used directly rather than
  a wrapper (e.g. `eofs` package) for full transparency and portability.
- `N_EOFS = 3` and `N_SCREE = 10` are module-level constants, easily
  adjusted at the top of the script.
- Script is self-contained: loads `flux_3d.nc` directly, no dependency
  on `05_plot_results.py` or `06_fay_analysis.py`.

---

## [1.3.2] — 2026-07-29

### Fixed — HPC compatibility: graceful degradation when compute nodes lack internet

Leonardo (and similar HPC systems) run job steps on compute nodes with no
outbound internet access. v1.3.1 called `requests.get()` unconditionally
inside `load_climate_indices()`, which crashed the entire pipeline if the
ONI or SAM file was not already cached — even when the run had nothing to
do with climate indices.

#### `05_plot_results.py` and `06_fay_analysis.py`

**`load_climate_indices()` return type changed `dict` → `dict | None`.**
- Every network call is now wrapped in a `try/except Exception` block.
- On failure (connection error, timeout, HTTP error), the function prints a
  clearly formatted `[climate]` warning block with exact `wget` commands for
  manual download from a login node, then returns `None` instead of raising.
- An empty-parse guard added: if either file exists but parses to zero rows
  (e.g. corrupt or wrong-format file), the function also returns `None` with
  a diagnostic message rather than crashing downstream.

**`main()` updated in both scripts.**
- `idx = load_climate_indices()` result is checked with `if idx is not None`
  before calling `fig07_climate_index_regression()` (in 05) or
  `fig_fay_climate_regression()` (in 06).
- When `idx is None`, a `[skip]` message is printed and the run continues
  normally, producing all other figures without interruption.
- Final `[done]` message reports the correct figure count (6 or 7 depending
  on index availability).

**Manual download instructions** (run once from a login node with internet):
```bash
wget -O <DATA_DIR>/oni_monthly.txt \
    https://psl.noaa.gov/data/correlation/oni.data

wget -O <DATA_DIR>/sam_monthly.txt \
    https://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/aao/monthly.aao.index.b79.current.ascii.table
```
Once both files are present in `data/`, resubmit the job — the scripts detect
the files, skip the download block, and produce fig07 / fig_fay_climate_regression
without any further intervention.

---

## [1.3.1] — 2026-07-28

### Added — ENSO & SAM climate index regression (Bellacicco et al. 2025)

Implements the climate-mode correlation analysis from Bellacicco et al. (2025,
*Nature Communications* 16, 7100), applied here to the surface air-sea CO₂
flux J_net rather than the Physical Injection Pump export studied in that paper.
The core question: does the same ENSO/SAM teleconnection that modulates
subsurface carbon export also manifest at the air-sea interface?

#### `05_plot_results.py`

**`load_climate_indices()`** — new function (note: HPC robustness added in v1.3.2).
- Downloads ONI (Oceanic Niño Index) from NOAA PSL and SAM (Southern Annular
  Mode) from NOAA CPC on first run; caches both in `data/` as
  `oni_monthly.txt` and `sam_monthly.txt`.
- Parses space-delimited fixed-width formats. Missing ONI values (−99.9)
  dropped. Both series aligned to the 1993–2023 project window and returned
  as monthly and annual-mean `pandas.Series`.
- Scientific sources: ONI = ERSSTv5 Niño-3.4 3-month running mean (NOAA PSL);
  SAM = Marshall (2003) station-based index, NOAA CPC update.

**`fig07_climate_regression_global.png`** — new 3-panel figure.
- *Panel a*: bar chart of annual J_net (left axis) overlaid with z-score
  normalised ONI (red) and SAM (green) on the right axis.
- *Panel b*: scatter J_net vs ONI (circles) and vs SAM (triangles), separate
  OLS trend lines, Spearman ρ and p-value annotated per index.
- *Panel c*: scatter J_net vs OLS multilinear combination α·ONI + β·SAM,
  points coloured by year (`plasma`). Spearman ρ and combined-predictor
  coefficients annotated. Follows Bellacicco et al. (2025) Fig. 3 structure.

**`main()`** updated to call `load_climate_indices()` and pass result to
`fig07_climate_index_regression()`. (Skip logic added in v1.3.2.)

#### `06_fay_analysis.py`

**`load_climate_indices()`** — self-contained copy so script runs independently;
reuses the same cached files in `data/`.

**`fig_fay_climate_regression.png`** — new 17-panel figure.
- Each panel: three bars — ρ(J_net, ONI) red, ρ(J_net, SAM) green,
  ρ(J_net, α·ONI+β·SAM) blue — per Fay (2014) domain.
- Filled bars = significant at α = 0.05; hatched (///) = not significant.
  ✱ annotation on significant bars. Reference lines at ±0.4.
- Allows direct per-province comparison between surface flux ENSO/SAM
  response and the subsurface PIP signal reported by Bellacicco et al.

---

## [1.3.0] — 2026-07-28

### Added — `05_plot_results.py` major expansion + `06_fay_analysis.py` new script

This release adds five new figures to `05_plot_results.py`, significantly
upgrades two existing ones, and splits all Fay domain analysis into a
dedicated `06_fay_analysis.py` for performance reasons.

#### `05_plot_results.py`

**Global time cutoff applied in `load_all()`.**
All `xr.DataArray` and `pandas.Series` objects are sliced to
`1993-01-01 → 2023-12-31` at load time. Every downstream figure uses
pre-cropped data automatically — no per-figure masking.

**`fig01_flux_vs_co2.png`** — merged and upgraded (was simple twin-axis).
- Raw monthly J_net (faint fill + line) + LOESS smooth (bold, frac=0.15).
- PELT breakpoint detection; only B1 and B4 kept, defining 3 segments.
- Per-segment OLS trend annotated **above** the trend line [Pg C yr⁻²].
- NOAA CO₂ annual mean on right axis; CO₂ slope annotated **below** its line.
- B1/B4 labels on top twin x-axis. Fallback `_fig01_simple()` if
  `statsmodels` not installed.
- y-axis fixed at −0.2 → 3.0 Pg C yr⁻¹.

**`fig03_trend_significance_map.png`** — upgraded with significance overlay.
- Per-pixel Sen's slope (annual means) with Mann-Kendall significance at
  α = 0.01 (tightened from 0.05).
- Hatching (///) applied to ocean-only non-significant pixels
  (`insig_mask &= ocean_2d` — land never hatched).
- Gray `contour` outline at the significant/non-significant boundary.
- Colorbar fixed at ±0.5 mol C m⁻² yr⁻¹ per decade.
- `cartopy` land feature drawn last (on top of hatching).

**`fig02_mean_flux_map.png`** — range fixed at ±2.5 mol C m⁻² yr⁻¹; white-
centred diverging colourmap (`cmocean.balance` / `RdBu_r` fallback); title
simplified.

**`fig04_delta_pco2_map.png`** — colorbar range fixed at ±100 µatm.

**`fig06_sink_saturation.png`** (renumbered from fig07).
- Annual J_net vs annual atmospheric CO₂ scatter, coloured by year.
- OLS regression with β₁ [Pg C yr⁻¹ ppm⁻¹], Pearson r, p-value annotated.
- Data pre-cropped by `load_all()`.

**`mann_kendall_pvalue()`** helper: two-sided MK test via
`scipy.stats.kendalltau` on a monotone integer index.

**Removed from 05:** `fig08–fig10` (Fay figures) → moved to `06_fay_analysis.py`.
Old standalone `fig06_loess_trend()` merged into `fig01`.

#### `06_fay_analysis.py` — new script

All three Fay domain figures separated from `05_plot_results.py` for
computational performance (17-domain area-integration loop is heavy).

- **`fig_fay_ts.png`**: 17-panel monthly timeseries [Pg C yr⁻¹] per Fay domain,
  coloured per domain (jet colourmap), biome map in 18th panel.
- **`fig_fay_clim.png`**: 17-panel climatological seasonal cycle ±1σ per domain.
- **`fig_fay_trends.png`**: 17-panel Sen's slope bars with Mann-Kendall
  significance at α = 0.01; filled = significant, hatched = not.
- `load_fay_biomes()`, `domain_flux_timeseries()`, `_fay_global_biome_map()`,
  `compute_grid_cell_area()` helpers — self-contained.
- Applies same `T_START / T_END` cutoff as 05.

---

## [1.2.3] — 2026-07-28

### Changed — `04_validate.py`

**Global colour updates:**
- `C_REC` changed from `#0072B2` to `#0558AD` (deeper blue).
- `C_OBS` changed from `#CC0000` to `#CD2859` (deep rose/red).
- Segment colours for LOESS breakpoints replaced with 8 visually distinct
  hardcoded values — warm amber (`#E69F00`), sea green (`#2E8B57`), vivid
  violet (`#9400D3`), burnt orange (`#FF6F20`), cerulean (`#1ABDE8`), dark
  red (`#8B0000`), charcoal (`#3D3D3D`), dark goldenrod (`#B8860B`) — none
  overlapping with the two main curve colours.

**`plot_timeseries`:**
- Rolling Pearson r line changed from dash-dot to solid.
- First subplot legend `bbox_to_anchor` changed from `(0.5, -0.12)` to
  `(0.5, -0.07)`.
- `hspace` reduced from 0.45 to 0.30.

**`plot_validation_map`:**
- `hspace` reduced from 0.40 to 0.20.
- Colorbar left edge moved to tighter position.

**`plot_loess`:**
- Annotation `xy` anchored to trend line midpoint rather than LOESS curve.
- Arrow linewidth increased from 0.7 to 1.8. Label `fontsize` to 10,
  `fontweight="bold"` added. Trend line `lw` increased to 2.2.

**`plot_fay_ts` / `plot_fay_clim`:**
- `hspace=0.32, wspace=0.25`; `figsize=(20, 24)`.
- Subplot titles `fontsize=12, fontweight="bold"`. Axis labels `fontsize=11`.
- Domain inset maps removed.
- Reconstruction curve coloured per domain (jet colourmap, 1–17).
- MULTIOBS curve fixed at `C_OBS` across all domains.

---

## [1.2.2] — 2026-07-27

### Changed — `04_validate.py`

- `C_OBS` changed from `#D55E00` (orange on screen) to `#CC0000` (pure red).
- Panel labels increased to `fontsize=18, fontweight="bold"`.
- Both subplot legends moved below axes; `fig.subplots_adjust(bottom=0.15)`.
- Colorbar positioning tightened in `plot_validation_map`.
- Breakpoint vertical lines changed to `ls="-."`, `color="#444444"`, `lw=1.4`.
- B1–Bn labels added on top x-axis via `twiny()`.
- Breakpoint dates printed to stdout.
- Same breakpoints applied to both REC and OBS LOESS curves; same segment
  colour index per segment; REC labels below, OBS labels above trend lines.
- `plot_spectra`: Nyquist rightmost tick; `ax.set_xlim(right=nyq_yr * 0.85)`.
- `plot_cps`: reference lines `color="red", lw=2.0`; `sharex=False`.

---

## [1.2.1] — 2026-07-27

### Changed — `04_validate.py` — full styling and analysis pass

- `plot_style.py` integrated via `sns.set_theme` (paper context,
  `font_scale=1.2`). Full rcParams: 4-sided spines, bidirectional inout ticks,
  `axes.linewidth=1.5`, `#CCCCCC` grid, bordered legends, 300 dpi.
- Data hard cut at `CUTOFF = np.datetime64("2024-01-01")` in `load_data()`.
- CVD-safe colours: `C_REC=#0072B2`, `C_OBS=#D55E00`, `C_RMSD=#CC79A7`,
  `C_CORR=#009E73`.
- `_clean()`, `_label_panel()`, `_add_colorbar()` helpers added.
- `plot_timeseries` → 2-row: monthly timeseries + 12-month rolling RMSD/r.
- `plot_loess` → PELT breakpoint detection with per-segment OLS trends.
- `plot_validation_map` → `contourf`, Robinson projection, fixed colour limits.
- New: `plot_fay_ts()`, `plot_fay_clim()`, `plot_spectra()`, `plot_cps()`.

---

## [1.2.0] — 2026-07-25

### Summary
Dropped MULTIOBS reconstruction. `03_compute_flux.py` produces single
reconstruction only (`GLOBAL_MULTIYEAR_BGC_001_029` pCO₂ + Wanninkhof k
with optional wind variance correction).

### Changed — `03_compute_flux.py`
- `compute_flux()` label parameter removed.
- `main()` produces one flux (`fgco2`). Removed `fgco2_improved` and
  `J_net_improved_PgC` from output files.

### Changed — `04_validate.py`
- Full rewrite: reconstruction vs MULTIOBS only. Two-curve timeseries.
  Validation map reverts to 1×2 layout. Taylor diagram: single model point.

---

## [1.1.0] — 2026-07-24

### Added — wind variance correction (σ²_u)

### Changed — `01_download_data.py`
- `download_era5_wind_daily()` added: ERA5 daily 10 m wind in yearly chunks
  via CDS API, merged into `era5_wind10m_daily.nc`.

### Changed — `02_preprocess.py`
- `load_spco2_ocean_obs()` added.
- `load_wind_variance()` added: σ²_u from ERA5 daily via `resample("1ME").var()`.

### Changed — `03_compute_flux.py`
- `gas_transfer_velocity()` extended with optional `wind_variance` argument.
  When provided: k ∝ (⟨u⟩² + σ²_u).

### Changed — `04_validate.py`
- Extended to 3 curves (original, improved, MULTIOBS).

---

## [1.0.1] — 2026-07-23

### Fixed — `01_download_data.py`
- `download_era5_wind()` body uncommented and made callable.

### Fixed — `02_preprocess.py`
- `valid_time` → `time` rename for ERA5 files.
- MULTIOBS flux regridded to BGC 0.25° grid before time alignment.

### Fixed — `03_compute_flux.py`
- Sign convention corrected: `delta_pco2 = spco2_atm - spco2_ocean`.
- `global_integral()` return type annotation fixed.

### Fixed — `04_validate.py`
- `fgco2_obs` sign flip added.
- `compute_grid_cell_area()` helper added (was missing; caused NameError).
- `global_integral_obs()` sign corrected.

---

## [1.0.0] — 2026-07-23

### Summary
First complete implementation of Stage 1: surface net air-sea CO₂ flux
reconstruction, global integration, cross-validation, and figure suite.

### Added
- `config.py` — central configuration module.
- `01_download_data.py` — data acquisition from CMEMS, NOAA, ERA5.
- `02_preprocess.py` — regrid, unit convert, harmonize, merge datasets.
- `03_compute_flux.py` — Schmidt number, K0, k, flux, global integral.
- `04_validate.py` — reconstruction vs MULTIOBS skill metrics and figures.
- `05_plot_results.py` — 5 primary result figures (global only).
- `README.md`, `CHANGELOG.md`, `requirements.txt`.
