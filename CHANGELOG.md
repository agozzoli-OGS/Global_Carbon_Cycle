# Changelog

All notable changes to this project are documented here.  
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).  
Versioning follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [1.2.3] — 2026-07-28

### Changed — `04_validate.py`

**Global colour updates:**
- `C_REC` changed from `#0072B2` to `#0558AD` (deeper blue).
- `C_OBS` changed from `#CC0000` to `#CD2859` (deep rose/red).
- Segment colours for LOESS breakpoints replaced with 8 visually distinct hardcoded values — warm amber (`#E69F00`), sea green (`#2E8B57`), vivid violet (`#9400D3`), burnt orange (`#FF6F20`), cerulean (`#1ABDE8`), dark red (`#8B0000`), charcoal (`#3D3D3D`), dark goldenrod (`#B8860B`) — none overlapping with the two main curve colours.

**`plot_timeseries`:**
- Rolling Pearson r line changed from dash-dot (`ls="-."`) to solid (`ls="-"`).
- First subplot legend `bbox_to_anchor` changed from `(0.5, -0.12)` to `(0.5, -0.07)` — less padding between plot and legend.
- `hspace` in `GridSpec` reduced from 0.45 to 0.30 — less vertical gap between the two rows.

**`plot_validation_map`:**
- `hspace` reduced from 0.40 to 0.20 — less gap between RMSD and bias panels.
- `right` in `GridSpec` increased from 0.91 to 0.93; colorbar left edge moved from 0.923 to 0.945 — tighter spacing between map panels and colorbars.

**`plot_loess`:**
- Annotation `xy` now anchored to the midpoint of the **trend line** (`y_fit[mid_idx]`) rather than the LOESS smooth curve (`y_seg[mid_idx]`). Arrow thus originates from the fitted line, not the data.
- Arrow linewidth increased from 0.7 to 1.8.
- Label `fontsize` increased from 8 to 10; `fontweight="bold"` added.
- Trend line `lw` increased from 2.0 to 2.2.

**`plot_fay_ts` / `plot_fay_clim`:**
- `hspace=0.32, wspace=0.25` (was 0.42/0.30) — tighter subplot spacing.
- `figsize=(20, 24)` (was `(20, 26)`) — more compact overall figure.
- Subplot titles `fontsize=12, fontweight="bold"` (was `fontsize=10`).
- Axis labels `fontsize=11` (was `fontsize=9`). Tick labels `fontsize=10` (was `fontsize=8/9`).
- Suptitle `fontsize=14` (was 12). Figure legend `fontsize=11` (was 10).
- Domain inset maps removed entirely from both figures.
- Reconstruction curve now coloured **per domain** using `plt.cm.jet` scaled over domains 1–17 — matching the global biome map colormap so each domain's curve colour corresponds to its region on the map.
- MULTIOBS curve fixed at `C_OBS` (`#CD2859`) across all domains — consistent identification regardless of domain.
- Figure legend updated: reconstruction entry labelled "domain colour" to clarify the per-domain colouring scheme.

---

## [1.2.2] — 2026-07-27

### Changed — `04_validate.py`

**Global colour fix:**
- `C_OBS` changed from `#D55E00` (Wong vermilion, displays orange on screen) to `#CC0000` (pure red). Applied globally to all figures.

**`plot_timeseries`:**
- Panel labels `a` and `b` increased to `fontsize=18, fontweight="bold"`.
- Both subplot legends moved below their respective axes using `bbox_to_anchor=(0.5, -0.07)` (row 1) and `(0.5, -0.22)` (row 2), `ncol=2`, `frameon=True`, `edgecolor="#555555"`.
- `fig.subplots_adjust(bottom=0.15)` added to give legends space below the figure.

**`plot_validation_map`:**
- Panel titles reformatted as multiline strings: product name on its own line, no units in the title (units only in the colorbar label).
- Colorbar positioning tightened: `right=0.93` (was 0.91), colorbar `left=0.945` (was 0.923) — gap between map edge and colorbar reduced from ~0.033 to ~0.015 in figure-fraction units.

**`plot_loess`:**
- Breakpoint vertical lines changed from `ls="--"` to `ls="-."`, `color="#444444"`, `lw=1.4` — darker and dash-dot style.
- B1–Bn labels added on the top x-axis via a `twiny()` axis with ticks placed at breakpoint dates; spine hidden, only tick marks and labels shown.
- Breakpoint dates printed to stdout (year-month format) for reference.
- Same set of breakpoints (detected on the reconstruction LOESS smooth) applied to both the reconstruction and MULTIOBS trend lines — same integer indices used for both curves.
- Same segment colour index used for both curves within each segment — visual correspondence between reconstruction and MULTIOBS trend lines per segment.
- Reconstruction trend labels placed **below** the trend line (`xytext=(0, -16)`); MULTIOBS trend labels placed **above** (`xytext=(0, 16)`).
- Trend line `lw` increased to 2.0. Arrow `lw=0.7` (was 0.7 — unchanged here, improved further in v1.2.3).
- Segment colour palette updated to skip `C_REC`/`C_OBS`: Wong orange, sky blue, yellow, black, mid-gray, saddlebrown, purple, dark turquoise.

**`plot_spectra`:**
- Rightmost x-axis tick set to actual Nyquist period: $2/f_s = 2/12$ yr = 2 months. `_set_period_ticks()` updated to accept `min_period_yr` argument and prepend the Nyquist tick when it is shorter than the standard 30 d tick.
- `ax.set_xlim(right=nyq_yr * 0.85)` added to prevent the x-axis extending beyond the resolvable period.

**`plot_cps`:**
- All threshold/reference lines changed to `color="red", lw=2.0` (gain threshold at 1, phase reference at 0, coherence threshold at 0.5) — more visible.
- Same Nyquist rightmost tick fix as `plot_spectra`.
- `sharex=False` between panels — each panel (gain, phase, coherence) has its own independent x-limit control.

---

## [1.2.1] — 2026-07-27

### Changed — `04_validate.py` — full styling and analysis pass

**Global:**
- `plot_style.py` integrated throughout via `sns.set_theme` (paper context, `font_scale=1.2`). Full rcParams: 4-sided spines, bidirectional inout ticks, `axes.linewidth=1.5`, `#CCCCCC` dashed grid behind data, bordered legends with `#AAAAAA` edge, `savefig.dpi=300`.
- Data hard cut at `CUTOFF = np.datetime64("2024-01-01")` applied in `load_data()` — all figures automatically use the trimmed record.
- Hardcoded CVD-safe colours introduced: `C_REC=#0072B2` (Wong blue), `C_OBS=#D55E00` (Wong vermilion), `C_RMSD=#CC79A7` (reddish-purple), `C_CORR=#009E73` (bluish-green).
- Helper functions `_clean()`, `_label_panel()`, `_add_colorbar()` applied consistently across all figure functions.

**`plot_timeseries` → converted to 2-row figure:**
- Row 1: monthly global time series (unchanged from v1.2.0).
- Row 2 (new): `rolling_rmsd_corr()` computes 12-month rolling RMSD and Pearson r on raw monthly integrals. RMSD plotted on left y-axis (`C_RMSD`), Pearson r on right y-axis (`C_CORR`). Panel labels `a`/`b` added.

**`plot_loess` → added PELT breakpoint detection:**
- LOESS smooth unchanged (frac=0.15, `statsmodels.lowess`).
- If `ruptures` is installed: PELT algorithm with RBF cost (`pen=3`) detects breakpoints in the reconstruction LOESS smooth. For each segment: piecewise linear trend (dotted, colour-coded per segment from `plt.cm.tab10`) plotted; slope in Pg C yr⁻² annotated with an arrowhead. Breakpoint vertical dashed lines added. Graceful fallback (no annotation) if `ruptures` not installed.

**`plot_fay_domains` → split into `plot_fay_ts` and `plot_fay_clim`:**
- `fig_fay_ts.png`: 17 panels of monthly time series only.
- `fig_fay_clim.png`: 17 panels of climatological seasonal cycle with ±1σ inter-annual shading. Month labels rotated 45°.
- Both: domain inset maps added (cartopy PlateCarree, `set_global()`), global biome map in panel 18 via `gridspec.GridSpec[5,2]`, figure legend with correct `bbox_to_anchor`, `gridspec` `hspace/wspace` explicitly set.

**`plot_validation_map`:**
- 2×2 layout (original vs improved) replaced with 1×2 (RMSD top, bias bottom) — improved reconstruction dropped.
- `cmocean.cm.amp` for RMSD, `cmocean.cm.balance` for bias (with `matplotlib` fallbacks).
- `contourf` with 21 levels, fixed colourbar extents (`VAL_RMSD_MAX=5`, `VAL_BIAS_MIN=-5`, `VAL_BIAS_MAX=5` from `config.py`). Robinson projection.

**`plot_spectra` (new):**
- Welch PSD (`nperseg = max(n//3, 24)`). Log-log. Period in years, x-axis inverted (long periods left). Custom period ticks: 30d, 90d, 180d, 1yr, 5yr, 10yr. Major grid `#999999/lw=0.8`, minor grid `#CCCCCC/lw=0.4`. No −5/3 slope reference.

**`plot_cps` (new):**
- **Coherence bug fixed:** replaced manual `|Pxy|²/(Pxx·Pyy)` (gives coherence ≡ 1 without segment averaging) with `scipy.signal.coherence()`. `nperseg=n//4` ensures multiple averaging segments. `nfft=2^(ceil(log2(n))+3)` for smooth output.
- Three panels: gain (blue), phase in days (orange), magnitude-squared coherence (green).
- Annual/semi-annual marker lines removed. Title split into two lines with `—` divider; "(reference: ...)" parenthetical removed.
- Custom period ticks applied to all three panels.

**`validation_metrics.csv`:**
- Now includes per-domain statistics for all 17 Fay biomes (`domain_1` through `domain_17`) in addition to global scalar metrics, written as a wide-format CSV with `pd.DataFrame(all_metrics).T`.

---

## [1.2.0] — 2026-07-27

### Summary
Project settles on a single, clean reconstruction using exclusively `GLOBAL_MULTIYEAR_BGC_001_029` as the pCO₂ source. The MULTIOBS pCO₂ reconstruction introduced in v1.1.0 is dropped. The MULTIOBS product is retained as a validation reference only. Validation pipeline substantially expanded with Fay domain analysis, LOESS variability, and spectral diagnostics.

### Changed — `03_compute_flux.py`
- All code related to the MULTIOBS pCO₂ improved reconstruction removed. Single reconstruction only.
- `compute_flux()` `label` parameter removed — output always named `fgco2`.
- `gas_transfer_velocity()` retains the `wind_variance` optional argument (applied equally to the single reconstruction when ERA5 daily winds are available).
- Output files simplified: `flux_3d.nc` contains `fgco2`, `k`, `K0`, `Sc` only. `global_flux.nc` contains `J_net_PgC` and `J_net_PgC_annual` only.
- Product label in all NetCDF attributes updated to `GLOBAL_MULTIYEAR_BGC_001_029`.

### Changed — `04_validate.py` — full rewrite

**Sign fix:** `fgco2_obs = ds_surf["fgco2_obs"]` — no negation. Previous versions applied a sign flip (`-ds_surf["fgco2_obs"]`) that was product-version-dependent and had been causing inconsistency.

**Dropped:**
- Taylor diagram (`fig_taylor.png`).
- MULTIOBS pCO₂ reconstruction (green) curve — two curves only (reconstruction vs MULTIOBS flux).
- 2×2 map layout from v1.1.x.

**Added:**
- `fig_validation_loess.png` — LOESS multi-year variability (frac=0.15, `statsmodels`).
- `fig_fay_domains.png` — 17-panel Fay (2014) biome domain validation. Loads `data/Time_Varying_Biomes.cmems.nc`. Area-weighted domain mean time series computed via `domain_mean_timeseries()`. Per-panel: monthly time series + climatological seasonal cycle on twin y-axis + inset domain map.
- `fig_spectra.png` — Welch PSD, log-log, period in years.
- `fig_cps.png` — cross-power spectrum: gain, phase [days], magnitude-squared coherence. Uses `scipy.signal.coherence()`. MULTIOBS as reference, reconstruction as target.
- `validation_metrics.csv` extended with per-domain Fay biome statistics.
- All product labels in figure text updated to `GLOBAL_MULTIYEAR_BGC_001_029`.

---

## [1.1.2] — 2026-07-27

### Changed — `04_validate.py`
- All figure labels updated from "Original"/"Improved" nomenclature to explicit pCO₂ source names.
- Time series legend: `"Original reconstruction (PISCES pCO₂)"` → `"pCO₂: GLOBAL_MULTIYEAR_BGC_001_029"`. `"Improved reconstruction (MULTIOBS pCO₂)"` → `"pCO₂: MULTIOBS (SOCAT-NN)"`.
- 2×2 map panel titles: `"RMSD — Original"` / `"RMSD — Improved"` → `"RMSD — pCO₂: GLOBAL_MULTIYEAR_BGC_001_029"` / `"RMSD — pCO₂: MULTIOBS"`.
- Taylor diagram markers: `"Original (r = ...)"` / `"Improved (r = ...)"` → `"pCO₂: GLOBAL_MULTIYEAR_BGC_001_029 (r = ...)"` / `"pCO₂: MULTIOBS (r = ...)"`.
- All docstrings updated to match.

---

## [1.1.1] — 2026-07-27

### Added — `01_download_data.py`
- `download_era5_wind_daily()` function: downloads ERA5 daily 10 m wind components (`u10`, `v10`) from the CDS API year-by-year (1993–2025) to avoid CDS timeout and size limits. Each year saved as `era5_wind10m_daily_YYYY.nc`, then merged into `era5_wind10m_daily.nc` via `xr.open_mfdataset`. The function lists yearly chunk files after merging so they can be manually deleted to recover disk space (~40–50 GB total).
- `download_era5_wind_daily()` called from `main()` after the existing `download_era5_wind()` call.

> **Note:** The daily download was absent in v1.1.0 despite `ERA5_WIND_DAILY_FILE` being defined in `config.py` and `load_wind_variance()` being implemented in `02_preprocess.py`. This omission meant the wind variance correction never activated in v1.1.0 because the source file was never generated.

---

## [1.1.0] — 2026-07-23

### Summary
Introduced a second improved reconstruction alongside the original, addressing the two dominant error sources identified during validation (r = 0.445, normalised σ = 0.22 in v1.0.0): (1) the free-running PISCES pCO₂ field lacks regional variance; (2) monthly-mean winds underestimate the gas transfer velocity. Both reconstructions preserved in all outputs and figures. This approach was subsequently reconsidered in v1.2.0 which dropped the MULTIOBS reconstruction.

### Changed — `config.py`
- `ERA5_WIND_DAILY_FILE` path added for the ERA5 daily wind download.
- `VAL_RMSD_MAX = 5.0`, `VAL_BIAS_MIN = -5.0`, `VAL_BIAS_MAX = 5.0` added — fixed colourbar extents for validation comparison maps, ensuring both reconstruction maps share the same scale.
- `FIGURE_DPI` raised from 150 to 300.

### Changed — `02_preprocess.py`
- `load_spco2_ocean_obs()` added: loads MULTIOBS `spco2` variable, auto-detects units (Pa or µatm), converts to atm, regrids to the BGC 0.25° grid via bilinear interpolation, saves as `spco2_ocean_obs` in `processed_surface.nc`.
- `load_wind_variance()` added: opens ERA5 daily wind file, computes daily scalar speed $\sqrt{u_{10}^2 + v_{10}^2}$, resamples to monthly variance via `resample("1ME").var()`, regrids to 0.25°, saves as `wind_variance` in `processed_surface.nc`. Returns `None` gracefully if `ERA5_WIND_DAILY_FILE` does not exist — downstream scripts handle `None` without error.
- Both new loaders are optional — if source files are absent, `processed_surface.nc` is produced with the v1.0.x variable set and all downstream scripts fall back to the original reconstruction only.

### Changed — `03_compute_flux.py`
- `gas_transfer_velocity()` extended with optional `wind_variance: xr.DataArray | None` argument. When provided: $k \propto (\langle u \rangle^2 + \sigma_u^2)$ instead of $\langle u \rangle^2$ alone. The correction accounts for the calibration of the Wanninkhof (2014) coefficient against the full wind speed distribution rather than monthly means.
- `compute_flux()` generalised with a `label: str` parameter so the same function produces both `fgco2` and `fgco2_improved` without code duplication.
- `main()` computes and saves both reconstructions. `flux_3d.nc` gains `fgco2_improved` and `k_improved`. `global_flux.nc` gains `J_net_improved_PgC` and its annual resample. Gracefully skips the improved reconstruction if `spco2_ocean_obs` is absent.

### Changed — `04_validate.py`
- Time series figure: extended from 2 to 3 curves — original (blue), improved (green), MULTIOBS reference (red dashed).
- Validation maps: 1×2 layout (RMSD, bias) replaced by 2×2 — left column: original reconstruction; right column: improved reconstruction. Both rows use the same fixed colourbar extents from `config.py`. `contourf` with 21 levels, Robinson projection.
- Taylor diagram: extended from 1 to 2 model markers — ● original (blue circle), ▲ improved (green triangle). Angle = arccos(r), radius = σ_rec/σ_obs.
- `validation_metrics.csv` now contains metrics for both reconstructions with `orig_*` and `imp_*` prefixes.

---

## [1.0.1] — 2026-07-23

### Fixed — `01_download_data.py`
- `download_era5_wind()` body uncommented and made callable. In v1.0.0 the CDS API call was left as a commented template, preventing ERA5 wind download.

### Fixed — `02_preprocess.py`
- Added `valid_time` → `time` dimension rename for ERA5 files: CDS sometimes returns `valid_time` instead of `time` as the temporal coordinate, causing alignment failures downstream.
- MULTIOBS flux regridded to the BGC 0.25° grid via `interp(latitude=..., longitude=..., method="linear")` before time alignment. In v1.0.0, time alignment was attempted on the original MULTIOBS grid (potentially different resolution), which silently produced misaligned arrays.

### Fixed — `03_compute_flux.py`
- **Sign convention corrected throughout.** `delta_pco2 = spco2_atm - spco2_ocean` (positive = uptake). v1.0.0 had `delta_pco2 = spco2_ocean - spco2_atm` (positive = outgassing), inconsistent with the stated convention and with MULTIOBS.
- `global_integral()` return type annotation fixed to `xr.Dataset` (was incorrectly typed as `xr.DataArray`).
- `J_net_PgC` sign: with the corrected flux convention, `J_net_PgC = F_global_PgC` (no additional sign flip needed). v1.0.0 applied an extra `-` sign.

### Fixed — `04_validate.py`
- `fgco2_obs` sign flip added (`fgco2_obs = -ds_surf["fgco2_obs"]`) to align the MULTIOBS sign convention with the (then-current) reconstruction sign convention.
- `compute_grid_cell_area()` helper duplicated into `04_validate.py` — it was missing, causing a `NameError` when `global_integral_obs()` was called.
- `global_integral_obs()` sign corrected to produce positive = uptake, consistent with `global_integral()` in `03_compute_flux.py`.

---

## [1.0.0] — 2026-07-23

### Summary
First complete implementation of Stage 1: surface net air-sea CO₂ flux reconstruction, global integration, cross-validation, and figure suite.

### Added — `config.py`
- Central configuration module: all project-wide constants, CMEMS product/dataset/variable identifiers, Wanninkhof (2014) coefficient `a=0.251` and Schmidt number polynomial coefficients, Weiss (1974) solubility coefficients A1–A3/B1–B3, unit conversion factors (`PA_TO_ATM`, `UATM_TO_ATM`, `CMHR_TO_MS`, `S_TO_YR`, `MOL_C_TO_G`, `MOL_C_TO_PG`), `EARTH_RADIUS_M`, `TOTAL_OCEAN_AREA_M2`, plotting defaults (`FIGURE_DPI=150`, colormaps), `VAL_RMSD_MAX/BIAS_MIN/BIAS_MAX` for validation maps.

### Added — `01_download_data.py`
- `download_bgc_hindcast()`: CMEMS `copernicusmarine.subset` for `spco2` from `GLOBAL_MULTIYEAR_BGC_001_029`, surface level only, full record 1993–2026, monthly.
- `download_physical_reanalysis()`: `thetao` and `so` from `GLOBAL_MULTIYEAR_PHY_001_030` (GLORYS12V1), surface level, full record.
- `download_multiobs_surface()`: `fgco2` and `spco2` from `MULTIOBS_GLO_BGC_CARBON_SURFACE_MYNRT_015_008`, for validation.
- `download_noaa_co2()`: NOAA GML `co2_mm_gl.csv` via HTTP `requests.get`.
- `download_era5_wind()`: ERA5 monthly wind `u10`, `v10` template (CDS API call — left commented in v1.0.0, fixed in v1.0.1).
- `--test` CLI flag: restricts download to 2010 only for credential verification.
- All download functions are idempotent: skip if output file already exists.

### Added — `02_preprocess.py`
- `load_spco2_ocean()`: opens `bgc_hindcast_spco2.nc`, squeezes depth dimension, converts Pa → atm (`×1/101325`).
- `load_sst_sss()`: opens `phy_reanalysis_sst_sss.nc`, squeezes to surface, bilinearly regrids 1/12° GLORYS12 → 0.25° BGC grid via `xr.DataArray.interp`.
- `load_wind_speed()`: opens ERA5 monthly file, computes $\sqrt{u_{10}^2 + v_{10}^2}$, regrids to 0.25°. Returns `None` gracefully if file absent.
- `load_atmospheric_co2()`: reads NOAA CSV with `comment="#"`, strips whitespace from column names, coerces numeric types, filters missing values, constructs datetime index, converts ppm → atm (`×1e-6`), returns 1-D `xr.DataArray`.
- `load_multiobs_flux()`: auto-detects `fgco2` units (kg C m⁻² s⁻¹ or mol C m⁻² s⁻¹) and converts to mol C m⁻² yr⁻¹.
- `harmonize_time()`: converts cftime objects to numpy datetime64[ns] for `xr.merge` compatibility.
- `build_ocean_mask()`: derives 2-D static ocean mask from PISCES NaN pattern — `notnull().any(dim="time")`.
- `broadcast_atm_co2()`: interpolates 1-D NOAA time series to CMEMS monthly time axis, expands to (time, lat, lon) for pixel-wise arithmetic.
- `main()`: assembles all components, aligns to BGC hindcast time axis, merges into `processed_surface.nc`. Prints dimension summary for sanity check.

### Added — `03_compute_flux.py`
- `schmidt_number_co2(sst_degC)`: Wanninkhof (2014) Table 1 polynomial, 4th order in SST [°C].
- `gas_transfer_velocity(wind_speed, Sc)`: Wanninkhof (2014) quadratic wind parameterisation, converts cm hr⁻¹ → m s⁻¹.
- `co2_solubility_K0(sst_degC, sss)`: Weiss (1974) natural-log formulation in T [K] and S [PSU], converts mol L⁻¹ atm⁻¹ → mol m⁻³ atm⁻¹ (×1000).
- `compute_flux(k, K0, spco2_atm, spco2_ocean, ocean_mask)`: bulk formula, converts mol m⁻² s⁻¹ → mol m⁻² yr⁻¹, applies ocean mask.
- `compute_grid_cell_area(lat, lon)`: spherical geometry $A = R^2 \cos\varphi \,\Delta\varphi\,\Delta\lambda$, returns (lat, lon) DataArray [m²].
- `global_integral(flux, cell_area, ocean_mask)`: area-weighted sum over lat/lon at each time step, converts mol C yr⁻¹ → Pg C yr⁻¹.
- `main()`: computes Sc → K0 → k → F → global integral; saves `data/flux_3d.nc` (fgco2, k, K0, Sc) and `output/global_flux.nc` (J_net_PgC, J_net_PgC_annual).

### Added — `04_validate.py`
- `pixel_rmsd()`, `pixel_bias()`: per-pixel RMSD and bias over time dimension.
- `scalar_metrics()`: global Pearson r, RMSD, bias, σ_rec, σ_obs over all finite ocean pixels.
- `global_integral_obs()`: area-weighted MULTIOBS global integral for time series comparison.
- `plot_timeseries()`: two-curve time series (reconstruction vs MULTIOBS), both in Pg C yr⁻¹.
- `plot_rmsd_map()`: 2-panel RMSD and bias map with cartopy Robinson projection. Fallback to `imshow` without cartopy.
- `plot_taylor_diagram()`: polar Taylor diagram — angle = arccos(r), radius = σ_rec/σ_obs. MULTIOBS reference at (0, 1).
- `validation_metrics.csv` saved with global scalar statistics.

### Added — `05_plot_results.py`
- `fig01_flux_vs_co2.png`: twin-axis: J_net [Pg C yr⁻¹] (left) vs atmospheric CO₂ [ppm] (right) over 1993–2026. Monthly values faint fill; annual means solid lines.
- `fig02_annual_flux_map.png`: time-mean air-sea flux map, contourf, Robinson projection.
- `fig03_flux_trend_map.png`: per-pixel Sen's slope of annual flux, Theil-Sen estimator via `scipy.stats.theilslopes` in a pixel loop. Units: mol C m⁻² yr⁻¹ per decade.
- `fig04_delta_pco2_map.png`: time-mean ΔpCO₂ = pCO₂(ocean) − pCO₂(atm) [µatm].
- `fig05_monthly_seasonal_cycle.png`: climatological monthly mean J_net (bar) and atmospheric CO₂ (line) on twin axes.

### Added — `README.md`, `CHANGELOG.md`, `requirements.txt`
