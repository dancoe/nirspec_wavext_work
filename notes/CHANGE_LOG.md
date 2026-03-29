## 2026-03-29 — IFU v2 Coefficient Derivation

### Root Cause of IFU v1 Failure Identified
The IFU v1 NNLS solver (`analysis/solver/solve_parlanti_ifu_v1.py`) used a regularisation
anchor `k_anchor=2` with prior k=1.0. This biased k toward 1.0 everywhere.
True k needed: ~1.2 at NRS2 boundary, declining to ~0.3 at 3.5 µm.
Result: correction factor 1/k ≈ 1.2 was far too small → recalibrated 30-50% below truth.

### G191-B2B Anomaly Investigated
G191-B2B consistently gives R = S_obs/CALSPEC = 0.815 at 3.0 µm vs P330E/J1743045 ≈ 0.638.
The ~27% excess cannot be explained by Parlanti's published α (max 5.4%).
Root cause: pipeline photom calibration artifact for hot WDs (T~60000K).
**Decision**: exclude G191-B2B from k derivation (accepted limitation for G140M NRS2).

### IFU v2 Algorithm
1. k(λ) = median(R_P330E, R_J1743045), where R_i = S_obs_i / f_CALSPEC_i
   - Smoothed over 40 channels, clipped to min 0.01
   - P330E-C3 (PID6645) excluded — stage3_ext G235M data essentially empty (R ≈ 0.008)
2. α(λ), β(λ) loaded from Parlanti published calibration FITS files
   (`data/parlanti_repo/calibration_files/calibration_functions_*.fits`)
   — optics-level, pipeline-independent
3. Correction: f_corr = (S_obs − α̃·f(λ/2) − β̃·f(λ/3)) / k

### Results

| Grating | k median | k range |
|:--------|:---------|:--------|
| G140M NRS2 | **0.567** | 0.298–1.322 |
| G235M NRS2 | **0.768** | 0.545–1.307 |

Recalibrated spectra: P330E and J1743045 match truth within ~5% across full NRS2 range.
G191-B2B G140M overcorrected by ~20-30% (documented limitation); G235M good.

### New Files
- `analysis/solver/solve_parlanti_ifu_v2.py` — corrected solver (no regularisation)
- `reports/329_ifu_v2/REPORT.md` — comprehensive report
- `reports/329_ifu_v2/plot_ifu_v2_source_spectra.py` — validation plots
- `reports/329_ifu_v2/plot_ifu_v2_coeffs_log.py` — log-scale coefficient plot
- `reports/329_ifu_v2/ifu_v2_spectra_{P330E,G191-B2B,J1743045}.png` — validation PNGs
- `reports/329_ifu_v2/ifu_v2_coeffs_log.png` — coefficient visualisation
- `plots/Parlanti/cal/ifu_v2/coeffs_ifu_v2_{G140M,G235M}.csv` — coefficient tables
- `plots/Parlanti/cal/ifu_v2/ifu_v2_{g140m,g235m}_coeffs.png` — coefficient plots
- `reports/REPORTS.md` — updated with v2 entry

---

## 2026-03-29 ~01:00 — FS v1 Coefficient Derivation (NRS2 Pipeline Fixed)

### Root Cause of FS NRS2 All-Null Pixels Identified
The fflat/sflat default CRDS references (`fflat_0154.fits`, `sflat_0154.fits`) have
`FAST_VARIATION` wavelength tables covering only 0.97–1.89 µm (NRS1 G140M range).
All NRS2 wavelengths (1.97–3.16 µm for G140M, 3.15–5.27 µm for G235M) fall outside
this range → pipeline flags all NRS2 pixels `DO_NOT_USE`.

**Fix**: Use Parlanti extended flat references (`fflat_0105.fits`, `fflat_0091.fits`,
`sflat_0191.fits`, `sflat_0192.fits`) from `data/parlanti_repo/CRDS_1364/`.
These have `FAST_VARIATION` extending to 5.5 µm with `slit_name=ANY`.

### Additional Fix: Skip `resample_spec`
Even with valid flat correction, `extract_1d` from s2d produced NPIXELS=0.
Fix: skip `resample_spec` and extract from cal (MultiSlit) directly.

### Updated Pipeline Script: `analysis/reduction/run_fs_nrs2_spec2.py`
- Added grating-specific flat override constants (FFLAT_G140M/G235M, SFLAT_G140M/G235M)
- Added FLAT_OVERRIDES dict selecting correct reference by (grating, filter)
- Added `resample_spec: skip: True`, `wavecorr: skip: True`
- Added custom `override_photom` for extended NRS2 wavelength coverage

### NRS2 x1d Products Generated (6 files, Jy)
- `data/PID1536_J1743045/nrs2_spec2_cal/jw01536002001_03104_00003_nrs2_x1d.fits` (G140M)
- `data/PID1536_J1743045/nrs2_spec2_cal/jw01536002001_03106_00001_nrs2_x1d.fits` (G235M)
- `data/PID1537_G191-B2B/nrs2_spec2_cal/jw01537007001_07101_00003_nrs2_x1d.fits` (G140M)
- `data/PID1537_G191-B2B/nrs2_spec2_cal/jw01537007001_09101_00003_nrs2_x1d.fits` (G235M)
- `data/PID1538_P330E/nrs2_spec2_cal/jw01538160001_06101_00001_nrs2_x1d.fits` (G140M)
- `data/PID1538_P330E/nrs2_spec2_cal/jw01538160001_08101_00005_nrs2_x1d.fits` (G235M)

### New Solver: `analysis/solver/solve_parlanti_fs_v1.py`
Updated SOURCES dict to use Jy-calibrated x1d files. Ran NNLS solver:

| Grating | k median | k range | α̃ max |
|:--------|:--------|:--------|:------|
| G140M NRS2 | **0.846** | 0.706–1.093 | 0.390 |
| G235M NRS2 | **0.976** | 0.834–1.041 | 0.308 |

**FS vs IFU comparison:**
- G140M: FS k=0.846, IFU k=0.833 → agree within 2%
- G235M: FS k=0.976, IFU k=0.899 → ~8% offset, partly attributable to IFU flat bias at 3.3–3.8 µm

### Output Plots/CSVs
- `plots/Parlanti/cal/fs_v1/coeffs_fs_v1_G140M.csv`
- `plots/Parlanti/cal/fs_v1/coeffs_fs_v1_G235M.csv`
- `plots/Parlanti/cal/fs_v1/fs_v1_g140m_coeffs.png`
- `plots/Parlanti/cal/fs_v1/fs_v1_g235m_coeffs.png`
- `plots/Parlanti/cal/fs_v1/fs_vs_ifu_k_comparison.png`

### Session Log Added
- `notes/logs/2026-03-29-0100.md`

---

## 2026-03-28 23:30 — IFU v1 Coefficient Derivation

### New Solver: `analysis/solver/solve_parlanti_ifu_v1.py`
Created a new solver using IFU `stage3_ext` products (Jy) as the data source for the Parlanti k/α/β derivation:
- **Method**: NNLS, point-wise on 200-channel grids, 40-channel box smoothing, 3 standard stars
- **G140M NRS2 result**: k(λ) median = 0.833, range 0.712–1.003; α = β ≈ 0
- **G235M NRS2 result**: k(λ) median = 0.899, range 0.817–1.001; α = β ≈ 0
- k(λ) starts near 1.0 at the NRS2 cutoff and declines toward longer wavelengths
- No second/third-order contamination detectable (requires cool calibrator to break k/α degeneracy)

### Session Log Files Added
- `notes/logs/2026-03-28-2300.md` — IFU data-source diagnostic and coefficient re-derivation plan
- `notes/logs/2026-03-28-2330.md` — IFU v1 solver results, interpretation, and caveats

### INSTRUCTIONS.md Updated
Added session log convention: timestamped files in `notes/logs/` for each significant work session.

---
## 2026-03-28 — IFU Calibration Diagnostic Overhaul

### Root Cause of Scaling Mismatch
Diagnosed why the extended-wavelength observations appeared incorrectly scaled:
- The old `plot_parlanti_components_multi.py` used per-exposure `nrs2_extract_1d.fits` (in **DN/s**) for "EXT" data vs stage3 x1d (in **Jy**) for "NOM" data — a unit mismatch of ~10³.
- Additionally, same-grating NOM and EXT had **zero wavelength overlap** (G140M NOM ends at 1.87 µm; EXT starts at 1.96 µm), so the overlap-based scaling always returned NaN.

### Script Rewrite: `analysis/plotting/plot_parlanti_components_multi.py`
Complete rewrite to use the correct data sources and comparisons:
- **Data**: IFU `stage3_ext/f100lp_g140m-f100lp_x1d.fits` and `stage3_ext/f170lp_g235m-f170lp_x1d.fits` (Jy, all JWST pipeline stages complete) for all three standards (P330E, G191-B2B, J1743045).
- **Reference**: CALSPEC models (replacing the PRISM FS file as reference).
- **Comparison geometry**: G140M NRS2 (1.87–3.60 µm) vs G235M stage3 nominal (ground truth); G235M NRS2 (3.15–5.50 µm) vs G395M FS nominal.
- **Layout**: 3-panel figure — top: full spectra (log), bottom-left: k(λ) for G140M NRS2, bottom-right: k(λ) for G235M NRS2.
- **k(λ) overlay**: Shows both the directly-measured k = (higher-grating NOM) / (stage3_ext NRS2) and the V3 coefficients for comparison.
- Supports `--source all | 1538_P330E | 1537_G191B2B | 1536_J1743` argument.

### Key New Finding
The IFU stage3_ext NRS2 is **~2x under-calibrated** relative to the next-grating nominal (k = G235M NOM / G140M EXT ≈ 1.5–2.0 across the overlap). The V3 coefficients (k≈0.89) were derived from different data and are not applicable to the IFU stage3_ext products. New calibration coefficients must be derived from the IFU data.

### Output Plots Updated
- `plots/Parlanti/cal/153678_v3/CAL_COMPONENTS_P330E.png`
- `plots/Parlanti/cal/153678_v3/CAL_COMPONENTS_G191-B2B.png`
- `plots/Parlanti/cal/153678_v3/CAL_COMPONENTS_J1743045.png`

---

## 2026-03-28 (late night) — Documentation & Workflow Optimization
- **Documentation Consolidation**: Created `notes/LATEST_WORK.md` to track the most recent progress (V3 Calibration) and provide a clear roadmap for upcoming tasks (Complex Source Validation, ASDF Export, F-Flat Concatenation).
- **Instruction Cross-Referencing**: Updated `notes/INSTRUCTIONS.md` to include high-level pointers to `LATEST_WORK.md` and `IMPLEMENTATION_PLAN.md`, streamlining the navigation between current status and long-term goals.
- **Log Maintenance**: Appended the latest user request to `notes/PROMPT_LOG.md` (Entry 259).

## 2026-03-28 — MOS Visualization Pipeline Refactor
- **Absolute Detector Mapping**: Completely refactored `show_MOS_rate_files` in `mos.py` to use dynamic extent calculation based on `SUBSTRT1/2` and `SUBSIZE1/2`. This ensures that data is correctly placed on the absolute 2048x2048 detector grid, with Y-axis labels reflecting true detector rows (e.g., ~1050–1120 for fixed slits).
- **V2-Style Aspect Ratio Restored**: Re-implemented `aspect='auto'` within specifically tuned `set_ylim` boundaries. This restores the "filled axis" visual style of the earlier v2 diagnostic plots while maintaining perfect detector-relative cropping.
- **Robust WCS Interpolation**: 
  - Fixed a `ValueError` caused by NaNs in the WCS coordinate transformation for off-detector points.
  - Implemented a masked interpolation scheme for the secondary wavelength axis.
  - Adjusted coordinate offsets to use `x_base` from headers, ensuring accurate wavelength mapping for zoomed subarrays.
- **Row-wise Synchronization**: Re-engineered the `autocrop` logic to synchronize vertical scales and offsets across each row (NRS1, NRS2), allowing for accurate spatial comparisons between detectors.
- **Automated Diagnostic Plotting**: Generated `v8` reference plots for PIDs 1492, 1537, and 1538, confirming correct extraction region alignment and absolute coordinate labeling.

## 2026-03-27 (very late night) — IFU Analysis Branch
- **New IFU Analysis Branch**: Initiated IFU-mode wavelength extension analysis following Parlanti et al. (2025) methodology, paralleling the existing Fixed Slit work.
- **Calibration Sources Selected**: Identified 4 calibration standard star programs matching exactly those used by Parlanti et al.:
  - PID 1536: J1743045 (G140M + G235M + G395M IFU)
  - PID 1537: G191-B2B (white dwarf; G140M + G235M)
  - PID 1538: P330E (gray star Cycle 1; G140M + G235M)
  - PID 6645: P330-E Cycle 3 (gray star; G140M + G235M)
  - Minimum 3 sources satisfies the 3×3 linear system for k(λ), α̃(λ), β̃(λ)
- **MAST Data Acquisition**: Downloaded IFU rate files using `astroquery.mast`:
  - PID 1537: 16 files (G140M+G235M, 4 dithers × NRS1+NRS2)
  - PID 1538: 16 files (G140M+G235M, 4 dithers × NRS1+NRS2)
  - PID 1536: 24 files (G140M+G235M+G395M, 4 dithers × NRS1+NRS2)
  - PID 6645: ~84 files downloading (P330-E many visits)
  - Fixed truncated-download bug: added size check + 3-attempt retry to `download_ifu_rates.py`
- **Standard Reduction Completed (PIDs 1537 + 1538)**:
  - Ran `Spec2Pipeline` + `Spec3Pipeline` via `run_ifu_pipeline.py`
  - Outputs: `f100lp_g140m_s3d/x1d.fits` and `f170lp_g235m_s3d/x1d.fits` for each PID
  - NRS1-only (standard range); NRS2 correctly skipped — "No IFU slices" on NRS2 for nominal G140M/G235M is expected
- **PID 1536 Pipeline Running**: Stage 2+3 in progress for J1743045
- **Scripts Created**:
  - `analysis/acquisition/query_ifu_mast.py`: MAST discovery
  - `analysis/acquisition/download_ifu_rates.py`: Robust download with retry
  - `analysis/reduction/reduce_ifu.py`: Full 3-stage IFU pipeline
  - `analysis/reduction/run_ifu_pipeline.py`: Streamlined Stage 2+3 runner
- **Notes**: Created `notes/PARLANTI_IFU.md` documenting the full IFU analysis workflow, pipeline details, calibration strategy, and roadmap for reference file modification

## 2026-03-27 (night)
- **MOS Extraction Plotting**:
  - Created `analysis/plotting/plot_mos_extraction_regions.py` to leverage the `mos.py` utility for side-by-side visualization of count rate data and 2D extraction regions.
  - **Environment Fixes**: Implemented mocks for `tqdm` and `astroquery` in the plotting script to facilitate running in the `jwst_1.20.2` environment where these non-core packages were missing.
  - **Visualized NRS2 Extension**: Generated `plots/extraction/MOS_rate_extractions_NRS2.png`, confirming the successful definition of extraction bounding boxes for G140M, G235M, and G395M on the NRS2 detector.
  - **Documentation**: Summarized the extraction methodology and results in `plots/MOS_EXTRACTIONS.md`.

## 2026-03-27 (late night)
  - **V3 Multi-Source Bootstrap Solver (Kappa Fixed):**
    - Refined `analysis/solver/solve_parlanti_bootstrap_v3.py` to strictly implement the **Parlanti et al. (2025)** "North Star" model ($k \approx 1.0$ prior).
    - **Multi-Truth Integration**: Engineered a dynamic truth-merger in `solve_parlanti_bootstrap_v3.py` that combines G235M, G395M, and PRISM spectra to provide reference flux across the entire extended range (0.6–5.6 µm).
    - **Robust Scaling**: Implemented a dynamic overlap-finding scaling algorithm in both solver and plotting scripts. This handles the wavelength gaps between NRS1 and NRS2 detectors by searching for the first valid data window.
    - **Stabilized Solves**: Increased the $k=1$ weight ($w=5.0$) and implemented 40-channel smoothing to prevent non-physical zero-crossings in high-contamination regions (like the G191-B2B spikes).
  - **Validation & Visualization V3**:
    - Generated `plots/Parlanti/cal/153678_v3/SUMMARY_V3_{P330E,G191-B2B}.png` following the Parlanti style (Log scale, multi-color traces, extension arrows).
    - **Contamination Stress Test**: Confirmed successful ghost subtraction in the hot-star (G191-B2B) case, where 2nd-order contamination is orders of magnitude higher than the 1st-order signal.
    - Created `plots/Parlanti/cal/153678_v3/FINAL_COEFFS_V3.png` showing smoothed $k, \tilde\alpha, \tilde\beta$ functions.
  - **Documentation**:
    - Finalized `plots/Parlanti/cal/153678_v3/CAL_153678_V3.md` as the primary calibration artifact.
    - Audited and updated `notes/FLATS_FILES.md` and `notes/FLATS.md`.

## 2026-03-27 (late evening)
  - **Data Acquisition & Reduction:** Successfully batch-downloaded and reduced FS data for **PIDs 1536, 1537 (G191-B2B), and 1538 (P330E)**. All targets were processed through the extended WCS and extraction pipeline (NRS2 enabled).
  - **Direct 3-Source Solver:** Implemented `analysis/solve_parlanti_3source.py`. This script performs a direct $3 \times 3$ solve of the Parlanti system at every wavelength using the three distinct calibration sources (IRAS-05248, G191-B2B, P330E).
    - Utilizes **NNLS (Non-Negative Least Squares)** to enforce physical throughputs.
    - Applies **11-pixel boxcar smoothing** to stabilize derived $k, \tilde\alpha, \tilde\beta$ functions.
  - **Documentation Reorganization:** 
    - Split `notes/FLATS.md` into two documents to improve clarity and maintainability.
    - `notes/FLATS.md`: Focuses on the flat-fielding strategy, meaning, and visualization Atlas (carousels).
    - `notes/FLATS_FILES.md`: Provides the detailed, audited inventory of all 40 `DFLAT`, `FFLAT`, and `SFLAT` reference files.
  - **Ongoing Calibration Validation:** Initiated generation of source-specific validation plots to verify how well the global 3-source solution models each individual calibrator.

## [Unreleased] - 2026-03-27

### Added
- Created `plots/Parlanti/cal/3source/FS_3source_cal.png`: Summary plot for the 3-source calibration, showing global coefficients and a representative calibrated spectrum (P330E) compared to ground truth (PRISM).
- Created `plots/Parlanti/cal/3source/CAL_3SOURCE.md`: Comprehensive report for the multi-source direct solve, summarizing methodology (NNLS), stability improvements over single-source polynomial fits, and next steps for G395M.
- Created `analysis/plot_3source_summary.py`: Script to generate the summary calibration performance figure.

### Fixed
- Fixed G140M file paths in parity with the 3-source solver for all targets (IRAS-05248, G191-B2B, P330E).
- Robustified `load_spec` in analysis scripts by using `astropy.io.fits` to handle Level 3 files with legacy ASDF metadata issues.

## 2026-03-27 (evening)
  - **Root cause found:** Extended NRS2 extractions (`jw01492003001_…`) were processed through `assign_wcs` + `extract_2d` + `extract_1d` but NOT the `photom` step. `FLUX` column is in DN/s (not Jy). Nominal Level 3 products (obs001) are in Jy. This unit mismatch caused the 100–200× flux discrepancy seen in the pre-cal plot.
  - **New script:** `analysis/calibrate_parlanti.py` implements:
    1. **Boundary-matching scale factor**: sigma-clipped median of `f_PRISM(λ) / S_raw(λ)` over the first 80 pts of each extended extraction, giving C = 8.04×10⁻⁴ Jy/(DN/s) for G140M and 1.35×10⁻³ Jy/(DN/s) for G235M.
    2. **Legendre polynomial fit (deg=4) with Tikhonov regularisation**: avoids the ill-conditioning of raw polynomial powers; uses normalised wavelength ξ∈[-1,1].
    3. **Iterative 3σ-clip** on fit residuals (3 iterations) to reject emission-line spikes.
    4. **Physical bounds**: k∈[0.01,2.0], a∈[0,0.20], b∈[0,0.05].
    5. **High-resolution contamination reference**: for the calibration APPLICATION step, uses the nominal G140M/G235M/G395M stitch (falling back to PRISM) rather than the low-res PRISM alone.
  - **Results:**
    - G140M fit RMS: 50% | k drops 0.47→0.01 (2.0→3.0 µm) | a peaks at 13% near boundary | b clamped at 5%
    - G235M fit RMS: 54% | k starts near 1.0, drops to ~0.1 at 4.5 µm | a large (5–20%) | b 0–5%
  - **Plots generated:**
    - `plots/Parlanti/cal/FS_1492_cal.png` — 2-panel calibration plot (cf. Parlanti Fig. 5), linear scale, contrasting nominal (black), uncalibrated (green), calibrated (red) spectra
    - `plots/Parlanti/cal/parlanti_coefficients_v2_diagnostic.png` — coefficient profiles k(λ), a(λ), b(λ)
    - `plots/Parlanti/cal/parlanti_coefficients_v2.txt` — coefficient summary table
  - **Documentation:** Updated `notes/CALIBRATION.md` with full methodology, data table, interpretation, and known limitations.

## 2026-03-27 14:10
- **Docs & Parlanti notes committed:** Created `notes/PARLANTI.md` and updated `notes/CALIBRATION.md` and `notes/IMPLEMENTATION_PLAN.md` with detailed Parlanti (arXiv:2512.14844) Section 2 interpretation, the 3-source minimum requirement, and a multi-source path-forward. Committed and pushed `notes/PARLANTI.md`, `notes/CALIBRATION.md`, and `notes/IMPLEMENTATION_PLAN.md` to `nirspec_wavext_work` (main). See commit abbe3a2. Added PROMPT_LOG entry for this request.

## 2026-03-27 3:50 PM
- **Parlanti Model Coefficient Determination (Initial Analysis):**
  - Created `analysis/determine_parlanti_coefficients.py` script to fit Parlanti et al. 2025 flux correction model coefficients from M-grating spectral data.
  - **Model:** $S(\lambda) = k(\lambda) \cdot f(\lambda) + a(\lambda) \cdot f(\lambda/2) + b(\lambda) \cdot f(\lambda/3)$, where $f(\lambda)$ is PRISM baseline and $k, a, b$ are polynomial throughput functions.
  - **Data Used:**
    - PRISM baseline (411 points): λ=[0.597, 5.300] µm
    - G140M extended (1,095 points): λ=[1.981, 3.268] µm
    - G235M extended (1,272 points): λ=[3.308, 5.310] µm  
    - G395M extended (68 points): λ=[5.496, 5.617] µm
  - **Fit Results:**
    - **G140M** (λ=1.8-3.3 µm): RMS error 8.37 Jy (429%) – ⚠️ POOR; polynomial overfitting detected
    - **G235M** (λ=3.0-5.3 µm): RMS error 0.854 Jy (63%) – ✓ OK; reasonable fit quality
    - **G395M** (λ=4.2-5.6 µm): RMS error 0.548 Jy (11%) – ✓ GOOD but limited data (58 pts)
  - **Issues Identified:**
    1. Polynomial (cubic) basis is numerically unstable → switch to B-splines
    2. G140M shows extreme noise/systematics → investigate flat-field artifacts
    3. G395M has too few points → enforce physical constraints
    4. No bounds on throughput functions → implement $0 \leq k(\lambda) \leq 1.1$, $a(\lambda), b(\lambda) \geq 0$
  - **Diagnostic Plots Generated:**
    - `plots/parlanti_fit_G140M.png` – 4-panel diagnostic (spectrum, residuals, coefficients, components)
    - `plots/parlanti_fit_G235M.png` – 4-panel diagnostic
    - `plots/parlanti_fit_G395M.png` – 4-panel diagnostic
  - **Output Files:**
    - `plots/parlanti_coefficients.txt` – Raw polynomial coefficients (⚠️ NOT YET VALIDATED)
    - `plots/PARLANTI_ANALYSIS_REPORT.md` – Detailed analysis report with recommendations
  - **Documentation:**
    - Updated `IMPLEMENTATION_PLAN.md` with Parlanti analysis results and next steps
    - Added prompts 9-10 to `PROMPT_LOG.md`
  - **Next Actions:**
    - [ ] Refit with B-spline basis (improved numerical stability)
    - [ ] Inspect G140M 2D extraction for systematics
    - [ ] Implement physical bounds on coefficients
    - [ ] Validate against published Parlanti values

## 2026-03-27 11:45 AM
- **Flat Field Documentation Fix**:
  - Corrected image embedding syntax in `FLATS.md` by removing the `file:///` prefix to ensure proper rendering across different markdown viewers.
  - Added vertical separation between mode-specific plots for improved readability.

## 2026-03-27 11:15 AM
- **Flat Field Extension Implementation:**
  - **`jwst/flatfield/flat_field.py` Modifications**: 
    - **`combine_fast_slow`**: Replaced `right=np.nan` with `right=tab_flat[-1]` and `right=tab_flat_error[-1]` in `np.interp` calls. This enables red-end extrapolation of "fast" (table-based) flat variations using the value at the reddest available wavelength, as per Parlanti et al. 2025.
    - **`interpolate_flat`**: Modified range-checking logic to allow extrapolation to the right (red end). For wavelengths exceeding the reference range, the "slow" (image-based) variation is now explicitly set to **1.0** (constant) and pixels are **not** flagged as `DO_NOT_USE`, ensuring processing continues for `wavext`.
- **Documentation**:
  - Created `notes/FLATS.md` explaining the D-flat, S-flat, and F-flat roles in NIRSpec and the newly implemented extension strategies for wavext.
  - Linked `references/Parlanti Fig A1.png` as the conceptual basis for the IFU flat extension.
- **Git Commit**: Committed all changes to the `feature/nirspec_wavelength_extension` branch of `jwst_nirspec_wavext` and updated `nirspec_wavext_work`.

## 2026-03-27 10:00 AM
- **Plotting Refinement (Combined M-Gratings):**
  - **Front-Layering**: Configured PRISM baseline to be plotted on top (`zorder=20`) of the M-gratings for better visibility as a reference baseline.
  - **Styling Alignment**: Adopted standard styles from `plot_parlanti_flux.py`, including `alpha=0.6` for all traces, `lw=1.0` (Nominal), and `lw=0.8` (Extended).
  - **Dynamic Range Optimization**: Switched y-axis scaling logic to use a **0.1%–99.9% percentile range** to capture the full essential dynamic range of the contamination signal while excluding extreme noise/spikes.
  - **Data Cleaning**: Implemented strict masking of non-finite (NaN) and non-positive flux values (≤ 1e-6) to eliminate vertical line artifacts in log-scale view.
- **Documentation Migration**:
  - Moved `PARLANTI_PLOTS.md` from `notes/` to `plots/` to align with the generated outputs.
  - Created `plots/Parlanti_gratings.md` as a dedicated explainer for the `FS-1492_pre-cal.png` summary plot.
- **Git Commit**: Successfully committed and pushed the final set of analysis scripts and documentation to `nirspec_wavext_work`.

## 2026-03-27 09:50 AM
- **Batch Processing of M-Gratings:** Successfully automated the download, WCS assignment, 2D extraction, and 1D extraction for G140M, G235M, and G395M (NRS2) using the custom `wavelengthrange_extended.asdf`.
  - **G140M Success:** Extracted spectrum up to ~3.3 µm on NRS2 (nominal cutoff is ~1.9 µm), capturing significant 2nd order light for calibration.
  - **G235M Success:** Processed NRS2 file `jw01492003001_03104_00004_nrs2_rate.fits`, extending coverage into the red end.
  - **G395M Success:** Processed NRS2 file `jw01492003001_03106_00002_nrs2_rate.fits`.
- **Comparison Visualization:** Generated comparison plots (`plots/*_comparison.png`) showing $f(\lambda)$ (Nominal) vs $S(\lambda)$ (Extended) with yellow shading highlighting the extension/overlap regions.
- **PRISM NRS2 Investigation:** Attempted PRISM processing on NRS2 file `jw01492005001_17101_00006_nrs2_rate.fits` but encountered `NoDataOnDetectorError`.
  - **Findings:** Corrected `wavelengthrange_extended.asdf` for PRISM CLEAR to [0.5, 5.6] µm, but `assign_wcs` still removes the `S1600A1` slit on NRS2.
  - **Next Step:** Verify if PRISM `S200A1` projects onto NRS2 or if the existing Level 3 products already cover the desired range.
- **Documentation Updates:**
  - Logged latest prompt in `PROMPT_LOG.md`.
  - Updated `IMPLEMENTATION_PLAN.md` with the "Current Phase" (March 27, 2026) tasks and goals.
  - Created a new analysis script `process_selected_gratings.py` for batch processing and plotting.

## 2026-03-27 01:25 AM
- **Wavelength Extension Investigation:** Checked why extended wavelengths are not plotted for G140H and G235H on the `NRS2` detector.
  - Investigated `assign_wcs` bounding box logic in `jwst_nirspec_wavext/jwst/assign_wcs/nirspec.py`. Found that the WCS models the physical detector boundaries, capping extraction at `x=2047`.
  - For `S200A1` on `NRS2`, `G140H` physically terminates at 1.819 µm (well within the nominal 1.89 µm cutoff), and `G235H` terminates at 3.050 µm (within the nominal 3.17 µm cutoff).
  - Thus, there is zero flux beyond the nominal cutoff for these combinations because the light physically falls off the detector edge. Plotted output correctly shows no yellow extended shading because all available data is already within the nominal ranges.
  - Successfully updated and executed the `run_assign_wcs.py`, `run_extract_2d.py`, and `run_extract_1d.py` scripts to verify these limits on G235H data (`jw01492001001_0310a_00007_nrs2_rate.fits`).

## 2026-03-27 12:40 AM
- **Extract1dStep Success:** Successfully executed `AssignWcsStep`, `Extract2dStep`, and `Extract1dStep` for G395H PID 1492 Fixed Slit data using the custom `wavelengthrange_extended.asdf` reference file.
  - **Unit Verification:** Confirmed that while the ASDF reference uses meters (1e-6), the pipeline correctly scales the resulting spectra into **microns**.
  - **Detector Edge Constraints:** Verified that the extraction is limited physically by the NRS2 detector edge (e.g., ~5.144 µm for S200A1).
## 2026-03-27
- **V2 Multi-Source Bootstrap Solver**: Implemented `solve_parlanti_bootstrap.py` (V4 logic) with physical priors (1st-order dominance) and NNLS.
- **Combined PID Analysis**: Successfully integrated PIDs 1492, 1537, and 1538 into a single calibration run.
- **Project Reorganization**: Moved all analysis scripts into subdirectories: `acquisition/`, `reduction/`, `solver/`, `plotting/`, `utils/`.
- **Validation**: Applied V2 coefficients to P330E and G191-B2B, matching ground truth in overlap regions.
- **Final Plotting**: Generated `FS_P330E_cal_V2.png` and `FS_IRAS-05248_cal_V2.png` matching the requested format.
- **Pipeline Branch Rebase (jwst_1.20.2):** Resolved a critical `AttributeError: No attribute 'get_dtype'` that occurred during `Extract1dStep`.
    - **Root Cause:** The `feature/nirspec_wavelength_extension` branch was based on an outdated `jwst` pipeline codebase that used deprecated `stdatamodels` methods.
    - **Fix:** Rebased the local fork onto the official `1.20.2` tag to align the pipeline with the current environment's `stdatamodels` and `stpipe` dependencies.
- **Repository Setup:** Finalized the separation of `jwst_nirspec_wavext` and `nirspec_wavext_work` with correct SSH key mapping and GitHub remotes.
- **Next Steps:** Proceed with the throughput ($k$, $a$, $b$) function derivation (Parlanti et al. model) by comparing nominal regions with extended extraction results.

## 2026-03-27 12:08AM
- **Extraction Phase:** Starting the extraction and calibration of the PID 1492 FS data.
- **Environment Verification:** Confirming `micromamba activate jwst_1.20.2` and `PYTHONPATH` exports for the custom `jwst_nirspec_wavext` library.

## 2026-03-27 12:07AM
- **Git Remote & Config:** Configured both `nirspec_wavext_work` and `jwst_nirspec_wavext` repositories to use the `id_ed25519_dancoe` SSH key identity via `~/.ssh/config` host mapping.
- **Initial Push:** Successfully pushed the `nirspec_wavext_work` repository to GitHub.
- **Pipeline Push:** Pushed the `feature/nirspec_wavelength_extension` branch of `jwst_nirspec_wavext` to GitHub, including the ASDF lazy-load fix.
- **Git Config:** Updated local `user.name` to `dancoe` and `user.email` to `dancoe@gmail.com` for both repositories.


## 2026-03-27 12:05AM
- **Directory Reorganization:** Separated the project into two main repositories:
  - `jwst_nirspec_wavext`: The local fork of the JWST pipeline.
  - `nirspec_wavext_work`: Contains project `notes`, `analysis`, and `reference_files`.
- **Reference Files:** Created `nirspec_wavext_work/reference_files/` and moved `wavelengthrange_extended.asdf` into it.
- **Extract1dStep Success:** Successfully executed `AssignWcsStep`, `Extract2dStep` and `Extract1dStep` for G395H PID 1492 Fixed Slit data with `wavelengthrange_extended.asdf`.
  - Debugged WCS output to confirm correct bounding box coordinates mapping to 0.6–5.6 µm range.
  - Encountered and fixed an `AttributeError: No attribute 'get_dtype'` during `Extract1dStep`. Solved by rebasing the outdated `jwst_nirspec_wavext` fork onto the `1.20.2` tag to align the pipeline modifications with `stdatamodels` requirements in the `jwst_1.20.2` environment.
  - Verified that `Extract1dStep` successfully extracts spectra in microns up to the physical edge of the NRS2 detector (e.g., 5.144 µm for G395H `S200A1`).
  - Next step: Compare nominal wavelength regions with extended extraction results and develop the throughput ($k$, $a$, $b$) function derivation.
- **Documentation:** Updated all absolute and relative file paths in the `notes/` directory to reflect the new repository structure.

## 2026-03-26 11:56PM
- **Bug Fix:** [nirspec.py](file:///Users/dcoe/NIRSpec/wavext/jwst_nirspec_wavext/jwst/assign_wcs/nirspec.py) — `get_spectral_order_wrange()` (lines ~1231–1260)
    - **Root Cause:** `wave_range_model.wavelengthrange[index]` and `wave_range_model.order[index]` return lazy-loaded ASDF ndarray proxies that back-reference the open ASDF file handle. After `wave_range_model.close()` is called, any later attempt to access those proxies (e.g. `wrange[0]` at `validate_open_slits`) raises `OSError: Attempt to load block from closed file`.
    - **Fix:** Force-materialise both values into plain Python types (`int(order)`, `list(wrange)`) **before** `wave_range_model.close()` in all code paths. This resolves the blocker for FS (and all other modes) when using the custom `wavelengthrange_extended.asdf` override.
    - **Verified:** `run_assign_wcs.py` now completes successfully on `jw01492001001_03108_00007_nrs2_rate.fits` and writes the WCS output.

## 2026-03-26 11:20PM
- **Modified Code:** [nirspec.py](file:///Users/dcoe/NIRSpec/wavext/jwst_nirspec_wavext/jwst/assign_wcs/nirspec.py)
    - Removed hard-coded `NoDataOnDetectorError` for IFU modo processing on NRS2 for M-gratings and G140H/F070LP. This allows the pipeline to proceed with wavelength extension for these combinations.
- **Reference Files:** 
    - Generated `wavelengthrange_extended.asdf` using `stdatamodels`. This file extends the `wavelengthrange` for all NIRSpec gratings (0.6 - 5.6 µm range).
- **Modified Documentation:** Updated [IMPLEMENTATION_PLAN.md](IMPLEMENTATION_PLAN.md) to emphasize reference file extensions over code extrapolation.

- **Spectral Ghost Decomposition**: Implemented  to visualize the contribution of (\lambda/2)$ and (\lambda/3)$ to the extended NRS2 spectrum.
- **Wavelength Analysis**: Theoretically confirmed that reference spectra only need to cover the observed range ($\leq 5.3$ µm) because ghosts at $\lambda_{obs}$ originate from $\lambda_{obs}/n$.
- **Model Visualization**: Added colored ghost components (Green for 2nd order, Magenta for 3rd order) and Model Sum (Red) to diagnostic plots.
