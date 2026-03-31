# Analysis Plan v7 — Full Flux Calibration & Cross-Validation

**Date:** March 30, 2026  
**Project:** NIRSpec Wavelength Extension Calibration  
**Objective:** Transition to native Jy-unit extractions for NRS2 extended regions using custom `photom` reference files, and implement systematic hold-out cross-validation to assess solver reliability.

---

## 1. Motivation (The "DN/s Problem")

The v6 analysis demonstrated that the all-source joint solve (CALSPEC + Science Targets) provides stable coefficients $k, \alpha, \beta$ for both FS and IFU modes. However, a significant technical hurdle remains:

- **NRS2 Units:** Extended extractions currently remain in **DN/s** (counts per second) because the standard pipeline `photom` step doesn't cover the extended wavelengths.
- **Normalization Sensitivity:** In v6, we bridged this gap by scaling DN/s to Jy using a median scale factor derived from wavelength overlap. This introduces a potential source of error if the overlap region is noisy or very small.
- **Goal for v7:** Use **custom extended `photom` reference files** to re-process all NRS2 data. This will output spectra directly in **Jy**, removing the arbitrary boundary-matching scale factor and making the pipeline more self-consistent.

---

## 2. Key Objectives (v7)

### 1. Enable Full Flux Calibration (Jy)
- Use `analysis/create_photom_fs_ext.py` to generate a custom `photom` reference file that includes sensitivity data for the extended wavelength ranges (~1.8–3.5 µm for G140M, ~3.0–5.5 µm for G235M).
- Reprocess PID 1492 and standard star data using these reference file overrides.
- Confirm that the output `x1d` spectra are in Jy and match nominal spectra in overlap regions without manual scaling.

### 2. Systematic Hold-out Cross-Validation
- In v6, all sources were used to solve for a single set of coefficients.
- In v7, we will implement a "leave-one-out" loop:
    - Solve for coefficients using $(N-1)$ sources.
    - Test the solution on the $N$-th (held-out) source.
    - Quantify the residual between the corrected held-out spectrum and its "truth" (CALSPEC or adjacent-grating nominal).
- This will identify if any specific source (e.g., the cool star NGC2506-G31 or the science target PID 1492) is biasing the solution.

### 3. Solver Refinement (Parlanti v7)
- Update the `solve_parlanti_fs_v6.py` script (v7 version) to remove the normalization logic.
- Directly read the new Jy-unit files.
- Compare v7 coefficients to v6 to quantify the impact of the DN/s scale factor approximation.

---

## 3. Reference Files & Reprocessing

The v7 work relies on the following custom reference files (located in `reference_files/`):

| File Type | Filename | Purpose |
| :--- | :--- | :--- |
| **Wavelength Range** | `wavelengthrange_extended.asdf` | Extends detection bounding boxes to NRS2 edges. |
| **Photom** | `jwst_nirspec_photom_fs_nrs2_ext.fits` | Enables flux calibration for extended wavelengths. |

Reprocessing will be done via:
```bash
# Example Step 2 override
micromamba run -n jwst_1.20.2 python scripts/run_fs_nrs2_spec2.py \
    --override_photom reference_files/jwst_nirspec_photom_fs_nrs2_ext.fits
```

---

## 4. Expected Results

1.  **Unit Consistency:** All spectra used in the NNLS loop will be in Jy natively.
2.  **Solver Stability:** Hold-out residuals should be $< 10\%$ for FS M-gratings across all types (WD to G1V).
3.  **Coefficient Convergence:** $k(\lambda)$ should more naturally approach 1.0 (or the expected throughput ratio) without manual normalization.

---

## 5. Deliverables

- `reports/v7/fs_v7/REPORT_fs_v7.md`: Final v7 validation report.
- `results/v7/calib_v7_fs_*.fits`: Refined calibration coefficients.
- `reports/v7/fs_v7/plots/fs_v7_crossval_residuals.png`: New diagnostic showing hold-out performance.

---

*Created: 2026-03-30*
*Author: Antigravity*
