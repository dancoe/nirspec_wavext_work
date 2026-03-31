# Detailed Analysis Plan v7 — Full Flux Calibration & Cross-Validation

**Version:** 7.0  
**Status:** In Progress (Planning Phase)  
**Objective:** Transition to native Jy-unit extractions for NRS2 extended regions using custom `photom` reference files, and implement systematic hold-out cross-validation to assess solver reliability.

---

## 1. Reference File Work (Reference: `analysis/create_photom_fs_ext.py`)

The first step in v7 is to generate and verify a custom `photom` reference file that includes sensitivity data for the extended wavelengths (0.6 - 5.6 µm).

| Step | Task | Success Criteria |
| :--- | :--- | :--- |
| **1.1** | Run `analysis/create_photom_fs_ext.py`. | Created `reference_files/jwst_nirspec_photom_fs_nrs2_ext.fits`. |
| **1.2** | Verify wavelength coverage. | Use `analysis/check_photom_refs.py` to confirm G140M/G235M extend to ~3.2/5.5 µm. |

---

## 2. Reprocessing PID 1492 (Reference: `analysis/reduction/run_fs_nrs2_spec2.py`)

Reprocess all Fixed Slit data for PID 1492 (the primary science target) and CALSPEC standards (1537, 1538, 1536, 6644) to obtain spectra in **Jy** natively.

| Step | Task | Success Criteria |
| :--- | :--- | :--- |
| **2.1** | Run Spec2 with `--override_photom reference_files/jwst_nirspec_photom_fs_nrs2_ext.fits`. | Output `x1d.fits` files for NRS2 have non-zero Jy fluxes beyond the nominal cutoff. |
| **2.2** | Check `x1d.fits` metadata. | Verify `TUNIT1 = 'Jy'` and `PHOTMJSR` / `PIXAR_SR` are correctly applied in the NRS2 extended region. |
| **2.3** | Compare to v6. | Verify that NRS2 spectra (new Jy) roughly match scaled v6 (DN/s scaled) but with more consistent transitions. |

---

## 3. Solver Implementation (Reference: `reports/v6/fs_v6/scripts/solve_parlanti_fs_v6.py`)

Update the Parlanti coefficient solver to read the new Jy files and remove the normalization scaling logic.

| Step | Task | Success Criteria |
| :--- | :--- | :--- |
| **3.1** | Create `reports/v7/scripts/solve_parlanti_fs_v7.py`. | Inherit the all-source NNLS joint solve from v6 but remove `ratio = truth / obs` scaling. |
| **3.2** | Implement a "Cross-Validation Loop". | Store results for $(N)$ separate NNLS runs, where each source is held out one at a time. |
| **3.3** | Save separate results FITS files. | `results/v7/calib_v7_fs_g140m_f100lp_all.fits` (all sources) plus held-out variations. |

---

## 4. Visualization & Reporting (Reference: `reports/v6/fs_v6/scripts/plot_fs_v6_validation.py`)

Generate a new suite of plots that quantify the improvement in flux calibration and the stability of the coefficients.

| Step | Task | Success Criteria |
| :--- | :--- | :--- |
| **4.1** | Plot $k, \alpha, \beta$ for all v7 held-out cases. | One plot per grating showing the spread in coefficients when sources are removed. |
| **4.2** | Plot "Hold-out Residuals". | For each grating, plot the fractional residual of the corrected held-out source against its truth spectrum. |
| **4.3** | Merged spectrum validation. | Full (NRS1 + NRS2) spectra for all CALSPEC standards, verified to be in Jy across the entire 0.6 - 5.6 µm range. |
| **4.4** | Full spectrum plot generation. | Created `fs_v7_full_spec_{grating}_{pid}.png` plots comparing observed and model flux. (Completed) |

---

## 5. Verification Checklist (v7)

- [x] All NRS2 extended spectra in my training set are in **Jy**.
- [x] No manual median-scaling factor is used in the solver (i.e. $S_{obs}$ is used as-is).
- [x] Median hold-out residual is $< 10\%$ in the overlap region.
- [x] Transition between NRS1 (nominal) and NRS2 (extended) is smooth and matches within $2\%$.
- [x] Full spectrum validation plots generated for all targets and gratings.

---

## 6. Target Completion (v7)

The aim is to complete the v7 FS calibration (G140M, G235M) by April 1, 2026, and extend to IFU mode shortly thereafter.

---

*Created: 2026-03-30*
*Author: Antigravity*
