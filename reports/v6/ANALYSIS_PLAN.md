# Analysis Plan v6 — All-Source Joint Solve

**Date:** March 30, 2026  
**Project:** NIRSpec Wavelength Extension Calibration  
**Objective:** Incorporate all available science-target data as calibration constraints alongside CALSPEC standards, eliminating the arbitrary distinction between "calibration" and "science" observations.

---

## 1. Motivation

The v5 analysis demonstrated that the 4-source NNLS for FS (WD + G2V + A8III + G1V) successfully resolved the k/α degeneracy. However, the IFU v5 corrected spectra showed poor agreement with the observed data in the extended NRS2 region (Section 8 of `ifu_v5/REPORT_ifu_v5.md`). Key problems:

- Corrected spectra did not match observed data in extended NRS2 region
- IFU mode lacks a G1V cool star — no IFU PID 6644 analog
- 3-star degenerate system (WD + A + G2V) leaves k/α partially unconstrained

**v6 solution:** Any target with a *known* or *derivable* reference spectrum at the extended wavelengths can serve as a calibration constraint — regardless of whether it is a CALSPEC standard or a science target. For a science target, the reference spectrum is derived from a cross-grating nominal pipeline product.

---

## 2. Updated Dataset (v6)

### IFU Mode

| Source | PID | Spectral Type | G140M constraint | G235M constraint |
|:-------|:----|:-------------|:----------------|:----------------|
| P330E | 1538 | G2V | CALSPEC truth | CALSPEC truth |
| G191-B2B | 1537 | WD-DA.8 | CALSPEC truth | CALSPEC truth |
| J1743045 | 1536 | A8III | CALSPEC truth | CALSPEC truth |
| UGC-5101 | 2186 | ULIRG (AGN+SB) | — | **Cross-grating: G395M nominal** |

UGC-5101 is added only for G235M, because the G395M nominal spectrum overlaps the G235M NRS2 extended region (3.15–5.27 µm).

### Fixed Slit (FS) Mode

| Source | PID | Spectral Type | G140M constraint | G235M constraint |
|:-------|:----|:-------------|:----------------|:----------------|
| G191-B2B | 1537 | WD-DA.8 | CALSPEC truth | CALSPEC truth |
| P330E | 1538 | G2V | CALSPEC truth | CALSPEC truth |
| J1743045 | 1536 | A8III | CALSPEC truth | CALSPEC truth |
| NGC2506-G31 | 6644 | G1V | CALSPEC truth | CALSPEC truth |
| PID 1492 | 1492 | Unknown/science | **Cross-grating: G235M MAST** | **Cross-grating: G395M MAST** |

No new data downloads are required. All sources use data already on disk.

---

## 3. Cross-Grating Truth Method

For a science target with no CALSPEC model:

1. Load the **NRS2 extended** spectrum (our pipeline, DN/s or Jy units).
2. Load the **MAST nominal** spectrum for the adjacent grating (Jy units).
3. Find the wavelength overlap between the two.
4. Compute `scale = median(MAST_Jy / NRS2_obs)` over the overlap region.
5. Apply scale: `NRS2_scaled = NRS2_obs × scale` → now in Jy.
6. Stitch: use `NRS2_scaled` as the observed spectrum and `MAST_Jy` (stitched with nominal NRS1) as the truth interpolator.

This converts any science target into a self-consistent observed/truth pair for the NNLS system, without requiring a stellar atmosphere model.

---

## 4. Implementation

### IFU v6 Solver
- Script: `reports/v6/ifu_v6/scripts/solve_parlanti_ifu_v6.py`
- Output: `results/v6/calib_v6_ifu_g140m_f100lp.fits`, `calib_v6_ifu_g235m_f170lp.fits`
- Run command:
  ```bash
  micromamba run -n jwst_1.20.2 python reports/v6/ifu_v6/scripts/solve_parlanti_ifu_v6.py
  ```

### FS v6 Solver
- Script: `reports/v6/fs_v6/scripts/solve_parlanti_fs_v6.py`
- Output: `results/v6/calib_v6_fs_g140m_f100lp.fits`, `calib_v6_fs_g235m_f170lp.fits`
- Run command:
  ```bash
  micromamba run -n jwst_1.20.2 python reports/v6/fs_v6/scripts/solve_parlanti_fs_v6.py
  ```

### IFU v6 Validation Plots
- Script: `reports/v6/ifu_v6/scripts/plot_ifu_v6_validation.py`
- Plots: `reports/v6/ifu_v6/plots/`

### FS v6 Validation Plots
- Script: `reports/v6/fs_v6/scripts/plot_fs_v6_validation.py`
- Plots: `reports/v6/fs_v6/plots/`

---

## 5. Expected Results

| Mode | Grating | Expected k | Basis |
|:-----|:--------|:----------|:------|
| IFU v6 | G140M | ~0.5 (unchanged from v4/v5) | No new constraint for G140M IFU |
| IFU v6 | G235M | ~0.6 (higher than v5 ~0.47) | UGC-5101 4th source constrains G235M |
| FS v6 | G140M | ~0.7–0.9 | PID 1492 adds 5th G140M source |
| FS v6 | G235M | ~0.9 (near v5 0.96) | PID 1492 adds 5th G235M source |

---

## 6. Validation Objectives

1. **Full spectrum flush:** Do corrected NRS1+NRS2 merged spectra match CALSPEC truth for all standards?
2. **Cross-grating consistency:** Does the corrected PID 1492 or UGC-5101 spectrum match the adjacent-grating nominal?
3. **Coefficient comparison:** Do v6 coefficients improve on v5 in the diagnostic sense (lower residuals, more physical α)?
4. **Condition number:** Is the NNLS system better conditioned (lower median condition number) with the added source?

---

## 7. Successor Plans

- **v7:** Hold out one science target from the NNLS training set; use it for true cross-validation. Design a systematic grid of held-out sources.
- **v8 (future):** Incorporate PRISM or H-grating data for additional wavelength coverage and constraint at λ > 5.0 µm.

---

*Created: 2026-03-30*
