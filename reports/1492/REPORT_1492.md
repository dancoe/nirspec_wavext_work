# PID 1492 (IRAS-05248-7007) — v7 Calibration Report

**Date:** March 30, 2026  
**Project:** NIRSpec Wavelength Extension (v7)  
**Target:** IRAS-05248-7007 (LIRG, Science Target)

---

## 1. Overview

Program ID 1492 provides high-SNR NIRSpec Fixed Slit (FS) observations across all Medium-Resolution gratings (G140M, G235M, G395M) and the PRISM. Because it covers overlapping wavelength regions, it serves as a critical "stitching" target for cross-validating the Parlanti spectral recovery model where CALSPEC standards are unavailable or overlap is limited.

## 2. Data Processing Workflow (v7)

The v7 processing of PID 1492 data differs from previous versions (v5/v6) by adopting native **Jy-unit extractions** at the `Spec2Pipeline` stage.

### Stage 2 Reprocessing
- **Reference Files:**
    - `wavelengthrange`: `jwst_nirspec_wavelengthrange_0008.asdf` (Extends WCS into NRS2 gaps).
    - `photom`: `jwst_nirspec_photom_fs_nrs2_ext.fits` (Custom sensitivity curves for NRS2 wavelengths).
    - `fflat`/`sflat`: Parlanti extended flats (covering 0.6–6.0 µm).
- **Execution:** Run via `analysis/reduction/run_fs_nrs2_spec2.py`.
- **Note on units:** The `photom` step is applied during `Spec2Pipeline`, producing `cal.fits` files in Jy.

### Stage 3 Combination
- **Execution:** `analysis/reduction/run_fs_nrs2_spec3.py` combines multiple `cal.fits` exposures into Level-3 products.
- **PHOTOM in Spec3?**: While the forked `Spec3Pipeline` *can* call `photom` (to avoid Stage 2 errors), for the v7 FS workflow, **photom is run during Stage 2**. Stage 3 receives Jy-unit inputs, performs `resample_spec` (optional), and executes `extract_1d` to produce the final `x1d.fits` files.
- **Critical Logic:** In v7, we verify that `S_PHOTOM` is `COMPLETE` in the input `cal` files. This eliminates the need for the manual DN/s-to-Jy median scaling used in v6.

## 3. Validation Results

For PID 1492, the "Truth" spectrum is defined by the flux-calibrated overlap with an adjacent grating (e.g., G140M "truth" is derived from the G235M nominal region).

- **NRS2 Extended Recovery:** v7 successfully restores flux beyond the nominal NIRSpec detector cutoffs (~3.1 µm for G140M, ~5.2 µm for G235M).
- **Transition Continuity:** Native Jy calibration ensures a smoother transition between detectors (NRS1 vs NRS2) compared to manual v6 scaling.

## 4. Diagnostic Plots

The following plots verify the accuracy of the spectral recovery:

- **1492 NRS2 Validation (G140M):** [1492_nrs2_validation_g140m.png](plots/1492_nrs2_validation_g140m.png) — *NRS2 observed flux compared to v7 ghost-correction model.*
- **1492 NRS2 Validation (G235M):** [1492_nrs2_validation_g235m.png](plots/1492_nrs2_validation_g235m.png)

---

## 5. Next Steps
- Implement v7 coefficients in the automated `solve_parlanti_fs_v7.py` script.
- Extend the native Jy workflow to IFU targets (PIDs 2654, 2186).
- Update the collective [DATA.md](../notes/DATA.md) and [REPORT_fs_v7.md](../v7/REPORT_fs_v7.md).

*Author: Antigravity*
