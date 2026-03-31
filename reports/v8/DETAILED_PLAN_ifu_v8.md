# Detailed Analysis Plan ifu_v8 — Science Target Extraction & Final Validation

**Version:** 8.0  
**Status:** Planning Phase  
**Objective:** Execute custom-aperture extractions for IFU science targets and generate the final v8 validation suite featuring full-spectrum (0.6 - 5.6 µm) plots with rigorous gap preservation.

---

## 1. Scope & Data
This phase analyzes all Integral Field Unit (IFU) targets identified in `OBSERVATIONS.md`, assuming the **v7 natively flux-calibrated (Jy)** cubes are complete.

| PID  | Target | Mode | Notes |
|:-----|:-------|:-----|:------|
| 1536 | J1743045   | IFU | Standard Star |
| 1537 | G191-B2B   | IFU | Standard Star |
| 1538 | P330E      | IFU | Standard Star |
| 6645 | P330E      | IFU | Cycle 3 Repeat |
| 2186 | UGC-5101   | IFU | ULIRG Science Target |
| 2654 | SDSSJ0841  | IFU | AGN Science Target |

---

## 2. Parlanti-Style Science Extractions (Custom Aperture)

To ensure consistency with the Parlanti et al. (2025) methodology for science validation, the extractions for PIDs 2186 and 2654 will use a fixed circular aperture.

| Source | Target | Extraction Strategy |
|:-------|:-------|:---------------------|
| **AGN SDSSJ0841** | PID 2654 | **0.5" radius circular aperture** centered on the brightest pixel. |
| **ULRIG UGC-5101** | PID 2186 | **0.5" radius circular aperture** centered on the brightest pixel. |

### Task 2.1: Custom Extraction Script
Create `reports/v8/scripts/extract_ifu_science_v8.py`:
- Load `_s3d.fits` (IFU cubes).
- Locate the peak pixel in the NRS1 nominal range (as a spatial reference).
- Define a 0.5" circular mask (translated to pixels using the WCS or `CDELT` parameters).
- Sum fluxes across the mask for all wavelengths (NRS1 + NRS2).
- Output new `v8_extract_0.5.x1d.fits` files for analysis.

---

## 3. Gap-Aware Validation Plots

Identical to the FS requirement, the IFU plotting routines must leave spectral gaps as-is.

| Grating | Nominal Gap (NRS1/NRS2) |
|:--------|:-------------------------|
| **G140M** | ~2.17 – 2.28 µm |
| **G235M** | ~3.66 – 3.82 µm |

### Task 3.1: Generate Full Spectrum Plots
Run `reports/v8/scripts/plot_ifu_v8_full_spectra.py`:
- Plot the 0.6 – 5.6 µm range for all 6 targets.
- Ensure the **G140M** spectrum at **2.2 µm** and **G235M** spectrum at **3.7 µm** show intentional gaps (no connecting lines).
- Highlight the ghost-subtracted NRS2 regions.

---

## 4. Verification Checklist (ifu_v8)

- [ ] All 6 IFU targets (1536, 1537, 1538, 2186, 2654, 6645) processed.
- [ ] Circular 0.5" extractions verified for 2186 and 2654.
- [ ] Gap preservation confirmed in final `.png` plots (break observed at 2.17 - 2.28 µm and 3.66 - 3.82 µm).
- [ ] Residuals between corrected science targets and nominal-range benchmarks minimized.

---

*Created: 2026-03-31*  
*Author: Antigravity*
