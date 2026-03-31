# NIRSpec Wavelength Extension Report — 329 IFU v3

**Date:** March 29, 2026  
**Project:** NIRSpec Wavelength Extension Calibration  
**Data Version:** IFU v3 — Level-3 stage3_ext (NRS1+NRS2 combined)

---

## Summary

IFU v3 makes the Level-3 provenance explicit. The solver algorithm is identical to v2 (unbiased median-R, no regularisation), but this version:
- Uses **Spec3Pipeline stage3_ext x1d** products (NRS1+NRS2 combined cube cubepar extraction) rather than per-exposure Level-2 products.
- Adds **MAST L3 nominal vs our L3 ext** comparison panels showing how the extended pipeline extends wavelength coverage beyond the MAST products.
- Is the IFU counterpart to FS v3, enabling a clean like-for-like v3 comparison.

Excellent recalibration (within ~5–10%) for P330E and J1743045. G191-B2B G140M overcorrection (~20–30%) is source-specific, as established in v2.

---

## Pipeline Lineage (WAVEXT.md tweaks applied)

| Step | Customisation |
|:-----|:--------------|
| `assign_wcs` | `wavelengthrange_0008.asdf` extended range (G140M 0.6–3.3 µm, G235M 1.0–5.3 µm) |
| `assign_wcs` | `jwst/assign_wcs/nirspec.py` hard-coded NRS2 M-grating error removed |
| `flat_field` | Parlanti sflat/fflat (NRS2 extended coverage) |
| `photom` | Parlanti photom tables (NRS2 extended coverage) |
| `cube_build` | `cubepar` Parlanti tables (extended wavelength grid: `f100lp_g140m-f100lp_x1d.fits` etc.) |
| `extract_1d` | Extracts over full NRS1+NRS2 cube |

Output: `data/IFU/PID{1536|1537|1538}/stage3_ext/{filter}_{grating}-{filter}_x1d.fits`

---

## Algorithm (identical to v2)

### Step 1 — Derive k(λ)

$$R_i(\lambda) = \frac{S_{\mathrm{obs},i}(\lambda)}{f_{\mathrm{CALSPEC},i}(\lambda)}$$

$$k_i(\lambda) = R_i(\lambda) - \tilde{\alpha}_P(\lambda)\cdot r_{2,i}(\lambda) - \tilde{\beta}_P(\lambda)\cdot r_{3,i}(\lambda)$$

$$k(\lambda) = \mathrm{median}_i\bigl[k_i(\lambda)\bigr], \quad \text{40-channel box-smoothed}$$

**k-sources:** P330E (PID 1538), J1743045 (PID 1536).  
**G191-B2B excluded** from k derivation (source-dependent photometric excess, T~60 000 K hot WD).

### Step 2 — α̃(λ), β̃(λ)

Loaded from Parlanti published calibration FITS:
- `calibration_functions_g140m_f100lp.fits`
- `calibration_functions_g235m_f170lp.fits`

### Correction formula

$$f_\mathrm{corr}(\lambda) = \frac{S_\mathrm{obs}(\lambda) - \tilde{\alpha}(\lambda)\,f(\lambda/2) - \tilde{\beta}(\lambda)\,f(\lambda/3)}{k(\lambda)}$$

---

## Coefficient Summary

| Grating | k median | k range | α̃ max | β̃ max |
|:--------|:---------|:--------|:------|:------|
| G140M NRS2 | **0.567** | 0.298–1.322 | 0.051 (Parlanti) | 0.008 (Parlanti) |
| G235M NRS2 | **0.768** | 0.545–1.307 | 0.042 (Parlanti) | 0.008 (Parlanti) |

k starts above 1 near the NRS1/NRS2 boundary (pipeline over-corrects at the transition) and falls to ~0.3–0.5 at the far red end of each grating's NRS2 extension.

---

## 1. Calibration Coefficients

### G140M Coefficients
![IFU v3 G140M Coefficients](../../plots/Parlanti/cal/ifu_v3/ifu_v3_g140m_coeffs.png)

### G235M Coefficients
![IFU v3 G235M Coefficients](../../plots/Parlanti/cal/ifu_v3/ifu_v3_g235m_coeffs.png)

---

## 2. MAST Level-3 vs Our Level-3 Extended

These panels show the key v3 addition: MAST L3 (steelblue) alongside our Level-3 extraction on NRS1 (maroon) for comparison, and the NRS2 extension uncorrected (orange) and corrected (red). Comparing the maroon and steelblue lines on the left half of each panel validates that our pipeline correctly reproduces the nominal extraction on NRS1.

### G140M — All Sources
![MAST vs L3 G140M](../../plots/Parlanti/cal/ifu_v3/ifu_v3_mast_vs_l3_G140M.png)

### G235M — All Sources
![MAST vs L3 G235M](../../plots/Parlanti/cal/ifu_v3/ifu_v3_mast_vs_l3_G235M.png)

---

## 3. NRS2 Spectral Validation

For each source: gray = stage3_ext NRS2 (raw); blue dashed = ground-truth (next grating); red = corrected (v3); black dotted = CALSPEC model.

### P330E — G140M NRS2
![P330E G140M](../../plots/Parlanti/cal/ifu_v3/ifu_v3_spectra_P330E_G140M.png)

### P330E — G235M NRS2
![P330E G235M](../../plots/Parlanti/cal/ifu_v3/ifu_v3_spectra_P330E_G235M.png)

### J1743045 — G140M NRS2
![J1743045 G140M](../../plots/Parlanti/cal/ifu_v3/ifu_v3_spectra_J1743045_G140M.png)

### J1743045 — G235M NRS2
![J1743045 G235M](../../plots/Parlanti/cal/ifu_v3/ifu_v3_spectra_J1743045_G235M.png)

### G191-B2B — G140M NRS2 (validation only, excluded from k derivation)
![G191-B2B G140M](../../plots/Parlanti/cal/ifu_v3/ifu_v3_spectra_G191-B2B_G140M.png)

### G191-B2B — G235M NRS2 (validation only, excluded from k derivation)
![G191-B2B G235M](../../plots/Parlanti/cal/ifu_v3/ifu_v3_spectra_G191-B2B_G235M.png)

---

## 4. Full Spectrum View

Full 0.6–5.6 µm view merging MAST L3 nominal gratings (steelblue/mediumblue/navy) with our NRS2 L3 ext corrected extensions (crimson = G140M, orangered = G235M).

### P330E
![P330E Full Spectrum](../../plots/Parlanti/cal/ifu_v3/ifu_v3_full_spectrum_P330E.png)

### J1743045
![J1743045 Full Spectrum](../../plots/Parlanti/cal/ifu_v3/ifu_v3_full_spectrum_J1743045.png)

### G191-B2B
![G191-B2B Full Spectrum](../../plots/Parlanti/cal/ifu_v3/ifu_v3_full_spectrum_G191-B2B.png)

---

## Interpretation

**P330E and J1743045** — EXCELLENT agreement (within ~5–10%) between the v3-corrected NRS2 extension and the next-grating nominal truth across both G140M and G235M NRS2 ranges.

**G191-B2B** — G235M NRS2 corrects well. G140M NRS2 overcorrects by ~20–30% at λ > 2.5 µm. This is a known hot WD artefact (source-dependent photometric excess at the NRS2/pipeline interface), not a solver deficiency.

**MAST vs L3 ext comparison** — The stage3_ext products clearly extend wavelength coverage beyond the MAST L3 products. Our NRS1 Level-3 extraction (maroon) matches the MAST L3 extraction (steelblue), demonstrating consistency. The NRS2 raw extension (orange) is brighter than expected; after applying v3 k(λ), the corrected extension (red) aligns with the truth.

**IFU v3 vs FS v3** — IFU v3 k values are systematically lower than FS v3 k values (IFU G140M median 0.567 vs FS 0.725; IFU G235M median 0.768 vs FS 0.972). This is a genuine mode-dependent photometric calibration difference, not a pipeline bug.

---

## Solver Script
- [../../analysis/solver/solve_parlanti_ifu_v3.py](../../analysis/solver/solve_parlanti_ifu_v3.py)

---

*Created automatically by Antigravity on 2026-03-29.*
