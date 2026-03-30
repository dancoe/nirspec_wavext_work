# Analysis Plan v5 — The Degeneracy Breaker

**Date:** March 29, 2026  
**Project:** NIRSpec Wavelength Extension Calibration  
**Objective:** Incorporate **NGC2506-G31** (G1V; PID 6644) into the simultaneous $\kappa, \alpha, \beta$ derivation to break the hot-star degeneracy and stabilize the instrument calibration.

---

## 1. Goal
As identified in v4, solving for $\kappa, \alpha, \text{and } \beta$ with only hot stars (WDs, A-stars) leads to a solution where $\alpha$ is over-estimated (2–6× Parlanti values) and $\kappa$ is under-estimated. 

The primary goal of v5 is to add a **cool star (G1V)**. Since these stars have significant molecular absorption features at $\lambda/2$ that map onto $\lambda$ in the second order (but not the first), the solver can uniquely identify the second-order component.

---

## 2. Updated Dataset (v5)

### Fixed Slit (FS)
*   **P330-E** (G2V; Solar-analog baseline)
*   **G191-B2B** (WD; Hot continuum baseline)
*   **J1743045** (A8III; Warm continuum baseline)
*   **NGC2506-G31 (G1V; PID 6644) — CRITICAL ADDITION**

### IFU (Conditional)
*   If PID 6645 (P330E centered) and archival IFU observations of cool stars exist, they will be prioritized here as well. Note: PID 6644 appears focused on FS.
*   **PID 2654**: SDSSJ0749, SDSSJ0841 (AGN; z~2 for G140M validation).
*   **PID 2186**: UGC-5101 and 2 other local ULIRGs (for G235M/F170LP validation).

---

## 3. Validation Objectives (v5-IFU)
Beyond the standard stars, v5 will utilize these science targets to duplicate the Parlanti et al. validation cases:
1.  **G140M/F100LP:** Recover H-alpha for the z~2.6 QSOs in PID 2654 and confirm it matches the nominal expectations.
2.  **G235M/F170LP:** Verify the extended continuum for the ULIRGs in PID 2186 against the nominal G395M overlaps.

---

## 4. Implementation Steps

1.  **Acquisition:**
    *   Download PID 6644 FS rate files (S1600A1).
    *   Download PID 2654 and 2186 IFU rate files for the validation science targets.
2.  **Processing:**
    *   Run `calwebb_spec2` with standard detector-1 limits.
    *   Run `calwebb_spec2` with the wavelength-extended NRS2 model.
    *   Extract spectra for all 4 gratings (G140M, G235M, G395M) and PRISM.
3.  **Solver Execution:**
    *   Initialize the NNLS solver with all 4 FS sources.
    *   Compare the resulting $\kappa(\lambda)$ to the v3 (Parlanti-fixed) and v4 (Hot-star only) results.
4.  **Data Conditioning and Refinement:**
    *   **Data Gaps:** Address areas where coefficients drop to zero (e.g., in v4). Zero-filling creates artificial features when smoothed.
    *   **Interpolation vs. Gaps:** Instead of setting missing data to zero, v5 will evaluate interpolation or maintaining explicit gaps (as per Parlanti 2025). 
    *   **Smoothing Strategy:** Refine the smoothing kernel to respect these gaps to prevent edge artifacts.

---

## 5. Expected Results
*   A "re-inflation" of $\kappa(\lambda)$ back toward higher values consistent with standard sensitivity drops.
*   The $\alpha$ coefficient settling closer to the Parlanti et al. (2025) values (likely < 0.1).
*   Corrected spectra for all 4 temperature classes (WD, A, G-solar, G-cool) showing consistency with ground-truth nominal gratings.

---
*Created by Antigravity; v5 Roadmap pending execution.*
