# Analysis Plan v4 — Full Coefficient Derivation

**Date:** March 29, 2026  
**Project:** NIRSpec Wavelength Extension Calibration  
**Objective:** Calculate empirical $\kappa(\lambda)$, $\tilde{\alpha}(\lambda)$, and $\tilde{\beta}(\lambda)$ coefficients for both IFU and FS modes, moving beyond the Parlanti-fixed models in v2/v3.

---

## 1. Roadmap of Versions

### v4 — The Three-Star Baseline (Active)
*   **Targets:** P330-E (G2V), J1743045 (A8III), G191-B2B (WD)
*   **Goal:** Establish the mathematical stability of the 3-coefficient solve ($\kappa, \alpha, \beta$) with the primary calibration set.
*   **Current Verdict:** Preliminary v4 results show $\alpha$ values 2–6× higher than Parlanti, suggesting significant degeneracy with $\kappa$ when only using hot stars.

### v5 — The Degeneracy Breaker
*   **Targets:** v4 set + **NGC2506-G31** (G1V; PID 6644)
*   **Goal:** Use the unique molecular features of a cool star to definitively separate throughput ($\kappa$) from ghosting ($\alpha$).
*   **Status:** Awaiting reduction of PID 6644.

### v6 — Full Statistical Refinement
*   **Targets:** v5 set + **J1808347** (PID 1128) + **Ashby A-stars** (PID 1678)
*   **Goal:** Finalize the instrument-level calibration for public release with maximum statistical weight and bias reduction.

---

## 2. Problem Statement
The v2 and v3 analyses used the Parlanti et al. (2025) published values for second and third-order contamination ($\tilde{\alpha}$ and $\tilde{\beta}$) while only solving for the throughput ($\kappa$). For a truly independent calibration, we must solve for all three parameters simultaneously using our own dataset.

### The Degeneracy Constraint
Solving for three unknowns at each wavelength requires at least **three independent sources** with different spectral shapes. Using only hot stars (WDs and A-stars) leads to a $k/\alpha$ degeneracy because their SEDs look similar at these wavelengths.

---

## 3. Dataset Selection Plan

### A. IFU Calibration Set
*   **Target 1 (Solar):** P330-E (G2V; PIDs 1538, 6645, 6606) — *Primary Anchor*
*   **Target 2 (Hot):** G191-B2B (WD; PID 1537)
*   **Target 3 (Warm):** J1743045 (A8III; PID 1536)

### B. Fixed Slit (FS) Calibration Set
*   **Target 1 (Solar):** P330-E (G2V; PID 1538/6606)
*   **Target 2 (Hot):** G191-B2B (WD; PID 1537)
*   **Target 3 (Warm):** J1743045 (A8III; PID 1536)
*   **Target 4 (CRITICAL):** **NGC2506-G31** (G1V; PID 6644) — *Key for separating k from α via absorption line reflection.*

---

## 3. Implementation Workflow

1.  **Data Reduction:**
    *   Reduce NGC2506-G31 (PID 6644) FS rate files for NRS2.
    *   Verify S/N > 5 across the full 1.7–5.3 µm range.
2.  **Solver Implementation:**
    *   Implement a 3-component NNLS solver at each wavelength grid point.
    *   Apply a 40-channel smoothing window (consistent with Parlanti 2025).
    *   Enable cross-validation by holding out one A-star (e.g., J1808347).
3.  **Visualization:**
    *   Compare our derived $\alpha$ and $\beta$ directly to Parlanti's.
    *   Validate by applying the new coefficients to recalibrate nominal-range datasets.

---

## 4. Expected Outcomes
*   Independent derivation of detector ghosting fractions for FS mode.
*   Quantification of mode-dependence in the throughput $\kappa$ between IFU and FS.
*   Validation against standard CALSPEC models with < 5% residuals.

---
*Created by Antigravity based on analysis of OBSERVATIONS.md*
