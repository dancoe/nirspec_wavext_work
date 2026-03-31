# Detailed Analysis Plan v8 — IFU Comprehensive Calibration

**Version:** 8.0  
**Status:** Planning  
**Objective:** Broaden the IFU calibration suite by including all available sources (standards and AGN) and refining gap handling.

---

## 1. Data Inventory & Preparation (IFU v7 Stage 3)
Ensure all IFU Stage 3 `spec3` products are in hand, unified by the v7 natively calibrated (Jy) pipeline.

**Targets to include:**
- **Standard Stars:** 1536 (J1743045), 1537 (G191-B2B), 1538 (P330E), 6645 (P330E Cycle 3).
- **AGN/ULIRGs (Validation):** 2186 (UGC-5101), 2654 (SDSSJ0749).

---

## 2. Solver Refinement: Gap Awareness (Reference: v8 FS)
- **Identify Gaps:** Systematic masking for gaps at ~2.2 µm (G140M) and ~3.7 µm (G235M).
- **No Interpolation:** Re-run the IFU NNLS solver ensuring we do not interpolate smooth coefficients across the regions of missing data.

---

## 3. Visualization: Full Spectrum Baseline (v5 Style)
- Create new plots for all 6 sources showing the final NRS1,2 spectra compared to CALSPEC (for stars) or self-calibration truth (for AGN).
- Ensure "NRS1,2" labeling is used across all gratings.
- Gap-Safe plotting: Do not draw lines across the ~2.2/3.7 µm gaps.

---

## 4. Verification Checklist (v8)
- [ ] Gaps are visible in the final plots and not smoothed over.
- [ ] All 3 Gordon standards + 2 AGN targets + Zeidler P330E are successfully analyzed.
- [ ] Median hold-out residuals for standard stars are documented.
- [ ] Final IFU $k, \alpha, \beta$ coefficients saved in `results/v8/`.

---

*Created: 2026-03-31*  
*Author: Antigravity*
