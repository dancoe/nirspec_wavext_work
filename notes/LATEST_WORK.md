# Latest Work: NIRSpec Wavelength Extension

## 2026-03-28 Update: V3 Calibration (Kappa Fixed)

We have successfully completed the **V3 Calibration** following the **Parlanti et al. (2025)** methodology. This version enforces a physical 1st-order throughput ($k \approx 1.0$) while solving for 2nd and 3rd order ghost contamination ($\tilde\alpha, \tilde\beta$) using a multi-source NNLS bootstrap solver.

### Key Achievements
- **Alignment**: The recalibrated extended spectra (NRS2) now precisely match nominal higher-dispersion extractions (ground truth) across overlapping regions.
- **Coefficient Analysis**: Tested alternative parameterizations (constant vs. polynomial) for ghost coefficients. Confirmed that point-wise bootstrap solvers are superior for capturing non-polynomial grating structures.
- **Spectral Decomposition**: Implemented high-visibility diagnostic plots ([CAL_COMPONENTS_P330E.png](../plots/Parlanti/cal/153678_v3/CAL_COMPONENTS_P330E.png)) that show raw, scaled, and ghost-corrected components. These plots revealed a significant unit mismatch ($~10^3$ factor) between the uncalibrated raw extension (ADU/s) and the flux-calibrated baseline (Jy).
- **Stability**: Derived coefficients are physically consistent, staying positive and avoiding non-physical corrections.

## Plan Going Forward

1.  **Complex Source Validation**: Apply the V3 solution to `IRAS-05248` (PID 1492) to verify accuracy in the presence of strong absorption features.
2.  **ASDF Integration**: Export the derived coefficients into the custom `wavelengthrange_extended.asdf` for use in the standard pipeline workflow.
3.  **F-Flat Normalization**: Complete the flux calibration pipeline by merging `FAST_VARIATION` tables for concatenated flat-field correction.
4.  **MOS Visualization**: Fix pixel-to-wavelength mapping in MOS diagnostic plots by ensuring consistent detector-specific WCS model usage.
---
For more details on the long-term project plan, see [IMPLEMENTATION_PLAN.md](IMPLEMENTATION_PLAN.md). To see general project instructions and setup, see [INSTRUCTIONS.md](INSTRUCTIONS.md).
