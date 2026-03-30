# NIRSpec Wavelength Extension Report — 330 FS v5

**Date:** March 29, 2026  
**Project:** NIRSpec Wavelength Extension Calibration  
**Data Version:** FS v5 — 4-Source joint-solve (The Degeneracy Breaker)

---

## 1. Summary

The FS v5 calibration integrates a G1V cool star (NGC2506-G31; PID 6644) as a fourth source alongside the hot white dwarf (G191-B2B), G2V solar-like (P330E), and A8III star (J1743045). The goal was to resolve the $k/\alpha$ degeneracy encountered in v4, where simultaneous NNLS solutions for $k$, $\alpha$, and $\beta$ over-estimated the second-order contamination and under-estimated the first-order throughput.

**Key results:**
- **Degeneracy Resolved:** The inclusion of the G1V spectrum, which has a distinct SED and spectral features, provides the necessary constraints to decouple throughput from ghosting.
- **Throughput Recovery:** $k(\lambda)$ has returned to $\sim$0.96, consistent with the expected baseline and moving away from the degenerate low $k$ (0.7) found in v4.
- **Physical Contamination:** The derived second-order contamination factor $\alpha$ is substantially lower and more physically plausible ($0.003–0.012$) than the v4 estimates ($0.03–0.04$).
- **Spectral Validation:** Residuals for all four standards are optimized, with the core regions ($1.0–3.0$ µm) matching ground truth to within $\sim$2–5%.

---

## 2. Coefficient Comparison (v5 vs. v4 vs. Parlanti)

The plots below compare the v5 derived coefficients against the degenerate v4 results and the Parlanti et al. (2025) reference values.

### G140M/F100LP
![G140M Coefficients](plots/coeffs_g140m_f100lp.png)

### G235M/F170LP
![G235M Coefficients](plots/coeffs_g235m_f170lp.png)

### FS v5 Coeffs (2-Panel Log View)
The plot below provides a side-by-side log-scale view of the G140M and G235M coefficients, emphasizing the separation of throughput (black) from the smaller contamination terms (orange, blue).

![v5 Log Coeffs](plots/fs_v5_coeffs_log_2panel.png)

---

## 3. Residuals Summary

Spectral validation panels for all four standard stars confirm that the v5 joint-solve provides a consistent fit across diverse spectral types (WD, G2V, A8III, and G1V).

### G140M/F100LP Residuals
![G140M Residuals](plots/residuals_g140m_f100lp.png)

### G235M/F170LP Residuals
![G235M Residuals](plots/residuals_g235m_f170lp.png)

---

## 4. Derived Median Coefficients

| Grating | k median | alpha median | beta median |
|:--------|:---------|:-------------|:------------|
| G140M FS | **0.959** | **0.0028** | **0.0003** |
| G235M FS | **0.957** | **0.0124** | **0.0011** |

*(v4 comparison: G140M k=0.701, G235M k=0.723; the v5 solve recovers ~25% more throughput.)*

---

## 5. Discussion

The v5 results demonstrate that the NIRSpec wavelength extension can be robustly calibrated using a small but diverse set of standard stars. The previous finding of very high $\alpha$ in v4 was definitively an artifact of the $k/\alpha$ degeneracy in a 3-source system dominated by hot stars. The G1V star, with its strong absorption features and cooler SED, is instrumental in breaking this degeneracy.

The high $\alpha$ values at the NRS1/NRS2 boundary (found in v4) have been resolved; the v5 solution is much more stable at the red-end limits of the gratings.

---

## 6. Next Steps

1. **IFU Integration:** Run the 4-source joint-solve for IFU once the PID 6644 IFU data is processed.
2. **Science Validation (z-agn):** Apply these v5 coefficients to the SDSS QSOs (PID 2654) and ULIRGs (PID 2186) to verify ghost recovery of H-alpha features on the extended range.
3. **Smooth and Release:** Perform final spline smoothing of the derived coefficients for implementation in the main project pipeline.

---

## 7. Full Spectrum Baseline Verification

The following plots show the final merged NIRSpec spectrum (NRS1 + Cor-NRS2) compared to the CALSPEC truth model. These confirm the recovery of the first-order SED after accurate ghost subtraction.

### P330E (G2V)
![P330E Full Spectrum](plots/full_spectrum_merged_v5_P330E.png)

### G191-B2B (WD)
![G191-B2B Full Spectrum](plots/full_spectrum_merged_v5_G191-B2B.png)

### J1743045 (A8III)
![J1743045 Full Spectrum](plots/full_spectrum_merged_v5_J1743045.png)

### NGC2506-G31 (G1V)
![NGC2506-G31 Full Spectrum](plots/full_spectrum_merged_v5_NGC2506-G31.png)

---
*Maintained by Antigravity*
