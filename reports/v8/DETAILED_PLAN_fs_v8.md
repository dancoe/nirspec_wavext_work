# Detailed Analysis Plan fs_v8 — Final Validation & Gap-Aware Reporting

**Version:** 8.0  
**Status:** Planning Phase  
**Objective:** Finalize the Fixed Slit (FS) NIRSpec wavelength extension calibration with gap-aware validation plots and full-spectrum (0.6 - 5.6 µm) consistency checks for all primary targets.

---

## 1. Scope & Data
This analysis assumes the **v7 natively flux-calibrated (Jy) extractions** are complete for all FS targets.

| PID  | Target | Mode | Status |
|:-----|:-------|:-----|:-------|
| 1536 | J1743045   | FS S1600A1 | Ready |
| 1537 | G191-B2B   | FS S1600A1 | Ready |
| 1538 | P330E      | FS S1600A1 | Ready |
| 6644 | NGC2506-G31| FS S1600A1 | Ready |
| 1492 | IRAS-05248 | FS (various slits) | Ready |

---

## 2. Gap-Aware Validation Plots (Reference: `fs_v5` style)

A key requirement for v8 is to handle spectral gaps correctly in visualization. We must ensure that the plotting routines do **not** interpolate or smooth over the physical gaps between detectors or filters.

| Grating | Nominal Gap (NRS1/NRS2) |
|:--------|:-------------------------|
| **G140M** | ~2.17 – 2.28 µm |
| **G235M** | ~3.66 – 3.82 µm |

### Task 2.1: Update Plotting Script
Modified version of `reports/v7/scripts/plot_v7_full_spectra.py` to:
- Detect wavelength jumps > 0.05 µm (or specific gap ranges).
- Use `ax.plot(..., linestyle='-', ...)` separately for each segment (NRS1, NRS2).
- **Rule:** Do not allow `interp1d` or `plt.plot` to connect the last pixel of NRS1 to the first pixel of NRS2 if they are separated by a physical gap.

---

## 3. Full Spectrum Generation (0.6 - 5.6 µm)

Generate a single, premium-quality validation plot per source that merges all available gratings (G140M, G235M, G395M, PRISM) into a unified spectral timeline.

| Step | Task | Success Criteria |
|:-----|:-----|:-----------------|
| **3.1** | Run `reports/v8/scripts/plot_fs_v8_full_spectra.py`. | Created `reports/v8/fs_v8/plots/full_spectrum_merged_v8_{TARGET}.png`. |
| **3.2** | Verify gap preservation. | Visual check: G140M spectrum breaks clearly at ~2.2 µm. |
| **3.3** | Merged model comparison. | Plot corrected NRS2 spectra alongside v7 coefficients ($k f_1 + \alpha f_2 + \beta f_3$) to confirm residuals are minimized. |

---

## 4. Verification Checklist (fs_v8)

- [ ] All plots use HSL-harmonized colors for gratings (similar to v7 style).
- [ ] Log-scale y-axis used for full spectrum plots (0.6 - 5.6 µm).
- [ ] No interpolation lines visible across the 2.2 µm (G140M) and 3.7 µm (G235M) gaps.
- [ ] Coefficients from v7 applied globally to PID 1492 (IRAS-05248) as the final science check.

---

*Created: 2026-03-31*  
*Author: Antigravity*
