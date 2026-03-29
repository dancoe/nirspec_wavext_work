# Latest Work: NIRSpec Wavelength Extension

## 2026-03-28 Update: IFU Calibration Diagnostic — Data-Source Correction

### Problem Identified
The `plot_parlanti_components_multi.py` diagnostic was comparing the wrong datasets:
- **Old (wrong)**: Used per-exposure `nrs2_extract_1d.fits` files (in **DN/s**) from `data/PID*/` versus nominal FS x1d products (in **Jy**). The G140M NOM and EXT had **zero wavelength overlap** (NOM: 0.97–1.87 µm, EXT: 1.96–3.16 µm), so the overlap-based scale factor always returned NaN.
- **New (correct)**: Uses **IFU stage3_ext x1d** products (in **Jy**, from `data/IFU/*/stage3_ext/`) which contain BOTH NRS1 and NRS2 in a single flux-calibrated file. These are compared against the IFU stage3 nominal of the *next* grating as ground truth.

### Key Findings from Revised Diagnostic
The new comparison (`CAL_COMPONENTS_*.png` for all three standards) reveals:
- **G140M stage3_ext NRS2** (1.87–3.60 µm): ~**2x too faint** relative to the G235M stage3 nominal in the overlap region (1.87–3.17 µm). The throughput factor k(λ) rises from 1.0 at 1.87 µm to ~2.0 at 3.0 µm — consistent across P330E, G191-B2B, and J1743045.
- **G235M stage3_ext NRS2** (3.15–5.50 µm): ~**1.5x too faint** relative to the G395M FS nominal in the overlap region. k(λ) is near 1.0–1.7 with wavelength dependence.
- **V3 coefficients** (k≈0.89 from FS DN/s solver) are **not applicable** to the IFU stage3_ext products — they were derived from per-exposure FS NRS2 extractions that needed their own unit conversion and boundary scaling.

### Comparison geometry (both in Jy, same calibration system):
| Extended data | Ground-truth reference | Overlap region |
|:---|:---|:---|
| G140M stage3_ext NRS2 | G235M stage3 nominal (IFU) | 1.87–3.17 µm |
| G235M stage3_ext NRS2 | G395M FS nominal | 3.15–5.14 µm |

### Next Steps
1. **Re-derive k/α/β from IFU stage3_ext data**: Update `solve_parlanti_bootstrap_v3.py` to use IFU stage3_ext products instead of per-exposure FS extractions. The new k(λ) for G140M NRS2 is ~1.5–2.0 (not 0.89 from V3).
2. **Complex Source Validation**: Apply the corrected solution to `IRAS-05248` (PID 1492).
3. **ASDF Integration**: Export the derived coefficients into `wavelengthrange_extended.asdf`.
4. **F-Flat Normalization**: Complete flux calibration by merging `FAST_VARIATION` tables.

---
For more details on the long-term project plan, see [IMPLEMENTATION_PLAN.md](IMPLEMENTATION_PLAN.md). To see general project instructions and setup, see [INSTRUCTIONS.md](INSTRUCTIONS.md).
