# Latest Work: NIRSpec Wavelength Extension

## 2026-03-28 Update: IFU v1 Coefficient Derivation

See session logs: [2026-03-28-2300](logs/2026-03-28-2300.md) · [2026-03-28-2330](logs/2026-03-28-2330.md)

### IFU v1 Solver Results (`solve_parlanti_ifu_v1.py`)

Using IFU `stage3_ext` products (Jy) and CALSPEC models as the reference:

| Grating | k(λ) median | k(λ) range | α | β |
|:--------|:-----------|:-----------|:--|:--|
| G140M NRS2 | **0.833** | 0.712–1.003 | ≈ 0 | ≈ 0 |
| G235M NRS2 | **0.899** | 0.817–1.001 | ≈ 0 | ≈ 0 |

- k(λ) starts near 1.0 at the NRS2 cutoff and decreases toward longer wavelengths
- No second/third-order contamination detected (α = β = 0)
- **Caveat**: k/α degeneracy cannot be broken with only 3 hot stars of similar spectral shape
  (f(λ/2)/f(λ) is source-independent for power-law continua; need a cool calibrator to separate)
- After applying k correction, ~15–25% residual offset from truth likely reflects G235M/G395M IFU
  nominal calibration accuracy at wavelengths > 2.2 µm, not uncorrected contamination

### Diagnostic (root cause, fixed in previous session)
The `plot_parlanti_components_multi.py` diagnostic was comparing wrong datasets:
- **Old (wrong)**: Per-exposure `nrs2_extract_1d.fits` in **DN/s** vs stage3 FS in **Jy** → 10³× offset
- **New (correct)**: IFU `stage3_ext` x1d in **Jy** vs IFU `stage3` of the *next* grating as truth

### Comparison geometry (both in Jy)
| Extended data | Truth reference | Overlap |
|:---|:---|:---|
| G140M stage3_ext NRS2 | G235M stage3 IFU nominal | 1.87–3.17 µm |
| G235M stage3_ext NRS2 | G395M FS nominal | 3.15–5.14 µm |

### Next Steps
1. **Break k/α degeneracy**: add a cool calibrator or adopt Parlanti+25 published α directly
2. **Validation**: Apply IFU v1 k(λ) to IRAS-05248 (PID 1492)
3. **ASDF Integration**: Export coefficients to `wavelengthrange_extended.asdf`
4. **Investigate G235M NOM calibration**: understand the 15–25% residual at λ > 2.2 µm
5. **F-Flat normalization**: merge `FAST_VARIATION` tables for full flux calibration

---
For more details on the long-term project plan, see [IMPLEMENTATION_PLAN.md](IMPLEMENTATION_PLAN.md). To see general project instructions and setup, see [INSTRUCTIONS.md](INSTRUCTIONS.md).
