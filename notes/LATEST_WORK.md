# Latest Work: NIRSpec Wavelength Extension

## 2026-03-29 Update: Fixed Slit v1 Coefficient Derivation

See session log: [2026-03-29 session](logs/)

### FS v1 Solver Results (`solve_parlanti_fs_v1.py`)

Using per-exposure FS NRS2 x1d in Jy (produced by `run_fs_nrs2_spec2.py` with Parlanti flat overrides)
and CALSPEC models as the reference. 3 sources: P330E (G2V), G191-B2B (WD), J1743045 (A8III).

| Grating | k(λ) median | k(λ) range | α̃ max | β̃ max |
|:--------|:-----------|:-----------|:------|:------|
| G140M NRS2 | **0.846** | 0.706–1.093 | 0.390 | 0.003 |
| G235M NRS2 | **0.976** | 0.834–1.041 | 0.308 | 0.029 |

**Comparison with IFU v1:** G140M agrees to <2% (0.846 vs 0.833); G235M differs by ~8% (0.976
vs 0.899). The k(λ) SHAPES track each other well at λ > 2.5 µm (G140M) and λ > 4.0 µm (G235M).

Notable: FS α̃ is larger than IFU at λ = 2.0–2.25 µm (G140M) and 3.3–3.7 µm (G235M), indicating
the Parlanti model attributes part of the flux excess to 2nd-order contamination.

### NRS2 Calibration Pipeline (root cause found and fixed this session)

The standard JWST pipeline fails to calibrate FS NRS2 for extended wavelengths because:
1. **fflat/sflat `FAST_VARIATION`** tables only cover the NRS1 wavelength range → all NRS2 pixels
   flagged as DO_NOT_USE
2. **photom**: uses `np.interp(..., left=np.nan, right=np.nan)` for out-of-range wavelengths →
   NaN conversion factor → all NRS2 pixels marked DO_NOT_USE

**Fix applied:**
- Use Parlanti (2025) extended flat references: `fflat_0105.fits` (G140M) / `fflat_0091.fits` (G235M)
  and `sflat_0191.fits` (G140M/NRS2) / `sflat_0192.fits` (G235M/NRS2) as overrides
- Use custom extended photom reference `jwst_nirspec_photom_fs_nrs2_ext.fits`
- Skip `wavecorr` (reference lacks NRS2 coverage) and `resample_spec` (s2d→x1d fails;
  extract directly from cal instead)

Key script: `analysis/reduction/run_fs_nrs2_spec2.py`
Produces: `nrs2_spec2_cal/*_nrs2_x1d.fits` in `data/PID1536`,`PID1537`,`PID1538`

---

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

### Next Steps (updated 2026-03-29)
1. **Break k/α degeneracy**: use PID 6644 cool star (NGC2506-G31, G1V) when available
2. **FS NRS2 flat calibration**: current uses IFU flats; investigate FS-specific flat for improved accuracy
3. **Validation**: Apply IFU/FS v1 k(λ) to IRAS-05248 (PID 1492)
4. **ASDF Integration**: Export coefficients to `wavelengthrange_extended.asdf`
5. **More exposures**: download remaining NRS2 rate files for all PIDs (currently 1 per grating/PID)

---
For more details on the long-term project plan, see [IMPLEMENTATION_PLAN.md](IMPLEMENTATION_PLAN.md). To see general project instructions and setup, see [INSTRUCTIONS.md](INSTRUCTIONS.md).

