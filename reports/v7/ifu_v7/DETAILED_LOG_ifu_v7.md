# NIRSpec IFU v7 Calibration — Execution Log (2026-03-30)

### 5a. IFU Photom Analysis
- **Finding**: `photom_0016.fits` (IFU) has `photmj=1.0` and `relresponse=1.0` for ALL entries (nominal + extended wavelengths).
- **Implication**: No extended photom reference file needed for IFU (unlike FS). Flux calibration in IFU comes from cube_build geometry, pathloss, and Parlanti sflat/fflat — NOT from the photom step.
- **IFU x1d units**: Jy (already, from all previous stage3_ext runs via cube_build+extract_1d).

### 5b. Data Availability Check
| PID | Mode | Rate Files | Status |
|:----|:-----|:-----------|:-------|
| 1536 | IFU | Re-downloaded from MAST → `data/IFU/PID1536_J1743045/stage1/` (24 files) | ✅ |
| 1537 | IFU | Re-downloaded from MAST → `data/IFU/PID1537_G191-B2B/stage1/` (24 files, obs IDs: 04101=G140M, 06101=G235M, 08101=G395H) | ✅ |
| 1538 | IFU | Re-downloaded from MAST → `data/IFU/PID1538_P330E/stage1/` (obs IDs: 03102=G140M, 03104=G235M — different from 1536/1537!) | ✅ |
| 2186 | IFU | Already local in `data/PID2186_UGC-5101/rate/` (8 files, obs ID 04101=G235M for UGC5101) | ✅ |

### 5c. Pipeline Re-Reductions
| PID | Stage2 | Stage3 | Total time | Output |
|:----|:-------|:-------|:-----------|:-------|
| 1536 | ✅ 16 cal (Mar 30 18:50–19:02) | ✅ 2 S3D + 2 X1D (Mar 30 19:08–19:14) | 24.3 min | G140M: 0.970–3.600µm; G235M: 1.661–5.500µm |
| 1537 | ✅ 16 cal (Mar 30 19:00–19:13) | ✅ 2 S3D + 2 X1D (Mar 30 19:19–19:25) | 25.4 min | G140M: 0.970–3.600µm; G235M: 1.661–5.500µm |
| 1538 | ✅ 16 cal (Mar 30 19:15–19:30) | ✅ 2 S3D + 2 X1D (Mar 30 19:36–19:42) | 27.2 min | G140M: 0.970–3.600µm; G235M: 1.661–5.500µm |
| 2186 | — (used existing cal files) | ✅ G235M only: 1.661–5.500µm (Mar 30 18:57) | 9.7 min | UGC5101 extended |

- **Note on 1538 obs IDs**: PID 1538 uses `03102`=G140M and `03104`=G235M (not `04101`/`06101` like 1536/1537). The download watcher  
  was updated with this info for future re-runs.

### 5d. IFU v7 Solver
- **Script**: `reports/v7/scripts/solve_parlanti_ifu_v7.py`
- **Sources**: G191-B2B (PID 1537), P330E (PID 1538), J1743045 (PID 1536)
- **Note**: NGC2506-G31 (G1V) not in IFU dataset → cannot fully break k/α degeneracy (only G2V + hot stars)

#### Final Results (all 3 freshly re-reduced sources):
| Grating | Sources | Hold-out G191-B2B | Hold-out P330E | Hold-out J1743045 |
|:--------|:--------|:------------------|:---------------|:------------------|
| G140M | 3 | median 2.4% | median 5.7% | median -6.1% |
| G235M | 3 | median -2.6% | median 3.5% | median -1.6% |

_Preliminary (old data): G140M [2.5%, 6.0%, -6.3%]; G235M [-3.0%, 3.7%, -1.9%] — consistent within ~0.4%_
_Note: PID 1538 re-reduction yielded identical x1d to old data (photom=1.0 makes stage2+3 deterministic with same Parlanti reffiles)_

- All median residuals < 10% ✅ (v7 target criterion met)
- IFU k(λ) consistently lower than FS k(λ): IFU ≈ 0.4–1.0 vs FS ≈ 1.0–1.6 (G140M), IFU ≈ 0.3–1.0 vs FS ≈ 0.7–2.1 (G235M)
- IFU α(λ) peak at 2.0µm ≈ 0.18 (G140M), consistent with first-harmonic order leakage; α shows second peak at 4.5µm (G235M)
- G235M β spike at >5.05µm — artifact from near-zero SNR at long-wavelength edge; calibration invalid >4.8µm
- G235M reliable range: 3.15–4.8µm; G140M reliable range: 1.87–2.9µm (>2.9µm affected by 3µm H₂O absorption)

### 5e. Plots Generated (intermediate, in `reports/v7/ifu_v7/plots/`)
- `ifu_v7_coefficients_g140m.png` — k/α/β for all sources + hold-outs ✅
- `ifu_v7_coefficients_g235m.png` ✅
- `ifu_v7_holdout_residuals_g140m.png` — per-source fractional residuals ✅
- `ifu_v7_holdout_residuals_g235m.png` ✅
- `ifu_v7_full_spectrum_g140m.png` — NRS1+corrected NRS2 vs CALSPEC ✅
- `ifu_v7_full_spectrum_g235m.png` ✅
- `ifu_v7_vs_fs_g140m.png` — IFU vs FS coefficient comparison ✅
- `ifu_v7_vs_fs_g235m.png` ✅

### 5f. Next Steps (as of Mar 30, 2026)
- ✅ 1538 stage2+3 COMPLETE (27.2 min, Mar 30 19:15–19:42)
- ✅ Final solver re-run on all 3 fresh sources (results unchanged vs intermediate — see 5d)
- ✅ OBSERVATIONS.md: all 4 PIDs marked ✅
- [ ] Add UGC-5101 (PID 2186) G235M validation diagnostic (NRS1–NRS2 join quality check at 3.15µm)
- [ ] Consider adding IFU data from Cycle 2/3 programs (PID 4498, 6606) for additional G2V coverage
- [ ] Write up IFU v7 calibration as a paper section/appendix (analogous to FS v7)
