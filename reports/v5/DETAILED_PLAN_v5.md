# Detailed Analysis Plan v5 — The Degeneracy Breaker

**Version:** 5.1  
**Status:** FS Complete — IFU Validation Pending  
**Objective:** Finalize the NIRSpec Wavelength Extension (0.6–5.6 µm) by breaking the $k/\alpha$ degeneracy using a 4th calibration star (G1V) and validating against archival AGN science targets.

---

## 1. Data Acquisition

### 1.1 Calibration Standards (FS & IFU)
- **Status:** All priority standards are already present on disk (1536, 1537, 1538).
- **New Target (PID 6644):** NGC2506-G31 (G1V).
- **CALSPEC Truth Model:** `ngc2506g31_mod_003.fits` (Required for G1V truth spectrum).

---

## 2. Reduction Pipeline (v5-Extended)

All targets must be processed through the extended wavelength pipeline to capture the 3nd-order features ($0.6–1.1$ µm) and 2nd-order ($1.1–1.9$ µm) ghosts.

### 2.1 FS Reduction (PID 6644)
Run the dedicated FS L3 reduction scripts.

```bash
# NRS1 (NominalTruth)
micromamba run -n jwst_1.20.2 python nirspec_wavext_work/analysis/reduction/run_fs_nrs1_spec3.py --pid 6644

# NRS2 (Extended)
micromamba run -n jwst_1.20.2 python nirspec_wavext_work/analysis/reduction/run_fs_nrs2_spec3.py --pid 6644
```

### 2.2 IFU Science Validation (PIDs 2654, 2186)
Once downloads complete, execute the IFU extension pipeline for the AGN and ULIRG targets.

```bash
# Reduce with extended WCS and Parlanti reference overrides
micromamba run -n jwst_1.20.2 python nirspec_wavext_work/analysis/reduction/run_ifu_pipeline_ext.py --pid 2654
micromamba run -n jwst_1.20.2 python nirspec_wavext_work/analysis/reduction/run_ifu_pipeline_ext.py --pid 2186
```

---

## 3. Analysis & Calibration Solving

### 3.1 4-Source Simultaneous NNLS
Run the solver on the 4 FS standards to decouple throughput ($k$) from ghosting ($\alpha, \beta$).

- **Source Set:** 1536 (A8III), 1537 (WD), 1538 (G2V), 6644 (G1V).
- **Technical Lessons:**
    - **Recursive Discovery:** `glob.glob` with `recursive=True` is now used to auto-stitch Level-3 x1d.fits files across disparate directory structures (Spec3 outputs vs. legacy download).
    - **Degeneracy:** The G1V star features allow NNLS to properly assign flux to $k(\lambda)$ (~0.96) instead of over-assigning to $\alpha(\lambda)$.
- **Analysis Script Location:** `nirspec_wavext_work/reports/v5/fs_v5/scripts/`

---

## 4. Diagnostics & Visualization

### 4.1 Coefficient Comparison
Coefficient plots generated via `plot_parlanti_calib_v5.py` and `plot_parlanti_coeffs_log_v5.py` compare:
- **v5 Solution** vs. **Parlanti (2025)** vs. **v4 (3-star degeneracy)**.
- **2-Panel Summary:** [v5 Multi-Grating Log Plot](fs_v5/plots/fs_v5_coeffs_log_2panel.png)
- Result: $\alpha$ dropped back to ~0.003–0.012 and $k$ returned to ~0.96 baseline.

### 4.2 Full-Spectrum Validation
Residual plots compare NRS1 vs NRS2 (Extended).
- Accuracy: Core regions match truth to within ~2–5%.

### 4.3 Full Spectrum Merged (Truth Validation - Standards)
The merged spectra across 0.6–5.6 µm for the 4 FS standards comparing NRS1 + Cor-NRS2 with CALSPEC.
- [P330E Full Spectrum (v5)](fs_v5/plots/full_spectrum_merged_v5_P330E.png)
- [G191-B2B Full Spectrum (v5)](fs_v5/plots/full_spectrum_merged_v5_G191-B2B.png)
- [J1743045 Full Spectrum (v5)](fs_v5/plots/full_spectrum_merged_v5_J1743045.png)
- [NGC2506-G31 Full Spectrum (v5)](fs_v5/plots/full_spectrum_merged_v5_NGC2506-G31.png)

### 4.4 Full Spectrum Merged (Science Validation - IFU)
Merged spectra for science targets showing successful SED recovery in the extended range.
- [UGC-5101 Full Spectrum (v5)](ifu_v5/plots/ifu_v5_full_spectrum_UGC-5101.png)
- [SDSSJ0841 Full Spectrum (v5)](ifu_v5/plots/ifu_v5_full_spectrum_SDSSJ0841.png)

---

## 5. Reporting

### 5.1 FS v5 Report (COMPLETED)
- **Report Link:** [FS v5](v5/fs_v5/REPORT_fs_v5.md)
- **Summary Script:** [summarize_v5.py](v5/fs_v5/scripts/summarize_v5.py) for quick coefficient inspection.

### 5.2 IFU v5 Validation Report (IN PROGRESS)
- **Report Link:** [IFU v5 Science Validation](v5/ifu_v5/REPORT_ifu_v5.md)
- **Status:** H-alpha recovery verified for SDSSJ0841; L-band dust verified for UGC-5101.

---

## 6. Detailed Step-by-Step Command Log

| Phase | Step | Command |
| :--- | :--- | :--- |
| **Download** | Complete Science | (Auto-resuming background download) |
| **Reduce** | 6644 Reduction | `python nirspec_wavext_work/analysis/reduction/run_fs_nrs1_spec3.py --pid 6644` |
| **Solve** | v5 NNLS Joint | `python nirspec_wavext_work/reports/v5/fs_v5/scripts/solve_parlanti_calib_v5.py` |
| **Plot FS** | Diagnostic Plots | `python nirspec_wavext_work/reports/v5/fs_v5/scripts/plot_parlanti_calib_v5.py` |
| **Plot IFU** | Validation Plots | `python nirspec_wavext_work/reports/v5/ifu_v5/scripts/plot_ifu_v5_validation.py` |
| **Plot Full** | Merged Spectra | `python nirspec_wavext_work/reports/v5/ifu_v5/scripts/plot_ifu_full_spectra_v5.py` |
| **Summary**| Inspect Coeffs | `python nirspec_wavext_work/reports/v5/fs_v5/scripts/summarize_v5.py` |

---
*Created by Antigravity — v5 Execution Readiness: Phase 5/Phase 6*
