# 🚀 NIRSpec FS v7 — Full Flux Calibration & Validation Report
**Natively calibrated Jy Spectra and Hold-out Cross-Validation**

## 1. Summary of Achievements
In v7, we successfully transitioned the NIRSpec wavelength extension calibration to a natively flux-calibrated workflow in **Jy** units. 
- **Native Jy Extractions**: NRS2 spectra produced by `Spec2Pipeline` are now natively in Jy, eliminating the need for manual scaling (DN/s → Jy) used in v6.
- **Improved Source Data**: Leveraged all exposures and multiple dithers across PIDs 1492, 1536, 1537, 1538, and 6644.
- **Robust Model Validation**: Verified generality via systematic **leave-one-out cross-validation**.

## 2. Full Spectrum Baseline Verification (v5 Style)
The following plots show the final NIRSpec NRS1,2 spectrum compared to the CALSPEC truth model across the full 0.6–5.6 µm range. These confirm the recovery of the first-order SED after accurate ghost subtraction using the v7 coefficients.

### P330E (G2V) — PID 1538
![P330E Full Spectrum](plots/full_spectrum_merged_v7_P330E.png)

### G191-B2B (WD) — PID 1537
![G191-B2B Full Spectrum](plots/full_spectrum_merged_v7_G191-B2B.png)

### J1743045 (A8III) — PID 1536
![J1743045 Full Spectrum](plots/full_spectrum_merged_v7_J1743045.png)

### NGC2506-G31 (G1V) — PID 6644
![NGC2506-G31 Full Spectrum](plots/full_spectrum_merged_v7_NGC2506-G31.png)

### IRAS-05248 (AGN) — PID 1492
![IRAS-05248 Full Spectrum](plots/full_spectrum_merged_v7_IRAS-05248.png)

## 3. Key Results: v7 vs v6 Consistency
The coefficients derived in v7 match the previous v6 results to high precision, confirming consistency while streamlining the data processing.

**G140M**  
![v7 vs v6 Comparison G140M](plots/fs_v7_v6_comparison_g140m.png) 

**G235M**  
![v7 vs v6 Comparison G235M](plots/fs_v7_v6_comparison_g235m.png) 

## 4. Generality & Hold-out Validation
The hold-out validation proves that the ghost model correctly predicts spectral contamination even for sources not included in the fit.

### G140M Hold-out Performance
![G140M Holdouts](plots/fs_v7_holdout_validation_g140m.png)

### G235M Hold-out Performance
![G235M Holdouts](plots/fs_v7_holdout_validation_g235m.png)

## 5. Full Spectrum Validation (NRS2 only)
For each source, the observed NRS2 flux (including ghosts) is compared to the model $k f_1 + \alpha f_2 + \beta f_3$ using the derived coefficients.

### G140M Full Spectra
![G191-B2B](plots/fs_v7_full_spec_g140m_1537.png) (PID 1537)
![P330E](plots/fs_v7_full_spec_g140m_1538.png) (PID 1538)
![J1743045](plots/fs_v7_full_spec_g140m_1536.png) (PID 1536)
![NGC2506-G31](plots/fs_v7_full_spec_g140m_6644.png)  (PID 6644)
![IRAS-05248](plots/fs_v7_full_spec_g140m_1492.png) (PID 1492)

### G235M Full Spectra
![G191-B2B](plots/fs_v7_full_spec_g235m_1537.png) (PID 1537)
![P330E](plots/fs_v7_full_spec_g235m_1538.png) (PID 1538)
![J1743045](plots/fs_v7_full_spec_g235m_1536.png) (PID 1536)
![NGC2506-G31](plots/fs_v7_full_spec_g235m_6644.png)  (PID 6644)
![IRAS-05248](plots/fs_v7_full_spec_g235m_1492.png) (PID 1492)

## 6. Derived Products
Final calibration coefficients are saved in:
- `results/v7/calib_v7_fs_g140m_all.fits`
- `results/v7/calib_v7_fs_g235m_all.fits`

The next phase will apply these coefficients as a final global correction to all NIRSpec science targets for the project.
