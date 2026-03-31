# 🚀 NIRSpec IFU v7 — Full Flux Calibration & Validation Report
**Natively calibrated Jy Spectra and Hold-out Cross-Validation**

## 1. Summary of Achievements
In v7, we completed the full re-reduction and calibration of the NIRSpec IFU dataset (PIDs 1536, 1537, 1538, and 2186).
- **Native Jy Extractions**: IFU data are natively in Jy units (photom = 1.0). Flux calibration is handled via `cube_build` and Parlanti S-flat / F-flat extensions.
- **Improved Source Data**: Freshly re-processed all calibration sources through Stage 2 and Stage 3 with these updated reference files.
- **Robust Model Validation**: Verified generality via systematic **leave-one-out cross-validation**.

## 2. Full Spectrum Baseline Verification (v5 Style)
The following plots show the final NIRSpec NRS1,2 spectrum compared to the CALSPEC truth model across the full 0.6–5.6 µm range.

### G191-B2B (WD)
![G191-B2B Full Spectrum](plots/full_spectrum_merged_v7_G191-B2B.png)

### P330E (G2V)
![P330E Full Spectrum](plots/full_spectrum_merged_v7_P330E.png)

### J1743045 (A8III)
![J1743045 Full Spectrum](plots/full_spectrum_merged_v7_J1743045.png)

## 3. Key Results: IFU vs FS Comparison
The coefficients derived for IFU show distinct differences from the Fixed Slit (FS) coefficients, particularly in the $k(\lambda)$ throughput term, likely due to different point-source extraction geometries and pathloss corrections.

**G140M**  
![IFU vs FS Comparison G140M](plots/ifu_v7_vs_fs_g140m.png) 

**G235M**  
![IFU vs FS Comparison G235M](plots/ifu_v7_vs_fs_g235m.png) 

## 4. Generality & Hold-out Validation
The joint NNLS solve and hold-out residuals demonstrate the reliability of the v7 IFU model across different standard stars.

### G140M Hold-out Performance
![G140M Holdouts](plots/ifu_v7_holdout_residuals_g140m.png)

### G235M Hold-out Performance
![G235M Holdouts](plots/ifu_v7_holdout_residuals_g235m.png)

## 5. Full Spectrum Validation (IFU Components)
Comparison of NRS1 (nominal) and corrected NRS2 (ghost-subtracted) components against CALSPEC models.

### G140M Components
![G140M Full Spectrum](plots/ifu_v7_full_spectrum_g140m.png)

### G235M Components
![G235M Full Spectrum](plots/ifu_v7_full_spectrum_g235m.png)

## 6. Derived Products
Final IFU calibration coefficients are saved in:
- `results/v7/calib_v7_ifu_g140m_all.fits`
- `results/v7/calib_v7_ifu_g235m_all.fits`

The next phase will incorporate these into the global NIRSpec wavelength extension correction suite.
