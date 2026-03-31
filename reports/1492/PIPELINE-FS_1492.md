# NIRSpec PID 1492 Full Fixed Slit Reduction (v7)

This document records the full re-reduction of NIRSpec PID 1492 data, starting from the Level 1 `rate` files and incorporating the v7 wavelength extension for NRS2.

## 1. Data Inventory
- **Location:** `data/PID1492/MAST/rate/` (Original MAST copies)
- **Local Workspace:** `data/PID1492/`
- **Total Exposures:**
    - **NRS1:** ~48 rate files (Nominal coverage)
    - **NRS2:** ~7 rate files (Extended coverage validation)

## 2. Reduction Strategy: v7 Native Jy Workflow
Following the `JWPipeNB-NIRSpec-FS.ipynb` notebook and the `WAVEXT.md` plan:

1.  **Stage 2 (Spec2Pipeline)**
    - **All Exposures (NRS1 + NRS2):** Run through standard `Spec2Pipeline`.
    - **NRS2 Custom Overrides:** Apply the following to enable wavelength extension:
        - `--override_wavelengthrange wavelengthrange_0008.asdf`
        - `--override_photom jwst_nirspec_photom_fs_nrs2_ext.fits` (Native Jy flux beyond nominal cutoffs)
        - `--override_fflat`: Parlanti-extended fibre-flats.
        - `--override_sflat`: Parlanti-extended slit-flats.
    - **Output:** Produced in `data/PID1492/spec2_full_cal/`.

2.  **Stage 3 (Spec3Pipeline)**
    - Group by Grating/Filter combination.
    - Combine `cal.fits` exposures into final drizzled products.
    - **Wavelength Merging:** Merge NRS1 and NRS2 into a single unified spectrum.

## 3. Scripts Used
- `reduction/run_fs_full_1492.py`: Custom wrapper for full PID 1492 reduction.
- `reports/1492/scripts/plot_1492_extractions.py`: Plotting script for diagnostic 2D extractions.

## 4. Current Status
| Step | Status | Notes |
| :--- | :--- | :--- |
| **Stage 1 (Rate)** | DONE | Downloaded from MAST (in `MAST/rate/`) |
| **Stage 2 (Spec2)** | DONE | Processed 54 rate files. v7 Photom/Flats applied. |
| **Stage 3 (Spec3)** | DONE | Dithered exposures combined (45 products in `spec3_full_cal/`). |
| **Validation Plots** | DONE | `rate_1492_G140M_nrs1_nrs2.png` and combined `x1d` ready. |

## Final Diagnostic Plots

### All Rate Configurations (Representative)
This plot shows one representative exposure for every unique Grating/Filter/Detector configuration in PID 1492, utilizing the pixel-tick styling for improved alignment verification.

![All Rate Configs](file:///Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/reports/1492/plots/all_1492_rate_configs.png)

### Combined Level 3 Spectra (NRS1 + NRS2 Extended)
The following plot shows the final Level 3 extracted spectra for the three primary dithered observations. The **NRS2 (Red)** spectra demonstrate the successful recovery of extended wavelength coverage using the v7 wavelength calibration.

![Level 3 Combined Spectra](file:///Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/reports/1492/plots/spec3_1492_nrs1_nrs2_sidebyside.png)

## Final Level 3 Products
The following combined products are available in `data/FS/PID1492/spec3_full_cal/`:

| Configuration | Product Type | Filename Snippet |
| :--- | :--- | :--- |
| G140M/F100LP | x1d (Combined) | `jw01492_g140m_f100lp_s200a1_nrs1_..._x1d.fits` |
| G140M/F100LP | s2d (Combined) | `jw01492_g140m_f100lp_s200a1_nrs1_..._s2d.fits` |
| G235M/F170LP | x1d (Combined) | `jw01492_g235m_f170lp_s200a1_nrs1_..._x1d.fits` |
| G395M/F290LP | x1d (Combined) | `jw01492_g395m_f290lp_s200a1_nrs1_..._x1d.fits` |
| ... | ... | ... |
