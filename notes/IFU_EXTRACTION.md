# NIRSpec IFU Extraction Methodologies

This document compares the standard JWST pipeline `extract_1d` methodology for IFU data with the custom **v8** extraction used in this project.

## 1. Standard Pipeline Default (`extract_1d`)

For NIRSpec IFU data processing through the standard `calwebb_spec3` pipeline, the `extract_1d` step follows these rules:

### A. Point-Source Extraction (`SRCTYPE = POINT`)
- **Aperture:** Uses a **circular aperture** centered on the source (based on header coordinates or `ifu_autocen` centroid).
- **Radius:** The aperture radius is **wavelength-dependent**. It scales with wavelength to follow the changing Point Spread Function (PSF) / Airy disk across the spectral range.
- **Reference File:** Radii are defined in the `EXTRACT1D` reference file, specifically **`jwst_nirspec_extract1d_0002.asdf`** for IFU data in the current context.
- **Background:** Traditionally applies a **circular background annulus** with wavelength-dependent inner and outer radii.
- **Correction:** Applies a wavelength-dependent **aperture correction** (throughput correction) to the extracted flux to account for light lost outside the aperture, bringing the final flux to the "infinite aperture" level.

### B. Extended-Source Extraction (`SRCTYPE = EXTENDED`)
- **Aperture:** Sums the **entire spatial area** of the IFU cube for each wavelength slice (rectangular extraction).
- **Background:** No background subtraction is performed.
- **Correction:** No aperture correction is applied.

---

## 2. Project v8 Custom Extraction

For the v8 science validation of targets (e.g., **UGC-5101**, **SDSSJ0841**), a custom aperture extraction is implemented in the script `extracted_ifu_v8.py`.

### v8 Parameters:
- **Aperture:** Fixed circular aperture.
- **Radius:** **0.5" radius** (non-wavelength dependent). 
- **Centering:** The aperture is centered on the **peak pixel** located in a "white-light" image (the mean of the `_s3d.fits` cube along the spectral axis). This ensures robust alignment regardless of header coordinate offsets.
- **Background:** No background annulus subtraction is performed (direct summation of the 0.5" circle).
- **Correction:** **No aperture correction** is applied. Because the 0.5" radius is significantly larger than the standard pipeline default (which tracks the inner PSF core), it captures a much static fraction of the total flux across all wavelengths.

### Rationale for v8 Customization:
Standard pipeline extractions often provide smaller apertures that are sensitive to the accuracy of the ASDF-based PSF model. The larger 0.5" fixed radius is used for these complex science targets to capture more extended emission and provide a more stable comparison between different pipeline/reference file versions (e.g., v7 vs. v8).

---

## Technical Summary Comparison

| Feature | Standard Pipeline (Point Source) | Project v8 Custom (Science) |
| :--- | :--- | :--- |
| **Aperture Shape** | Circular | Circular |
| **Aperture Size** | **Wavelength-dependent (variable)** | **Fixed 0.5" Radius** |
| **Background** | Subtracted using Annulus | No background subtraction |
| **Correction** | Aperture Correction applied | No correction applied |
| **Centering** | Header / Auto-centroid | **Peak of White-Light Image** |
| **Source File** | `*_s3d.fits` Level 3 Data Cube | `*_s3d.fits` Level 3 Data Cube |

*Reference Implementation: [extracted_ifu_v8.py](../reports/v8/scripts/extracted_ifu_v8.py)*
