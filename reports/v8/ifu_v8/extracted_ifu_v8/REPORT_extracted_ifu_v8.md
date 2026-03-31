# IFU v8 Extracted Spectra Diagnostics

Detailed visualizations and comparison between the v8 custom 0.5" fixed circular aperture vs. the default pipeline extraction.

## 1. Pipeline Reference Parameters
The default pipeline extraction radii were obtained from the CRDS reference file:
- **Reference File**: `jwst_nirspec_extract1d_0002.asdf`
- **POINT Source Radius**: 0.45" (Fixed across wavelength in this version)
- **EXTENDED Source**: Sums entire spatial window (full frame)

Full spectrum radii table: [pipeline_radii.csv](pipeline_radii.csv)

## UGC-5101 (PID 2186)

### Spectral Slices and Extraction Regions
- **Red Solid Circle**: v8 extraction (r=0.5")
- **Cyan Dashed**: Default Pipeline (Fixed 0.45" for POINT or Full-Frame for EXTENDED)

![UGC-5101 Slices](slices_2186.png)

### Spectrum and Ratio Comparison
![UGC-5101 Extraction](extraction_2186.png)

---

## SDSSJ0841 (PID 2654)

### Spectral Slices and Extraction Regions
- **Red Solid Circle**: v8 extraction (r=0.5")
- **Cyan Dashed**: Default Pipeline (Fixed 0.45" for POINT or Full-Frame for EXTENDED)

![SDSSJ0841 Slices](slices_2654.png)

### Spectrum and Ratio Comparison
![SDSSJ0841 Extraction](extraction_2654.png)

---

