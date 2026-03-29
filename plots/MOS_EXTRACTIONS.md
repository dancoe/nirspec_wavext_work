# NIRSpec Wavelength Extension: MOS-Style Extraction Regions

This document visualizes the count rate data and the 2D extraction regions used for NIRSpec wavelength extension on the NRS2 detector.

## Extraction regions on Detector (NRS2)

### Diagnostic Plot Grids (Absolute Detector Coordinates)
The following grids show the absolute detector coordinates and synchronized row-wise cropping for NRS1 and NRS2. The red dashed line denotes the 5.3 µm calibration limit.

````carousel
![PID 1492 v8 Grid](/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/plots/extraction/PID_1492_All_Gratings_NRS1_NRS2_v8.png)
<!-- slide -->
![PID 1537 v8 Grid](/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/plots/extraction/PID_1537_All_Gratings_NRS1_NRS2_v8.png)
<!-- slide -->
![PID 1538 v8 Grid](/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/plots/extraction/PID_1538_All_Gratings_NRS1_NRS2_v8.png)
````

### Key Features:
- **Extended Coverage**: The red extraction regions cover the full physical footprint of the spectra on the NRS2 detector, extending beyond the nominal cutoffs (e.g., reaching up to ~3.3 µm for G140M, where it typically stops at ~1.9 µm).
- **Contamination Modeling**: By extracting these extended regions, we capture the 2nd and 3rd order light necessary to solve the Parlanti et al. (2025) throughput model.
- **Detector Projects**: For these slit/grating combinations, we maximize the information from the available detector pixels, right up to the edges.

## Methodology
The plots were generated using `mos.py` which interfaces with the JWST datamodels to accurately overlay `MultiSlitModel` extraction parameters on the raw detector data.

The extraction was enabled by:
1.  **Reference File Override**: Using `wavelengthrange_extended.asdf` during `assign_wcs`.
2.  **Pipeline Patching**: Fixing lazy-load ASDF issues in `jwst` to prevent premature file closure.
3.  **Spectral Order Mapping**: Correctly mapping wavelength to pixel coordinates for the NRS2 detector in gratings that nominally only project onto NRS1.
