# NIRSpec Wavelength Extension: "wavext"

We're adding NIRSpec a wavelength extension feature to the JWST pipeline.

## Plan

We'll start with Fixed Slit (FS) mode.

### Data

We'll start with JWST PID 1492:  
* https://www.stsci.edu/jwst/science-execution/program-information?id=1492  
* https://www.stsci.edu/jwst-program-info/download/jwst/pdf/1492/  

## Resources

### Parlanti et al. (IFU)

Parlanti et al.explains how to do wavext for the IFU:  
* https://arxiv.org/abs/2512.14844  
* https://github.com/eleonoraparlanti/nirspecIFU-extended  

1. pipeline / ref file tweaks (simple?)
2. S-flat and F-flat extensions (simple way okay?)
3. CALIBRATION a bit more involved
    * They solve for 2nd and 3rd order spectra coefficients with 3 observations of the same source in both nominal and extended wavelengths

Data: CAL 1536, 1537, 1538, 6645 and GO 2564 & 2186. 

Accuracy: ~5% G235M; ~10% G140M.

### Valentino et al.: DJA msaexp v4 (MOS)


* Valentino et al.: https://arxiv.org/html/2503.01990
* DJA: https://dawn-cph.github.io/dja/spectroscopy/nirspec/
* msaexp wavext: https://github.com/gbrammer/msaexp/blob/main/msaexp/pipeline_extended.py

extending the `wavelengthrange` pipeline reference file

## Technical Implementation Detail

### Code Changes
- **IFU NRS2 Fix:** In `jwst/assign_wcs/nirspec.py`, remove the hard-coded `NoDataOnDetectorError` for NRS2 M-gratings (lines 284-292).

### Reference File Overrides
- **wavelengthrange:** Create a custom ASDF file with extended ranges (e.g., [0.6, 3.3] µm for G140).
- **Usage:** Override in pipeline via `--override_wavelengthrange my_wavrange.asdf`.

### Mode Handling
- **FS/MOS/BOTS:** Works via `wavelengthrange` override as long as WCS maps to detector pixels.
- **IFU:** Requires both the code fix and the `wavelengthrange` override.
