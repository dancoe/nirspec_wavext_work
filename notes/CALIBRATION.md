# CALIBRATION.md

## NIRSpec Wavelength Extension Calibration Theory

This theory is adopted from Parlanti et al. (arXiv:2512.14844), which accounts for 2nd and 3rd order spectral contamination during extended wavelength extraction.

### Order Overlap Equation

The observed spectrum $S(\lambda)$ at extended wavelengths is modeled as a sum of the 1st, 2nd, and 3rd order contributions:

$$S(\lambda) \approx k(\lambda) f(\lambda) + a(\lambda) f(\lambda/2) + b(\lambda) f(\lambda/3)$$

Where:
*   $f(\lambda)$ = Intrinsic spectrum of the source.
*   $S(\lambda)$ = Measured (observed) spectrum from the extended extraction.
*   $k(\lambda)$ = 1st-order throughput contribution (drops to ~30% at the longest extended wavelengths).
*   $a(\lambda)$ = 2nd-order contamination contribution (< ~5%).
*   $b(\lambda)$ = 3rd-order contamination contribution (< ~1%).

### Calibration Strategy

To determine the unknown throughput functions $k$, $a$, and $b$, we require a minimum of three sources with known intrinsic spectra:

1.  **Requirement:** At least 3 target sources where we have:
    *   **$f(\lambda)$:** Intrinsic flux obtained from the standard nominal extraction (where contamination is negligible or well-calibrated).
    *   **$S(\lambda)$:** Measured flux from the extended extraction (including contaminated wavelengths).
2.  **Step 1 (Preprocessing):** Extract the data using the custom `wavelengthrange` reference file to extended wavelengths.
3.  **Step 2 (Analysis):** Measure the discrepancies between the nominal and extended regions to isolate $S(\lambda)$.
4.  **Step 3 (Solving):** Solve the system of equations across the target samples to calculate $k(\lambda)$, $a(\lambda)$, and $b(\lambda)$ for each grating/filter combination.

### Priority
The initial calibration efforts will focus on **Fixed Slit (FS)** data from Program ID 1492.
