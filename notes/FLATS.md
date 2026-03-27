# NIRSpec Flats for Wavelength Extension (wavext)

This document explains the various flat field components for NIRSpec and how they have been extended or modified to support the **wavext** feature (wavelength extension to 0.6 – 5.6 µm).

## Overview of NIRSpec Flat Field Components

NIRSpec flat fielding is composite, consisting of three primary components:

1.  **D-flat (Detector Flat)**: Corrects for pixel-to-pixel and large-scale variations in the detector response.
2.  **S-flat (Spectrograph Flat)**: Corrects for throughput variations from the slit/aperture through the spectrograph optics to the disperser.
3.  **F-flat (Fore-optics Flat)**: Corrects for throughput variations from the telescope and NIRSpec fore-optics.

Each component typically consists of two parts:
- **Slow variation**: A 2D image (or 3D image with image planes for different wavelengths) capturing spatial variations.
- **Fast variation**: A 1D table capturing high-frequency wavelength-dependent variations (the "fast" variation).

## Extensions for wavext

To support wavelength extension beyond the nominal limits (e.g., NRS2 for M-gratings), we follow the strategy outlined in **Parlanti et al. (2025)**.

### S-flats (Spectroscopic Flat)
- **Component**: Aperture –> Disperser.
- **Extension Strategy**:
    - **Slow Variation**: Set to **1.0** (constant) for all wavelengths in the extended region.
    - **Fast Variation**: Uses the **value at the reddest wavelength** from the nominal flat field table to extrapolate into the red.
- **Implementation**: The pipeline's `CombineFastSlow` and `InterpolateFlat` functions in `jwst.flatfield.flat_field` were modified to support this. Instead of flagging regions outside the reference file range as `DO_NOT_USE`, the code now extrapolates using the reddest available value and sets the slow variation to 1.0.

### F-flats (Fore-optics Flat)
- **Component**: Telescope + NIRSpec fore-optics.
- **Extension Strategy**:
    - **Table-based concatenation**: Fore-optics throughput is broader than the gratings. Parlanti et al. concatenated the fast variations across different gratings (e.g., G140M + G235M + G395M) to build a unified F-flat model.
    - **Red-end extension**: Similar to S-flats, for any region remaining outside the concatenated range, the value at the reddest available wavelength is used.

### D-flats (Detector Flat)
- **Component**: Detector-only variations.
- **Extension Strategy**: Rely on the last available detector plane for the red-end extension, assuming the detector response is relatively well-behaved at these longer wavelengths.

## Mode-Specific Handling

| Mode | S-flat Extension | F-flat Extension |
| :--- | :--- | :--- |
| **Fixed Slit (FS)** | Reddest value extension in `flat_field.py`. | Concatenated M-grating table in reference files. |
| **MOS** | Same as FS; handled slitlet-by-slitlet. | Concatenated reference files. |
| **IFU** | **Critical:** Slow variation = 1.0 beyond reference range. | Concatenated gratings as per Fig A.1. |
| **BOTS** | Follows FS strategy (S1600A1). | Same as FS. |

## Visualization of Extensions (Parlanti Strategy)

The Parlanti et al. strategy for IFU S-flats is illustrated in `references/Parlanti Fig A1.png`. By setting the slow variation to 1.0 and keeping the fast variation constant at the reddest wavelength, the pipeline can provide a stable (if approximate) flat field correction for the extended 2nd and 3rd order light, which is then refined during the throughput derivation phase.

![Parlanti Fig A1](/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/references/Parlanti Fig A1.png)
*Fig. A.1: Illustration of the IFU S-flat extension strategy (Parlanti et al. 2025).*

## Pipeline Modifications applied to `jwst/flatfield/flat_field.py`

The following logic was implemented:
- **`combine_fast_slow`**: Replaced `right=np.nan` with `right=tab_flat[-1]` to allow red-end extrapolation using the last available table value.
- **`interpolate_flat`**: Modified to avoid flagging the red-end extended region as `DO_NOT_USE`. For wavelengths longer than the reference range, the slow variation (image component) is now explicitly set to `1.0`.
