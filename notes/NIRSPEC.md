# JWST NIRSpec Technical Specifications

This document provides foundational technical details for the JWST Near-Infrared Spectrograph (NIRSpec) instrument, as used in the Wavelength Extension (WAVEXT) project.

## General Instrument Parameters

| Parameter | Specification | Notes |
| :--- | :--- | :--- |
| **Wavelength Range** | 0.6 – 5.3 μm | Extended to 0.5 – 5.6 μm in this project. |
| **Pixel Scale** | ~0.1 arcsec/pixel | SAM-sampled on the detectors. |
| **Detectors** | 2 H2RG arrays (NRS1, NRS2) | 2048 x 2048 pixels each. |

---

## Observing Modes

### Integral Field Unit (IFU)
- **Field of View (FOV)**: 3.0" x 3.0"
- **Design**: Dissected into 30 slices (0.1" wide and 3.0" long) to provide 3D data cubes.

### Fixed Slit (FS)
NIRSpec has five fixed slits located on the focal plane, primarily for point-source and high-contrast spectroscopy:
- **S200A1**: 0.2" wide
- **S200A2**: 0.2" wide
- **S400A1**: 0.4" wide
- **S1600A1**: 1.6" x 1.6" (Square aperture, used for high throughput and Bright Object Time-Series / BOTS mode)
- **S200B1**: 0.2" wide (Primarily for redundancy)

### Multi-Object Spectroscopy (MOS / MSA)
- **Field of View (FOV)**: 3.6' x 3.4' (Over the full 4 quadrants)
- **Micro-Shutter Assembly (MSA)**:
  - **Quadrants**: 4 quadrants arranged in a 2x2 grid.
  - **Shutters per Quadrant**: 365 (in the spectral dispersion direction) x 171 (in the spatial direction).
  - **Total Shutters**: 62,415 shutters per quadrant (~250k total across the array).
- **Physical Specifications**:
  - **Clear Aperture**: 0.20" x 0.46" (per shutter)
  - **Shutter Pitch**: 0.27" x 0.53" (center-to-center)
  - **MOS Bars**: Each bar between shutters is 0.07" wide.

---

## Gratings and Resolving Powers (R)

NIRSpec achieves spectroscopy at three spectral resolutions ($R = \lambda / \Delta \lambda$):

| Type | Grating / Dispersion Element | Resolving Power ($R$) |
| :--- | :--- | :--- |
| **Low** | CLEAR / PRISM | ~100 |
| **Medium** | G140M, G235M, G395M | ~1000 |
| **High** | G140H, G235H, G395H | ~2700 |

### Grating-Filter Combinations & Wavelength Ranges
Note: "Extended Range" refers to the scientific potential of the detector real-estate captured during the WAVEXT effort (leveraging NRS2 coverage).

| Grating | Filter | Nominal Range (μm) | Extended Range (μm) |
| :--- | :--- | :--- | :--- |
| **PRISM** | CLEAR | 0.6 – 5.3 | 0.50 – 5.60 |
| **G140M/H** | F070LP | 0.70 – 1.27 | 0.60 – 3.30 |
| **G140M/H** | F100LP | 0.97 – 1.89 | 0.90 – 3.30 |
| **G235M/H** | F170LP | 1.66 – 3.17 | 1.50 – 5.30 |
| **G395M/H** | F290LP | 2.87 – 5.27 | 2.60 – 5.60 |

---

## Mode-Specific Notes
- **Detector Gaps**: Due to the physical gap between NRS1 and NRS2, there is a small loss in spectral coverage (wavelength "gap") that varies by grating and filter. Dithering is typically used to recover these wavelengths.
- **WAVEXT Impact**: The project specifically targets the "extra" spectral area on NRS2 that is not currently part of the standard pipeline calibration, extending the science utility beyond 5.3 μm or into overlapping grating regions.
