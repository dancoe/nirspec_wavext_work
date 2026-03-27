# NIRSpec Flats for Wavelength Extension (wavext)

This document explains the various flat field components for NIRSpec and how they are extended to support the **wavext** feature (wavelength extension to 0.6 – 5.6 µm).

> [!IMPORTANT]
> To support wavelength extension, we follow the strategy of **extending reference files** directly. I have discovered that the latest CRDS NIRSpec S-flat library (02XX series, activated 2025-06) **already incorporates many of these extensions**, including MOS PRISM spatial maps that now extend to **5.55 µm**.

## Overview of NIRSpec Flat Field Components

NIRSpec flat fielding is composite, consisting of three primary components:

1.  **D-flat (Detector Flat)**: Pixel-to-pixel and large-scale variations in the detector response.
2.  **F-flat (Fore-optics Flat)**: Throughput variations from the telescope and NIRSpec fore-optics.
3.  **S-flat (Spectrograph Flat)**: Throughput variations from the slit/aperture through the spectrograph optics to the disperser.

---

## 2D Spatial Sensitivity Atlas (Side-by-Side Plots)

The maps below use a diverging **seismic** colormap (0.9 – 1.1) to highlight spatial artifacts. All plots include explicit x-column/y-row labels and colorbars.

### Wide-Field Detector Flat (D-flat)
The D-flat captures pixel-level noise, quadrant artifacts, and large-scale gradients.

````carousel
![D-flat: 0.80000 µm](../plots/flats/d-flats/d-flat_0.80000um.png)
<!-- slide -->
![D-flat: 1.00000 µm](../plots/flats/d-flats/d-flat_1.00000um.png)
<!-- slide -->
![D-flat: 5.30000 µm](../plots/flats/d-flats/d-flat_5.30000um.png)
````

### MOS PRISM S-flat (12 Planes)
These plots, located in `plots/flats/s-flats/mos/prism/`, capture the full spectral evolution of spatial sensitivity for MOS PRISM.
*(Note: The latest downloaded PRISM S-flat plots currently only have data on NRS1, when there should also be data on NRS2.)*

````carousel
![PRISM: 0.36393 µm](../plots/flats/s-flats/mos/prism/mos_prism_s-flat_0.36393um.png)
<!-- slide -->
![PRISM: 1.30774 µm](../plots/flats/s-flats/mos/prism/mos_prism_s-flat_1.30774um.png)
<!-- slide -->
![PRISM: 2.72346 µm](../plots/flats/s-flats/mos/prism/mos_prism_s-flat_2.72346um.png)
<!-- slide -->
![PRISM: 4.13918 µm](../plots/flats/s-flats/mos/prism/mos_prism_s-flat_4.13918um.png)
<!-- slide -->
![PRISM: 5.55489 µm](../plots/flats/s-flats/mos/prism/mos_prism_s-flat_5.55489um.png)
````

### MOS M-Gratings S-flats
We prioritize spatial extension plots for the Medium resolution (M) gratings.

*(Additional visual maps can be found under `plots/flats/s-flats/mos/g140m`, etc.)*

### Fixed Slit (FS) and IFU S-flats
These usually consist of a single spatial footprint characterizing the slit geometry.

**IFU S-flat Example (G140M F100LP):** Dual-detector comparison showing spatial sensitivity variations across NRS1 and NRS2.

![IFU S-flat](../plots/flats/dual_detector_ifu_sflat.png)

![FS S-flat](../plots/flats/s-flats/fs/fs_s-flat.png)

---

## Reference File Extension Strategy (wavext)

Following **Parlanti et al. (2025)**, we extend any reference files that do not already reach 5.61 µm:

- **Slow Variation (x, y, λ)**: Append a plane initialized to **1.0** (unity) at 5.61 µm.
- **Fast Variation (λ)**: Append a **top-hat** extension from the reddest nominal wavelength to 5.61 µm.

Reproduction of the strategy from the literature:

![Parlanti Repro](../plots/Parlanti/FS-1492_pre-cal.png)

---

---

## Flat Reference File Inventory

For a comprehensive, audited list of all `DFLAT`, `FFLAT`, and `SFLAT` files used in this project, including wavelength sampling details and detector coverage, please refer to:

👉 **[FLATS_FILES.md](FLATS_FILES.md)**

### Inventory Summary
- **D-flats:** 2 files (One per detector, 0.6 – 5.3 µm)
- **F-flats:** 4 files (Shared across modes, supporting multi-grating concatenation)
- **S-flats:** 34 unique files across MOS, FS, and IFU (BOTS reuses FS)
- **Total:** 40 flat-field reference files in the `wavext` CRDS context.

---

## Future Plan: F-flat Concatenation
As identified in [IMPLEMENTATION_PLAN.md](IMPLEMENTATION_PLAN.md), the next critical step for **wavext** is the concatenation of the F-flat `FAST_VARIATION` tables across consecutive gratings. This will allow the pipeline to apply a consistent fore-optics throughput correction across the entire 0.6–5.6 µm range, regardless of the primary grating in use.
