# NIRSpec Flat Reference Files (Inventory)

This document provides an audited inventory of the `DFLAT`, `FFLAT`, and `SFLAT` reference files used for the **wavext** (wavelength extension) project. 

For the strategy, meaning, and visualization of these flats, see [FLATS.md](FLATS.md).

---

## 2D Spatial Sensitivity Atlas Context
The files listed below quantify the spatial sensitivity variations across the detectors. For visual maps and spectral evolution carousels, refer to the [Atlas in FLATS.md](FLATS.md#2d-spatial-sensitivity-atlas-side-by-side-plots).

---

## Complete Reference File Inventory (CRDS Context: jwst_1464.pmap)

### D-flats (Detector Flats)

`jwst_nirspec_dflat_00##.fits`

Both NRS1 and NRS2 detectors share the **same wavelength sampling** (39 planes total):

**Wavelengths (µm):** 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.8, 2.0, 2.3, 2.6, 2.9, 3.2, 3.5, 3.8, 4.1, 4.4, 4.6, 4.625, 4.65, 4.675, 4.7, 4.725, 4.75, 4.775, 4.8, 4.825, 4.85, 4.875, 4.9, 4.925, 4.95, 4.975, 4.995, 5.0, 5.005, 5.1, 5.2, 5.3

- **Detector NRS1:** `jwst_nirspec_dflat_0001.fits`
- **Detector NRS2:** `jwst_nirspec_dflat_0002.fits`

**Coverage:** **0.6 – 5.3 µm** (39 planes). Fine sampling around 4.6–5.0 µm captures detector responsivity region.

---

### F-flats (Fore-optics Flats)

`jwst_nirspec_fflat_0###.fits`

All observational modes (MOS, FS, IFU, BOTS) share the same fore-optics flat field definitions. Fore-optics throughput is handled coherently across all grating/filter combinations.

| Grating | Filter | NRS1 | NRS2 | Notes |
|---------|--------|------|------|-------|
| **PRISM** | CLEAR | ✓ | ✓ | Lowest resolution; broadest spectral coverage |
| **G140M** | F070LP | ✓ | ✓ | M-grating: medium resolution, extended red |
| **G140M** | F100LP | ✓ | ✓ |  |
| **G140H** | F070LP | ✓ | ✓ | H-grating: high resolution |
| **G140H** | F100LP | ✓ | ✓ |  |
| **G235M** | F170LP | ✓ | ✓ | Mid-IR M-grating |
| **G235H** | F170LP | ✓ | ✓ | Mid-IR H-grating |
| **G395M** | F290LP | ✓ | ✓ | Longest wavelength M-grating |
| **G395H** | F290LP | ✓ | ✓ | Longest wavelength H-grating |

**Total:** 4 F-flat reference files (not all enumerated individually; concatenated across modes).

**The 4 F-flat Files:**

1. **`jwst_nirspec_fflat_0146.fits`** — Mid-to-Long Wavelength (G395M/G395H Region)
   - 4 SCI HDUs with spatial sensitivity maps
   - FAST_VARIATION: 1,356 wavelength points covering **2.848–5.286 µm**
   - Used for IFU mode with long-wavelength grating observations

2. **`jwst_nirspec_fflat_0152.fits`** — Broad Wavelength Coverage (Concatenated Multi-Grating)
   - 4 SCI HDUs covering all channels
   - FAST_VARIATION: 431 wavelength points covering **0.553–5.373 µm** (broadest range)
   - Concatenated F-flat combining G140H + G235M + G395H coverage per Parlanti et al. (2025)
   - Used for IFU and MOS modes requiring unified wavelength range

3. **`jwst_nirspec_fflat_0154.fits`** — Fixed Slit S200A1 Aperture
   - Single SCI HDU (no spatial variations for single slit)
   - FAST_VARIATION: Slit identifier **`S200A1`** (200×200 mas aperture), 1,379 wavelength points covering **0.553–1.840 µm**
   - Aperture-specific fore-optics correction for Fixed Slit mode

4. **`jwst_nirspec_fflat_0173.fits`** — Short-to-Mid Wavelength (G140M/G140H Region)
   - Single SCI HDU
   - FAST_VARIATION: Generic slit (`ANY`), 1,447 wavelength points covering **0.970–1.890 µm**
   - Short-wavelength-optimized fore-optics correction for MOS/IFU with G140M/G140H

**Wavelength Sampling in FAST_VARIATION Table:**
- **G140H + F100LP:** ~431 wavelength points  
  Range: 0.55322–0.80000+ µm (first 20 points: 0.55322, 0.56443, 0.57564, 0.58685, 0.59805, 0.60926, 0.62047, 0.63168, 0.64289, 0.65410, 0.66531, 0.67652, 0.68773, 0.69894, 0.71015, 0.72135, 0.73256, 0.74377, 0.75498, 0.76619 µm)
  
- **G395H + F290LP:** ~1,356 wavelength points  
  Range: 2.84786–5.507+ µm (first 20 points: 2.84786, 2.84966, 2.85146, 2.85325, 2.85505, 2.85685, 2.85865, 2.86045, 2.86225, 2.86405, 2.86585, 2.86765, 2.86945, 2.87125, 2.87305, 2.87484, 2.87664, 2.87844, 2.88024, 2.88204 µm)  
  
  Fine sampling (~0.002 µm steps) in the G395H band supports extended wavelength extrapolation.

**Note:** Parlanti et al. (2025) concatenated M-grating F-flats to build a unified model spanning 0.6–5.6 µm. Red-end regions beyond reference file limits use the reddest available wavelength value (top-hat extension).

---

### S-flats (Spectrograph Flats)

`jwst_nirspec_sflat_0###.fits`

#### MOS S-flats (Multi-Object Spectroscopy)

Grating | Filter | NRS1 | NRS2  
------- | ------ | ---- | ----  
**PRISM** | CLEAR | `0222` | `0233`  
**G140M** | F070LP | `0171` | `0200`  
**G140M** | F100LP | `0231` | `0226`  
**G140H** | F070LP | `0172` | `0199`  
**G140H** | F100LP | `0221` | `0228`  
**G235M** | F170LP | `0229` | `0224`  
**G235H** | F170LP | `0220" | `0230`  
**G395M** | F290LP | `0225` | `0232`  
**G395H** | F290LP | `0227` | `0223`  

MOS S-flat total: **18 reference files** (9 grating/filter combinations × 2 detectors).

**Wavelength Coverage Examples:**

***PRISM + CLEAR (12 planes):*** 0.36393 – 5.55489 µm
- 0.36393, 0.83584, 1.30774, 1.77965, 2.25156, 2.72346, 3.19537, 3.66727, 4.13918, 4.61108, 5.08299, 5.55489

***G140M + F070LP (5 planes):*** 0.70000 – 1.27000 µm
- 0.70000, 0.84250, 0.98500, 1.12750, 1.27000

***G140M + F100LP (12 planes):*** 0.91558 – 1.94275 µm
- 0.91558, 1.00896, 1.10234, 1.19572, 1.28910, 1.38248, 1.47586, 1.56924, 1.66261, 1.75599, 1.84937, 1.94275

***G235M + F170LP (12 planes):*** 1.57039 – 3.25682 µm
- 1.57039, 1.72370, 1.87701, 2.03032, 2.18364, 2.33695, 2.49026, 2.64357, 2.79689, 2.95020, 3.10351, 3.25682

***G395M + F290LP (12 planes):*** 2.72476 – 5.44049 µm
- 2.72476, 2.97165, 3.21853, 3.46542, 3.71230, 3.95919, 4.20607, 4.45296, 4.69984, 4.94672, 5.19361, 5.44049

---

#### FS S-flats (Fixed Slit)

Grating | Filter | NRS1 | NRS2
------- | ------ | ---- | ----
**PRISM** | CLEAR | `0153` | `0155`
**G140M** | F070LP | `0161` | `0150`
**G140M** | F100LP | `0147` | `0154`
**G140H** | F070LP | `0160` | `0151`
**G140H** | F100LP | `0146` | `0152`
**G235M** | F170LP | `0149` | `0159`
**G235H** | F170LP | `0162` | `0156`
**G395M** | F290LP | `0158` | `0148`
**G395H** | F290LP | `0157` | `0145`

FS S-flat total: **18 reference files** (9 grating/filter combinations × 2 detectors).

**Wavelength Coverage Examples:**

***PRISM + CLEAR (12 planes):*** 0.36393 – 5.55489 µm (same as MOS PRISM)
- 0.36393, 0.83584, 1.30774, 1.77965, 2.25156, 2.72346, 3.19537, 3.66727, 4.13918, 4.61108, 5.08299, 5.55489

***G140M + F100LP (12 planes):*** 0.91558 – 1.94275 µm (same as MOS G140M F100LP)
- 0.91558, 1.00896, 1.10234, 1.19572, 1.28910, 1.38248, 1.47586, 1.56924, 1.66261, 1.75599, 1.84937, 1.94275

**Note:** FS S-flats typically have narrower spatial footprints but identical wavelength sampling to MOS equivalents.

---

#### IFU S-flats (Integral Field Spectroscopy)

Grating | Filter | NRS1 | NRS2
------- | ------ | ---- | ----
**PRISM** | CLEAR | `0214` | `0196`
**G140M** | F070LP | `0210` | `0184`
**G140M** | F100LP | `0208` | `0191`
**G140H** | F070LP | `0218` | `0188`
**G140H** | F100LP | `0216` | `0212`
**G235M** | F170LP | `0211` | `0192`
**G235H** | F170LP | `0217` | `0213`
**G395M** | F290LP | `0209` | `0195`
**G395H** | F290LP | `0219` | `0215`

IFU S-flat total: **18 reference files** (9 grating/filter combinations × 2 detectors).

**Wavelength Coverage Examples:**

***PRISM + CLEAR (12 planes):*** 0.36393 – 5.55489 µm (same as MOS/FS PRISM)
- 0.36393, 0.83584, 1.30774, 1.77965, 2.25156, 2.72346, 3.19537, 3.66727, 4.13918, 4.61108, 5.08299, 5.55489

***G140M + F100LP (12 planes):*** 0.91558 – 1.94275 µm (same as MOS/FS G140M F100LP)
- 0.91558, 1.00896, 1.10234, 1.19572, 1.28910, 1.38248, 1.47586, 1.56924, 1.66261, 1.75599, 1.84937, 1.94275

**Note:** IFU S-flats use identical wavelength sampling to MOS and FS, but cover the full 2D IFU field of view.

---

### Summary Statistics

| Component | Type | Total Count | Notes |
|-----------|------|-----------|-------|
| **D-flats** | Detector | 2 | One per detector (NRS1, NRS2) |
| **F-flats** | Fore-optics | 4 | Shared across modes |
| **S-flats (MOS)** | Spectrograph | 18 | 9 grating/filter combos × 2 detectors |
| **S-flats (FS)** | Spectrograph | 18 | 9 grating/filter combos × 2 detectors |
| **S-flats (IFU)** | Spectrograph | 18 | 9 grating/filter combos × 2 detectors |
| **S-flats (BOTS)** | Spectrograph | *(same as FS)* | Bright object time-series uses FS S-flats |
| **Total Unique S-flats** | — | **34** | Across MOS, FS, IFU (BOTS reuses FS) |
| **All Flat Files** | — | **40** | 2 D + 4 F + 34 S |
