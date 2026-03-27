# Instructions for NIRSpec Flat-Field Extension Plots (wavext)

This document describes how to generate and interpret the flat-field extension plots for the NIRSpec wavelength extension project.

## Requirements
- Python 3.x
- `astropy`, `matplotlib`, `numpy`
- Local `crds_cache` containing the relevant reference files.

## Plotting Scripts

| Script | Purpose | Output Location |
| :--- | :--- | :--- |
| `analysis/plot_flats_batch.py` | Generate 1D S & F-flat table extension plots for a single mode. | `plots/flats/<mode>_<type>_<config>.png` |
| `analysis/plot_flats_all_modes.py` | Batch generate 1D extension plots for FS, IFU, and MOS modes. | `plots/flats/*.png` |
| `analysis/plot_dual_detector_flats.py` | Generate 2D side-by-side NRS1/NRS2 spatial maps (SCI images). | `plots/flats/dual_detector_*.png` |
| `analysis/plot_parlanti_reproduction.py` | Create two-panel conceptual reproductions (top: Fast Var, bottom: Slow Var). | `plots/flats/parlanti_repro_*.png` |

---

## 1D Spectral Extension Plots (Fast Variation)
The 1D plots show the `FAST_VARIATION` table data from the S-flat or F-flat.
- **X-axis**: Wavelength in microns (µm).
- **Y-axis**: Relative throughput value.
- **Blue Line**: Nominal calibration range (available in current CRDS).
- **Red Dashed Line**: Injected `wavext` extension (holding the reddest available value constant).
- **Highlight (Orange/Yellow)**: Region where extrapolation is performed (redwards of the nominal limit).

**Usage**:
```bash
python3 analysis/plot_flats_all_modes.py
```

---

## 2D Spatial Sensitivity Plots (Slow Variation)
The 2D plots show the `SCI` image extensions for NRS1 and NRS2 side-by-side. 
- **Scale**: Linear scale from **0.9 to 1.1** (to capture typical spatial variations and artifacts).
- **NRS1/NRS2 Mapping**: NRS1 is on the left, NRS2 on the right.
- **Pixel Range**: Full detector frame (2048 x 2048).
- **Interpretation**: In the `wavext` strategy, the extended planes (at 5.61 µm) are forced to 1.0 (unity), meaning they will appear as dark flat fields in a 0.9–1.1 color scale if the data isn't present, or perfectly uniform if initialized to 1.0.

**Usage**:
```bash
python3 analysis/plot_dual_detector_flats.py
```
> [!TIP]
> To match the Parlanti Fig A1 format exactly, ensure the scale is locked to [0.9, 1.1] and use a perceptual colormap like `viridis`.

---

## Highlighting Mode-Specific Extensions
When adding new modes, ensure the following keywords are correctly identified in the header to find the correct reference files:
- `EXP_TYPE` (e.g., `NRS_FIXEDSLIT`, `NRS_MSASPEC`, `NRS_IFU`)
- `GRATING` and `FILTER`
- `DETECTOR` (`NRS1` or `NRS2`)

For IFU modes, F-flats are ideally concatenated across M-gratings to provide the full extended baseline.
