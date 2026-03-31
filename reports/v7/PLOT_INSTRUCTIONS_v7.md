# 📊 NIRSpec v7 Plotting Instructions

This document provides instructions for generating all diagnostic and validation plots for the NIRSpec v7 report.

## 1. Environment Requirements
All scripts are designed to run with `python3` and require:
- `numpy`, `matplotlib`, `astropy`, `scipy`

## 2. Core Plotting Scripts

### 2.1 Full Spectrum Merged Plots (v5 Style)
Generates the comprehensive 0.6–5.6 µm spectra showing NRS1/2 overlap, ghost components, and CALSPEC truth.
- **Script**: `scripts/plot_v7_full_spectra_v5_style.py`
- **Output**: `fs_v7/plots/full_spectrum_merged_v7_{Source}.png`
    - **Observed Data**: Black dots (`alpha=0.6`) with a faint line (`alpha=0.2`).
    - **PRISM Data**: Green thick line (`lw=5`, `alpha=0.2`).
    - **CALSPEC Model**: Black thick line.
    - **Labels**: Large font sizes (14pt) for all components; "NRS1,2" suffix removed.
    - **PRISM Label**: Green at full alpha.
    - **2nd Order Ghost**: Black pixels (`,`, `alpha=0.5`) + faint dashed line.
    - **3rd Order Ghost**: Black pixels (`,`, `alpha=0.4`) + faint dotted line.
    - **Labels**: "2nd order (λ/2)" and "3rd order (λ/3)".
    - **Y-Axis**: Tightened to `1.4x` CALSPEC maximum.

### 2.2 Per-Source Validation (NRS2 only)
Compares observed NRS2 flux against the ghost model (k*f1 + α*f2 + β*f3).
- **Script**: `scripts/plot_v7_full_spectra.py`
- **Output**: `fs_v7/plots/fs_v7_full_spec_{grating}_{pid}.png`
- **Plotting Style**:
    - Includes individual plots for 2nd order and 3rd order components.
    - Ghost Marker: Black pixels (`,`).
    - Ghost Alphas: 0.5 (2nd), 0.4 (3rd).
    - Labels: "2nd Order (α*f2) λ/2" and "3rd Order (β*f3) λ/3".
    - **Grating Colors**: G235M is `goldenrod`, PRISM is `m` (magenta).

### 2.3 Hold-out Cross-Validation
Validates the generality of the ghost model coefficients by plotting predictions for sources excluded from the fit.
- **Script**: `scripts/plot_v7_holdouts.py`
- **Output**: `fs_v7/plots/fs_v7_holdout_validation_{grating}.png`

### 2.4 v7 vs v6 Comparison
Compares the coefficients (k, α, β) derived in v7 (natively flux-calibrated) against the v6 results.
- **Script**: `scripts/plot_v7_v6_comparison.py`
- **Output**: `fs_v7/plots/fs_v7_v6_comparison_{grating}.png`

## 3. IFU Validation Plots
Similar scripts tailored for the IFU v7 report.
- **Script**: `scripts/plot_ifu_v7_full_spectra_v5_style.py`
- **Script**: `scripts/plot_ifu_v7_validation.py`

## 4. Execution Workflow
To refresh all v7 FS plots, run the following command from the `v7/scripts` directory:

```bash
python3 plot_v7_full_spectra_v5_style.py
python3 plot_v7_full_spectra.py
python3 plot_v7_holdouts.py
python3 plot_v7_v6_comparison.py
```

## 5. Design Decisions (Maintenance)
- **Ghost Components**: Always use high `zorder` (20+) for labels to keep them visible.
- **Scaling**: Use `np.nanmax() * 1.5` for consistent headroom across diverse sources.
- **Transparency**: Layering dots and lines provides a "premium" feel while maintaining data density.
