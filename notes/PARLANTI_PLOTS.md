# NIRSpec Spectral Overlap Visualization Guide (Parlanti et al. 2025)

This document documents the logic and color schemes used for the spectral overlap plots.

### Color Scheme
- **G140M**: Blue (Nominal) / Cornflower Blue (Extended)
- **G235M**: Dark Goldenrod (Nominal) / Goldenrod (Extended)
- **G395M**: Red (Nominal) / Light Coral (Extended)
- **PRISM**: Black (Reference f(λ)) with alpha=0.3
- **Ghost Reference (f(λ/2))**: Green (#2ecc71), Dashed
- **Ghost Reference (f(λ/3))**: Magenta (#e056fd), Dashed
- **Ghost Component (α·f(λ/2))**: Solid Green (#2ecc71)
- **Ghost Component (β·f(λ/3))**: Solid Magenta (#e056fd)
- **Model Sum**: Red Solid

### Scripts
1. `analysis/plot_combined_gratings.py`: Generates the multicomponent log-scale plot (`FS-1492_pre-cal.png`) with robust outlier rejection, NaN/non-positive masking, and **front-layered PRISM** baseline.
2. `analysis/process_selected_gratings.py`: Handles batch processing and per-grating comparison plots.
3. `plots/Parlanti_gratings.md`: Detailed explainer for the combined overlay results.

### Plot Features
- **Log Scale**: Better visualization of lower-flux spectral extensions.
- **Robust Limits**: 0.1%–99.9% percentile clipping to capture full dynamic range while excluding extreme outliers/spikes.
- **Layering**: PRISM baseline is plotted on top (`zorder=20`, `alpha=0.6`) for maximum visibility against the gratings.
- **Micro-Styling**: Removed internal grid lines for a cleaner, premium presentation.

### How to Reproduce
Run the following script to generate the combined summary plot in `plots/`:
```bash
python3 analysis/plot_combined_gratings.py
```
This script requires the nominal `x1d` files and the extended `extract_1d.fits` (from `process_selected_gratings.py`) in `data/PID1492/`.
