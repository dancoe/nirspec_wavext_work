# 📔 NIRSpec wavext v8 — Detailed Analysis Log

## 2026-03-31: v8 Calibration Characterization & Extraction

### Actions Performed:
1.  **IFU Custom Extraction (PID 2186/2654)**: 
    - Implemented a 0.5" radius circular aperture extraction centered on the brightest pixel.
    - Verified `PIXAR_SR` and `BUNIT=MJy/sr` handling to ensure results are in Jy.
    - Compared results to default `extract_1d` pipeline products:
        - UGC-5101 (2186): Ratio ~0.75
        - SDSSJ0841 (2654): Ratio ~0.89
    - Documented cases in `EXTRACTED_ifu_v8.md`.

2.  **IFU v8 Baseline Comparisons**:
    - Script: `extractions_all_ifu_v8.py`.
    - Implemented **Ratio Subplots**: Visualized all custom extractions as ratios of the pipeline `extract1d` spectrum to pinpoint flux offsets.
        - Ratio y-axis: Log-scale [0.5, 2.0].
    - **Spatial Diagnostic Overlays**: Added `slices_*.png` for all IFU targets.
        - Overlaid circular extraction apertures (0.5" and 0.45") on 2D spectral cube slices.
        - Confirmed peak-pixel vs. nominal-pointing offsets visually.
    - Updated report hierarchy in `REPORT_all_extractions.md` with object descriptions and clear sectioning.

3.  **Visualization Finalization**:
    - Tightened top-panel spectral plot scaling (percentile-based) for better clarity in the overlap and gap regions.
    - Merged targets from PIDs 1536, 1537, 1538, 2186, 2654, 6645.
    - Added brief target labels (e.g. ULIRG, Solar Analog, White Dwarf) inline in final reports.

4.  **Documentation Updates**:
    - Updated `PARLANTI_IFU.md` with explicit details on the 0.5" aperture extraction methodologies.
    - Confirmed consistency between FS and IFU mode processing in v8.

### 2026-03-31: Parlanti Coefficient Reporting & Diagnostics
1.  **Coefficient Visualization Suite**:
    - Developed `plot_v8_parlanti_coefficients.py` to generate a standardized reporting suite for the ghost contamination model.
    - **Visualization Modes**:
        - **This work (v8)**: Clean results from the unified v8 solver.
        - **v8 vs. Literature Overlay**: Direct comparison using solid v8 lines and faint dashed Parlanti et al. (2025) lines.
        - **Breakout Plots**: Separate comparison panels for each coefficient ($k$, $\tilde{\alpha}$, $\tilde{\beta}$) to allow granular inspection of offsets.
    - **Advanced Styling**:
        - Log-scale $y$-axis capped at 10.0 to focus on the dynamic range of the $\tilde{\alpha}$ and $\tilde{\beta}$ coefficients.
        - Added a thin, faint black reference line at $y=1$ to highlight neutrality.
        - Applied `fill_between` (alpha=0.08) beneath literature curves in individual plots to distinguish the reference baseline.
        - Color-coded by variable: Violet ($k$), Orange ($\tilde{\alpha}$), Blue ($\tilde{\beta}$).
    - **X-axis Standard**: 0.85 – 3.7 µm (G140M), 1.5 – 5.5 µm (G235M) to encompass the full spectral reach of the model.

2.  **Organization & Report Integration**:
    - Centralized all coefficient-related figures in `reports/v8/parlanti/` to serve both FS and IFU reports.
    - Updated `REPORT_fs_v8.md` and `REPORT_ifu_v8.md` with structured sections:
        - Unified coefficient overview.
        - Individual variable breakout comparisons.
        - Standalone literature reference.
    - Applied the "This work (v8)" labeling convention across all diagnostics for clarity.

3.  **Gap Handling & Detector Geometry**:
    - Evaluated wavelength gaps in the full-spectrum merged plots.
    - **Gap Wavelengths**: NRS1/NRS2 gaps are approximately 1.3–1.45 µm (G140M) and 2.35–2.55 µm (G235M) for centered FS targets.
    - **Spatio-Spectral Variability**: Documented that while gap wavelengths are stable for fixed slits, they shift smoothly across IFU slices and vary by slit ID (e.g., S200A1 vs. S200A2).
    - Preservation of gaps was maintained in v8 validation to ensure "honest" representation of raw detector coverage.

### Environment:
- Python 3.11 / JWST 1.20.2
- `PYTHONPATH`: `/Users/dcoe/NIRSpec/wavext/jwst_nirspec_wavext` (patched)
- CRDS Context: `jwst_1464.pmap`

### Open Issues:
- None for v8. Final validation plots and coefficient derivations look excellent across all configurations.
