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

### Environment:
- Python 3.11 / JWST 1.20.2
- `PYTHONPATH`: `/Users/dcoe/NIRSpec/wavext/jwst_nirspec_wavext` (patched)
- CRDS Context: `jwst_1464.pmap`

### Open Issues:
- None for v8. Final validation plots look excellent across all grating configurations.
- Gaps are intentional and reflect the detector physical geometry.
