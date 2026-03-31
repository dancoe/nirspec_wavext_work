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

2.  **IFU v8 Analysis**:
    - Generated merged full-spectrum (0.6–5.6 µm) plots for all 6 IFU targets (1536, 1537, 1538, 2186, 2654, 6645).
    - Preserved physical gaps at 2.2 µm and 3.7 µm by suppressing interpolation in plotting.
    - Final report: `REPORT_ifu_v8.md`.

3.  **FS v8 Analysis**:
    - Generated merged full-spectrum plots for FS targets (1536, 1537, 1538, 6644, 1492).
    - Consistent gap handling with IFU.
    - Final report: `REPORT_fs_v8.md`.

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
