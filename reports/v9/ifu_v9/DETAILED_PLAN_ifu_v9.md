# Detailed Analysis Plan v9 — IFU Uniform 0.5" Extraction

**Version:** 9.0  
**Status:** Planning  
**Objective:** Standardize the IFU calibration validation using a uniform $r=0.5"$ fixed-radius aperture for all sources, ensuring consistency across standard stars and science targets.

---

## 1. Target Selection & Data Inventory
Following specific v9 requirements, we are excluding sources with nearby contamination:
- **Included Standard Stars:** 1536 (J1743045), 1537 (G191-B2B), 1538 (P330E), 6645 (P330E).
- **Included Science Target:** 2186 (UGC-5101).
- **Excluded:** 2654 (SDSSJ0841) — due to nearby object contamination.

---

## 2. Methodology: Uniform Aperture Extraction
- **Aperture Tweak:** Use $r=0.5"$ circular aperture for **ALL** observations.
- **Centering:** Centered on the brightest pixel (peak-centered) to capture maximum flux regardless of WCS offsets.
- **Background:** No background subtraction applied, consistent with the v8 custom extraction baseline.
- **Unit Conversion:** Full MJy/sr to Jy conversion using `PIXAR_SR` from `spec3` headers.

---

## 3. Visualization & Reporting
Generate a full suite of diagnostics in `reports/v9/ifu_f9/`:
- **Full Spectrum Merged Plots**: Compare the $r=0.5"$ extraction to reference (CALSPEC for stars).
- **Spatial Slices**: Document the centering and aperture footprint for all 5 targets.
- **Standalone Slice Report**: Create `REPORT_ifu_slices_v9.md` for spatial-only verification.
- **Comprehensive Validation**: Update `REPORT_ifu_v9.md` with integrated spectra and ratios.

---

## 4. Verification Checklist
- [ ] 1536, 1537, 1538, 2186, 6645 are all processed.
- [ ] $r=0.5"$ aperture is consistently used across all targets.
- [ ] PID 2654 is successfully excluded from the suite.
- [ ] Gaps at ~1.3-1.45 µm and ~2.35-2.55 µm are correctly handled in plotting.

---

*Created: 2026-03-31*  
*Author: Antigravity*
