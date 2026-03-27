# NIRSpec Wavelength Extension ("wavext") Investigation & Plan

This plan outlines the strategy for extending the NIRSpec wavelength coverage in the JWST pipeline for FS, MOS, IFU, and BOTS modes.

## Goal
Extract full spectra (0.6 - 5.3 µm) across all NIRSpec modes by removing hard-coded limitations and overriding reference files.

## Proposed Changes

### 1. Code Changes

#### [MODIFY] [nirspec.py](file:///Users/dcoe/NIRSpec/wavext/jwst_nirspec_wavext/jwst/assign_wcs/nirspec.py)
- **IFU Mode:** Remove the hard-coded `NoDataOnDetectorError` for NRS2 M-gratings and G140H/F070LP (lines 283-292).
- **Graceful Fallback:** Implement a data-driven check (e.g., checking if any valid wavelengths fall on the detector) instead of hard-coded grating/filter checks. This ensures processing continues if the user provides an extended wavelength range that *does* fall on NRS2.

### 2. Reference File Changes

We will prioritize **extending reference files** as needed over implementing extrapolation logic in the pipeline code. This ensures the pipeline relies on verified (even if approximate) calibration data.

#### Environment Setup
Refer to [PIPELINE_RUNNING.md](PIPELINE_RUNNING.md) for CRDS server URLs, local cache paths (`~/crds_cache`), and official JWST pipeline notebooks.

#### Wavelength Range Reference File (`wavelengthrange`)
- **Action:** Create a custom ASDF reference file that extends the `wavelengthrange` for all NIRSpec filter/grating combinations.
- **FS/MOS/BOTS:** Update the ranges to reflect the full detectable extent on both NRS1 and NRS2.
- **Reference Values (Adopted from `msaexp`/Parlanti):**
    - `F070LP_G140M/H`: [0.6, 3.3] µm
    - `F100LP_G140M/H`: [0.9, 3.3] µm
    - `F170LP_G235M/H`: [1.5, 5.3] µm
    - `F290LP_G395M/H`: [2.6, 5.6] µm
    - `CLEAR_PRISM`: [0.5, 5.6] µm

#### Flat Field Reference Files (`sflat`, `fflat`)
- **Action:** Extend existing FFLAT and SFLAT reference files to cover the full range (0.6 - 5.3 µm). For extended regions without measured flats, we will implement placeholders (e.g., constant or nearest-neighbor extensions) within the reference file itself so the pipeline's native interpolation can handle it.

### 3. Mode-Specific Extraction Strategies

| Mode | Extraction Strategy | Key Challenges |
| :--- | :--- | :--- |
| **FS** | Relies on `extract_2d` using WCS bounding box. | Ensures all 5 fixed slits are handled correctly when NRS2 is active. |
| **MOS** | Similar to FS; uses `extract_2d` for each slitlet. | Managing overlap on NRS2 and ensuring `extract_2d` doesn't skip NRS2 slits. |
| **IFU** | Uses `assign_wcs` then `cube_build`. | Removing the NRS2 error in `assign_wcs` is the primary blocker. |
| **BOTS** | Handled as `S1600A1` fixed slit. | Ensure TSO-specific logic in `extract_2d` respects the extended range. |

## Verification Plan

### Current State (March 2026)
- **Patch Done:** NRS2 IFU WCS patch applied in [nirspec.py](file:///Users/dcoe/NIRSpec/wavext/jwst_nirspec_wavext/jwst/assign_wcs/nirspec.py).
- **Reference Files:** [wavelengthrange_extended.asdf](file:///Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/reference_files/wavelengthrange_extended.asdf) created (0.6 - 5.6 µm).
- **Environment:** `jwst_1.20.2` with `asdf<5.0` is the required baseline.
- **Blocker:** ~~Encountering `OSError: Attempt to load block from closed file` during `AssignWcsStep`.~~ **RESOLVED.** The root cause was that `get_spectral_order_wrange()` in `nirspec.py` returned lazy-loaded ASDF ndarray proxy objects for `order` and `wrange` *after* calling `wave_range_model.close()`. These proxies attempted to load their data from the already-closed ASDF file when later accessed in `validate_open_slits()`. Fix: force-materialise both values with `int(order)` and `list(wrange)` before `close()` in all code paths.
- **`assign_wcs` verified working** on `jw01492001001_03108_00007_nrs2_rate.fits` (NRS2 FS, G140M) with [wavelengthrange_extended.asdf](file:///Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/reference_files/wavelengthrange_extended.asdf) override.
- **Extraction Verification:** Run `extract_2d` and check the `wavelength` array in the output `MultiSlitModel` or `SlitModel`.
- **Existing Tests:**
    - `pytest jwst/assign_wcs/tests/test_nirspec.py`
    - `pytest jwst/extract_2d/tests/test_nirspec.py`

### Automated Tests
- **WCS Verification:** Run `assign_wcs` on NRS2 data with the override and verify the `wcs` object contains the extended range.

### Manual Verification
- **PID 1492 Data:** Process Fixed Slit data from PID 1492 through `assign_wcs`, `extract_2d`, and `extract_1d` using the custom [wavelengthrange_extended.asdf](file:///Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/reference_files/wavelengthrange_extended.asdf) file.
- **FS/MOS/BOTS:** Fully enabled via the `wavelengthrange` reference file override. These modes rely on `get_spectral_order_wrange` to define the extraction bounding box. By extending the ranges in the reference file, the WCS objects will automatically cover the extended detector regions.
- **IFU:** Requires both the code fix (to allow NRS2 processing) and the `wavelengthrange` override.

### Priority
The current priority is **Fixed Slit (FS)** data. Once FS is calibrated and verified, we will proceed to MOS, IFU, and finally BOTS.

### Technical Hurdles & Investigation
- **ASDF Version Mismatch:** The default environment contained `asdf 5.2.0`, while `jwst 1.20.2` (released late 2025) is more compatible with `asdf 4.x`.
- **`Resolved to relative URL` Error:** Encountered when opening data models in older environments (e.g. `jwst 1.19.1`). Resolved by switching to `jwst 1.20.2` and pinning `asdf==4.5.0`.
    - **`Attempt to load block from closed file` Error:** ~~(Persistent Blocker)~~ **FIXED (2026-03-26).** This occurred when `get_spectral_order_wrange()` returned ASDF lazy-load proxies for `order` and `wrange` after calling `wave_range_model.close()`. The proxy then tried to load data from the closed file when `validate_open_slits()` accessed `wrange[0]`.
        - **Fix Applied:** Force-materialise with `int(order)` and `list(wrange)` in `get_spectral_order_wrange()`, before `wave_range_model.close()`, for all code paths.
        - **Attempted Fixes (Failed, Superseded):**
            - Using `memmap=False`: Fails because `memmap` no longer affects current `asdf` schema loading in these versions.
            - Using `with datamodels.open(...)`: The handle is still closed prematurely within the `AssignWcsStep` internal deep-copy logic.
    - **Future Workaround:** Investigate a "force-load" of all SLIT/IMAGE meta arrays into memory before passing to the step.
### 4. Calibration & Overlap Modeling (Parlanti et al. 2025)

The wavelength extension strategy relies on modeling the spectral overlap between 1st, 2nd, and 3rd order light. Using the model from **Parlanti et al. 2025**, we solve for the throughput functions $k, a, b$ where:
$S(\lambda) \approx k(\lambda) f(\lambda) + a(\lambda) f(\lambda/2) + b(\lambda) f(\lambda/3)$

- **Strategy:** Use PRISM data as the intrinsic baseline $f(\lambda)$.
- **Status (March 2026):**
    - **Medium (M) Gratings (Priority):** G140M/G235M are the primary targets for calibration. These stay on the NRS2 detector long enough to capture significant 2nd-order overlap.
    - **High (H) Gratings (Finding):** Due to high dispersion, wavelengths above ~1.8 µm (G140H) or ~3.0 µm (G235H) physically roll off the detector edge (`x=2047`) for the central fixed slit (`S200A1`). Consequently, H-gratings do not experience the same level of spectral contamination as M-gratings in FS mode.
    - **Verification:** Custom 1D extractions for G140M show successful mapping up to ~3.3 µm on NRS2, compared to the ~1.9 µm nominal cutoff.

### Current Phase (March 27, 2026)
- **Goal:** Download, reduce, analyze, compare, and plot data for G140M, G235M, G395M, and PRISM.
- **Tasks:**
    - ✓ Download all remaining M-grating and PRISM Level 2 (`_rate.fits`) and Level 3 (`_s2d.fits` / `_x1d.fits`) data for PID 1492.
    - Run extended WCS assignment on all Level 2 files.
    - Run extended 2D and 1D extraction.
    - Compare extended extractions with nominal extractions to verify wavelength alignment and overlap.
    - ✓ **Implement Parlanti flux correction model coefficient fitting** (Initial analysis complete)
    - Generate comparison plots showing nominal vs. extended spectra with overlap regions highlighted.

### Parlanti Model Coefficient Determination (March 27, 2026)

**Completed v1:** Initial polynomial fit (cubic order) — numerically unstable.

**Completed v2 (March 27, 2026 evening):** Improved fit using Legendre polynomials + Tikhonov regularisation.

**Key Finding (Units):** Extended NRS2 extractions are in **DN/s** (no photom step applied), while nominal Level 3 products are in **Jy**. This 100–200× discrepancy is bridged using a boundary-matching scale factor derived from comparing the extended spectrum to the PRISM at the nominal–extended transition (where k≈1 and contamination is minimal).

**Data Used:**
- PRISM baseline (obs001 Level 3): 411 points, λ=[0.597, 5.300] µm, Jy
- G140M extended (obs003 NRS2): 1,095 points, λ=[1.981, 3.268] µm, DN/s → scaled to Jy
- G235M extended (obs003 NRS2): 1,272 points, λ=[3.308, 5.310] µm, DN/s → scaled to Jy
- G395M extended (obs003 NRS2): 68 points, λ=[5.496, 5.617] µm — too few, PRISM doesn't extend here

**v2 Fit Results:**

| Grating | λ Range (µm) | RMS | k range | a range | b |
|---------|-------------|-----|---------|---------|---|
| G140M | 1.98–3.27 | 50% | 0.01–0.47 | 0–13% | 5% (clamped) |
| G235M | 3.31–5.31 | 54% | 0.10–0.96 | 5–20% | 0–5% |

**Calibrated Plot:** `plots/Parlanti/cal/FS_1492_cal.png` — 2-panel Parlanti Fig. 5 style

**Script:** `analysis/calibrate_parlanti.py`

**Next Actions:**

- [ ] Run photom step on NRS2 extended extractions to get proper Jy units (remove boundary-scale approximation)
- [ ] Increase b(λ) upper bound to ~10% — the 3rd-order contamination may be higher than 5%
- [ ] Use multiple sources at different redshifts (as Parlanti does) for instrument-level coefficient determination
- [ ] Compare with Parlanti et al. 2025 published values when available

**Outputs:**
- `plots/Parlanti/cal/FS_1492_cal.png` — calibrated spectra comparison plot
- `plots/Parlanti/cal/parlanti_coefficients_v2.txt` — coefficient summary
- `plots/Parlanti/cal/parlanti_coefficients_v2_diagnostic.png` — coefficient profile plots
- `plots/Parlanti/cal/PARLANTI_ANALYSIS_REPORT.md` — v1 analysis report (superseded by v2)
- Notes: `notes/CALIBRATION.md`

---
[IMPLEMENTATION_PLAN.md ends here]
