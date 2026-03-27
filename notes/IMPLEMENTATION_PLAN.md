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
- **`extract_1d` Failure (`AttributeError: No attribute 'get_dtype'`)**: Encountered due to a version mismatch where the `jwst_nirspec_wavext` fork was based on an older `jwst` version that called a deprecated `stdatamodels` method.
    - **Fix Applied:** Rebased the `jwst_nirspec_wavext` branch onto the `1.20.2` tag to align the pipeline code with the installed `jwst_1.20.2` environment. `extract_1d` now runs successfully and produces extended 1D spectra up to 5.144 µm for S200A1 (limited physically by the detector edge).
