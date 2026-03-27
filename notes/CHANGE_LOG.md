# CHANGE_LOG.md

## 2026-03-27 12:05AM
- **Directory Reorganization:** Separated the project into two main repositories:
  - `jwst_nirspec_wavext`: The local fork of the JWST pipeline.
  - `nirspec_wavext_work`: Contains project `notes`, `analysis`, and `reference_files`.
- **Reference Files:** Created `nirspec_wavext_work/reference_files/` and moved `wavelengthrange_extended.asdf` into it.
- **Documentation:** Updated all absolute and relative file paths in the `notes/` directory to reflect the new repository structure.

## 2026-03-26 11:56PM
- **Bug Fix:** [nirspec.py](file:///Users/dcoe/NIRSpec/wavext/jwst_nirspec_wavext/jwst/assign_wcs/nirspec.py) — `get_spectral_order_wrange()` (lines ~1231–1260)
    - **Root Cause:** `wave_range_model.wavelengthrange[index]` and `wave_range_model.order[index]` return lazy-loaded ASDF ndarray proxies that back-reference the open ASDF file handle. After `wave_range_model.close()` is called, any later attempt to access those proxies (e.g. `wrange[0]` at `validate_open_slits`) raises `OSError: Attempt to load block from closed file`.
    - **Fix:** Force-materialise both values into plain Python types (`int(order)`, `list(wrange)`) **before** `wave_range_model.close()` in all code paths. This resolves the blocker for FS (and all other modes) when using the custom `wavelengthrange_extended.asdf` override.
    - **Verified:** `run_assign_wcs.py` now completes successfully on `jw01492001001_03108_00007_nrs2_rate.fits` and writes the WCS output.

## 2026-03-26 11:20PM
- **Modified Code:** [nirspec.py](file:///Users/dcoe/NIRSpec/wavext/jwst_nirspec_wavext/jwst/assign_wcs/nirspec.py)
    - Removed hard-coded `NoDataOnDetectorError` for IFU modo processing on NRS2 for M-gratings and G140H/F070LP. This allows the pipeline to proceed with wavelength extension for these combinations.
- **Reference Files:** 
    - Generated `wavelengthrange_extended.asdf` using `stdatamodels`. This file extends the `wavelengthrange` for all NIRSpec gratings (0.6 - 5.6 µm range).
- **Modified Documentation:** Updated [IMPLEMENTATION_PLAN.md](IMPLEMENTATION_PLAN.md) to emphasize reference file extensions over code extrapolation.
