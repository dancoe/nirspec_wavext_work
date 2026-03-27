# CHANGE_LOG.md

## 2026-03-27 12:40 AM
- **Extract1dStep Success:** Successfully executed `AssignWcsStep`, `Extract2dStep`, and `Extract1dStep` for G395H PID 1492 Fixed Slit data using the custom `wavelengthrange_extended.asdf` reference file.
  - **Unit Verification:** Confirmed that while the ASDF reference uses meters (1e-6), the pipeline correctly scales the resulting spectra into **microns**.
  - **Detector Edge Constraints:** Verified that the extraction is limited physically by the NRS2 detector edge (e.g., ~5.144 µm for S200A1).
- **Pipeline Branch Rebase (jwst_1.20.2):** Resolved a critical `AttributeError: No attribute 'get_dtype'` that occurred during `Extract1dStep`.
    - **Root Cause:** The `feature/nirspec_wavelength_extension` branch was based on an outdated `jwst` pipeline codebase that used deprecated `stdatamodels` methods.
    - **Fix:** Rebased the local fork onto the official `1.20.2` tag to align the pipeline with the current environment's `stdatamodels` and `stpipe` dependencies.
- **Repository Setup:** Finalized the separation of `jwst_nirspec_wavext` and `nirspec_wavext_work` with correct SSH key mapping and GitHub remotes.
- **Next Steps:** Proceed with the throughput ($k$, $a$, $b$) function derivation (Parlanti et al. model) by comparing nominal regions with extended extraction results.


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
