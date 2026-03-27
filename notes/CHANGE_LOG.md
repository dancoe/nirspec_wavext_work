# CHANGE_LOG.md

Add new entries to the top (most recent first). Don't delete anything.

## 2026-03-27 01:25 AM
- **Wavelength Extension Investigation:** Checked why extended wavelengths are not plotted for G140H and G235H on the `NRS2` detector.
  - Investigated `assign_wcs` bounding box logic in `jwst_nirspec_wavext/jwst/assign_wcs/nirspec.py`. Found that the WCS models the physical detector boundaries, capping extraction at `x=2047`.
  - For `S200A1` on `NRS2`, `G140H` physically terminates at 1.819 µm (well within the nominal 1.89 µm cutoff), and `G235H` terminates at 3.050 µm (within the nominal 3.17 µm cutoff).
  - Thus, there is zero flux beyond the nominal cutoff for these combinations because the light physically falls off the detector edge. Plotted output correctly shows no yellow extended shading because all available data is already within the nominal ranges.
  - Successfully updated and executed the `run_assign_wcs.py`, `run_extract_2d.py`, and `run_extract_1d.py` scripts to verify these limits on G235H data (`jw01492001001_0310a_00007_nrs2_rate.fits`).

## 2026-03-27 12:40 AM
- **Extract1dStep Success:** Successfully executed `AssignWcsStep`, `Extract2dStep`, and `Extract1dStep` for G395H PID 1492 Fixed Slit data using the custom `wavelengthrange_extended.asdf` reference file.
  - **Unit Verification:** Confirmed that while the ASDF reference uses meters (1e-6), the pipeline correctly scales the resulting spectra into **microns**.
  - **Detector Edge Constraints:** Verified that the extraction is limited physically by the NRS2 detector edge (e.g., ~5.144 µm for S200A1).
- **Pipeline Branch Rebase (jwst_1.20.2):** Resolved a critical `AttributeError: No attribute 'get_dtype'` that occurred during `Extract1dStep`.
    - **Root Cause:** The `feature/nirspec_wavelength_extension` branch was based on an outdated `jwst` pipeline codebase that used deprecated `stdatamodels` methods.
    - **Fix:** Rebased the local fork onto the official `1.20.2` tag to align the pipeline with the current environment's `stdatamodels` and `stpipe` dependencies.
- **Repository Setup:** Finalized the separation of `jwst_nirspec_wavext` and `nirspec_wavext_work` with correct SSH key mapping and GitHub remotes.
- **Next Steps:** Proceed with the throughput ($k$, $a$, $b$) function derivation (Parlanti et al. model) by comparing nominal regions with extended extraction results.

## 2026-03-27 12:08AM
- **Extraction Phase:** Starting the extraction and calibration of the PID 1492 FS data.
- **Environment Verification:** Confirming `micromamba activate jwst_1.20.2` and `PYTHONPATH` exports for the custom `jwst_nirspec_wavext` library.

## 2026-03-27 12:07AM
- **Git Remote & Config:** Configured both `nirspec_wavext_work` and `jwst_nirspec_wavext` repositories to use the `id_ed25519_dancoe` SSH key identity via `~/.ssh/config` host mapping.
- **Initial Push:** Successfully pushed the `nirspec_wavext_work` repository to GitHub.
- **Pipeline Push:** Pushed the `feature/nirspec_wavelength_extension` branch of `jwst_nirspec_wavext` to GitHub, including the ASDF lazy-load fix.
- **Git Config:** Updated local `user.name` to `dancoe` and `user.email` to `dancoe@gmail.com` for both repositories.


## 2026-03-27 12:05AM
- **Directory Reorganization:** Separated the project into two main repositories:
  - `jwst_nirspec_wavext`: The local fork of the JWST pipeline.
  - `nirspec_wavext_work`: Contains project `notes`, `analysis`, and `reference_files`.
- **Reference Files:** Created `nirspec_wavext_work/reference_files/` and moved `wavelengthrange_extended.asdf` into it.
- **Extract1dStep Success:** Successfully executed `AssignWcsStep`, `Extract2dStep` and `Extract1dStep` for G395H PID 1492 Fixed Slit data with `wavelengthrange_extended.asdf`.
  - Debugged WCS output to confirm correct bounding box coordinates mapping to 0.6–5.6 µm range.
  - Encountered and fixed an `AttributeError: No attribute 'get_dtype'` during `Extract1dStep`. Solved by rebasing the outdated `jwst_nirspec_wavext` fork onto the `1.20.2` tag to align the pipeline modifications with `stdatamodels` requirements in the `jwst_1.20.2` environment.
  - Verified that `Extract1dStep` successfully extracts spectra in microns up to the physical edge of the NRS2 detector (e.g., 5.144 µm for G395H `S200A1`).
  - Next step: Compare nominal wavelength regions with extended extraction results and develop the throughput ($k$, $a$, $b$) function derivation.
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
