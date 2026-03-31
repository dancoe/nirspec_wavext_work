# NIRSpec FS v7 Calibration — Execution Log

## Step 1: Flux Calibration Strategy
- **Goal**: Enable native **Jy** units for NRS2 extended extractions directly within the `Spec2Pipeline`.
- **Action**: Created custom `photom` reference file with extended `relresponse` arrays (equal to 1.0 beyond nominal cutoff) using `jwst_nirspec_photom_0014.fits` as the source to maintain `slit` column support.
- **Script**: `analysis/create_photom_fs_ext.py` (v7-robust version uses 500-element column buffers).

## Step 2: Pipeline Reprocessing (v7)
- **Goal**: Process all Fixed Slit (FS) datasets through Level 2 Spec2 with custom overrides.
- **Targets**: PIDs 1492, 1536, 1537, 1538, 6644.
- **Overrides**:
  - `assign_wcs`: `wavelengthrange_extended.asdf`
  - `flat_field`: Parlanti-extended `fflat` and `sflat`
  - `photom`: `jwst_nirspec_photom_fs_nrs2_ext.fits` (v7)
- **Status**: Completed for all targets. Verified `_x1d.fits` files now have `TUNIT2 = Jy` and valid flux coverage for > 2.0 um.

## Step 3: Solver Implementation (v7)
- **Goal**: Perform joint NNLS solve for $k, \alpha, \beta$ and assess stability via leave-one-out cross-validation (CV).
- **Solver**: `reports/v7/scripts/solve_parlanti_fs_v7.py`.
- **Inputs**: Natively flux-calibrated Jy spectra from `nrs2_spec2_cal/` directories. No manual scaling applied.
- **Sources**: G191-B2B, P330E, J1743045, NGC2506-G31, and PID 1492 (stitched truth).
- **Status**: Successfully generated `calib_v7_fs_{grating}_all.fits` and per-source hold-out fits.

## Step 4: Validation & Plotting
- **v7 vs v6 Comparison**: Verified that the transition to native Jy-units produces the same coefficients as the v6 manual-scaling approach.
- **Hold-out Cross-Validation**: Generated plots for each source, showing the predicted observed spectrum using coefficients derived from all other sources.
- **Result**: High generality across sources; model effectively reconstructs observed spectra despite source exclusion.
