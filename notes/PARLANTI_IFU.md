# PARLANTI_IFU.md

## Overview

This document records the IFU-mode analysis branch of the NIRSpec wavelength-extension
project, initiated 2026-03-27. It parallels the existing Fixed Slit (FS) work but applies
the Parlanti et al. (2025) methodology directly to **NIRSpec IFU (NRS_IFU)** calibration
star data.

**Reference:** Parlanti et al. 2025, arXiv:2512.14844  
*"Doubling NIRSpec/IFS capability to calibrate the single epoch black hole mass relation at high redshift"*  
Code repo: https://github.com/eleonoraparlanti/nirspecIFU-extended

---

## Goal

Reproduce and extend the Parlanti et al. calibration of the IFU wavelength extension
using the same calibration programs they used (PIDs 1537 and 1538), running the reduction
with our own pipeline version (1.20.2) and deriving the k(λ), α̃(λ), β̃(λ) curves
independently for comparison.

---

## Targets

| PID  | Target    | Obs ID | Gratings                              | Notes                                    |
|------|-----------|--------|---------------------------------------|------------------------------------------|
| 1536 | J1743045  | 003    | G140M/F100LP, G235M/F170LP, G395M/F290LP | Standard star — **3rd source**        |
| 1537 | G191-B2B  | 008    | G140M/F100LP, G235M/F170LP            | White dwarf standard star                |
| 1538 | P330E     | 062    | G140M/F100LP, G235M/F170LP            | Gray star (CALSPEC standard) Cycle 1     |
| 6645 | P330-E    | —      | G140M/F100LP, G235M/F170LP            | Same gray star, Cycle 3 (**4th source**) |

These are the same four programs Parlanti et al. used as their primary calibration sources.
With calibration stars whose CALSPEC SEDs are well known, each PID provides one independent
equation per wavelength. Three sources are the minimum needed to solve for k(λ), α̃(λ), β̃(λ)
simultaneously; the fourth (PID 6645) provides redundancy and residual error estimation.

---

## Gratings of Interest

| Grating / Filter  | Nominal Range  | Extended Range | Gap (NRS1/NRS2 gap) | Requires for f(λ)   |
|-------------------|---------------|----------------|---------------------|---------------------|
| G140M / F100LP    | 0.97–1.88 µm  | 0.97–3.55 µm   | 2.17–2.28 µm        | CALSPEC or G235M    |
| G235M / F170LP    | 1.70–3.15 µm  | 1.66–5.27 µm   | 3.66–3.82 µm        | CALSPEC or G395M    |
| G395M / F290LP    | 2.87–5.27 µm  | (not primary)  | —                   | Provides f(λ) ref   |

For standard stars with CALSPEC models, the intrinsic spectrum f(λ) is known at all wavelengths,
so no additional grating observations are required. G395M is still downloaded for PID 1536
(which has it) to allow cross-checks using the overlap region.

---

## Data Acquisition

### Discovery

Queried MAST with `astroquery.mast` for IFU-mode observations from all four calibration programs.

- PID 1536 has 18 IFU observations; **Obs 003** = J1743045 with G140M + G235M + G395M (24 rate files)
- PID 1537 has 18 IFU observations total; **Obs 008** = G140M/F100LP + G235M/F170LP
- PID 1538 has 23 IFU observations total; **Obs 062** = G140M/F100LP + G235M/F170LP
- PID 6645 has 90 IFU observations; P330-E with G140M + G235M (multiple visits)

Query script: `analysis/acquisition/query_ifu_mast.py`

### Download

Downloaded calibration Level-2 rate files (`_rate.fits`), skipping TA (target
acquisition) images (< 50 MB).

Download script: `analysis/acquisition/download_ifu_rates.py`

#### Issues Encountered

**URI zero-padding bug:** The MAST API returns program IDs zero-padded to 5 digits
in URIs (`jw01538…`) while the program metadata returns unpadded (`1538`). The fix:
do not filter URIs by PID string — the observation-level query already restricts to
the correct program.

**`MaskedColumn` iteration:** MAST `productFilename` is a `MaskedColumn`; must call
`str(uri)` to iterate correctly, avoiding type errors.

**Truncated downloads (FITS header stubs):** The first download run produced 64KB
stub files for most G235M files (PID 1537) and all PID 1538 files. These stubs are
exactly one FITS block (65536 bytes = 128 records × 512 bytes) — the MAST connection
timed out after writing the FITS primary header but before the data extensions. The
symptom in pipeline processing: `KeyError: "Extension 'ASDF' not found."` because the
ASDF metadata extension was not written. Fixed by:
1. Adding a post-download size check in `download_ifu_rates.py`: verify downloaded
   file is > 50MB before accepting it
2. Deleting any existing stub before re-downloading
3. Added 3-attempt retry loop per file

The root cause was sequential `Observations.download_file()` calls timing out on
sequential large (84MB) files. The fix ensures partial downloads are detected and retried.

### Downloaded Files

**PID 1537 — G191-B2B** (`data/IFU/PID1537_G191-B2B/stage1/`): 16 files
- G140M/F100LP: 4 dithers × NRS1+NRS2 = 8 files
- G235M/F170LP: 4 dithers × NRS1+NRS2 = 8 files

**PID 1538 — P330E** (`data/IFU/PID1538_P330E/stage1/`): 16 files
- G140M/F100LP: 4 dithers × NRS1+NRS2 = 8 files
- G235M/F170LP: 4 dithers × NRS1+NRS2 = 8 files

**PID 1536 — J1743045** (`data/IFU/PID1536_J1743045/stage1/`): ~24 files
- G140M/F100LP: 4 dithers × NRS1+NRS2 = 8 files
- G235M/F170LP: 4 dithers × NRS1+NRS2 = 8 files
- G395M/F290LP: 4 dithers × NRS1+NRS2 = 8 files

**PID 6645 — P330-E Cycle 3** (`data/IFU/PID6645_P330E-C3/stage1/`): TBD
- G140M/F100LP + G235M/F170LP (multiple visits)

### Header Verification

Spot-checked FITS headers confirmed:
- `EXP_TYPE = NRS_IFU`  (not FS mode)
- `FILTER = F100LP` (G140M set) or `F170LP` (G235M set)
- `GRATING = G140M` or `G235M`
- `DETECTOR = NRS1` or `NRS2`

---

## Directory Layout

```
data/IFU/
  PID1536_J1743045/
    stage1/         # downloaded _rate.fits (24 files: G140M+G235M+G395M)
    stage2/         # Spec2Pipeline output: _cal.fits
    stage3/         # Spec3Pipeline output: _s3d.fits, _x1d.fits
    asn/
  PID1537_G191-B2B/
    stage1/         # downloaded _rate.fits (16 files: G140M+G235M) ✅
    stage2/         # _cal.fits (8 NRS1 files) ✅
    stage3/         # f100lp_g140m_s3d/x1d + f170lp_g235m_s3d/x1d ✅
    asn/
  PID1538_P330E/
    stage1/         # downloaded _rate.fits (16 files) ✅
    stage2/         # _cal.fits (8 NRS1 files) ✅
    stage3/         # f100lp_g140m_s3d/x1d + f170lp_g235m_s3d/x1d ✅
    asn/
  PID6645_P330E-C3/
    stage1/         # downloaded _rate.fits (TBD)
    stage2/
    stage3/
    asn/
```

Note: Stage 1 (Detector1Pipeline: `_uncal.fits` → `_rate.fits`) was **skipped**
because MAST provides pre-calibrated rate files. Stage 1 would only need to be
re-run if we needed custom jump detection settings.

---

## Pipeline Reduction

### Environment

- JWST pipeline version: **1.20.2** (micromamba env `jwst_1.20.2`)
- Build: 12.1.1
- CRDS context: `jwst_1464.pmap` (auto-selected at runtime)
- `CRDS_PATH`: `~/crds_cache`
- `CRDS_SERVER_URL`: `https://jwst-crds.stsci.edu`
- `PYTHONPATH`: `/Users/dcoe/NIRSpec/wavext/jwst_nirspec_wavext` (forked pipeline)

Reference notebook: [JWPipeNB-NIRSpec-IFU.ipynb](https://github.com/spacetelescope/jwst-pipeline-notebooks/blob/main/notebooks/NIRSPEC/IFU/JWPipeNB-NIRSpec-IFU.ipynb) (pipeline v1.20.2, Build 12.1.1)

### Scripts

| Script | Purpose |
|--------|---------|
| `analysis/acquisition/query_ifu_mast.py` | Discover IFU observations on MAST |
| `analysis/acquisition/download_ifu_rates.py` | Download rate files for PIDs 1536, 1537, 1538, 6645 |
| `analysis/reduction/reduce_ifu.py` | Full 3-stage IFU pipeline (original comprehensive script) |
| `analysis/reduction/run_ifu_pipeline.py` | Streamlined Stage 2+3 execution for all 4 PIDs |

### Stage 2 — Spec2Pipeline (`calwebb_spec2`)

Runs individually for each rate file. Produces a `_cal.fits` per exposure.

**ASN construction (`writel2asn`):**
- Creates a Level 2 ASN with the science file as primary member
- Other same-grating/detector exposures added as `selfcal` members (for `badpix_selfcal`)
- If the dither pattern includes NOD, non-current dithers are added as `background` members
- Steps **disabled**: `bkg_subtract`, `cube_build`, `extract_1d`
  (cube building done in Stage 3 where all dithers are combined)

**CRDS reference files fetched for first exposure (PID 1537, G140M/F100LP):**
- `area`, `fflat` (nirspec_fflat_0173.fits), `sflat` (nirspec_sflat_0208.fits)
- `photom` (nirspec_photom_0016.fits), `pathloss` (nirspec_pathloss_0003.fits)
- IFU geometry: `ifufore_0015.asdf`, `ifupost_0016.asdf`, `ifuslicer_0015.asdf`
- `wavelengthrange_0008.asdf`
- `msa_0011.asdf`, `msaoper_0011.json`, `camera_0010.asdf`, `collimator_0010.asdf`
- `fore_0064.asdf`, `fpa_0011.asdf`, `disperser_0079.asdf`

**Steps observed running:**
1. `assign_wcs` — builds the NIRSPEC nrs_ifu WCS pipeline; applies barycentric velocity correction
2. `msa_flagging` — flags 23 failed-open shutters
3. `badpix_selfcal` — skipped (no imprint/background)
4. `nsclean` — skipped
5. `imprint_subtract` — skipped (no L/R microshutter array)
6. `srctype`, `wavecorr`, `flat_field`, `pathloss`, `photom`, `pixel_replace`

### Stage 3 — Spec3Pipeline (`calwebb_spec3`)

Combines all dithered `_cal.fits` files for a given grating/filter into 3D data cubes
and extracts 1D spectra.

**ASN construction (`writel3asn`):**
- Groups `_cal.fits` files by `(FILTER, GRATING)` header keywords
- Produces one ASN per group, named e.g. `f100lp_g140m_l3asn.json`

**Output products:**
- `*_s3d.fits` — 3D IFU data cube (RA, Dec, wavelength)
- `*_x1d.fits` — 1D extracted spectrum (aperture-integrated over the PSF)

**Steps configured:**
- `master_background`: skip
- `outlier_detection`: kernel_size `3 3`
- `cube_build`, `extract_1d`: run (default settings)

---

## Relation to Parlanti et al. Calibration

The Parlanti method requires three pieces for each standard star:

1. **Extended spectrum** $S_\lambda(\lambda)$: from the modified pipeline (modified S-flat
   and F-flat reference files that allow extraction beyond the nominal range)
2. **Nominal spectrum** $f_\lambda(\lambda)$: from the standard pipeline (or from the
   nominal-range portion of the same data cube)
3. **Known intrinsic SED**: for standard stars, a CALSPEC model spectrum provides
   $f_\lambda(\lambda)$ at all wavelengths

The calibration solves:
$$R(\lambda) = k(\lambda) + \tilde\alpha(\lambda)\,r_2(\lambda) + \tilde\beta(\lambda)\,r_3(\lambda)$$

at each wavelength λ > λ_max_nom, where $R = S_\lambda / f_\lambda$ is measured from the data
and $r_2 = f_\lambda(\lambda/2)/f_\lambda(\lambda)$, $r_3 = f_\lambda(\lambda/3)/f_\lambda(\lambda)$
are known from the intrinsic SED.

### Our Pipeline Version vs. Parlanti

| Item | Parlanti (2025) | This work |
|------|----------------|-----------|
| Pipeline version | 1.18.1 | 1.20.2 |
| CRDS context | jwst_1364.pmap | jwst_1464.pmap (auto) |
| Reference modifications | Modified S-flat + F-flat to extend range | Standard reffiles for now (first reduction) |
| Stage 2 output | _cal.fits with extended λ | _cal.fits (nominal range, standard) |

**Current reduction uses unmodified reference files.** This means the first Stage 2/3
pass will produce standard (nominal range) cubes. The wavelength extension requires
applying the modified S-flat and F-flat — that step comes next.

---

## What Happens Next

### Step 1: Standard reduction (in progress)

`run_ifu_pipeline.py --stage 23` produces nominal-range cubes as a baseline:
- Verify pipeline runs cleanly on our data
- Inspect 1D spectra to confirm correct flux levels in nominal range
- This gives us $f_\lambda(\lambda)$ from the nominal data

### Step 2: Modify reference files (following Parlanti Sec. 2.1)

For each grating/detector:

**S-flat modification** (`analysis/reffiles/modify_sflat.py` — to be created):
- NRS1: extend SCI trace (set extended pixels to 1) and clear DQ flags beyond λ_max_nom
- NRS1: extend `FAST_VARIATION` binary table with constant last value
- NRS2: clear all DQ flags; extend `FAST_VARIATION`

**F-flat modification** (`analysis/reffiles/modify_fflat.py` — to be created):
- G140M/F100LP: concatenate `FAST_VARIATION` from G140M + G235M + G395M tables
- G235M/F170LP: concatenate `FAST_VARIATION` from G235M + G395M tables

**Wavelength range modification** (`analysis/reffiles/modify_wavelengthrange.py` — to be created):
- Extend the `wavelengthrange` reference file entries for G140M and G235M to cover the full extended range

**Pipeline code patch:**
- Parlanti also modified `calwebb_spec2` to not error when NRS2 files are missing.
  In our forked pipeline (`jwst_nirspec_wavext/`), this patch can be applied to
  `jwst/pipeline/calwebb_spec2.py`.

### Step 3: Re-run Stage 2+3 with modified reffiles

Re-run Spec2Pipeline + Spec3Pipeline pointing to the modified reference files via
`override_*` parameters (or via a local CRDS override file). Output will be extended-range cubes.

### Step 4: Extract and compare spectra

For each standard star and each grating:
- Extract 1D spectrum from the extended cube ($S_\lambda(\lambda)$)
- Compare to known CALSPEC SED ($f_\lambda(\lambda)$) over the extended range
- Compute $R(\lambda) = S_\lambda / f_\lambda$

### Step 5: Solve for k(λ), α̃(λ), β̃(λ)

Using `analysis/solver/` machinery (analog of existing FS solver):
- With two standard stars (G191-B2B + P330E), we have 2 equations per λ
- Minimum is 3 sources; with 2 we need to fix β̃ or use prior from Parlanti's published curves
- Add a third source (e.g., J1743045 from PID 1537 if available in IFU mode) for full solution

### Step 6: Compare to Parlanti's published calibration curves

Parlanti et al. publish k(λ), α̃(λ), β̃(λ) in their GitHub repo. Direct comparison will
show whether pipeline version 1.20.2 and CRDS context differences change the calibration.

---

## Current Status

| Step | Status | Notes |
|------|--------|-------|
| Data discovery (MAST query) | ✅ Done | All 4 PIDs queried |
| Download PIDs 1537 + 1538 | ✅ Done | 32 files, 16 per PID |
| Download PID 1536 (J1743045) | 🔄 Running | ~24 files, G140M+G235M+G395M |
| Download PID 6645 (P330-E C3) | ⏳ Pending | queued after 1536 |
| Stage 2 — PID 1537 G191-B2B | ✅ Done | 8 NRS1 cal files |
| Stage 3 — PID 1537 G191-B2B | ✅ Done | G140M + G235M s3d + x1d |
| Stage 2 — PID 1538 P330E | ✅ Done | 8 NRS1 cal files |
| Stage 3 — PID 1538 P330E | ✅ Done | G140M + G235M s3d + x1d |
| Stage 2+3 — PID 1536 J1743045 | ⏳ Pending | after download |
| Stage 2+3 — PID 6645 P330E-C3 | ⏳ Pending | after download |
| Modify S-flat reference files | ⏳ Pending | following Parlanti Sec. 2.1.1 |
| Modify F-flat reference files | ⏳ Pending | following Parlanti Sec. 2.1.2 |
| Modify wavelengthrange reffile | ⏳ Pending | |
| Re-run pipeline with mod. reffiles | ⏳ Pending | |
| Solve for k, α̃, β̃ (IFU) | ⏳ Pending | |
| Compare to Parlanti results | ⏳ Pending | |

---

## Notes and Observations

- **IFU vs. FS:** The IFU mode places each slice on the same set of detector pixels
  at a given wavelength (slit-like stripes), in contrast to FS where the spatial
  dimension is along the slit. This means the NRS2 flat calibration is conceptually
  simpler for IFU: the entire extended-range region just needs DQ flags cleared and
  the flat extended with ones.

- **4-point dither pattern:** Both PIDs use a 4-point dither (PATT_NUM 1–4). For
  optimal cube reconstruction, all 4 dithers should be included in Stage 3. The
  `MAX_EXPOSURES_PER_GRATING = 4` setting in `run_ifu_pipeline.py` keeps all of them.

- **NRS2 in standard pipeline:** The standard `calwebb_spec3` discards NRS2 data for
  IFU medium-resolution gratings. To include it, either the pipeline must be patched
  (Parlanti's approach) or we run NRS1 and NRS2 through Spec2/3 separately and
  merge the resulting x1d spectra.

- **Pipeline v1.18.1 vs. v1.20.2:** The Parlanti paper used v1.18.1. Our work uses
  v1.20.2. The CRDS context has changed from `jwst_1364.pmap` to a more recent one.
  Reference file indices have incremented (e.g., sflat_0208 vs. earlier versions).
  Minor differences in calibration are expected; any systematic offsets in the solved
  k(λ) curves would be scientifically interesting.

- **Spatial independence confirmed by Parlanti:** The k(λ), α̃(λ), β̃(λ) curves do
  not depend on where in the IFU FOV the point source falls (verified using P330-E
  at different dither positions and different cycles). This means the same correction
  curves derived from our standard-star data apply to science targets anywhere in the FOV.

- **Two stars vs. three for the calibration solve:** With 3+ independent standard star
  observations (PIDs 1536 + 1537 + 1538, plus PID 6645 as 4th), we can directly solve
  the 3×3 linear system at each wavelength for k(λ), α̃(λ), β̃(λ) exactly as Parlanti does.
  Each star provides one equation per wavelength:
  $R_i(\lambda) = k(\lambda) + \tilde\alpha(\lambda)\,r_{2,i}(\lambda) + \tilde\beta(\lambda)\,r_{3,i}(\lambda)$
  where $r_{2,i}$ and $r_{3,i}$ are computed from each star's known CALSPEC SED.
  With 4 sources we can solve all combinations of 3 and average, matching Parlanti's approach.

- **G395M for PID 1536:** J1743045 was observed with G395M/F290LP in addition to G140M and G235M.
  This lets us independently verify the G235M extended range using the G395M nominal coverage
  (2.87–5.27 µm) as a cross-check on the CALSPEC model over 3–5 µm.

- **Flux Unit Discrepancy (Jy vs. ADU/s):** The raw NRS2 extension outputs (pre-Stage 3) often
  lack the `BUNIT` and `PHOTMJSR` keywords in their headers. This results in a flux scale
  that is approximately ~800 times higher (e.g., 25 Jy median vs. 0.03 Jy truth for P330E) 
  than the calibrated reference. For the Parlanti solve, a manual anchoring factor 
  $k_{raw} \approx 0.0013$ for P330E is required to bring these into physical alignment before 
  ghost subtraction. Artifacts/spikes in the raw data can reach $10^4$ or higher.
