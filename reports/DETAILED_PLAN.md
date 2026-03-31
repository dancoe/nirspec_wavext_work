# NIRSpec Wavelength Extension: Comprehensive Analysis Plan & Guide

This document serves as the master manual and roadmap for the NIRSpec Wavelength Extension project (0.6–5.6 µm). It consolidates the "wisdom" gained through v1–v5 iterations, providing links to all critical scripts, reference files, and workflow steps required to reproduce the full analysis.

---

## 1. Environment & Setup

The analysis requires a specific forked version of the JWST pipeline and a pinned environment to ensure compatibility with custom ASDF schemas and coordinate transforms.

### 1.1 Repositories
- **Project Workspace:** `nirspec_wavext_work/` (Analysis scripts, reports, and data management).
- **Forked Pipeline:** `jwst_nirspec_wavext/` (Custom patches in `assign_wcs`, `extract_2d`, etc.).

### 1.2 Environment Activation
```bash
# Activate the pinned environment
micromamba activate jwst_1.20.2

# Point to the custom patched pipeline
export PYTHONPATH=/Users/dcoe/NIRSpec/wavext/jwst_nirspec_wavext
```
> [!IMPORTANT]
> Always verify the environment using `python -c "import jwst; print(jwst.__version__)"`. It must be `1.20.2`.

---

## 2. Data Acquisition

We target specific calibration standards and science validation targets to build and verify the ghost contamination model.

### 2.1 Target List
- **Calibration Standards:** PIDs 1536 (A8III), 1537 (WD), 1538 (G2V), 6644 (G1V), 6645.
- **Science Validation:** PIDs 2654 (AGN), 2186 (ULIRG), 1492 (IRAS-05248).

### 2.2 Download & Audit
1. **Download Scripts:** 
   - [download_v5_data.py](download_v5_data.py) (Main)
   - [download_ifu_rates.py](download_ifu_rates.py) (Standard stars)
2. **Technical Gotchas (Download):**
   - **URI Zero-Padding:** MAST URIs are zero-padded to 5 digits (e.g., `jw01538...`) while metadata is unpadded (`1538`).
   - **FITS Header Stubs:** Watch for files exactly 64KB (one FITS block); these are stalled downloads that lack ASDF extensions.
   - **Size Check:** Ensure `_rate.fits` files are > 50MB to avoid target acquisition (TA) images or corrupted stubs.
3. **Verification:** Run the audit tool to ensure all `_rate.fits` files are present.
   - [audit_downloads.py](audit_downloads.py)
   - See [DATA_DOWNLOADS.md](DATA_DOWNLOADS.md) for strategy details.

---

## 3. Reduction Pipeline (Extended WCS)

Standard pipeline steps are modified to bypass nominal wavelength cutoffs (1.9 µm for G140M, 3.1 µm for G235M).

### 3.1 Key Reference Overrides
- **Wavelength Range:** [wavelengthrange_extended.asdf](wavelengthrange_extended.asdf) (Defines 0.6–5.6 µm range).
- **Flat Fields:** Parlanti-extended flats (`fflat_0105.fits`, `sflat_0191.fits`, etc.).
- **Photometry:** Extended photom files (`jwst_nirspec_photom_fs_nrs2_ext.fits`).

### 3.2 Reduction Execution
- **Fixed Slit (FS):**
  ```bash
  python nirspec_wavext_work/analysis/reduction/run_fs_nrs1_spec3.py --pid 6644
  python nirspec_wavext_work/analysis/reduction/run_fs_nrs2_spec3.py --pid 6644
  ```
- **IFU Mode:**
  - **4-Point Dithers:** Essential for optimal cube reconstruction; ensure all dithers are in the Stage 3 ASN.
  - **Spatial Independence:** Calibration ($k, \alpha, \beta$) is independent of the source position within the IFU FOV.
  ```bash
  python nirspec_wavext_work/analysis/reduction/run_ifu_pipeline_ext.py --pid 2654
  ```
- See [PIPELINE_RUNNING.md](PIPELINE_RUNNING.md) for command details.

---

## 4. Calibration Analysis (The Parlanti Model)

We utilize the model from **Parlanti et al. (2025)** to correct for throughput loss ($k$) and ghost contamination ($\alpha$ for 2nd order, $\beta$ for 3rd order).

### 4.1 The Model
$$S_{obs}(\lambda) \approx k(\lambda) f_{true}(\lambda) + \alpha(\lambda) f_{true}(\lambda/2) + \beta(\lambda) f_{true}(\lambda/3)$$

### 4.2 The v5 Joint Solve (NNLS)
To break the $k/\alpha$ degeneracy, we perform a simultaneous Non-Negative Least Squares (NNLS) fit using 4 standards (Hot stars + G1V).
- **Primary Solver:** [solve_parlanti_calib_v5.py](solve_parlanti_calib_v5.py)
- **Wisdom:** The G1V star (PID 6644) is essential because its unique spectral shape allows the solver to distinguish between a drop in $k$ and an increase in $\alpha$.

---

## 5. Diagnostics & Visualization

Validation is performed by comparing NRS1 (Truth) with corrected NRS2 (Extended) data.

### 5.1 Coefficient Plots
Compare derived coefficients against Parlanti (2025) and previous iterations.
- [plot_parlanti_calib_v5.py](plot_parlanti_calib_v5.py)
- [plot_parlanti_coeffs_log_v5.py](plot_parlanti_coeffs_log_v5.py)

### 5.2 Spectral Validation
Diagnostic plots showing full-spectrum recovery (0.6–5.6 µm).
- **Merged Spectrum:** [generate_final_plots.py](generate_final_plots.py)
- **Component Decomposition:** [plot_parlanti_components_multi.py](plot_parlanti_components_multi.py)

---

## 6. Reporting & Wisdom

All major results and technical breakthroughs are documented in versioned reports.

### 6.1 Report Index
- [REPORTS.md](REPORTS.md) — Master index of all project reports.
- [LATEST_WORK.md](LATEST_WORK.md) — Current status and immediate hurdles.

### 6.2 Key Wisdom (The "Gotchas")
- **Recursive Discovery:** Use `glob.glob` with `recursive=True` to find `x1d.fits` files, as the pipeline often scatters then in `target/` or `level3/` folders.
- **Unit Mismatch:** Raw extractions from `extract_2d` are in **DN/s**. You MUST apply the `photom` step or a boundary-matched scale factor to compare with CALSPEC **Jy** models.
- **Flux Anchoring (IFU):** If `BUNIT` or `PHOTMJSR` keywords are missing in NRS2 extensions, apply an anchoring factor (e.g., $k_{raw} \approx 0.0013$ for P330E) to align scales.
- **Lazy Loading:** ASDF objects in the pipeline can lead to "closed file" errors if pointers are accessed after the file handle closes. Always materialise lists (`list(wrange)`) before closing.
- **IFU-Specific Detail:** See [PARLANTI_IFU.md](PARLANTI_IFU.md) for deep-dives into the IFU-mode reduction branch and MAST query specifics.

---

## 7. Project Navigation

| Resource | Link |
| :--- | :--- |
| **Project Index** | [INDEX.md](INDEX.md) |
| **Instructions** | [INSTRUCTIONS.md](INSTRUCTIONS.md) |
| **Observation Log** | [OBSERVATIONS.md](OBSERVATIONS.md) |
| **Change Log** | [CHANGE_LOG.md](CHANGE_LOG.md) |

---
*Created by Antigravity — v5 Execution Readiness: Phase 6/Phase 6*
