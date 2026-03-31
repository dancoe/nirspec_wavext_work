# Detailed Analysis Plan v8 — FS Final Calibration & Gap Handling

**Version:** 8.0  
**Status:** Planning  
**Objective:** Finalize the Fixed Slit (FS) wavelength extension coefficients with refined gap handling and full-spectrum validation.

---

## 1. Data Inventory & Preparation
- **PIDs:** 1492, 1536, 1537, 1538, 6644.
- **Status:** All v7 native Jy-calibrated products are local.
- **Task:** Verify that all NRS2 extended spectra (0.6 - 5.6 µm) are present and correctly formatted in `data/FS/`.

---

## 2. Solver Refinement: Gap Awareness
Systematically handle the physical gaps in the NIRSpec detector layout to avoid artificial smoothing or interpolation across unobserved regions.

| Step | Task | Success Criteria |
| :--- | :--- | :--- |
| **2.1** | Identify precise gap regions. | Confirm gaps at ~2.2 µm (G140M) and ~3.7 µm (G235M). |
| **2.2** | Update `solve_parlanti_fs_v8.py`. | Ensure the spline/polynomial fit for $k, \alpha, \beta$ does not rely on "bridging" these gaps if data is missing. |
| **2.3** | Masking. | Explicitly mask out these wavelength regions in the training set to prevent them from biasing the fit. |

---

## 3. Visualization: "Gap-Safe" Plotting
Update plotting logic to respect physical gaps in the merged spectra.

| Step | Task | Success Criteria |
| :--- | :--- | :--- |
| **3.1** | Create `plot_v8_full_spectra_v5_style.py`. | Inherit v7 logic but split the spectrum into separate segments at gap boundaries. |
| **3.2** | No Interpolation. | Ensure `ax.plot` does not draw lines between the end of NRS1 and the start of NRS2 if there is a gap. |
| **3.3** | Merged Labeling. | Use the "NRS1,2" labeling convention established in v7. |

---

## 4. Verification Checklist (v8)
- [ ] Gaps at ~2.2 µm and ~3.7 µm are visible and not "smoothed over".
- [ ] All 5 primary FS targets (including cool star 6644) are included in the final solve.
- [ ] Hold-out residuals remain stable compared to v7.
- [ ] Final $k(\lambda), \alpha(\lambda), \beta(\lambda)$ coefficients saved to `results/v8/`.

---

*Created: 2026-03-31*  
*Author: Antigravity*
