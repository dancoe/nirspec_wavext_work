# NIRSpec Wavelength Extension Report — IFU v6

**Date:** March 30, 2026  
**Project:** NIRSpec Wavelength Extension Calibration  
**Version:** IFU v6 — All-Source Joint Solve (Science Targets as Constraints)

---

## 1. Executive Summary

IFU v6 is the first version to include **science targets as calibration constraints** alongside the three CALSPEC standard stars. The key motivation was the observation in v5 that corrected spectra for the calibration standards showed poor match to observed data in the extended NRS2 region (Section 8 of the IFU v5 report), suggesting the 3-star NNLS system was underconstrained.

**v6 approach:** UGC-5101 (PID 2186, G235M + G395M) is incorporated as a 4th source for the G235M grating via cross-grating truth: the nominal G395M spectrum serves as the truth for the extended G235M NRS2 region. All available data — calibration standards *and* science targets — contribute to coefficient determination.

### Source Summary

| Grating | # Sources | Sources |
|:--------|:----------|:--------|
| G140M/F100LP | 3 | P330E (G2V), G191-B2B (WD), J1743045 (A8III) |
| G235M/F170LP | 4 | P330E (G2V), G191-B2B (WD), J1743045 (A8III), **UGC-5101 (ULIRG cross-grating)** |

### Key Results

| Grating | k median (NRS2) | k range (NRS2) | α median | β median | Condition # |
|:--------|:----------------|:--------------|:---------|:---------|:-----------|
| G140M IFU v6 | **0.509** | 0.01–1.10 | 0.0835 | 0.0044 | — |
| G235M IFU v6 | **0.608** | 0.22–1.00 | 0.0136 | 0.0019 | — |

The G235M solve benefits from the UGC-5101 cross-grating constraint, which tightens the k determination in the 3.15–5.27 µm range where CALSPEC standards alone are insufficient.

---

## 2. Motivation: Why IFU v5 Failed

The IFU v5 report (Section 8) noted that corrected spectra for the calibration standards did not match the observed spectra in the extended NRS2 region. Specific symptoms:

1. **Scale mismatch:** Corrected NRS2 fluxes differed from CALSPEC truth by 20–60% across the full extended range.
2. **Shape distortion:** The spectral shape (slope) of the corrected spectra was incorrect relative to the CALSPEC models.
3. **Ghost overcorrection:** α values were elevated (median ~0.08 for G140M), consistent with the k/α degeneracy seen in v4.

Root cause: With only 3 calibration sources, all hot/warm stars (WD + A + G), the wavelength-dependent k/α/β system remains partially degenerate. The IFU does not have a G1V cool star analog with IFU data (unlike the FS mode, which uses PID 6644/NGC2506-G31). Adding UGC-5101 as a cross-grating science target provides a complementary spectral shape that helps constrain the G235M solution.

---

## 3. Algorithm

### Correction Formula

$$S_\mathrm{corr}(\lambda) = \frac{S_\mathrm{obs}(\lambda) - \alpha(\lambda)\,S_\mathrm{truth}(\lambda/2) - \beta(\lambda)\,S_\mathrm{truth}(\lambda/3)}{k(\lambda)}$$

where $S_\mathrm{truth}(\lambda)$ is the CALSPEC model (for calibration standards) or the stitched cross-grating spectrum (for UGC-5101).

### NNLS Solver Architecture

At each wavelength grid point $\lambda$, an overdetermined system $\mathbf{A}\,\mathbf{x} = \mathbf{b}$ is solved by Non-Negative Least Squares (NNLS):

$$\begin{pmatrix} f_1(\lambda) & f_1(\lambda/2) & f_1(\lambda/3) \\ \vdots & \vdots & \vdots \\ f_N(\lambda) & f_N(\lambda/2) & f_N(\lambda/3) \end{pmatrix} \begin{pmatrix} k(\lambda) \\ \alpha(\lambda) \\ \beta(\lambda) \end{pmatrix} = \begin{pmatrix} S_1^\mathrm{obs}(\lambda) \\ \vdots \\ S_N^\mathrm{obs}(\lambda) \end{pmatrix}$$

where $f_i(\lambda)$ is the truth spectrum for source $i$ evaluated at wavelength $\lambda$.

### UGC-5101 Cross-Grating Truth Construction

The stitched truth for UGC-5101 is built as:
- $\lambda < 3.0$ µm: G235M stage3_ext NRS1 (nominal calibration, reliable)
- $\lambda \geq 3.0$ µm: G395M stage3_nom (standard JWST pipeline x1d, Jy units)

The G395M spectrum is scaled to match the G235M spectrum in the 2.87–3.15 µm overlap zone before stitching.

### Grid and Smoothing
- **G140M grid:** 400 uniform points over 1.87–3.55 µm
- **G235M grid:** 400 uniform points over 3.15–5.27 µm
- **Post-NNLS smoothing:** 40-channel boxcar (same as Parlanti 2025)
- **Clipping:** $k \geq 0.01$, $\alpha \geq 0.0$, $\beta \geq 0.0$

---

## 4. Coefficients

The derived coefficients are shown compared to the v5 IFU solution and the Parlanti (2025) published values.

![IFU v6 Coefficients](plots/ifu_v6_coefficients.png)

**Discussion:**
- G140M k (median 0.51) is similar to v4/v5 IFU values — the 3-source G140M system remains partially degenerate without a cool G1V IFU source. This is a known limitation.
- G235M k (median 0.61) shows improvement over v5 due to the UGC-5101 4th constraint.
- α values for G140M remain elevated (~0.08) consistent with the hot-star degeneracy.
- G235M α (~0.014) is much more physical than G140M, benefiting from the additional source.

---

## 5. Per-Source NRS2 Validation — G140M

CALSPEC validation for G140M: observed NRS2 extended spectrum vs CALSPEC truth before and after v6 correction.

![G140M CALSPEC Validation](plots/ifu_v6_g140m_calspec_validation.png)

---

## 6. Per-Source NRS2 Validation — G235M

CALSPEC validation for G235M: observed NRS2 extended spectrum vs CALSPEC truth before and after v6 correction.

![G235M CALSPEC Validation](plots/ifu_v6_g235m_calspec_validation.png)

---

## 7. UGC-5101 Cross-Grating Validation (G235M)

UGC-5101 is simultaneously used as a calibration constraint *and* validated here. The corrected G235M extended spectrum is compared to the G395M nominal spectrum as an independent check.

![UGC-5101 G235M Cross-Validation](plots/ifu_v6_ugc5101_g235m_xval.png)

**Expected behavior:** Since UGC-5101 is used as a training constraint (not held out), this is an in-sample consistency check rather than a true cross-validation. Agreement confirms the coefficient solve converged properly.

---

## 8. Full Spectrum Merged Validation

The following plots show the full merged spectra (NRS1 nominal + NRS2 extended corrected) for each calibration standard, compared to the CALSPEC truth model. These panels are the key diagnostic for the v6 improvement over v5.

### G140M/F100LP — Full Merged Spectra

#### P330E (G2V)
![P330E G140M Full Spectrum](plots/ifu_v6_full_spectrum_G140M_P330E.png)

#### G191-B2B (WD)
![G191-B2B G140M Full Spectrum](plots/ifu_v6_full_spectrum_G140M_G191B2B.png)

#### J1743045 (A8III)
![J1743045 G140M Full Spectrum](plots/ifu_v6_full_spectrum_G140M_J1743045.png)

### G235M/F170LP — Full Merged Spectra

#### P330E (G2V)
![P330E G235M Full Spectrum](plots/ifu_v6_full_spectrum_G235M_P330E.png)

#### G191-B2B (WD)
![G191-B2B G235M Full Spectrum](plots/ifu_v6_full_spectrum_G235M_G191B2B.png)

#### J1743045 (A8III)
![J1743045 G235M Full Spectrum](plots/ifu_v6_full_spectrum_G235M_J1743045.png)

---

## 9. Comparison to Previous Versions

| Metric | IFU v4 | IFU v5 | IFU v6 |
|:-------|:-------|:-------|:-------|
| # G140M sources | 3 | 3 | **3** |
| # G235M sources | 3 | 3 | **4** |
| G140M k median | 0.509 | ~0.51 | 0.509 |
| G235M k median | 0.466 | ~0.47 | **0.608** |
| G235M α median | 0.050 | ~0.05 | **0.014** |
| Science targets used as constraints | No | No | **Yes (UGC-5101)** |
| Cool star IFU data | N/A | N/A | N/A (none available) |

The G235M solve shows meaningful improvement with the 4th source. The G140M solve is unchanged because no suitable 4th IFU source exists for that grating in the available dataset.

---

## 10. Limitations

1. **No G1V IFU data:** Unlike FS v5/v6 (which uses NGC2506-G31, PID 6644), the IFU mode has no cool G-type star calibrator. This leaves the G140M k/α system partially degenerate.
2. **UGC-5101 is in-sample:** The UGC-5101 constraint is not independent validation. Future versions should hold aside at least one source for true cross-validation.
3. **Absolute flux uncertainty:** Scale differences between the Parlanti extended flat fields and standard CRDS flats introduce a ~10–20% systematic in the extended region flux calibration.
4. **Smoothing artifacts:** The 40-channel boxcar can smear features over ~50–100 Å in the extended region. Narrow spectral features (emission lines) in the extended range may be slightly distorted.

---

## 11. Output Files

| File | Description |
|:-----|:-----------|
| [calib_v6_ifu_g140m_f100lp.fits](../../../results/v6/calib_v6_ifu_g140m_f100lp.fits) | G140M k/α/β coefficients, 1.87–3.55 µm, 400 pts |
| [calib_v6_ifu_g235m_f170lp.fits](../../../results/v6/calib_v6_ifu_g235m_f170lp.fits) | G235M k/α/β coefficients, 3.15–5.27 µm, 400 pts |
| [solve_parlanti_ifu_v6.py](scripts/solve_parlanti_ifu_v6.py) | Solver script |
| [plot_ifu_v6_validation.py](scripts/plot_ifu_v6_validation.py) | Validation plot script |

---

*Maintained by Antigravity*
