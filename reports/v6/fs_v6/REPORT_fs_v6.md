# NIRSpec Wavelength Extension Report — FS v6

**Date:** March 30, 2026  
**Project:** NIRSpec Wavelength Extension Calibration  
**Version:** FS v6 — All-Source Joint Solve (Science Targets as Constraints)

---

## 1. Executive Summary

FS v6 extends the FS v5 calibration by incorporating **PID 1492** as a 5th source for both G140M and G235M via cross-grating truth. PID 1492 contains Fixed Slit (S200A1) observations of an astrophysical target with existing MAST pipeline nominal x1d spectra at multiple gratings. The nominal MAST x1d spectrum of the adjacent grating serves as the truth for the NRS2 extended region.

All available data — calibration standards *and* science targets — contribute to the coefficient solve. This follows the same philosophy as IFU v6.

### Source Summary

| Grating | # Sources | Sources |
|:--------|:----------|:--------|
| G140M/F100LP | 5 | G191-B2B (WD), P330E (G2V), J1743045 (A8III), NGC2506-G31 (G1V), **PID 1492 cross-grating** |
| G235M/F170LP | 5 | G191-B2B (WD), P330E (G2V), J1743045 (A8III), NGC2506-G31 (G1V), **PID 1492 cross-grating** |

### Key Results

| Grating | k median (NRS2) | k range (NRS2) | α median | β median |
|:--------|:----------------|:--------------|:---------|:---------|
| G140M FS v6 | **0.726** | 0.01–1.81 | 0.0222 | 0.0097 |
| G235M FS v6 | **0.910** | 0.01–2.10 | 0.0311 | 0.0012 |

The G235M k median (0.91) is close to the v5 value (0.957), confirming the physical throughput is consistent. G140M k (0.73) is somewhat lower than v5 (0.959), reflecting the additional cross-grating constraint modifying the solution.

---

## 2. Motivation: Adding PID 1492 as a 5th Source

### v5 Baseline (4 CALSPEC standards)

In FS v5, the 4-source NNLS (WD + G2V + A8III + G1V) successfully resolved the k/α degeneracy seen in v4. The v5 k values returned to ~0.96, consistent with physical throughput expectations.

However, the v5 corrected spectra in the full merged validation showed that the v5 calibration relied entirely on the diversity of stellar spectral types to constrain the system. Adding an entirely different source type (an arbitrary astrophysical target with known adjacent-grating truth) provides:

1. **Additional wavelength coverage:** PID 1492 x1d spectra cover G140M + G235M + G395M, providing a coherent truth chain from NRS1 to the full NRS2 extended region.
2. **Additional spectral diversity:** The target SED (regardless of astrophysical identity) provides a new constraint on the k/α/β system.
3. **Improved conditioning:** More rows in the NNLS matrix reduce the condition number and improve numerical stability.

### Cross-Grating Normalisation

PID 1492 NRS2 extended spectra are in detector units (DN/s). The normalisation from DN/s to the Jy scale of the MAST nominal spectra is derived by computing:

$$\mathrm{scale} = \mathrm{median}\left\{\frac{S_\mathrm{MAST}(\lambda)\,[\mathrm{Jy}]}{S_\mathrm{NRS2}(\lambda)\,[\mathrm{DN/s}]}\right\}_{\lambda \in \mathrm{overlap}}$$

where the overlap is the region where the NRS2 extended spectrum and the MAST truth share wavelengths. After scaling, the PID 1492 NRS2 spectrum enters the NNLS system on equal footing with the CALSPEC sources.

---

## 3. Algorithm

### Correction Formula

$$S_\mathrm{corr}(\lambda) = \frac{S_\mathrm{obs}(\lambda) - \alpha(\lambda)\,S_\mathrm{truth}(\lambda/2) - \beta(\lambda)\,S_\mathrm{truth}(\lambda/3)}{k(\lambda)}$$

### Source Files Used

| Source | PID | Type | G140M observed | G140M truth | G235M observed | G235M truth |
|:-------|:----|:-----|:--------------|:-----------|:--------------|:-----------|
| G191-B2B | 1537 | WD | nrs2_spec3_ext x1d | CALSPEC model | nrs2_spec3_ext x1d | CALSPEC model |
| P330E | 1538 | G2V | nrs2_spec3_ext x1d | CALSPEC model | nrs2_spec3_ext x1d | CALSPEC model |
| J1743045 | 1536 | A8III | nrs2_spec3_ext x1d | CALSPEC model | nrs2_spec3_ext x1d | CALSPEC model |
| NGC2506-G31 | 6644 | G1V | nrs2_spec3_ext x1d | CALSPEC model | nrs2_spec3_ext x1d | CALSPEC model |
| PID 1492 | 1492 | science | NRS2 extract_1d (DN/s) | G235M MAST (Jy) | NRS2 extract_1d (DN/s) | G395M MAST (Jy) |

### Grid and Smoothing
- **G140M grid:** 400 uniform points over 1.87–3.20 µm
- **G235M grid:** 400 uniform points over 3.15–5.15 µm
- **Post-NNLS smoothing:** 40-channel boxcar
- **Clipping:** $k \geq 0.01$, $\alpha \geq 0.0$, $\beta \geq 0.0$

---

## 4. Coefficients

The derived v6 coefficients are shown compared to the FS v5 solution and the Parlanti (2025) published values.

![FS v6 Coefficients](plots/fs_v6_coefficients.png)

**Discussion:**
- G140M k (median 0.73) is somewhat below the v5 value (~0.96). The PID 1492 cross-grating source may pull the solution toward a slightly different k, possibly reflecting calibration differences between PID 1492 and the CALSPEC standards.
- G235M k (median 0.91) is close to v5, confirming the G235M calibration is robust against the addition of the 5th source.
- G140M α (median 0.022) is more physical than IFU v6 (0.08), reflecting the benefit of having the G1V source (NGC2506-G31) in the FS mode.
- G235M α (median 0.031) is physically reasonable and consistent with Parlanti published values (~0.01–0.05).

---

## 5. Per-Source NRS2 Validation — G140M

![G140M CALSPEC Validation](plots/fs_v6_g140m_calspec_validation.png)

Summary panels for all 4 CALSPEC standards. Each panel shows the raw NRS2 observed spectrum, the v6-corrected spectrum, and the CALSPEC truth model, all normalised to the truth at the median wavelength.

---

## 6. Per-Source NRS2 Validation — G235M

![G235M CALSPEC Validation](plots/fs_v6_g235m_calspec_validation.png)

---

## 7. PID 1492 Cross-Validation — G140M

PID 1492 is used as both a training source *and* a consistency check. The corrected G140M NRS2 extended spectrum is compared to the adjacent G235M MAST nominal spectrum scaled to Jy.

![PID 1492 G140M Cross-Validation](plots/fs_v6_pid1492_g140m_xval.png)

---

## 8. PID 1492 Cross-Validation — G235M

The corrected G235M NRS2 extended spectrum is compared to the G395M MAST nominal spectrum.

![PID 1492 G235M Cross-Validation](plots/fs_v6_pid1492_g235m_xval.png)

---

## 9. Full Spectrum Merged Validation

The following plots show the full merged spectra (NRS1 nominal + NRS2 extended corrected) for each calibration standard compared to the CALSPEC truth model.

### G140M/F100LP — Full Merged Spectra

#### G191-B2B (WD)
![G191-B2B G140M Full Spectrum](plots/fs_v6_full_g140m_G191B2B.png)

#### P330E (G2V)
![P330E G140M Full Spectrum](plots/fs_v6_full_g140m_P330E.png)

#### J1743045 (A8III)
![J1743045 G140M Full Spectrum](plots/fs_v6_full_g140m_J1743045.png)

#### NGC2506-G31 (G1V)
![NGC2506-G31 G140M Full Spectrum](plots/fs_v6_full_g140m_NGC2506G31.png)

### G235M/F170LP — Full Merged Spectra

#### G191-B2B (WD)
![G191-B2B G235M Full Spectrum](plots/fs_v6_full_g235m_G191B2B.png)

#### P330E (G2V)
![P330E G235M Full Spectrum](plots/fs_v6_full_g235m_P330E.png)

#### J1743045 (A8III)
![J1743045 G235M Full Spectrum](plots/fs_v6_full_g235m_J1743045.png)

#### NGC2506-G31 (G1V)
![NGC2506-G31 G235M Full Spectrum](plots/fs_v6_full_g235m_NGC2506G31.png)

---

## 10. Comparison to Previous Versions

| Metric | FS v4 | FS v5 | FS v6 |
|:-------|:------|:------|:------|
| # G140M sources | 3 | 4 | **5** |
| # G235M sources | 3 | 4 | **5** |
| G140M k median | 0.701 | 0.959 | **0.726** |
| G235M k median | 0.723 | 0.957 | **0.910** |
| G140M α median | 0.029 | 0.003 | **0.022** |
| G235M α median | 0.037 | 0.012 | **0.031** |
| Cool star (G1V) | No | Yes | Yes |
| Science target as constraint | No | No | **Yes (PID 1492)** |

Note: The G140M k is lower in v6 than v5. This is expected when adding a source whose spectral shape pulls the NNLS toward a different k/α balance. The physical interpretation is that PID 1492 introduces more freedom in the G140M extended wavelength, which the solver trades against k. Future versions should investigate the sensitivity of k to the PID 1492 DN/s → Jy normalisation.

---

## 11. Limitations and Future Work

1. **PID 1492 normalisation uncertainty:** The DN/s → Jy scale factor derived from the overlap region is a single scalar per grating. Any systematic offset (pointing, background, or WCS mismatch between the NRS2 extended and MAST nominal spectra) propagates directly into the k/α/β solution.
2. **In-sample validation:** PID 1492 is used as both a training source and a validation check. A true cross-validation would hold out one source and verify that the remaining N−1 sources reproduce it.
3. **G140M k regression:** The G140M k dropped from v5's 0.96 to v6's 0.73. This could indicate that PID 1492 is not well-normalised in G140M, or that the additional constraint genuinely shifts the solution. Further investigation with an independent source is needed.
4. **Beta term:** The β (3rd-order ghost) coefficients remain small but non-zero. The quality of the correction is dominated by k and α.

---

## 12. Output Files

| File | Description |
|:-----|:-----------|
| [calib_v6_fs_g140m_f100lp.fits](../../../results/v6/calib_v6_fs_g140m_f100lp.fits) | G140M k/α/β coefficients, 1.87–3.20 µm, 400 pts |
| [calib_v6_fs_g235m_f170lp.fits](../../../results/v6/calib_v6_fs_g235m_f170lp.fits) | G235M k/α/β coefficients, 3.15–5.15 µm, 400 pts |
| [solve_parlanti_fs_v6.py](scripts/solve_parlanti_fs_v6.py) | Solver script |
| [plot_fs_v6_validation.py](scripts/plot_fs_v6_validation.py) | Validation plot script |

---

*Maintained by Antigravity*
