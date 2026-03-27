# PARLANTI.md

## Reference

**Parlanti et al. 2025** (arXiv:2512.14844v1, submitted 16 Dec 2025)  
*"Doubling NIRSpec/IFS capability to calibrate the single epoch black hole mass relation at high redshift"*  
Authors: E. Parlanti, B. Trefoloni, S. Carniani, F. D'Eugenio, M. Perna, G. Tozzi, H. Übler, G. Venturi, S. Zamora

GitHub repository: https://github.com/eleonoraparlanti/nirspecIFU-extended

---

## Overview

Parlanti et al. develop and release a **modification of the JWST Science Calibration Pipeline** for **NIRSpec IFU** observations that more than doubles the wavelength coverage of the medium-resolution (R~1000) G140M/F100LP and G235M/F170LP configurations. The technique extends the spectra onto the second detector (NRS2), which is discarded in the standard reduction.

Key scientific application: recover H-alpha emission (from the extended G140M/F100LP) for 5 QSOs at z~2 that nominally only cover H-beta and [OIII]—enabling direct comparison of H-alpha and H-beta SE BH mass estimators against reverberation-mapping (RM) masses.

---

## Section 2 — Changes to the Standard Data Reduction (Reproduced in Detail)

### Preamble (start of Sec. 2)

> "Data from NIRSpec/IFS observed with medium-resolution (G140M/F070LP, G140M/F100LP, G235M/F170LP, G395M/F290LP), with resolving power R = λ/Δλ varying between ~700 and ~1300 and an average of R~1000, occupy just the first detector, namely 'NRS1', while the second detector is not used in the standard reduction process. However, we notice emission lines in the second detector and continuum emission until the red end of 'NRS2', as shown in Fig. 1, displaying the count-rate maps downloaded from the Mikulski Archive for Space Telescopes (MAST) of the G140M/F100LP observations from PID: 2057, which targeted a QSO at z~2.6. The count-rate maps reveal the presence of the [OIII] and H-beta emission lines in the nominal wavelength range, but also the H-alpha emission in the second detector beyond the nominal wavelength range for G140M/F100LP observations."

> "In this work, we extend the wavelength range of medium-resolution gratings to recover the emission in the second detector. Extending to the second detector basically doubles the wavelength coverage, allowing us to effectively obtain the equivalent of two grating/filter observations in a single observation."

> "The extension is made possible by the fact that the photon conversion efficiency does not steeply drop as a function of wavelength, in particular for the G140M/F100LP and the G235M/F170LP grating/filter configurations, which benefit most from the extension. For the G395M/F290LP, extending the wavelength range beyond the nominal limit of 5.27 µm is limited by the detector sensitivity, which drops significantly, reaching 0 at λ~5.5 µm."

**IFU note:** G140M/F070LP could in principle be extended, but with limited scientific gain and more effort (range 0.90–1.26 µm in IFU, second- and third-order contamination would require fourth- and fifth-order corrections). Not treated in this work.

### Table 1: Wavelength Coverage

| Config | Nominal Range | Extended Range | Gap |
|--------|--------------|----------------|-----|
| G140M/F100LP | 0.97–1.88 µm | 0.97–3.55 µm | 2.17–2.28 µm |
| G235M/F170LP | 1.70–3.15 µm | 1.66–5.27 µm | 3.66–3.82 µm |

---

### 2.1 Data Reduction

> "To extend the medium resolution gratings, we first create and modify the reference files to extract flux beyond the nominally calibrated wavelength range, and then we perform the flux calibration. For the initial tests and calibration steps, we exploited the observations of standard stars, observed in the calibration programs 6645, 1536, 1537, and 1538. The four programs observed three different stars with all the possible configurations of gratings and filters, with the goal of calibrating the flux and the flat fields for the IFU and other NIRSpec modes. We also used the observations of other targets (local galaxies, AGN, high-z galaxies) obtained with different gratings to verify the accuracy of the calibrations and calibrate our extension of the medium resolution gratings, which we then applied to other targets."

> "The complementary sources that help verify and calibrate our extension were selected to respect the following criteria:
> i) have observations with **at least two consecutive medium resolution gratings** (i.e., G140M/F100LP + G235M/F170LP or G235M/F170LP + G395M/F290LP)
> ii) being bright enough to have **continuum detection with S/N > 5** across all the wavelength range covered."

**Pipeline Version:** `calwebb_spec2` and `calwebb_spec3` version **1.18.1**, context **`jwst_1364.pmap`**.

**Reference file modifications made:**
1. Extrapolated each `S-flat` above the nominal wavelength range with ones (see Sec. 2.1.1)
2. Modified the `F-flat` files so each is extended up to 5.5 µm by concatenating the F-flat for every filter (see Sec. 2.1.2)
3. Modified the wavelength range for each grating in the `cube_build` and `nirspec_wavelength_range` reference files

**Pipeline code modification:** Removed the check for the presence of files from the second detector in the `calwebb_spec2` step (otherwise causes an error).

#### 2.1.1 Spectroscopic Flat (S-flat)

> "The S-flat files consist of one file for each detector and each filter/grating configuration, and account for all the losses in the optical path from the aperture plane until the disperser. Each file is constituted by three image extensions and 2 Binary tables. The small corrections, slowly varying with wavelength, are stored in the 'SCI' extension, while the Binary table extension labeled 'FAST_VARIATION' encapsulates the rapidly wavelength-dependent variations with larger relative amplitude."

> "To extend the wavelength range of the reduced data, we modified the 'SCI' and 'DQ' extensions, which contain the information of the detector response at each pixel and the data quality flag for each pixel. Since in the IFU, for each grating and filter, the various wavelengths fall on the same pixel in the detector, the pixels containing wavelengths outside the nominal range have a value of 'nan' in the SCI flat and are flagged as pixels not to be used in the DQ extension. Hence, for the S-flat in NRS1, we extended each trace until the end of the detector in the 'SCI' and 'DQ' images. Because the spectral traces exhibit a slight curvature and tilt on the detector, we fitted a linear polynomial to their edges and extended the resulting paths according to the best-fit model. The pixels in the extended 'SCI' version have been fixed at 1, while in 'DQ', we remove each flag, inserting the value 0."

> "Moreover also the 'FAST_VARIATION' table, which is only defined in the nominal wavelength range, must be modified. In particular, we extended it with a constant value after the last element defined until the maximum wavelength covered by the extended detectors."

**NRS2 S-flat (IFU only):** In IFU mode, NRS2 DQ is all flagged as "do not use" in the standard cal. The modification: change the "DQ" extension to remove the flag for each pixel, extend "FAST_VARIATION" as above. (NRS2 SCI is initialized with all ones in the standard cal, so no change needed there.)

#### 2.1.2 Fore-optics Flat (F-flat)

> "The F-flat accounts for the losses from all the reflections in the telescope and in the NIRSpec fore optics. The F-flat does not depend on the detector, but is defined just for each filter/grating combination. In order to extend the wavelength coverage, we concatenated the Binary table extension describing the 'FAST_VARIATION' of one filter/grating with the ones of the following one. For instance, the newly defined F-flat for the G140M/F100LP is now the concatenation with the ones of the G235M/F170LP until 3.15 µm and the one of the G395M/F290LP until 5.2 µm."

---

### 2.2 Flux Calibration

**Context:** After the modified pipeline runs, the output `calwebb_spec3` datacube has an extended wavelength range, but **only spectra within the nominal range are flux-calibrated**. The extended region requires an empirical correction.

**Key verification (spatial independence):** Using standard star P330-E (PIDs 1538 and 6645), observed at different positions in the FOV and in different cycles (cycle 1 and cycle 3): no significant differences found. **Conclusion: the empirical correction curve depends only on wavelength, not on the spatial position of the target in the FOV.**

#### The Model Equation

> "Assuming the contributions only up to the third order, the output of `calwebb_spec3` can be expressed as:"

$$
S_\lambda(\lambda) = \begin{cases}
f_\lambda(\lambda) & \lambda \le \lambda_\mathrm{max}^\mathrm{nom} \\
k(\lambda)\,f_\lambda(\lambda) + \tilde\alpha(\lambda)\,f_\lambda(\lambda/2) + \tilde\beta(\lambda)\,f_\lambda(\lambda/3) & \lambda > \lambda_\mathrm{max}^\mathrm{nom}
\end{cases}
$$

where:
- $\lambda_\mathrm{max}^\mathrm{nom}$ = maximum wavelength of the NIRSpec nominal range
- $f_\lambda(\lambda)$ = intrinsic source spectrum (1st order)
- $g_\lambda(\lambda) \propto f_\lambda(\lambda/2)$ = 2nd-order spectrum
- $h_\lambda(\lambda) \propto f_\lambda(\lambda/3)$ = 3rd-order spectrum
- $k(\lambda)$ = wavelength-dependent 1st-order transmission efficiency
- $\tilde\alpha(\lambda)$ = effective 2nd-order contamination coefficient (relative to $f_\lambda$)
- $\tilde\beta(\lambda)$ = effective 3rd-order contamination coefficient (relative to $f_\lambda$)

Dividing both sides by $f_\lambda(\lambda)$:

$$
R(\lambda) = \frac{S_\lambda(\lambda)}{f_\lambda(\lambda)} = k(\lambda) + \tilde\alpha(\lambda)\,r_2(\lambda) + \tilde\beta(\lambda)\,r_3(\lambda)
$$

where $r_2(\lambda) = f_\lambda(\lambda/2)/f_\lambda(\lambda)$ and $r_3(\lambda) = f_\lambda(\lambda/3)/f_\lambda(\lambda)$ are **known** from the intrinsic source spectrum.

#### How Many Sources Are Needed?

> "In conclusion, at a fixed wavelength $\bar\lambda$, there are **three unknown variables**: $k(\bar\lambda)$, $\tilde\alpha(\bar\lambda)$, and $\tilde\beta(\bar\lambda)$. Therefore, **if both $S_\lambda(\lambda)$ and $f_\lambda(\lambda)$ are known for at least three independent sources**, a system of equations at each wavelength can be solved for $k(\lambda)$, $\tilde\alpha(\lambda)$, and $\tilde\beta(\lambda)$."

**Answer: minimum 3 sources.** Each source provides one equation at each wavelength. With 3 sources, the 3×3 linear system is solvable. In practice Parlanti et al. use many more than 3 (all available calibration program and archival datasets), then average over all combinations of 3.

#### How They Solved It

> "The functions $k(\lambda)$, $\tilde\alpha(\lambda)$, and $\tilde\beta(\lambda)$ were recovered for the G140M/F100LP and G235M/F170LP by exploiting all the datasets for which we have both the extended spectrum $S_\lambda(\lambda)$ from our reduction and the intrinsic spectrum $f_\lambda(\lambda)$ obtained by reducing the three gratings within their nominal wavelength range."

> "For each combination of three datasets, we computed the three functions. To account for the presence of some bad pixels and spikes in these calibration functions, we averaged and smoothed them over a window of **40 channels** along the wavelength axis."

**Calibration sources used:**
- PIDs 6645, 1536, 1537, 1538 (standard stars: G-type (P330-E), A-type, white dwarf)
- PIDs 2564 and 2186 (bright QSOs and local galaxies)
- For G140M validation: AGN 1 (SDSSJ0749, PID 2654) and AGN 2 (SDSSJ0841, PID 2654)
- For G235M validation: three local ULIRGs (PID 2186, incl. UGC-5101)

#### Calibration Accuracy Results

> "Figure 4 highlights that the calibration described above achieves an accuracy of ~5% and ~10%, respectively, for the G235M/F170LP and the G140M/F100LP configurations. These values are well within the current NIRSpec flux calibration uncertainties and discrepancies between different gratings (of the order of ~20–25%)."

**Correction formula (Eq. 2):**

$$
f_\lambda(\lambda) = \frac{S_\lambda(\lambda) - \tilde\alpha(\lambda)\,f_\lambda(\lambda/2) - \tilde\beta(\lambda)\,f_\lambda(\lambda/3)}{k(\lambda)}
$$

with $\tilde\alpha(\lambda) = \tilde\beta(\lambda) = 0$ for $\lambda < 2\lambda_\mathrm{min}$ and $\lambda < 3\lambda_\mathrm{min}$ respectively; $k(\lambda) = 1$ for $\lambda \le \lambda_\mathrm{max}^\mathrm{nom}$.

#### Contamination Amplitude

> "We note that the $\tilde\alpha$ and $\tilde\beta$ functions...are of the order of 1–5% at their maximum value, that is the range between 2–2.3 µm for the G140M and 3.2–4 µm for the G235M."

---

### 2.3 Properties of the Extended Spectra

**RMS increase:** The RMS of the extended region is higher than the nominal range, increasing from 1× (at the overlap between consecutive gratings) up to **3×** at the reddest wavelengths of the extension. Causes:
1. Decrease in disperser transmission outside nominal range (G140M drops from 0.9 at ~1 µm to ~0.2 at 3.5 µm; G235M drops from peak 0.9 at ~2.1 µm to 0.3 at 5 µm)
2. Increase in spectral resolution (more resolution elements per same wavelength interval → lower S/N per channel)

**Spectral resolution:** Continues to increase with wavelength beyond the nominal range, following Jakobsen:2022 trend. At the reddest wavelengths, the extended G140M reaches R~2500 — comparable to the high-resolution G395H grating. Parlanti et al. adopt the Shajib:2025 in-flight resolution calibration (>10% better than pre-launch estimates).

**Wavelength calibration:** Verified by comparing centroids of emission lines detected in both extended and nominal grating. Discrepancy is < 1/4–1/10 of the spectral resolution (effectively zero systematic offset).

---

### 2.4 Different Applications

Best use cases:
- **Deep observations** where long exposures in a blue filter are needed, with bright lines in the red filter coming "for free"
- **Dust-obscured sources** where shorter wavelengths are heavily attenuated
- **Bright sources** (QSOs, nearby galaxies) where the RMS increase is not limiting

Less suitable:
- Detection of **faint emission lines** specifically in the extended region (better to use the nominal filter covering those wavelengths)

---

## Key Quantitative Summary for Calibration

| Parameter | Value / Note |
|-----------|-------------|
| Min. sources to solve for k, α̃, β̃ | **3** (one equation per source per wavelength) |
| Source criteria | ≥2 consecutive gratings observed; S/N > 5 in continuum |
| Consecutive grating pairs | G140M+G235M or G235M+G395M |
| Smoothing window | 40 wavelength channels |
| α̃, β̃ amplitude (G140M) | 1–5%, max near 2.0–2.3 µm |
| α̃, β̃ amplitude (G235M) | 1–5%, max near 3.2–4.0 µm |
| Calibration accuracy (G140M) | ~10% |
| Calibration accuracy (G235M) | ~5% |
| RMS increase in extended region | 1× at boundary → 3× at reddest λ |
| Pipeline version | calwebb 1.18.1, jwst_1364.pmap |
| IFU gap (G140M) | 2.17–2.28 µm |
| IFU gap (G235M) | 3.66–3.82 µm |

---

## Comparison with Our Approach (FS Mode, PID 1492)

### Mode Difference
Parlanti et al. work with **IFU** ("IFS") mode exclusively. We are working with **Fixed Slit (FS)** mode. The core model (the $S = k f + \tilde\alpha f(\lambda/2) + \tilde\beta f(\lambda/3)$ equation) is identical and mode-independent (it describes the optical physics of spectral orders overlapping on the detector). The reference file modifications are mode-specific, however.

### Critical Gap: Number of Sources

**Parlanti use:** Many calibration sources (standard stars from PIDs 6645, 1536, 1537, 1538 + archival bright sources) → brute-force solution of 3×3 system at each wavelength, averaged over all C(N,3) combinations of N sources.

**Our v2 approach:** Single source (IRAS 05248-7007). We compensated by treating k, α̃, β̃ as smooth Legendre polynomials (5-coefficient expansion each = 15 free parameters) and fitting to the residuals of one source. This approach is **fundamentally different and under-constrained for separating k from α̃, β̃** — the three functions are entangled unless multiple sources with different spectral shapes are available.

The physical reason: with one source, $r_2(\lambda) = f(\lambda/2)/f(\lambda)$ and $r_3(\lambda) = f(\lambda/3)/f(\lambda)$ are fixed. We have one equation with three unknowns at each λ. We can only solve if we impose strong priors (e.g., smoothness). With 3+ sources having different spectral shapes (different $r_2$, $r_3$), the system becomes fully determined at each wavelength independently — no smoothness prior required.

### What We Need for a Proper Calibration

For FS mode (analogous to Parlanti Appendix B / Sec. 2.2):

1. **At least 3 sources** observed with both G140M and G235M (for G140M extended calibration)
2. **At least 3 sources** observed with both G235M and G395M (for G235M extended calibration)
3. Each source must have S/N > 5 in continuum across the full wavelength range
4. The intrinsic spectrum $f(\lambda)$ is taken from the nominal-range grating covering the extended-wavelength region (e.g., G235M nominal covers 1.7–3.15 µm, which overlaps with G140M extended 1.88–3.55 µm)

With PID 1492 (IRAS 05248-7007), we have exactly 1 source. To build a proper FS calibration, we need to identify 2+ additional suitable archival FS observations with consecutive grating coverage and high S/N.

### F-flat Modification

Parlanti et al. explicitly **concatenate the F-flat FAST_VARIATION table** from consecutive gratings. This step has **not yet been implemented** in our FS pipeline work. The F-flat correction accounts for telescope + fore-optics losses, and without extending it, flux cal in the extended region is missing this correction.

For FS mode, the F-flat structure is the same as for IFU (F-flat does not depend on the observing mode). This should be straightforward to apply.

### S-flat Modification (IFU vs FS)

In IFU mode, the S-flat traces on the detector have slight curvature/tilt; Parlanti fit a linear polynomial to the trace edges to extend them. For FS mode, the fixed slits produce simpler, narrower traces on the detector — the same principle applies but the geometry is different. For FS, the simplest approach (extending the SCI pixels to 1 and removing DQ flags beyond the nominal range) is likely sufficient, as done in their NRS2 procedure.

### Reference File Status (Our Work vs Parlanti)

| Reference File | Parlanti (IFU) | Our FS Work |
|----------------|----------------|-------------|
| wavelengthrange | Extended for cube_build + nirspec_wavelength_range | ✅ `wavelengthrange_extended.asdf` created |
| S-flat NRS1 SCI/DQ | Trace-extended with linear polynomial | Partial (DQ flags removed; geometry not trace-fitted) |
| S-flat NRS2 DQ | Flags removed, FAST_VARIATION extended | Partially done |
| S-flat FAST_VARIATION | Extended with constant after nominal λmax | Not yet confirmed |
| F-flat FAST_VARIATION | Concatenated across consecutive gratings | **Not yet done** |
| nirspec_wavelength_range | Modified for cube_build | Not applicable (FS) |
| calwebb_spec2 NRS2 check | Removed | Applied for IFU; for FS check needed |

### Pipeline Code Modification Status

| Change | Parlanti | Our FS Work |
|--------|---------|------------|
| Remove NRS2 error in assign_wcs (IFU) | ✅ | N/A (IFU only) |
| Extended wavelength range | ✅ | ✅ |
| Remove calwebb_spec2 NRS2 file check | ✅ | To verify for FS |

---

## Summary: What Is Needed for Proper FS Calibration

Based on Parlanti et al., the minimum requirements to properly solve for k(λ), α̃(λ), β̃(λ) in FS mode are:

### For G140M extended calibration (1.88–3.55 µm):
- ≥3 FS sources, each observed with **both G140M/F100LP and G235M/F170LP**
- For each: know both $S_\lambda(\lambda)$ (G140M NRS2 extended) and $f_\lambda(\lambda)$ (G235M nominal covers 1.7–3.15 µm → used as intrinsic reference over the overlap region)
- Sources should have diverse spectral shapes (different $r_2$, $r_3$ ratios)

### For G235M extended calibration (3.15–5.27 µm):
- ≥3 FS sources, each observed with **both G235M/F170LP and G395M/F290LP**
- For each: know $S_\lambda(\lambda)$ (G235M NRS2 extended) and $f_\lambda(\lambda)$ (G395M nominal 2.87–5.09 µm)

### PID 1492 Status: 1 source, all 4 gratings observed
- IRAS 05248-7007 supplies **all 4** nominal gratings (G140M + G235M + G395M + PRISM) → we have the necessary $f(\lambda)$
- But single-source → must use polynomial/regularized fitting (priors on smoothness of k, α̃, β̃), not the direct 3-source linear solve
- This produces approximate calibration suitable for initial testing, but not the instrument-level calibration required for public release

### Practical Path Forward for Multi-Source Calibration

Identify archival JWST NIRSpec FS observations with:
1. Both G140M/F100LP **and** G235M/F170LP of the same target
2. OR: Both G235M/F170LP **and** G395M/F290LP of the same target
3. Bright continuum sources (S/N > 5)
4. Standard star programs (equivalents of PIDs 6645, 1536–1538 for FS mode) are ideal

The FS equivalent of Parlanti's standard star programs would be the flux calibration programs (e.g., PID 1537: white dwarf P177D; PID 1538: G-type P330-E) — these were observed across grating configurations and provide the ideal multi-source calibration dataset.
