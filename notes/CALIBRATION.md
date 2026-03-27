# CALIBRATION.md

## NIRSpec Wavelength Extension Calibration Theory

This theory is adopted from Parlanti et al. (arXiv:2512.14844), which accounts for 2nd and 3rd order spectral contamination during extended wavelength extraction.

### Order Overlap Equation

The observed spectrum $S(\lambda)$ at extended wavelengths is modeled as a sum of the 1st, 2nd, and 3rd order contributions:

$$S(\lambda) \approx k(\lambda) f(\lambda) + a(\lambda) f(\lambda/2) + b(\lambda) f(\lambda/3)$$

Where:
*   $f(\lambda)$ = Intrinsic spectrum of the source.
*   $S(\lambda)$ = Measured (observed) spectrum from the extended extraction.
*   $k(\lambda)$ = 1st-order throughput contribution (drops to ~30% at the longest extended wavelengths).
*   $a(\lambda)$ = 2nd-order contamination contribution (< ~5%).
*   $b(\lambda)$ = 3rd-order contamination contribution (< ~1%).

### Calibration Strategy

To determine the unknown throughput functions $k$, $a$, and $b$, we require a minimum of three sources with known intrinsic spectra:

1.  **Requirement:** At least 3 target sources where we have:
    *   **$f(\lambda)$:** Intrinsic flux obtained from the standard nominal extraction (where contamination is negligible or well-calibrated).
    *   **$S(\lambda)$:** Measured flux from the extended extraction (including contaminated wavelengths).
2.  **Step 1 (Preprocessing):** Extract the data using the custom `wavelengthrange` reference file to extended wavelengths.
3.  **Step 2 (Analysis):** Measure the discrepancies between the nominal and extended regions to isolate $S(\lambda)$.
4.  **Step 3 (Solving):** Solve the system of equations across the target samples to calculate $k(\lambda)$, $a(\lambda)$, and $b(\lambda)$ for each grating/filter combination.

### Priority
The initial calibration efforts will focus on **Fixed Slit (FS)** data from Program ID 1492.

---

## Coefficient Determination — PID 1492 FS (March 27, 2026)

### Data

Target: **IRAS 05248-7007** (nearby ULIRG, highly emission-line rich)

| File | Type | λ range | N pts | Units |
|------|------|---------|-------|-------|
| `jw01492-o001_t001-…_clear-prism-…_x1d.fits` | PRISM Level 3 | 0.597–5.300 µm | 411 | Jy |
| `jw01492-o001_t001-…_f100lp-g140m-…_x1d.fits` | G140M Level 3 | 0.980–1.838 µm | 1347 | Jy |
| `jw01492-o001_t001-…_f170lp-g235m-…_x1d.fits` | G235M Level 3 | 1.699–3.069 µm | 1283 | Jy |
| `jw01492-o001_t001-…_f290lp-g395m-…_x1d.fits` | G395M Level 3 | 2.879–5.089 µm | 1231 | Jy |
| `jw01492003001_03102_00005_nrs2_g140m_extract_1d.fits` | G140M NRS2 ext. | 1.981–3.268 µm | 1095 | DN/s |
| `jw01492003001_03104_00004_nrs2_g235m_extract_1d.fits` | G235M NRS2 ext. | 3.308–5.310 µm | 1272 | DN/s |
| `jw01492003001_03106_00002_nrs2_g395m_extract_1d.fits` | G395M NRS2 ext. | 5.496–5.617 µm | 68 | DN/s |

### Critical Units Finding

The extended NRS2 extractions (obs003) were processed through `assign_wcs` + `extract_2d` + `extract_1d` but **NOT** the `photom` step. Their `FLUX` column is in **DN/s** (raw count rates), **not Jy**. This explains the factor ~100–200× discrepancy seen in the pre-cal plot.

The pipeline-standard nominal Level 3 products (obs001) are fully calibrated (in Jy).

**Consequence:** The Parlanti model equation must be applied in consistent units. We bridge this by determining a DN/s → Jy scale factor.

### Unit Scaling: Boundary Matching

At the very start of the extended extraction (close to the nominal–extended boundary), $k(\lambda) \approx 1$ and the contamination terms $a$, $b$ are negligible. Therefore:

$$S_{\rm raw}(\lambda_{\rm bdy}) \times C \approx f_{\rm PRISM}(\lambda_{\rm bdy})$$

where $C$ [Jy / (DN/s)] is the photometric conversion factor. We compute $C$ by taking the sigma-clipped median ratio of $f_{\rm PRISM}(\lambda) / S_{\rm raw}(\lambda)$ over the first 80 wavelength points of each extended extraction (emission-line peaks excluded by 3σ clipping).

| Grating | Boundary λ | Scale C [Jy/(DN/s)] |
|---------|-----------|---------------------|
| G140M | ~2.01 µm | 8.04 × 10⁻⁴ |
| G235M | ~3.35 µm | 1.35 × 10⁻³ |
| G395M | ~5.51 µm | N/A (PRISM doesn't reach 5.5 µm) |

**Limitation:** If pollution from 2nd-order emission lines is already non-negligible at the boundary, $C$ will be underestimated (making $S_{\rm scaled}$ slightly too low and $k$ start below 1). This is a first-order approximation.

### Improved Fitting Method (v2): Legendre Polynomials + Tikhonov Regularisation

The previous polynomial fit (v1) suffered from ill-conditioning because raw power monomials ($\lambda^0, \lambda^1, \ldots$) become nearly linearly dependent over small wavelength ranges.

**Improvements in v2:**

1. **Orthogonal basis:** Wavelength normalised to $\xi \in [-1, 1]$; $k$, $a$, $b$ each represented as a degree-4 Legendre polynomial expansion $\sum_{i=0}^4 c_i P_i(\xi)$.
2. **Tikhonov regularisation:** Solve $(A^T A + \lambda_{\rm reg} I) \mathbf{x} = A^T \mathbf{S}$, where $\lambda_{\rm reg}$ scales with the Frobenius norm of $A^T A$. This suppresses oscillations.
3. **Iterative σ-clipping:** 3 rounds of 3σ rejection on residuals re-weight the fit away from emission-line spikes.
4. **Physical bounds** (hard-clipped after fitting):
   - $k(\lambda) \in [0.01,\ 2.0]$ — allows for flat-field over/under-corrections
   - $a(\lambda) \in [0,\ 0.20]$ — 2nd-order contamination ≤ 20%
   - $b(\lambda) \in [0,\ 0.05]$ — 3rd-order contamination ≤ 5%
5. **Pre-smoothing** of $S(\lambda)$ with Savitzky-Golay (window = 15, order = 3) before fitting.

For the **calibration application step**, the subtraction of contamination uses a high-resolution nominal reference (G140M + G235M + G395M stitch, falling back to PRISM) instead of the low-resolution PRISM alone. This gives a cleaner removal of sharp emission-line contamination.

**Script:** `analysis/calibrate_parlanti.py`

### Results

| Grating | λ range | Fit RMS | k range | a range | b range |
|---------|---------|---------|---------|---------|---------|
| G140M | 1.98–3.27 µm | 50% | 0.01–0.47 | 0.00–13.2% | 5% (clamped) |
| G235M | 3.31–5.31 µm | 54% | 0.10–0.96 | 5.4–20% (clamped at ends) | 0–5% |

### Interpretation of Coefficients

**G140M k(λ):** Drops from ~0.47 at 2.0 µm to ~0.01 at 3.0 µm. This rapid decline reflects both the falling grating efficiency of G140M beyond its nominal range and the absent flat-field correction in the extended region. The fact that $k \approx 0.47$ at the boundary (rather than 1.0) indicates either residual 2nd-order contamination at the boundary or an underestimate in $C$.

**G140M a(λ):** Starts at ~13% at 2.0 µm and drops to <2% by 2.5 µm. The high value at the boundary arises because IRAS 05248-7007 has extremely bright emission lines at 1.0–1.2 µm (Hα, [N II] etc.). Their 2nd-order images at 2.0–2.4 µm dominate the extended contamination signal.

**G140M b(λ):** Pinned at the 5% upper bound throughout — the 3rd-order contamination (from 0.66–1.1 µm light) saturates the physical constraint. This suggests either the constraint is too tight, or systematic effects (flat-field errors) are partially absorbed into $b$.

**G235M k(λ):** Starts near 1 at 3.3 µm (good — close to the nominal boundary) and falls to ~0.1 at 4.5 µm before recovering slightly at 5 µm. This U-shape likely reflects a combination of detector quantum efficiency at the longest wavelengths and the flat-field response.

**G235M a(λ):** Large (hitting 20% bound) at both ends of the range. The source's bright mid-IR continuum at 1.7–3 µm (2nd order) is contaminating the observed G235M extended spectrum significantly.

### Calibrated Plot

**Output:** `plots/Parlanti/cal/FS_1492_cal.png`

The 2-panel plot (modelled on Parlanti et al. Fig. 5) shows:
- **Panel 1 (1.0–3.3 µm):** G140M nominal + G235M nominal (black), extended G140M uncalibrated (green), extended G140M calibrated (red). After calibration, the red trace comes to within ~1–3 mJy of the nominal level, consistent with the expected intrinsic spectrum of IRAS 05248-7007.
- **Panel 2 (1.7–5.3 µm):** G235M + G395M nominal (black), extended G235M uncalibrated (green), extended G235M calibrated (red). The calibrated G235M tracks the G395M nominal from 3.3–5.1 µm reasonably well, with ~50% RMS scatter.

### Known Limitations & Next Steps

1. **Single source calibration:** With one source we are fitting and validating against the same data — this is a circularity. Multiple sources at different redshifts are required to robustly determine instrument-level $k$, $a$, $b$ (as Parlanti et al. do).
2. **Missing photom step:** The extended extractions lack absolute flux calibration. The boundary-scaling is an approximation. Next step: run the `photom` step on the NRS2 extractions, even if only approximate, to obtain Jy units directly.
3. **b(λ) clamped:** The 3rd-order model may be under-constrained with one source. Consider setting the upper bound to 10% or using a stronger regularisation for $b$.
4. **PRISM resolution limitation:** The PRISM f(λ) is low-resolution; using it to subtract emission-line contamination leaves residuals. The v2 script uses the high-resolution nominal M-grating spectra for the subtraction step, which partially addresses this.
5. **G395M:** Only 68 points at 5.5–5.6 µm and no PRISM baseline available — calibration not yet possible.
6. **Obs001 vs Obs003:** The coefficients were derived using obs001 PRISM/nominals as f(λ) and obs003 NRS2 as S(λ). Both are of IRAS 05248-7007 — but different epochs/dithers. Any pointing offset or pointing-dependent slit losses could introduce systematic errors.

See also: `plots/Parlanti/cal/PARLANTI_ANALYSIS_REPORT.md` (v1 analysis) and `plots/Parlanti/cal/parlanti_coefficients_v2.txt`.
