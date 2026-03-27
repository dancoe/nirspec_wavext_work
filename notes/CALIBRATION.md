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

To determine the unknown throughput functions $k$, $\tilde\alpha$, and $\tilde\beta$, we require a minimum of **three sources** with known intrinsic spectra (Parlanti et al. 2025, Sec. 2.2).

**Why 3 sources minimum?** At each wavelength $\bar\lambda$, we have exactly 3 unknowns: $k(\bar\lambda)$, $\tilde\alpha(\bar\lambda)$, $\tilde\beta(\bar\lambda)$. Each source provides one equation (since $r_2 = f(\lambda/2)/f(\lambda)$ and $r_3 = f(\lambda/3)/f(\lambda)$ differ between sources with different spectral shapes). Three sources → 3 independent equations → unique solution at each wavelength.

**Formal statement (Parlanti et al.):**
> "At a fixed wavelength $\bar\lambda$, there are three unknown variables: $k(\bar\lambda)$, $\tilde\alpha(\bar\lambda)$, and $\tilde\beta(\bar\lambda)$. Therefore, if both $S_\lambda(\lambda)$ and $f_\lambda(\lambda)$ are known for at least three independent sources, a system of equations at each wavelength can be solved for $k(\lambda)$, $\tilde\alpha(\lambda)$, and $\tilde\beta(\lambda)$."

**Source selection criteria (Parlanti et al.):**
1. Have observations with **at least two consecutive medium resolution gratings** (G140M+G235M or G235M+G395M)
2. Bright enough to have **continuum detection with S/N > 5** across all wavelength ranges covered

**Solution procedure (Parlanti et al.):**
1. Gather all N qualifying sources
2. For each combination of 3 sources, solve the 3×3 linear system at each wavelength
3. Average all C(N,3) solutions and **smooth over a 40-channel window** to reject bad pixels/spikes

**Pipeline steps:**
1. **Preprocessing:** Extract data using custom `wavelengthrange` reference file to extended wavelengths
2. **Reference file modification:** Extend S-flat (SCI→1, DQ flags removed) and F-flat (concatenate FAST_VARIATION tables from consecutive gratings)
3. **Analysis:** Assemble S(λ) from extended extraction and f(λ) from nominal grating at overlapping wavelengths
4. **Solving:** Solve 3×3 system (or use regularised polynomial fit if N < 3) for $k$, $\tilde\alpha$, $\tilde\beta$
5. **Applying:** Use Eq. 2 to recover calibrated flux in extended region

### Important Note: Single-Source Limitation

With **only 1 source** (as in our PID 1492 single-target work), the 3-unknown system is under-constrained. The v2 polynomial fitting approach (Legendre + Tikhonov regularisation) compensates by imposing smoothness priors on k, α̃, β̃. This produces an *approximate* calibration but **cannot independently separate** the three functions — the solution is non-unique and depends on the regularisation strength. This is why our v2 RMS is ~50%, far worse than Parlanti's ~5–10%.

### Priority
The initial calibration efforts will focus on **Fixed Slit (FS)** data from Program ID 1492, using the single-source approximate approach. A proper multi-source calibration requires identifying additional archival FS observations (see "Path Forward" below and PARLANTI.md).

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

1. **Single source calibration:** With one source we are fitting and validating against the same data — this is a circularity. **Parlanti et al. require ≥3 sources with diverse spectral shapes** to solve the 3-unknown system independently at each wavelength (see "Calibration Strategy" above). Our v2 polynomial approach is an approximation that works only because we impose strong smoothness priors — it cannot robustly separate k from α̃, β̃.
2. **Missing F-flat modification:** Parlanti et al. concatenate the F-flat FAST_VARIATION tables from consecutive gratings (G140M F-flat extended with G235M and G395M entries up to 5.5 µm). This step has **not yet been performed** for our FS work. Without it, the extended region is missing the fore-optics flat-field correction.
3. **Missing photom step:** The extended extractions lack absolute flux calibration. The boundary-scaling is an approximation. Next step: run the `photom` step on the NRS2 extractions, even if only approximate, to obtain Jy units directly.
4. **b(λ) clamped:** The 3rd-order model may be under-constrained with one source. Consider setting the upper bound to 10% or using a stronger regularisation for $b$.
5. **PRISM resolution limitation:** The PRISM f(λ) is low-resolution; using it to subtract emission-line contamination leaves residuals. The v2 script uses the high-resolution nominal M-grating spectra for the subtraction step, which partially addresses this.
6. **G395M:** Only 68 points at 5.5–5.6 µm and no PRISM baseline available — calibration not yet possible.
7. **Obs001 vs Obs003:** The coefficients were derived using obs001 PRISM/nominals as f(λ) and obs003 NRS2 as S(λ). Both are of IRAS 05248-7007 — but different epochs/dithers. Any pointing offset or pointing-dependent slit losses could introduce systematic errors.

### Path Forward: Multi-Source Calibration

To replicate Parlanti et al.'s approach for FS mode, we need to identify additional archival JWST/NIRSpec FS targets observed in consecutive grating pairs:

**For G140M extended cal (1.88–3.55 µm):**
- Need ≥2 additional FS sources observed with **both G140M/F100LP and G235M/F170LP**
- Ideal: JWST flux calibration standard star programs (PIDs 1537, 1538) which observed standard stars in all FS grating configs

**For G235M extended cal (3.15–5.27 µm):**
- Need ≥2 additional FS sources observed with **both G235M/F170LP and G395M/F290LP**
- Same standard star programs apply

**Source quality requirement:** S/N > 5 in continuum across all wavelengths covered (Parlanti et al. criterion)

See [PARLANTI.md](PARLANTI.md) for comprehensive analysis of their methodology and detailed comparison with our approach.

See also: `plots/Parlanti/cal/PARLANTI_ANALYSIS_REPORT.md` (v1 analysis) and `plots/Parlanti/cal/parlanti_coefficients_v2.txt`.
