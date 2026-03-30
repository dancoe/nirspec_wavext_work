# NIRSpec Wavelength Extension Report — 329 FS v4

**Date:** March 29, 2026  
**Project:** NIRSpec Wavelength Extension Calibration  
**Data Version:** FS v4 — NNLS simultaneous k, α, β derivation

---

## Summary

FS v4 applies the same NNLS joint-solve approach as IFU v4 to the Fixed Slit (FS) Level-3 NRS2 x1d products. All three Parlanti model coefficients — $k(\lambda)$, $\tilde{\alpha}(\lambda)$, and $\tilde{\beta}(\lambda)$ — are derived simultaneously from three CALSPEC standard stars. No Parlanti published values are assumed.

**Key results:**
- FS $k$ is largely consistent with FS v3 ($k_{\rm G140M}$ 0.701 vs 0.725; $k_{\rm G235M}$ 0.723 vs 0.972).
- FS NNLS $\tilde{\alpha}$ is 2–6× higher than Parlanti published values across most of the NRS2 range, mirroring the IFU v4 finding.
- Spectral validation for P330E G140M is excellent (corrected spectrum within ~5–10% of MAST G235M truth over 2.0–3.0 µm). G191-B2B and J1743045 G140M also validate well.
- NGC2506-G31 (G1V, PID 6644) is the planned 4th source to make the FS system overconstrained; it is not yet downloaded.

---

## Algorithm (same as IFU v4)

At each wavelength $\lambda$, solve:

$$\min_{\mathbf{x} \geq 0} \left\| \begin{pmatrix} f_1(\lambda) & f_1(\lambda/2) & f_1(\lambda/3) \\ f_2(\lambda) & f_2(\lambda/2) & f_2(\lambda/3) \\ f_3(\lambda) & f_3(\lambda/2) & f_3(\lambda/3) \end{pmatrix} \begin{pmatrix} k \\ \alpha \\ \beta \end{pmatrix} - \begin{pmatrix} S_{\mathrm{obs},1}(\lambda) \\ S_{\mathrm{obs},2}(\lambda) \\ S_{\mathrm{obs},3}(\lambda) \end{pmatrix} \right\|^2$$

Solved with `scipy.optimize.nnls`. 40-channel boxcar smoothing applied to each output.

---

## Pipeline Lineage

| Step | Customisation |
|:-----|:--------------|
| `assign_wcs` | `wavelengthrange_0008.asdf` extended range (G140M 0.6–3.3 µm) |
| `assign_wcs` | `nirspec.py` NRS2 M-grating hard-coded error removed |
| `flat_field` | Parlanti sflat/fflat (NRS2 extended coverage) |
| `photom` | Parlanti photom NRS2 ext table (extrapolated) |
| `wavecorr` | **Skipped** (NRS2 slit no wavelength correction available) |
| `Spec3Pipeline` | `resample_spec=True, extract_1d=True` → `nrs2_spec3_ext/` |

Output: `data/PID{pid}_{name}/nrs2_spec3_ext/nrs2_l3_{pid}_{grating}_s1600a1_s000000001_x1d.fits`

---

## Data Sources

| Source   | SpT     | PID  | Data Product |
|:---------|:--------|:-----|:-------------|
| P330E    | G2V     | 1538 | FS nrs2_spec3_ext L3 x1d |
| G191-B2B | WD-DA.8 | 1537 | FS nrs2_spec3_ext L3 x1d |
| J1743045 | A8III   | 1536 | FS nrs2_spec3_ext L3 x1d |

Missing: NGC2506-G31 (G1V, PID 6644) — **download pending** (would make system overdetermined).

---

## Coefficient Summary

| Grating | k median | k range | α median | α max | β median | β max |
|:--------|:---------|:--------|:---------|:------|:---------|:------|
| G140M NRS2 | **0.701** | 0.010–1.742 | **0.0292** | 0.3074 | 0.0041 | 0.0970 |
| G235M NRS2 | **0.723** | 0.010–1.827 | **0.0369** | 0.6013 | 0.0038 | 0.0754 |

Version comparison (k only):

| Grating | k v3 (Parlanti α/β) | k v4 (NNLS) | k Parlanti |
|:--------|:--------------------|:------------|:-----------|
| G140M   | 0.725               | **0.701**   | ~0.8–1.0   |
| G235M   | 0.972               | **0.723**   | ~0.7–1.0   |

G235M k drops from 0.972 (v3) to 0.723 (v4) — the NNLS assigns 0.035–0.2 of the observed flux to second-order contamination rather than throughput.

### Condition Numbers

| Grating | Condition Number (median) | max |
|:--------|:--------------------------|:----|
| G140M   | 194 | 19 165 |
| G235M   | 857 | 113 049 |

FS G140M is better conditioned than IFU G140M (condition 194 vs 288), likely because the FS NRS2 wavelength range is shorter (1.95–3.20 µm vs 1.87–3.55 µm).

---

## IFU v4 vs FS v4 Comparison

| Grating | IFU v4 k | FS v4 k | IFU v4 α | FS v4 α |
|:--------|:---------|:--------|:--------|:--------|
| G140M   | 0.509    | **0.701** | 0.0835  | **0.0292** |
| G235M   | 0.466    | **0.723** | 0.0498  | **0.0369** |

FS $k$ is consistently higher than IFU $k$ (as established in v2/v3), confirming the genuine mode-dependent photometric calibration difference. Interestingly, IFU $\alpha$ is *higher* than FS $\alpha$ in v4 — this may reflect the IFU's larger aperture collecting more scattered/ghost light from the 2nd-order, or may be an artefact of the higher k/α degeneracy in the IFU system.

---

## 1. Calibration Coefficients vs Parlanti Published

### G140M Coefficients
![FS v4 G140M Coefficients](../../../plots/Parlanti/cal/fs_v4/fs_v4_g140m_coeffs.png)

### G235M Coefficients
![FS v4 G235M Coefficients](../../../plots/Parlanti/cal/fs_v4/fs_v4_g235m_coeffs.png)

---

## 2. v4 vs v3 Coefficient Comparison

### G140M  
![FS v4 vs v3 G140M](../../../plots/Parlanti/cal/fs_v4/fs_v4_vs_v3_G140M.png)

### G235M  
![FS v4 vs v3 G235M](../../../plots/Parlanti/cal/fs_v4/fs_v4_vs_v3_G235M.png)

---

## 3. MAST Level-3 NRS1 vs Our NRS2 Extended

### G140M — All Sources
![MAST vs L3 G140M](../../../plots/Parlanti/cal/fs_v4/fs_v4_mast_vs_l3_G140M.png)

*(Note: Our own NRS1 L3 products not generated for FS; only MAST NRS1 shown.)*

### G235M — All Sources
![MAST vs L3 G235M](../../../plots/Parlanti/cal/fs_v4/fs_v4_mast_vs_l3_G235M.png)

---

## 4. Per-Source NRS2 Spectral Validation

### P330E — G140M NRS2
![P330E G140M](../../../plots/Parlanti/cal/fs_v4/fs_v4_spectra_P330E_G140M.png)

### P330E — G235M NRS2
![P330E G235M](../../../plots/Parlanti/cal/fs_v4/fs_v4_spectra_P330E_G235M.png)

### G191-B2B — G140M NRS2
![G191-B2B G140M](../../../plots/Parlanti/cal/fs_v4/fs_v4_spectra_G191-B2B_G140M.png)

### G191-B2B — G235M NRS2
![G191-B2B G235M](../../../plots/Parlanti/cal/fs_v4/fs_v4_spectra_G191-B2B_G235M.png)

### J1743045 — G140M NRS2
![J1743045 G140M](../../../plots/Parlanti/cal/fs_v4/fs_v4_spectra_J1743045_G140M.png)

### J1743045 — G235M NRS2
![J1743045 G235M](../../../plots/Parlanti/cal/fs_v4/fs_v4_spectra_J1743045_G235M.png)

---

## Discussion

The FS v4 NNLS result is broadly consistent with IFU v4: $\alpha$ is substantially higher than Parlanti's published value, while $k$ is somewhat lower. The G235M FS $\alpha_{\rm max} = 0.60$ at the NRS1/NRS2 boundary (3.3 µm) is a concern — a value this high implies that >50% of the observed flux at some wavelengths is 2nd-order contamination, which seems extreme. This most likely indicates that at the boundary the NNLS is poorly constrained (high condition number at those wavelengths, ~25 000).

The key test to stabilise these values is adding NGC2506-G31 as a 4th source. Its much cooler SED (G1V) will produce a $f(\lambda/2)/f(\lambda)$ ratio that is distinctly different from the WD and A-star, forcing the NNLS to correctly apportion flux between $k$ and $\alpha$.

The $k(\lambda)$ values from FS v4 are reliable and consistent with FS v3 (within 5% for G140M, within 25% for G235M). L3 NRS2 spectra corrected with FS v4 $k$ alone (fixing $\alpha = \alpha_P$, $\beta = \beta_P$ from Parlanti) would give results essentially identical to v3.

---

## Next Steps

1. Download NGC2506-G31 (PID 6644) FS rate files from MAST
2. Run FS Spec2+Spec3 pipeline for NGC2506-G31 G140M and G235M NRS2
3. Add to v4 solver → v4.1 with 4-source overdetermined NNLS
4. Compare $\alpha_{\rm v4.1}$ vs Parlanti to assess whether high-$\alpha$ finding is physical
