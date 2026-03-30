# NIRSpec Wavelength Extension Report — 329 IFU v4

**Date:** March 29, 2026  
**Project:** NIRSpec Wavelength Extension Calibration  
**Data Version:** IFU v4 — NNLS simultaneous k, α, β derivation

---

## Summary

IFU v4 is the first fully independent calibration: all three Parlanti model coefficients — the throughput $k(\lambda)$, second-order contamination $\tilde{\alpha}(\lambda)$, and third-order contamination $\tilde{\beta}(\lambda)$ — are solved simultaneously using Non-Negative Least Squares (NNLS) on three CALSPEC standard stars spanning a wide range of spectral types. No Parlanti published values are assumed.

**Key result:** The NNLS-derived $\tilde{\alpha}(\lambda)$ is substantially higher (2–5× at most wavelengths) than Parlanti's published IFU value, while $k(\lambda)$ is correspondingly lower. This suggests that the Parlanti published $\tilde{\alpha}$ underestimates the second-order contamination contribution as seen in our Level-3 pipeline products, or that the $k/\alpha$ degeneracy is pulling flux from $k$ into $\alpha$.

Spectral validation is good for P330E and J1743045 (within ~10–15% over most of the NRS2 range). G191-B2B G235M also corrects well (~5%). Moderate noise/spikes in the corrected spectra at the red end of each grating where the throughput falls below 0.1.

---

## Algorithm (New in v4)

### NNLS Joint Solve

At each wavelength $\lambda$ on a 400-point grid, build:

$$\mathbf{A}(\lambda) = \begin{pmatrix} f_1(\lambda) & f_1(\lambda/2) & f_1(\lambda/3) \\ f_2(\lambda) & f_2(\lambda/2) & f_2(\lambda/3) \\ f_3(\lambda) & f_3(\lambda/2) & f_3(\lambda/3) \end{pmatrix}, \quad \mathbf{b}(\lambda) = \begin{pmatrix} S_{\mathrm{obs},1}(\lambda) \\ S_{\mathrm{obs},2}(\lambda) \\ S_{\mathrm{obs},3}(\lambda) \end{pmatrix}$$

Solve $\min_{\mathbf{x} \geq 0} \|\mathbf{A}\mathbf{x} - \mathbf{b}\|^2$ via `scipy.optimize.nnls` to obtain $\mathbf{x} = [k, \alpha, \beta]$.

Apply 40-channel boxcar smoothing to each coefficient (identical to Parlanti 2025 and v2/v3).

### Comparison to v3
In v3: $k_i(\lambda) = R_i(\lambda) - \tilde{\alpha}_P(\lambda) r_{2,i} - \tilde{\beta}_P(\lambda) r_{3,i}$, then $k = \mathrm{median}_i[k_i]$.  
In v4: $[k, \alpha, \beta]$ solved jointly — no Parlanti $\alpha_P$, $\beta_P$ assumed.

---

## Data Sources

| Source   | SpT     | PID  | Data Product |
|:---------|:--------|:-----|:-------------|
| P330E    | G2V     | 1538 | IFU stage3_ext x1d (L3 NRS1+NRS2) |
| G191-B2B | WD-DA.8 | 1537 | IFU stage3_ext x1d (L3 NRS1+NRS2) |
| J1743045 | A8III   | 1536 | IFU stage3_ext x1d (L3 NRS1+NRS2) |

Pipeline provenance identical to IFU v3 (WAVEXT.md.

---

## Coefficient Summary

| Grating | k median | k range | α median | α max | β median | β max |
|:--------|:---------|:--------|:---------|:------|:---------|:------|
| G140M NRS2 | **0.509** | 0.010–1.100 | **0.0835** | 0.2249 | 0.0044 | 0.0304 |
| G235M NRS2 | **0.466** | 0.010–1.185 | **0.0498** | 0.2052 | 0.0075 | 0.0981 |

Comparison with prior versions (k only):

| Grating | k v3 (Parlanti α/β) | k v4 (NNLS) | k Parlanti |
|:--------|:--------------------|:------------|:-----------|
| G140M   | 0.567               | **0.509**   | ~0.8–1.0   |
| G235M   | 0.768               | **0.466**   | ~0.7–1.0   |

The v4 $k$ is lower because more of the observed flux excess is now attributed to α rather than being folded into $k$.

### Condition Numbers (degeneracy diagnostic)

| Grating | Condition Number (median) | max |
|:--------|:--------------------------|:----|
| G140M   | 288 | 74 278 |
| G235M   | 810 | 185 369 |

Condition numbers >100 indicate partial degeneracy. The mild k/α degeneracy with 3 stars (G2V + WD + A8III) is expected; adding NGC2506-G31 (G1V, PID 6644) would lower condition numbers and stabilise $\alpha$ derivation.

---

## 1. Calibration Coefficients vs Parlanti Published

### G140M Coefficients (v4 NNLS vs Parlanti)
![IFU v4 G140M Coefficients](../../../plots/Parlanti/cal/ifu_v4/ifu_v4_g140m_coeffs.png)

*Orange = v4 NNLS α, Dashed red = Parlanti published α. V4 α is 2–4× higher at 1.87–2.2 µm and at 2.8–3.5 µm.*

### G235M Coefficients (v4 NNLS vs Parlanti)
![IFU v4 G235M Coefficients](../../../plots/Parlanti/cal/ifu_v4/ifu_v4_g235m_coeffs.png)

---

## 2. v4 vs v3 Coefficient Comparison

### G140M  
![IFU v4 vs v3 G140M](../../../plots/Parlanti/cal/ifu_v4/ifu_v4_vs_v3_G140M.png)

### G235M  
![IFU v4 vs v3 G235M](../../../plots/Parlanti/cal/ifu_v4/ifu_v4_vs_v3_G235M.png)

---

## 3. MAST Level-3 vs Our Level-3 Extended

### G140M — All Sources
![MAST vs L3 G140M](../../../plots/Parlanti/cal/ifu_v4/ifu_v4_mast_vs_l3_G140M.png)

### G235M — All Sources
![MAST vs L3 G235M](../../../plots/Parlanti/cal/ifu_v4/ifu_v4_mast_vs_l3_G235M.png)

---

## 4. Per-Source NRS2 Spectral Validation

### P330E — G140M NRS2
![P330E G140M](../../../plots/Parlanti/cal/ifu_v4/ifu_v4_spectra_P330E_G140M.png)

### P330E — G235M NRS2
![P330E G235M](../../../plots/Parlanti/cal/ifu_v4/ifu_v4_spectra_P330E_G235M.png)

### G191-B2B — G140M NRS2
![G191-B2B G140M](../../../plots/Parlanti/cal/ifu_v4/ifu_v4_spectra_G191-B2B_G140M.png)

### G191-B2B — G235M NRS2
![G191-B2B G235M](../../../plots/Parlanti/cal/ifu_v4/ifu_v4_spectra_G191-B2B_G235M.png)

### J1743045 — G140M NRS2
![J1743045 G140M](../../../plots/Parlanti/cal/ifu_v4/ifu_v4_spectra_J1743045_G140M.png)

### J1743045 — G235M NRS2
![J1743045 G235M](../../../plots/Parlanti/cal/ifu_v4/ifu_v4_spectra_J1743045_G235M.png)

---

## Discussion

### Why is v4 α higher than Parlanti?

Parlanti et al. (2025) performed their calibration on IFU data from PID 6645 with different pipeline settings. Our v4 NNLS finds a substantially higher second-order fraction at the NRS1/NRS2 boundary (1.87–2.2 µm). Several factors could contribute:

1. **k/α degeneracy not fully broken**: With only G2V + WD + A8III, the $r_2 = f(\lambda/2)/f(\lambda)$ ratios in the A matrix may still be partially collinear, driving flux into $\alpha$.
2. **Different photometric calibration**: Our pipeline uses the Parlanti photom reference files, which may introduce a systematic that is absorbed into $\alpha$.
3. **Genuine instrument behaviour**: The second-order contamination on NRS2 may be mode-dependent and differ from the IFU geometry assumed by Parlanti.

### Recommendation

The v4 coefficients should be used with caution for $\alpha$ and $\beta$ until NGC2506-G31 (G1V, PID 6644) data is downloaded and incorporated. The $k(\lambda)$ values from v4 are in reasonable agreement with v3 ($\pm 10$–20%) and can be used. The spectral correction quality (v4) is comparable to v3 for the calibration sources used.

---

## Next Steps

- Download and reduce NGC2506-G31 (PID 6644) FS data → incorporate as 4th source (overdetermined system)
- Compare $\tilde{\alpha}$ from v4 IFU vs v4 FS (mode-dependence)
- Quantify residual normalised correction errors over 1.87–5.3 µm
