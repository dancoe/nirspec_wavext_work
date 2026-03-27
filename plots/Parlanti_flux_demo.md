# NIRSpec Spectral Overlap Analysis (Parlanti et al. 2025)

These plots demonstrate the spectral extraction performance for **PID 1492 FS** using the custom `wavelengthrange_extended.asdf` reference file.

![G140M Overlap](g140m_parlanti_overlap.png)

### G140M (Medium Resolution)
*   **Contamination Window**: The region from **1.9 µm to 3.3 µm** clearly shows the extended G140M spectrum ($S(\lambda)$) rising above the noise floor.
*   **Model Fit**: The extended red line (`alpha=0.6`, `lw=0.5`) tracking above the PRISM baseline indicates the captured 2nd order light, which is the primary focus for our throughput calibration.

![G235H Overlap](g235h_parlanti_overlap.png)

### G235H (High Resolution)
*   **Physical Cutoff**: Note that the extended spectrum (blue line) terminates abruptly and does not show significant contamination. This confirms that for high-resolution gratings, the first-order light at the red end physically misses the NRS2 detector, effectively gating the spectral overlap in this specific slit configuration.

## Equations for Throughput Derivation

Using the **Parlanti et al. 2025** model, we can solve for the throughput functions $k$, $a$, $b$ given $S(\lambda)$ (the extended signal) and $f(\lambda)$ (the intrinsic PRISM signal):

$S(\lambda) \approx k(\lambda) f(\lambda) + a(\lambda) f(\lambda/2) + b(\lambda) f(\lambda/3)$

The next steps will utilize the data plotted above to solve these coefficients for Medium resolution gratings.
