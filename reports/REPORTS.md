# NIRSpec Wavelength Extension Analysis Reports

This directory contains diagnostic reports for the NIRSpec wavelength extension calibration project.

## Index of Reports
- [329 FS v2](329_fs_v2/REPORT.md) — (2026-03-29) FS v2 coefficient derivation: same unbiased median-R algorithm as IFU v2. Comparison plots for IFU v2, FS v2, and Parlanti k(λ). Excellent correction for P330E and J1743045.
- [329 IFU v2](329_ifu_v2/REPORT.md) — (2026-03-29) IFU v2 coefficient derivation: correct k from median R (no regularisation), + Parlanti α/β. Recalibrated spectra match truth for P330E and J1743045.
- [329 Parlanti-comparison v1](329_parlanti-comparison_v1/REPORT.md) — (2026-03-29) Comparison of Parlanti original, FS v1, and IFU v1 coefficients.
- [329 FS v1](329_fs_v1/REPORT.md) — (2026-03-29) Fixed Slit v1 coefficient derivation and spectral validation.
- [329 IFU v1](329_ifu_v1/REPORT.md) — (2026-03-29) IFU v1 coefficient derivation and spectral validation.
- [Parlanti Original](Parlanti/REPORT.md) — Original coefficients for k, alpha, and beta from the original Parlanti et al. (2025) paper.

## Instructions for Creating Reports

1.  **Create a subdirectory**: Name it as `MMDD_Mode_Version` (e.g., `329_fs_v1`).
2.  **Generate a `REPORT.md`**: This file should summarize the analysis, data used, and results.
3.  **Include Plotting Scripts**: Any scripts used to generate the report should be included in the same directory.
4.  **Relative Paths**: **CRITICAL**: Use relative paths for all file links and images within `REPORT.md`.
    *   *Correct*: `![Plot](plot.png)` or `[Script](script.py)`
    *   *Incorrect*: `![Plot](/absolute/path/to/plot.png)`
5.  **Update this Index**: Add the new report to the "Index of Reports" section above, with the newest reports first.

---
*Maintained by Antigravity*
