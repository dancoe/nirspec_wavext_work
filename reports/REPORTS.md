# NIRSpec Wavelength Extension Analysis Reports

This directory contains diagnostic reports for the NIRSpec wavelength extension calibration project.

## 🚀 Master Analysis Roadmap
- **[Detailed Analysis Plan (Master)](DETAILED_PLAN.md)** — The definitive guide to data download, reduction, analysis, and reporting for the full 0.6–5.6 µm extension.

## Index of Reports
- [IFU v6 All-Source Joint Solve](v6/ifu_v6/REPORT_ifu_v6.md) — (2026-03-30) IFU v6 adds UGC-5101 (cross-grating G395M truth) as 4th G235M source. G235M k rises from IFU v5 ~0.47 to 0.61; G140M k unchanged (no IFU G1V available). Science targets used as calibration constraints for first time.
- [FS v6 All-Source Joint Solve](v6/fs_v6/REPORT_fs_v6.md) — (2026-03-30) FS v6 adds PID 1492 cross-grating as 5th source (DN/s NRS2 scaled to Jy via MAST truth overlap). G235M k=0.910; G140M k=0.726 (lower than v5 0.959 — PID 1492 normalisation may pull solution). 13 validation plots generated.
- [IFU v5 Science Validation](v5/ifu_v5/REPORT_ifu_v5.md) — (2026-03-29) Science validation of v5 extended pipeline on AGN (PID 2654: SDSSJ0749, SDSSJ0841 at z~2) and ULIRG (PID 2186: UGC-5101 at z=0.039). UGC-5101 G235M extended to 5.5 µm (vs 3.17 µm nominal); PAH 3.3 µm and L-band dust recovered. SDSSJ0841 H-α recovery in extended G140M at z~2.
- [Parlanti Comparison v5](v5/parlanti-comparison_v5/REPORT_parlanti-comparison_v5.md) — (2026-03-29) v5 FS coefficients vs Parlanti (2025) vs degenerate v4. k(λ) agreement to <5% in NRS1; degeneracy resolved by G1V star. All coefficients (k, α, β) compared across both gratings.
- [FS v5](v5/fs_v5/REPORT_fs_v5.md) — (2026-03-30) FS v5 final joint-solve coefficients from 4 stars (WD, G2V, A8III, G1V). Resolves k/alpha degeneracy: G140M k=0.959, G235M k=0.957. Contamination alpha returns to physical values (<0.01).
- [329 FS v4](v4/329_fs_v4/REPORT_329_fs_v4.md) — (2026-03-29) FS v4 NNLS simultaneous derivation of $k$, $\alpha$, $\beta$ from 3 CALSPEC stars. G140M k=0.701, G235M k=0.723. NNLS $\alpha$ is 2–6× Parlanti published values; G235M k drops from v3 0.97 to v4 0.72.
- [329 IFU v4](v4/329_ifu_v4/REPORT_329_ifu_v4.md) — (2026-03-29) IFU v4 NNLS simultaneous derivation of $k$, $\alpha$, $\beta$ from 3 CALSPEC stars. G140M k=0.509, G235M k=0.466. NNLS $\alpha$ is 2–5× Parlanti published values; NGC2506-G31 4th source pending to validate.
- [Analysis Plan v4](v4/ANALYSIS_PLAN_v4.md) — (2026-03-29) Plan for full $\kappa, \alpha, \beta$ derivation for both IFU and FS.
- [329 Parlanti-comparison v3](v3/329_parlanti-comparison_v3/REPORT_329_parlanti-comparison_v3.md) — (2026-03-29) v3 coefficient comparison: IFU v3 (Level-3 stage3_ext) vs FS v3 (Level-3 nrs2_spec3_ext) vs Parlanti (2025). IFU v3 tracks Parlanti ~1.1–1.3×; FS v3 is ~1.3–2× higher near NRS1/NRS2 boundary (mode-dependent photom).
- [329 FS v3](v3/329_fs_v3/REPORT_329_fs_v3.md) — (2026-03-29) FS v3 coefficient derivation using Level-3 NRS2 Spec3Pipeline products. MAST L3 NRS1 vs our L3 NRS2 ext comparison. Excellent correction for P330E and J1743045 (within ~5%).
- [329 IFU v3](v3/329_ifu_v3/REPORT_329_ifu_v3.md) — (2026-03-29) IFU v3 coefficient derivation using Level-3 stage3_ext x1d (Spec3+cubepar, NRS1+NRS2). Adds MAST L3 nominal vs our L3 ext comparison panels. Same algorithm as v2.
- [329 Parlanti-comparison v2](v2/329_parlanti-comparison_v2/REPORT_329_parlanti-comparison_v2.md) — (2026-03-29) Coefficient comparison for k, alpha, and beta across v2 FS/IFU and Parlanti (Thin Black, Red, Cyan).
- [329 FS v2](v2/329_fs_v2/REPORT_329_fs_v2.md) — (2026-03-29) FS v2 coefficient derivation: same unbiased median-R algorithm as IFU v2. Comparison plots for IFU v2, FS v2, and Parlanti k(λ). Excellent correction for P330E and J1743045.
- [329 IFU v2](v2/329_ifu_v2/REPORT_329_ifu_v2.md) — (2026-03-29) IFU v2 coefficient derivation: correct k from median R (no regularisation), + Parlanti α/β. Recalibrated spectra match truth for P330E and J1743045.
- [329 Parlanti-comparison v1](v1/329_parlanti-comparison_v1/REPORT_329_parlanti-comparison_v1.md) — (2026-03-29) Comparison of Parlanti original, FS v1, and IFU v1 coefficients.
- [329 FS v1](v1/329_fs_v1/REPORT_329_fs_v1.md) — (2026-03-29) Fixed Slit v1 coefficient derivation and spectral validation.
- [329 IFU v1](v1/329_ifu_v1/REPORT_329_ifu_v1.md) — (2026-03-29) IFU v1 coefficient derivation and spectral validation.
- [Parlanti Original](Parlanti/REPORT_Parlanti.md) — (2025) Original coefficients for k, alpha, and beta from the original Parlanti et al. (2025) paper.

## Instructions for Creating Reports

1. Create a new subdirectory: `MMDD_Mode_Version` (e.g., `329_fs_v1`).
2. Generate all diagnostic plots within that subdirectory.
3. Create a `REPORT_MMDD_Mode_Version.md` (e.g., `REPORT_329_fs_v1.md`) in the same subdirectory.
4. Use **relative paths** for all embedded images and script links.
5. Add the new report to the top of the **Index of Reports** in this file.

---
*Maintained by Antigravity*
