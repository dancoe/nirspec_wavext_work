# Detailed Analysis Plan v6 — All-Source Joint Solve

**Version:** 6.0  
**Status:** Complete  
**Objective:** Expand NIRSpec wavelength extension calibration to include science targets as calibration constraints, alongside the CALSPEC standard stars used in previous versions.

---

## 1. Data Inventory (no new downloads needed)

### 1.1 IFU Mode Data

| PID | Target | Dir | G140M stage3_ext x1d | G235M stage3_ext x1d | G395M stage3_nom x1d |
|:----|:-------|:----|:---------------------|:---------------------|:---------------------|
| 1538 | P330E | `data/IFU/PID1538_P330E/stage3_ext/` | `f100lp_g140m-f100lp_x1d.fits` | `f170lp_g235m-f170lp_x1d.fits` | — |
| 1537 | G191-B2B | `data/IFU/PID1537_G191-B2B/stage3_ext/` | same tag | same tag | — |
| 1536 | J1743045 | `data/IFU/PID1536_J1743045/stage3_ext/` | same tag | same tag | — |
| 2186 | UGC-5101 | `data/IFU/PID2186_UGC-5101/stage3_ext/` (`data/PID2186_UGC5101/stage3_nom/`) | — | `g235m_f170lp_g235m-f170lp_x1d.fits` | `g395m_f290lp_g395m-f290lp_x1d.fits` |

### 1.2 Fixed Slit Data

| PID | Target | G140M NRS2 ext | G235M NRS2 ext | G140M truth | G235M truth | G395M truth |
|:----|:-------|:--------------|:--------------|:-----------|:-----------|:-----------|
| 1537 | G191-B2B | `PID1537_G191-B2B/` level-3 | same | `g191b2b_mod_012.fits` | same | — |
| 1538 | P330E | `PID1538_P330E/` level-3 | same | `p330e_mod_008.fits` | same | — |
| 1536 | J1743045 | `PID1536_J1743045/` level-3 | same | `1743045_mod_007.fits` | same | — |
| 6644 | NGC2506-G31 | `PID6644_NGC2506G31/` level-3 | same | `ngc2506g31_mod_003.fits` | same | — |
| 1492 | Unknown | `PID1492/jw01492003001_03102_00005_nrs2_g140m_extract_1d.fits` | `PID1492/jw01492003001_03104_00004_nrs2_g235m_extract_1d.fits` | MAST G140M | MAST G235M | MAST G395M |

**PID 1492 MAST nominal files:**
- G140M: `jw01492-o001_t001-s000000001_nirspec_f100lp-g140m-s200a1-subs200a1_x1d.fits`
- G235M: `jw01492-o001_t001-s000000001_nirspec_f170lp-g235m-s200a1-subs200a1_x1d.fits`
- G395M: `jw01492-o001_t001-s000000001_nirspec_f290lp-g395m-s200a1-subs200a1_x1d.fits`

CALSPEC models in `data/CALSPEC/`.

---

## 2. Solver Implementation

### 2.1 NNLS System

At each wavelength grid point $\lambda$, solve by NNLS:

$$\begin{pmatrix} f_1(\lambda) & f_1(\lambda/2) & f_1(\lambda/3) \\ \vdots & \vdots & \vdots \\ f_N(\lambda) & f_N(\lambda/2) & f_N(\lambda/3) \end{pmatrix} \begin{pmatrix} k(\lambda) \\ \alpha(\lambda) \\ \beta(\lambda) \end{pmatrix} = \begin{pmatrix} S_1^\mathrm{obs}(\lambda) \\ \vdots \\ S_N^\mathrm{obs}(\lambda) \end{pmatrix}$$

where:
- $f_i(\lambda)$ = truth spectrum for source $i$ (CALSPEC or cross-grating stitched)
- $S_i^\mathrm{obs}(\lambda)$ = observed NRS2 extended spectrum (same units as truth)

### 2.2 Source Acceptance Criteria

For each source at wavelength $\lambda$, all four values must be finite and positive:
- $S_i^\mathrm{obs}(\lambda)$, $f_i(\lambda)$, $f_i(\lambda/2)$, $f_i(\lambda/3)$

If fewer than 2 sources pass at a given $\lambda$, that pixel is skipped and filled by linear interpolation from neighbouring valid pixels.

### 2.3 Post-Solve Processing

1. Clip: $k \geq 0.01$, $\alpha \geq 0$, $\beta \geq 0$
2. Fill isolated NaN gaps by linear interpolation
3. Apply 40-channel boxcar smooth (same as Parlanti 2025)
4. Clip again after smoothing

---

## 3. IFU v6 — Detailed Steps

### Step 1: Run the solver

From workspace root (`/Users/dcoe/NIRSpec/wavext`):

```bash
micromamba run -n jwst_1.20.2 python nirspec_wavext_work/reports/v6/ifu_v6/scripts/solve_parlanti_ifu_v6.py
```

**Expected output (stdout):**
```
IFU v6 NNLS — G140M
  P330E (G2V): NRS2 1.870–3.554 µm  R_mid=...
  G191-B2B (WD): NRS2 1.870–3.554 µm  R_mid=...
  J1743045 (A8III): NRS2 1.870–3.554 µm  R_mid=...
  Total sources: 3
  ...
IFU v6 NNLS — G235M
  ... (3 CALSPEC + 1 UGC5101)
  Total sources: 4
SAVED: results/v6/calib_v6_ifu_g140m_f100lp.fits
SAVED: results/v6/calib_v6_ifu_g235m_f170lp.fits
```

**Verify:** Check that both FITS files exist in `results/v6/` with column structure `WAVELENGTH, K, ALPHA, BETA`.

### Step 2: Run the validation plot script

```bash
micromamba run -n jwst_1.20.2 python nirspec_wavext_work/reports/v6/ifu_v6/scripts/plot_ifu_v6_validation.py
```

**Expected plots generated in `reports/v6/ifu_v6/plots/`:**
- `ifu_v6_coefficients.png` — k, α, β vs λ (3-panel, v6 + v5 + Parlanti comparison)
- `ifu_v6_g140m_calspec_validation.png` — 3 panels (one per CALSPEC source) for G140M
- `ifu_v6_g235m_calspec_validation.png` — 3 panels for G235M
- `ifu_v6_full_spectrum_G140M_{P330E,G191B2B,J1743045}.png` — merged spectra (6 files)
- `ifu_v6_full_spectrum_G235M_{P330E,G191B2B,J1743045}.png`
- `ifu_v6_ugc5101_g235m_xval.png` — UGC-5101 G235M corrected vs G395M nominal

---

## 4. FS v6 — Detailed Steps

### Step 1: Run the solver

```bash
micromamba run -n jwst_1.20.2 python nirspec_wavext_work/reports/v6/fs_v6/scripts/solve_parlanti_fs_v6.py
```

**Expected output (stdout):**
```
FS v6 NNLS — g140m
  G191-B2B (WD-DA.8): NRS2 1.870–3.200 µm  R_mid=...
  P330E (G2V): ...
  J1743045 (A8III): ...
  NGC2506-G31 (G1V): ...
  PID 1492 [cross-grating]: added as g140m source
  Total sources: 5
FS v6 NNLS — g235m
  ...
  Total sources: 5
SAVED: results/v6/calib_v6_fs_g140m_f100lp.fits
SAVED: results/v6/calib_v6_fs_g235m_f170lp.fits
```

### Step 2: Run the validation plot script

```bash
micromamba run -n jwst_1.20.2 python nirspec_wavext_work/reports/v6/fs_v6/scripts/plot_fs_v6_validation.py
```

**Expected plots generated in `reports/v6/fs_v6/plots/`:**
- `fs_v6_coefficients.png` — k, α, β comparison
- `fs_v6_g140m_calspec_validation.png` — 4-panel CALSPEC validation (G140M)
- `fs_v6_g235m_calspec_validation.png` — 4-panel CALSPEC validation (G235M)
- `fs_v6_full_g140m_{G191B2B,P330E,J1743045,NGC2506G31}.png` — 4 full-spectrum merged plots
- `fs_v6_full_g235m_{G191B2B,P330E,J1743045,NGC2506G31}.png` — 4 full-spectrum merged plots
- `fs_v6_pid1492_g140m_xval.png` — PID 1492 G140M cross-validation
- `fs_v6_pid1492_g235m_xval.png` — PID 1492 G235M cross-validation

---

## 5. Verification Checklist

**After running solvers:**
- [ ] `results/v6/calib_v6_ifu_g140m_f100lp.fits` exists and has 400 rows
- [ ] `results/v6/calib_v6_ifu_g235m_f170lp.fits` exists and has 400 rows
- [ ] `results/v6/calib_v6_fs_g140m_f100lp.fits` exists and has 400 rows
- [ ] `results/v6/calib_v6_fs_g235m_f170lp.fits` exists and has 400 rows

**Coefficient sanity checks:**
- [ ] IFU G235M: k median in 0.5–0.8 range (improvement over v4/v5 ~0.47)
- [ ] FS G235M: k median in 0.85–0.97 range (near v5 0.957)
- [ ] All α values ≥ 0 after clipping

**After running plot scripts:**
- [ ] 10 IFU plots generated
- [ ] 13 FS plots generated
- [ ] No matplotlib errors (check stderr for warnings)

---

## 6. Reports

- [IFU v6 Report](ifu_v6/REPORT_ifu_v6.md)
- [FS v6 Report](fs_v6/REPORT_fs_v6.md)

Both reports use **only relative-path image links** (e.g., `plots/ifu_v6_coefficients.png`). Do not embed image data inline in reports.

---

## 7. Notes for v7

- **True cross-validation needed:** Hold out one science source from the NNLS training set and use it as validation only. PID 1492 is a good candidate to hold out since it is also the primary new source in v6.
- **IFU G140M degeneracy:** Without a G1V IFU source, G140M k/α will remain degenerate. Check if any archival IFU data of cool stars is available in MAST.
- **PID 1492 normalisation:** The DN/s → Jy scale factor is derived from a partial wavelength overlap. Investigate how sensitive the solution is to this normalisation.
- **DETAILED_LOG.md** contains full execution trace for reference when implementing v7.

---

*Created: 2026-03-30*
