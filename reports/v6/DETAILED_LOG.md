# Detailed Execution Log — IFU v6 and FS v6

**Date:** March 30, 2026  
**Author:** Antigravity  
**Purpose:** Full step-by-step log of how IFU v6 and FS v6 were implemented. Intended to serve as a complete guide for a future session (e.g. v7) so that the AI can replicate or extend this work from scratch.

> **Important operational note:** Do NOT embed plots inline in the AI conversation response. Doing so overwhelms the model context window and causes early termination. Instead, run plot scripts, confirm they produce files, and refer the user to the report markdown files where plots are linked.

---

## 1. Context and Motivation

### Problem Statement (from v5 review)
The IFU v5 report, Section 8 ("Full Spectrum Merged Validation"), showed that corrected spectra for calibration standards did NOT match the observed spectra in the extended NRS2 wavelength region. The v5 IFU solver used 3 CALSPEC sources (P330E G2V, G191-B2B WD, J1743045 A8III), all hot/warm stars. Without a cool G1V star in IFU mode (unlike FS v5 which added NGC2506-G31 PID 6644), the NNLS system remained degenerate.

### v6 Solution
Use ALL available data — not just CALSPEC standards — as calibration constraints. Any target with a known or derivable reference spectrum at the extended wavelengths can serve as an NNLS source. For science targets, the reference spectrum is derived from a cross-grating nominal pipeline product:

- **IFU v6:** Add UGC-5101 (PID 2186) as a 4th source for G235M only. Truth = G395M nominal × scale factor derived from overlap.
- **FS v6:** Add PID 1492 as a 5th source for both G140M and G235M. Truth = MAST adjacent-grating nominal (G235M for G140M solve, G395M for G235M solve).

---

## 2. Directory Structure

All work is rooted at `/Users/dcoe/NIRSpec/wavext/` (the workspace BASE).

```
BASE/
├── data/
│   ├── CALSPEC/                           # CALSPEC stellar models (.fits)
│   │   ├── p330e_mod_008.fits
│   │   ├── g191b2b_mod_012.fits
│   │   ├── 1743045_mod_007.fits
│   │   └── ngc2506g31_mod_003.fits
│   ├── IFU/                               # IFU pipeline outputs
│   │   ├── PID1538_P330E/stage3_ext/      # Extended IFU cubes + x1d
│   │   ├── PID1537_G191-B2B/stage3_ext/
│   │   ├── PID1536_J1743045/stage3_ext/
│   │   ├── PID2186_UGC-5101/stage3_ext/   # UGC-5101 G235M extended x1d
│   │   └── PID2186_UGC5101/stage3_nom/    # UGC-5101 G395M nominal x1d (note: different dir for nom vs ext!)
│   ├── PID1492/                           # Fixed Slit PID 1492 raw + MAST files
│   │   ├── jw01492003001_03102_00005_nrs2_g140m_extract_1d.fits  # G140M NRS2 ext (DN/s)
│   │   ├── jw01492003001_03104_00004_nrs2_g235m_extract_1d.fits  # G235M NRS2 ext (DN/s)
│   │   ├── jw01492-o001_t001-s000000001_nirspec_f100lp-g140m-s200a1-subs200a1_x1d.fits  # MAST G140M
│   │   ├── jw01492-o001_t001-s000000001_nirspec_f170lp-g235m-s200a1-subs200a1_x1d.fits  # MAST G235M
│   │   └── jw01492-o001_t001-s000000001_nirspec_f290lp-g395m-s200a1-subs200a1_x1d.fits  # MAST G395M
│   ├── PID1536_J1743045/                  # FS standard star data
│   ├── PID1537_G191-B2B/
│   ├── PID1538_P330E/
│   ├── PID6644_NGC2506G31/                # FS G1V standard (also PID6644_NGC2506-G31)
│   └── parlanti_repo/
│       └── calibration_files/             # Parlanti (2025) published coefficients
│           ├── calibration_functions_g140m_f100lp.fits
│           └── calibration_functions_g235m_f170lp.fits
├── results/
│   └── v6/                                # Solver output FITS files
│       ├── calib_v6_ifu_g140m_f100lp.fits
│       ├── calib_v6_ifu_g235m_f170lp.fits
│       ├── calib_v6_fs_g140m_f100lp.fits
│       └── calib_v6_fs_g235m_f170lp.fits
└── nirspec_wavext_work/
    └── reports/v6/
        ├── ANALYSIS_PLAN.md
        ├── DETAILED_PLAN.md
        ├── DETAILED_LOG.md                # This file
        ├── ifu_v6/
        │   ├── REPORT_ifu_v6.md
        │   ├── scripts/
        │   │   ├── solve_parlanti_ifu_v6.py
        │   │   └── plot_ifu_v6_validation.py
        │   └── plots/                     # 10 PNG files
        └── fs_v6/
            ├── REPORT_fs_v6.md
            ├── scripts/
            │   ├── solve_parlanti_fs_v6.py
            │   └── plot_fs_v6_validation.py
            └── plots/                     # 13 PNG files
```

---

## 3. Data Prerequisites

Before running the v6 solvers, verify that these files exist:

```bash
# Check IFU stage3_ext x1d files
ls data/IFU/PID1538_P330E/stage3_ext/*g140m*x1d.fits
ls data/IFU/PID1537_G191-B2B/stage3_ext/*g140m*x1d.fits
ls data/IFU/PID1536_J1743045/stage3_ext/*g140m*x1d.fits
ls data/IFU/PID2186_UGC-5101/stage3_ext/*g235m*x1d.fits
ls data/IFU/PID2186_UGC5101/stage3_nom/*g395m*x1d.fits

# Check FS NRS2 ext files (level-3 x1d from extended pipeline)
# These are found via glob recursively under PID1537_G191-B2B/ etc.
# Filename pattern: *g140m*x1d.fits with 'nrs2' in the name

# Check PID 1492 files
ls data/PID1492/jw01492003001_03102_00005_nrs2_g140m_extract_1d.fits
ls data/PID1492/jw01492003001_03104_00004_nrs2_g235m_extract_1d.fits
ls data/PID1492/jw01492-o001_t001-s000000001_nirspec_f170lp-g235m-*_x1d.fits
ls data/PID1492/jw01492-o001_t001-s000000001_nirspec_f290lp-g395m-*_x1d.fits
```

The IFU stage3_ext x1d filename convention is:  
`{filter}_{grating}-{filter}_x1d.fits` e.g. `f100lp_g140m-f100lp_x1d.fits`

---

## 4. IFU v6 Solver Implementation

**Script:** `reports/v6/ifu_v6/scripts/solve_parlanti_ifu_v6.py`  
**Run command (from BASE):**
```bash
micromamba run -n jwst_1.20.2 python nirspec_wavext_work/reports/v6/ifu_v6/scripts/solve_parlanti_ifu_v6.py
```

### 4.1 Key Configuration Constants

```python
CALSPEC_SOURCES = [
    {'name': 'P330E',    'sptype': 'G2V',    'calspec': 'p330e_mod_008.fits',     'ifu_dir': '.../PID1538_P330E'},
    {'name': 'G191-B2B', 'sptype': 'WD-DA.8','calspec': 'g191b2b_mod_012.fits',   'ifu_dir': '.../PID1537_G191-B2B'},
    {'name': 'J1743045', 'sptype': 'A8III',  'calspec': '1743045_mod_007.fits',   'ifu_dir': '.../PID1536_J1743045'},
]
EXT_TAG = {'G140M': 'f100lp_g140m-f100lp_x1d.fits',
            'G235M': 'f170lp_g235m-f170lp_x1d.fits'}
NRS2_LO = {'G140M': 1.87, 'G235M': 3.15}
GRID    = {'G140M': np.linspace(1.87, 3.55, 400),
            'G235M': np.linspace(3.15, 5.27, 400)}
SMOOTH_BOX = 40
```

### 4.2 CALSPEC Unit Conversion

CALSPEC models store `FLUX` in $\text{erg/s/cm}^2/\text{Å}$ (fλ). Convert to fν (Jy):

```python
fnu_jy = flam * wl_angstrom**2 / C_ANG_S * 1e23
```

where `C_ANG_S = 2.99792458e18` (speed of light in Å/s).

### 4.3 UGC-5101 Cross-Grating Truth (G235M only)

```python
path_g235m = '.../PID2186_UGC-5101/stage3_ext/g235m_f170lp_g235m-f170lp_x1d.fits'
path_g395m = '.../PID2186_UGC5101/stage3_nom/g395m_f290lp_g395m-f290lp_x1d.fits'
```

**Stitch logic:**
1. Compute scale = median(fl235 / fl395) over 2.87–3.15 µm overlap
2. Combine: G235M (λ < 3.0) + scaled G395M (λ ≥ 3.0)

**Note:** The two UGC-5101 directories use slightly different names:
- `PID2186_UGC-5101` (hyphen) for the extended data
- `PID2186_UGC5101` (no hyphen) for the nominal G395M data

This is a historical naming inconsistency; both directories must exist for the solver to work.

### 4.4 NNLS Solve Loop

At each grid wavelength λ:
1. Evaluate f_obs(λ), f_ref(λ), f_ref(λ/2), f_ref(λ/3) for each source
2. Reject any source where any value is NaN, non-finite, or ≤ 0
3. Skip pixel if < 2 valid sources
4. Call `scipy.optimize.nnls(A, b)` where A is (N_sources × 3), b is (N_sources,)
5. Store raw k, α, β

Post-processing:
1. Fill NaN gaps by linear interpolation from neighboring valid pixels
2. Clip k ≥ 0.01, α ≥ 0, β ≥ 0
3. Apply 40-channel boxcar smooth
4. Clip again

### 4.5 Output Format

```python
col_wl  = fits.Column(name='WAVELENGTH', format='E', array=wav_grid)
col_k   = fits.Column(name='K',          format='E', array=ks)
col_a   = fits.Column(name='ALPHA',      format='E', array=als)
col_b   = fits.Column(name='BETA',       format='E', array=bts)
hdu = fits.BinTableHDU.from_columns([col_wl, col_k, col_a, col_b])
hdu.header['GRATING'] = grating
hdu.header['NPTS']    = len(wav_grid)
```

### 4.6 Observed IFU v6 Results

```
calib_v6_ifu_g140m_f100lp.fits:
  wl range: 1.87–3.55 µm  n=400
  NRS2 k: median=0.5087  min=0.01  max=1.10
  NRS2 α: median=0.08351
  NRS2 β: median=0.00443

calib_v6_ifu_g235m_f170lp.fits:
  wl range: 3.15–5.27 µm  n=400
  NRS2 k: median=0.6078  min=0.22  max=1.00
  NRS2 α: median=0.01357
  NRS2 β: median=0.00186
```

The G235M k (0.61) is higher than IFU v4/v5 (~0.47), confirming UGC-5101 adds useful constraint. G140M is unchanged (still degenerate without an IFU G1V source).

---

## 5. FS v6 Solver Implementation

**Script:** `reports/v6/fs_v6/scripts/solve_parlanti_fs_v6.py`  
**Run command (from BASE):**
```bash
micromamba run -n jwst_1.20.2 python nirspec_wavext_work/reports/v6/fs_v6/scripts/solve_parlanti_fs_v6.py
```

### 5.1 FS Standard Star Discovery

FS NRS2 extended x1d files are located by `glob.glob` recursively:
```python
pattern = os.path.join(root_path, '**', f'*{grating_name}*x1d.fits')
hits = [f for f in glob.glob(pattern, recursive=True)
        if 'nrs2' in os.path.basename(f).lower()]
```

PID_DIRS mapping (handles alternative directory names):
```python
PID_DIRS = {
    '1537': ['PID1537_G191-B2B'],
    '1538': ['PID1538_P330E'],
    '1536': ['PID1536_J1743045'],
    '6644': ['PID6644_NGC2506G31', 'PID6644_NGC2506-G31'],
}
```

### 5.2 PID 1492 Cross-Grating Integration

PID 1492 NRS2 extended spectra are in **DN/s** (not Jy). The solver normalises them:

```python
# Compute scale in wavelength overlap
ratio = truth_interp(wl_grid) / nrs2_interp(wl_grid)
scale = float(np.nanmedian(ratio[good]))
fl_nrs_jy = fl_nrs * scale
```

Stitched truth for ghost terms (λ/2, λ/3):
- λ < truth start: Use NRS1 MAST (lower grating nominal MAST x1d)
- λ ≥ truth start: Use MAST next-grating nominal

For G140M:
- G140M truth = G235M MAST x1d (1.70–3.07 µm)
- G140M NRS1 for ghost stitch = G140M MAST x1d
- Overlap for scale: ~1.87–3.07 µm

For G235M:
- G235M truth = G395M MAST x1d (2.88–5.09 µm)
- G235M NRS1 for ghost stitch = G235M MAST x1d
- Overlap for scale: ~3.15–5.09 µm

### 5.3 NNLS Grid

```python
GRID = {
    'g140m': np.linspace(1.87, 3.20, 400),
    'g235m': np.linspace(3.15, 5.15, 400),
}
```

Note: FS G140M grid extends only to 3.20 µm (not 3.55 µm as in IFU), reflecting the FS detector edge.

### 5.4 Observed FS v6 Results

```
calib_v6_fs_g140m_f100lp.fits:
  wl range: 1.87–3.20 µm  n=400
  NRS2 k: median=0.7258  min=0.01  max=1.81
  NRS2 α: median=0.02219
  NRS2 β: median=0.00972

calib_v6_fs_g235m_f170lp.fits:
  wl range: 3.15–5.15 µm  n=400
  NRS2 k: median=0.9101  min=0.01  max=2.10
  NRS2 α: median=0.03114
  NRS2 β: median=0.00117
```

G235M k (0.91) is near the v5 value (0.957). G140M k (0.73) is somewhat lower than v5 (0.959), likely because the PID 1492 cross-grating normalisation in G140M is sensitive to overlap window choice. This should be investigated in v7 (sensitivity analysis).

---

## 6. Validation Plot Generation

### 6.1 IFU v6 Plots

**Script:** `reports/v6/ifu_v6/scripts/plot_ifu_v6_validation.py`  
**Run command (from BASE):**
```bash
micromamba run -n jwst_1.20.2 python nirspec_wavext_work/reports/v6/ifu_v6/scripts/plot_ifu_v6_validation.py
```

Uses `matplotlib.use('Agg')` for non-interactive backend (no display required).

**Plots generated (in `reports/v6/ifu_v6/plots/`):**

| Filename | Content |
|:---------|:--------|
| `ifu_v6_coefficients.png` | 3-panel: k, α, β vs λ for IFU v6, IFU v5, Parlanti |
| `ifu_v6_g140m_calspec_validation.png` | 3-column panel: obs / corrected / CALSPEC truth for each G140M source |
| `ifu_v6_g235m_calspec_validation.png` | Same for G235M (4 panels including UGC-5101 check) |
| `ifu_v6_full_spectrum_G140M_P330E.png` | Full NRS1+NRS2 merged spectrum (G140M, P330E) |
| `ifu_v6_full_spectrum_G140M_G191B2B.png` | Full merged (G140M, G191-B2B) |
| `ifu_v6_full_spectrum_G140M_J1743045.png` | Full merged (G140M, J1743045) |
| `ifu_v6_full_spectrum_G235M_P330E.png` | Full merged (G235M, P330E) |
| `ifu_v6_full_spectrum_G235M_G191B2B.png` | Full merged (G235M, G191-B2B) |
| `ifu_v6_full_spectrum_G235M_J1743045.png` | Full merged (G235M, J1743045) |
| `ifu_v6_ugc5101_g235m_xval.png` | UGC-5101 corrected G235M vs G395M nominal cross-validation |

### 6.2 FS v6 Plots

**Script:** `reports/v6/fs_v6/scripts/plot_fs_v6_validation.py`  
**Run command (from BASE):**
```bash
micromamba run -n jwst_1.20.2 python nirspec_wavext_work/reports/v6/fs_v6/scripts/plot_fs_v6_validation.py
```

**Plots generated (in `reports/v6/fs_v6/plots/`):**

| Filename | Content |
|:---------|:--------|
| `fs_v6_coefficients.png` | 3-panel: k, α, β vs λ for FS v6, FS v5, Parlanti |
| `fs_v6_g140m_calspec_validation.png` | 4-source G140M validation |
| `fs_v6_g235m_calspec_validation.png` | 4-source G235M validation |
| `fs_v6_full_g140m_G191B2B.png` | Full merged spectrum (G140M, G191-B2B) |
| `fs_v6_full_g140m_P330E.png` | Full merged (G140M, P330E) |
| `fs_v6_full_g140m_J1743045.png` | Full merged (G140M, J1743045) |
| `fs_v6_full_g140m_NGC2506G31.png` | Full merged (G140M, NGC2506-G31 G1V) |
| `fs_v6_full_g235m_G191B2B.png` | Full merged (G235M, G191-B2B) |
| `fs_v6_full_g235m_P330E.png` | Full merged (G235M, P330E) |
| `fs_v6_full_g235m_J1743045.png` | Full merged (G235M, J1743045) |
| `fs_v6_full_g235m_NGC2506G31.png` | Full merged (G235M, NGC2506-G31) |
| `fs_v6_pid1492_g140m_xval.png` | PID 1492 corrected G140M vs G235M MAST cross-validation |
| `fs_v6_pid1492_g235m_xval.png` | PID 1492 corrected G235M vs G395M MAST cross-validation |

---

## 7. Report Writing

Reports are plain Markdown tables and text — no inline images in the source. Plot references use relative paths to the `plots/` subdirectory. The format follows the v5 reports.

**IFU v6 report:** `reports/v6/ifu_v6/REPORT_ifu_v6.md`  
**FS v6 report:** `reports/v6/fs_v6/REPORT_fs_v6.md`

After creating reports, update the master index `reports/REPORTS.md` by prepending new entries at the top of the "Index of Reports" section.

---

## 8. Git Workflow

All work in `nirspec_wavext_work/` is tracked in the `nirspec_wavext_work` git repository.

```bash
cd /Users/dcoe/NIRSpec/wavext/nirspec_wavext_work

# Stage all new/modified files
git add reports/v6/ notes/CHANGE_LOG.md notes/PROMPT_LOG.md reports/REPORTS.md

# Commit
git commit -m "Add IFU v6 and FS v6: all-source joint solve + reports + plan files"

# Push
git push
```

Results files (FITS) in `results/v6/` are NOT tracked by the `nirspec_wavext_work` repo (they live in the parent directory). They do not need to be committed.

---

## 9. Environment

All Python runs use:
```bash
micromamba run -n jwst_1.20.2 python <script>
```

Key packages: `numpy`, `scipy`, `astropy`, `matplotlib`. No special imports beyond standard scientific Python.

The conda environment `jwst_1.20.2` at `/Users/dcoe/mambaforge/envs/jwst_1.20.2/` is the production environment.

---

## 10. Common Pitfalls

### 10.1 Directory name inconsistency (UGC-5101)
The UGC-5101 directories use two different names. Always check both:
- Extended data: `data/IFU/PID2186_UGC-5101/` (hyphen)
- Nominal G395M: `data/PID2186_UGC5101/` (no hyphen, under `data/` not `data/IFU/`)

### 10.2 Heredoc/terminal state
When running Python commands in the terminal from the AI context, the terminal may be in a confused state from previous heredoc attempts. Use `create_file` to write Python scripts to `/tmp/`, then run them with `micromamba run -n jwst_1.20.2 python /tmp/script.py`.

### 10.3 PID 1492 DN/s units
The PID 1492 NRS2 extended spectra are in DN/s (raw extraction), not Jy. The solver MUST normalise them against the MAST Jy spectra before NNLS. If the normalisation scale comes out as ~1.0 (i.e., units are already Jy), that is suspicious — verify with the header.

### 10.4 Inline image prohibition
Do NOT ask the AI to display or describe inline plot images. This fills the context window. Instead, run the plot scripts, confirm the files were created (e.g., `ls reports/v6/ifu_v6/plots/`), and read the report for descriptions.

### 10.5 Condition number check
If `np.linalg.svd` gives a condition number `cond >> 1000` at many wavelength points, the solver is degenerate. Add more spectrally diverse sources or investigate which sources are causing degeneracy. Log the median condition number from the solver output.

---

## 11. Summary of Results

| File | Grating | k median | α median | β median | Notes |
|:-----|:--------|:---------|:---------|:---------|:------|
| `calib_v6_ifu_g140m_f100lp.fits` | IFU G140M | 0.509 | 0.0835 | 0.0044 | 3 sources; degenerate (no G1V IFU) |
| `calib_v6_ifu_g235m_f170lp.fits` | IFU G235M | 0.608 | 0.0136 | 0.0019 | 4 sources; UGC-5101 added value |
| `calib_v6_fs_g140m_f100lp.fits`  | FS G140M  | 0.726 | 0.0222 | 0.0097 | 5 sources; k lower than v5 (PID 1492 effect) |
| `calib_v6_fs_g235m_f170lp.fits`  | FS G235M  | 0.910 | 0.0311 | 0.0012 | 5 sources; near v5 0.957 |

---

*Last updated: 2026-03-30*
