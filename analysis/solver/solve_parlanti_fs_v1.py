"""
Parlanti et al. (2025) k/α/β solver — Fixed Slit (FS) version 1 (fs_v1)

Parallels the IFU v1 solver but uses FS NRS2 per-exposure x1d extractions
instead of IFU stage3_ext cubes. This allows direct comparison of FS vs IFU
calibration coefficients.

Data sources:
  - S_obs(λ):  Per-exposure NRS2 extract_1d.fits from FS extended pipeline
               (1 dither per source; the dither closest to slit center)
  - f(λ):      CALSPEC model SEDs (Jy)
  - S_truth(λ): MAST Level 3 x1d for the successive grating (same slit mode)
               * G140M NRS2 truth: G235M MAST L3 x1d (PID 1537/1538/1536)
               * G235M NRS2 truth: G395M MAST L3 x1d (PID 1537/1538/1536)

Three primary standard stars (same as IFU v1):
  - P330-E   (G2V,   PID 1538)
  - G191-B2B (WD,    PID 1537)
  - J1743045 (A8III, PID 1536)

Model (Parlanti Eq. 1):
    S_obs(λ) = k(λ)·f(λ) + α̃(λ)·f(λ/2) + β̃(λ)·f(λ/3)

Solver: NNLS bootstrap, point-wise, with 40-channel box smoothing.

Output:
  - CSV coefficient files for G140M and G235M
  - Diagnostic plots: raw and corrected ratios vs. truth
  - Comparison overlays with IFU v1 results
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import interp1d
from scipy.optimize import nnls
from scipy.ndimage import median_filter

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE        = '/Users/dcoe/NIRSpec/wavext'
FS_DIR      = f'{BASE}/data'
CALSPEC_DIR = f'{BASE}/data/CALSPEC'
IFU_CAL_DIR = f'{BASE}/nirspec_wavext_work/plots/Parlanti/cal/ifu_v1'
OUTPUT_DIR  = f'{BASE}/nirspec_wavext_work/plots/Parlanti/cal/fs_v1'
os.makedirs(OUTPUT_DIR, exist_ok=True)

C_ANG_S = 2.99792458e18  # Å/s

# ── Source table ──────────────────────────────────────────────────────────────
# For each source, we need:
#   nrs2_g140m: NRS2 FS extract_1d covering G140M extended range (~1.96–3.16 µm)
#   nrs2_g235m: NRS2 FS extract_1d covering G235M extended range (~3.28–5.31 µm)
#   l3_g140m:   MAST L3 x1d for G140M nominal (NRS1 truth for the G140M NRS2 analysis)
#   l3_g235m:   MAST L3 x1d for G235M nominal (NRS1 truth, and also the ground truth
#                for G140M NRS2 comparison at 1.87–3.15 µm overlap)
#   l3_g395m:   MAST L3 x1d for G395M nominal (truth for G235M NRS2 analysis)

SOURCES = {
    '1538': {
        'name':     'P330E',
        'calspec':  'p330e_mod_008.fits',
        # NRS2 FS extractions (Jy, produced by run_fs_nrs2_spec2.py with Parlanti flat overrides)
        'nrs2_g140m': f'{FS_DIR}/PID1538_P330E/nrs2_spec2_cal/jw01538160001_06101_00001_nrs2_x1d.fits',
        'nrs2_g235m': f'{FS_DIR}/PID1538_P330E/nrs2_spec2_cal/jw01538160001_08101_00005_nrs2_x1d.fits',
        # MAST L3 nominal x1ds (Obs 160)
        'l3_g140m': f'{FS_DIR}/PID1538_P330E/jw01538-o160_t002-s000000001_nirspec_f100lp-g140m-s1600a1-sub2048_x1d.fits',
        'l3_g235m': f'{FS_DIR}/PID1538_P330E/jw01538-o160_t002-s000000001_nirspec_f170lp-g235m-s1600a1-sub2048_x1d.fits',
        'l3_g395m': f'{FS_DIR}/PID1538_P330E/jw01538-o160_t002-s000000001_nirspec_f290lp-g395m-s1600a1-sub2048_x1d.fits',
    },
    '1537': {
        'name':     'G191-B2B',
        'calspec':  'g191b2b_mod_012.fits',
        'nrs2_g140m': f'{FS_DIR}/PID1537_G191-B2B/nrs2_spec2_cal/jw01537007001_07101_00003_nrs2_x1d.fits',
        'nrs2_g235m': f'{FS_DIR}/PID1537_G191-B2B/nrs2_spec2_cal/jw01537007001_09101_00003_nrs2_x1d.fits',
        'l3_g140m': f'{FS_DIR}/PID1537_G191-B2B/jw01537-o007_t001-s000000001_nirspec_f100lp-g140m-s1600a1-sub2048_x1d.fits',
        'l3_g235m': f'{FS_DIR}/PID1537_G191-B2B/jw01537-o007_t001-s000000001_nirspec_f170lp-g235m-s1600a1-sub2048_x1d.fits',
        'l3_g395m': f'{FS_DIR}/PID1537_G191-B2B/jw01537-o007_t001-s000000001_nirspec_f290lp-g395m-s1600a1-sub2048_x1d.fits',
    },
    '1536': {
        'name':     'J1743045',
        'calspec':  '1743045_mod_007.fits',
        # NRS2 FS extractions (Jy, run_fs_nrs2_spec2.py)
        'nrs2_g140m': f'{FS_DIR}/PID1536_J1743045/nrs2_spec2_cal/jw01536002001_03104_00003_nrs2_x1d.fits',
        'nrs2_g235m': f'{FS_DIR}/PID1536_J1743045/nrs2_spec2_cal/jw01536002001_03106_00001_nrs2_x1d.fits',
        'l3_g140m': f'{FS_DIR}/PID1536_J1743045/jw01536-o002_t004-s000000001_nirspec_f100lp-g140m-s1600a1-sub2048_x1d.fits',
        'l3_g235m': f'{FS_DIR}/PID1536_J1743045/jw01536-o002_t004-s000000001_nirspec_f170lp-g235m-s1600a1-sub2048_x1d.fits',
        'l3_g395m': f'{FS_DIR}/PID1536_J1743045/jw01536-o002_t004-s000000001_nirspec_f290lp-g395m-s1600a1-sub2048_x1d.fits',
    },
}

# Wavelength grids for solving (NRS2 region only, in µm)
GRID_G140M = np.linspace(1.95, 3.50, 200)   # G140M FS NRS2 range
GRID_G235M = np.linspace(3.30, 5.25, 200)   # G235M FS NRS2 range

# ── Helpers ────────────────────────────────────────────────────────────────────
def load_spec(path, label=''):
    """Load wavelength and flux from a _x1d or _extract_1d FITS file."""
    if not os.path.exists(path):
        print(f'  WARNING: missing {os.path.basename(path)}' +
              (f' [{label}]' if label else ''))
        return np.array([]), np.array([])
    with fits.open(path) as h:
        # Try EXTRACT1D first (per-exposure), then EXTRACT1D (L3)
        ext = None
        for name in ['EXTRACT1D', 'EXTRACT1D']:
            if name in [hh.name for hh in h]:
                ext = h[name].data
                break
        if ext is None:
            ext = h[1].data
        wl = ext['WAVELENGTH'].astype(float)
        fl = ext['FLUX'].astype(float)
        dq = ext['DQ'].astype(int) if 'DQ' in ext.dtype.names else np.zeros(len(wl), int)
    good = (wl > 0) & np.isfinite(fl) & (fl > 0) & (dq == 0)
    return wl[good], fl[good]


def load_calspec_jy(fname):
    """Load CALSPEC model, convert F_lambda (erg/s/cm²/Å) → F_nu (Jy)."""
    with fits.open(f'{CALSPEC_DIR}/{fname}') as h:
        d    = h[1].data
        wl_a = d['WAVELENGTH'].astype(float)
        flam = d['FLUX'].astype(float)
    fnu_jy = flam * (wl_a**2) / C_ANG_S * 1e23
    order  = np.argsort(wl_a)
    return wl_a[order] / 1e4, fnu_jy[order]   # µm, Jy


def box_smooth(y, box=40):
    """Box (top-hat) smoothing, preserving endpoints."""
    s = np.convolve(y, np.ones(box) / box, mode='same')
    s[:box//2]  = y[:box//2]
    s[-box//2:] = y[-box//2:]
    return s


def check_flux_units(path, label=''):
    """Quick sanity check: return typical flux in Jy (should be 1e-4 to 1e1 for standards)."""
    wl, fl = load_spec(path, label)
    if len(fl) == 0:
        return np.nan
    return np.nanmedian(fl)


# ── Core solver ────────────────────────────────────────────────────────────────
def solve_fs(grating, wav_grid):
    """
    Solve for k(λ), α̃(λ), β̃(λ) using FS NRS2 extractions vs next-grating truth.

    For G140M NRS2:
        S_obs   = per-exposure NRS2 extract_1d (extended pipeline, 1.96–3.16 µm)
        S_truth = MAST L3 G235M x1d (nominal range, 1.70–3.15 µm)
        Overlap: ~1.96–3.15 µm

    For G235M NRS2:
        S_obs   = per-exposure NRS2 extract_1d (extended pipeline, 3.28–5.31 µm)
        S_truth = MAST L3 G395M x1d (nominal range, 2.87–5.27 µm)
        Overlap: ~3.28–5.27 µm

    At each wavelength λ in wav_grid, NNLS solves across available sources:
        S_obs(λ) / f_calspec(λ) = k + α·[f(λ/2)/f(λ)] + β·[f(λ/3)/f(λ)]
    """
    print(f'\n--- Solving {grating} (FS v1) ---')

    all_data = []

    for pid, src in SOURCES.items():
        name = src['name']

        if grating == 'G140M':
            nrs2_path  = src['nrs2_g140m']
            truth_path = src['l3_g235m']
            nrs2_lo    = 1.93
        else:  # G235M
            nrs2_path  = src['nrs2_g235m']
            truth_path = src['l3_g395m']
            nrs2_lo    = 3.25

        # Load NRS2 per-exposure extract_1d
        w_nrs2, f_nrs2 = load_spec(nrs2_path, f'{name} NRS2')
        # Load truth (MAST L3 next grating)
        w_truth, f_truth = load_spec(truth_path, f'{name} truth')

        if len(w_nrs2) < 10:
            print(f'  Skipping {name}: no NRS2 data at {nrs2_path}')
            continue
        if len(w_truth) < 10:
            print(f'  Skipping {name}: no truth data')
            continue

        # Load CALSPEC as intrinsic reference
        wl_cs, fj_cs = load_calspec_jy(src['calspec'])
        f_ref = interp1d(wl_cs, fj_cs, bounds_error=False, fill_value=np.nan)

        # Quick unit / order-of-magnitude check
        med_nrs2  = np.nanmedian(f_nrs2)
        med_truth = np.nanmedian(f_truth[
            (w_truth > nrs2_lo) & (w_truth < min(w_truth.max(), 3.2 if grating=='G140M' else 5.3))
        ] if len(w_truth) > 0 else [np.nan])
        med_cs_at_overlap = np.nanmedian(f_ref(
            np.linspace(max(w_nrs2.min(), 2.0 if grating=='G140M' else 3.5),
                        min(w_nrs2.max(), 3.0 if grating=='G140M' else 5.0), 50)
        ))
        print(f'  {name}: NRS2 {w_nrs2.min():.3f}–{w_nrs2.max():.3f} µm ({len(w_nrs2)} px)')
        print(f'    median NRS2={med_nrs2:.4g} Jy, truth={med_truth:.4g} Jy, CALSPEC={med_cs_at_overlap:.4g} Jy')

        all_data.append({
            'name':    name,
            'nrs2':    interp1d(w_nrs2,  f_nrs2,  bounds_error=False, fill_value=np.nan),
            'truth':   interp1d(w_truth, f_truth, bounds_error=False, fill_value=np.nan),
            'calspec': f_ref,
        })

    if len(all_data) == 0:
        print(f'ERROR: No valid source data for {grating} FS analysis')
        return None, None, None, []

    n_src = len(all_data)
    print(f'  Solving with {n_src} source(s)')

    # ── Point-wise NNLS at each wavelength ────────────────────────────────────
    ks, alphas, betas = [], [], []

    for L in wav_grid:
        A_rows = []
        b_vec  = []

        for d in all_data:
            f_cs = d['calspec'](L)
            if not np.isfinite(f_cs) or f_cs <= 0:
                continue
            f_cs2 = d['calspec'](L / 2)
            f_cs3 = d['calspec'](L / 3)
            if not np.isfinite(f_cs2): f_cs2 = 0.0
            if not np.isfinite(f_cs3): f_cs3 = 0.0

            S_obs = d['nrs2'](L)
            if not np.isfinite(S_obs) or S_obs <= 0:
                continue

            # Normalise by CALSPEC: R = S_obs/f_cs = k + α*(f_cs2/f_cs) + β*(f_cs3/f_cs)
            norm = f_cs
            A_rows.append([1.0, f_cs2 / norm, f_cs3 / norm])
            b_vec.append(S_obs / norm)

        if len(A_rows) == 0:
            ks.append(np.nan); alphas.append(0.0); betas.append(0.0)
            continue

        # Regularized NNLS with weak k prior at 1.0 and ghost-order regularization
        k_prior  = 1.5    # allows k to be > 1 if data demands
        g_reg    = 0.01   # weak second/third-order suppression

        A_ext = np.vstack([
            np.array(A_rows),
            [k_prior, 0.0, 0.0],   # k prior
            [0.0, g_reg, 0.0],     # α regularization
            [0.0, 0.0, g_reg],     # β regularization
        ])
        b_ext = np.hstack([
            np.array(b_vec),
            [k_prior * 1.0],       # k prior target = 1.0
            [0.0],
            [0.0],
        ])

        sol, _ = nnls(A_ext, b_ext)
        ks.append(sol[0])
        alphas.append(sol[1])
        betas.append(sol[2])

    ks     = np.array(ks)
    alphas = np.array(alphas)
    betas  = np.array(betas)

    # Replace NaN k values with interpolated ones before smoothing
    valid = np.isfinite(ks)
    if valid.sum() < 5:
        print(f'  WARNING: only {valid.sum()} valid wavelength points for {grating}')
        return ks, alphas, betas, all_data

    if not valid.all():
        ks_fill = np.interp(wav_grid, wav_grid[valid], ks[valid])
        ks = ks_fill

    # ── Smooth coefficients (Parlanti: 40-channel box) ─────────────────────────
    ks_s     = box_smooth(ks,     40)
    alphas_s = box_smooth(alphas, 40)
    betas_s  = box_smooth(betas,  40)

    return ks_s, alphas_s, betas_s, all_data


# ── Plot results ───────────────────────────────────────────────────────────────
def plot_results(grating, wav_grid, ks, alphas, betas, all_data, nrs2_lo):
    if ks is None:
        return

    fig, axes = plt.subplots(3, 1, figsize=(13, 12), sharex=True)
    fig.suptitle(
        f'FS Calibration v1 — {grating}/F{"1" if "140" in grating else "1"}00LP NRS2\n'
        f'k/α/β coefficients ({len(all_data)} sources)',
        fontsize=13
    )

    colors = {'P330E': '#e67e22', 'G191-B2B': '#2980b9', 'J1743045': '#27ae60'}

    # ── Panel 1: k(λ) and comparison with IFU v1 ──────────────────────────────
    ax = axes[0]
    ax.axhline(1.0, color='0.6', lw=0.8, ls='--', label='k=1 reference')
    ax.plot(wav_grid, ks, color='k', lw=2.0, label=f'k(λ) FS v1 ({grating})')

    # Overlay IFU v1 k(λ) if available
    ifu_csv = os.path.join(IFU_CAL_DIR, f'coeffs_ifu_v1_{grating}.csv')
    if os.path.exists(ifu_csv):
        ifu_data = np.loadtxt(ifu_csv, delimiter=',', skiprows=1)
        ax.plot(ifu_data[:, 0], ifu_data[:, 1], color='navy', lw=1.5, ls='--',
                alpha=0.7, label=f'k(λ) IFU v1 ({grating})')

    ax.set_ylabel('k(λ)', fontsize=11)
    ax.set_ylim(0.0, 2.0)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.2)
    ax.set_title('First-order throughput correction  [FS vs IFU comparison]', fontsize=10)

    # ── Panel 2: α̃(λ) and β̃(λ) ────────────────────────────────────────────────
    ax = axes[1]
    ax.axhline(0.0, color='0.6', lw=0.8, ls='--')
    ax.plot(wav_grid, alphas, color='darkorange', lw=1.5, label='α̃(λ) 2nd-order')
    ax.plot(wav_grid, betas,  color='purple',     lw=1.5, label='β̃(λ) 3rd-order')

    # Overlay IFU v1 alpha
    if os.path.exists(ifu_csv):
        ax.plot(ifu_data[:, 0], ifu_data[:, 2], color='darkorange', lw=1.0, ls='--',
                alpha=0.6, label='α̃(λ) IFU v1')
        ax.plot(ifu_data[:, 0], ifu_data[:, 3], color='purple',     lw=1.0, ls='--',
                alpha=0.6, label='β̃(λ) IFU v1')

    ax.set_ylabel('α̃ / β̃', fontsize=11)
    ax.set_ylim(-0.05, 0.5)
    ax.legend(fontsize=9, ncol=2)
    ax.grid(True, alpha=0.2)
    ax.set_title('Higher-order contamination coefficients', fontsize=10)

    # ── Panel 3: Ratio to truth (raw vs. corrected) ────────────────────────────
    ax = axes[2]
    k_interp = interp1d(wav_grid, ks,     bounds_error=False, fill_value=np.nan)
    a_interp = interp1d(wav_grid, alphas, bounds_error=False, fill_value=np.nan)
    b_interp = interp1d(wav_grid, betas,  bounds_error=False, fill_value=np.nan)

    lo_wl = nrs2_lo + 0.05
    hi_wl = wav_grid.max() - 0.05
    w_eval = np.linspace(lo_wl, hi_wl, 400)

    for d in all_data:
        name  = d['name']
        color = colors.get(name, 'gray')

        f_obs   = d['nrs2'](w_eval)
        f_truth = d['truth'](w_eval)
        f_cs    = d['calspec'](w_eval)
        f_cs2   = d['calspec'](w_eval / 2)
        f_cs3   = d['calspec'](w_eval / 3)

        ki = k_interp(w_eval)
        ai = a_interp(w_eval)
        bi = b_interp(w_eval)

        ok = (np.isfinite(f_obs) & np.isfinite(f_truth) & (f_obs > 0) & (f_truth > 0)
              & np.isfinite(ki) & (ki > 0))
        if ok.sum() < 5:
            continue

        ratio_raw  = f_obs / f_truth
        f_corr     = (f_obs - ai * f_cs2 - bi * f_cs3) / ki
        ratio_corr = f_corr / f_truth

        sm = lambda y: median_filter(y, size=21)
        ax.plot(w_eval[ok], sm(ratio_raw)[ok],  color=color, lw=1.0, alpha=0.35,
                ls='--', label=f'{name} raw')
        ax.plot(w_eval[ok], sm(ratio_corr)[ok], color=color, lw=1.8,
                label=f'{name} corrected')

    ax.axhline(1.0, color='k', lw=1.0)
    ax.axhspan(0.90, 1.10, color='limegreen', alpha=0.08)
    ax.axhspan(0.80, 1.20, color='gold',      alpha=0.04)
    ax.set_ylabel('S_corrected / S_truth', fontsize=11)
    ax.set_xlabel('Wavelength (µm)', fontsize=11)
    ax.set_ylim(0.2, 2.5)
    ax.legend(fontsize=8, ncol=3)
    ax.grid(True, alpha=0.2)
    ax.set_title('Ratio to truth: raw (dashed) vs Parlanti-corrected (solid)', fontsize=10)

    plt.tight_layout()
    outpath = f'{OUTPUT_DIR}/fs_v1_{grating.lower()}_coeffs.png'
    plt.savefig(outpath, dpi=200, bbox_inches='tight')
    plt.close()
    print(f'  Saved: {outpath}')


def save_csv(grating, wav_grid, ks, alphas, betas):
    """Save coefficient arrays to CSV."""
    if ks is None:
        return
    outpath = f'{OUTPUT_DIR}/coeffs_fs_v1_{grating}.csv'
    header = 'wavelength_um,k,alpha,beta'
    np.savetxt(outpath, np.column_stack([wav_grid, ks, alphas, betas]),
               delimiter=',', header=header, comments='')
    print(f'  CSV saved: {outpath}')


def summary_comparison():
    """
    Produce a side-by-side summary plot comparing FS v1 vs IFU v1 k(λ) for both gratings.
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)
    fig.suptitle('FS v1 vs IFU v1 — k(λ) coefficient comparison\n'
                 '3 sources: P330E (G2V), G191-B2B (WD), J1743045 (A8III)', fontsize=13)

    for ax, grating, label in zip(axes, ['G140M', 'G235M'],
                                  ['G140M/F100LP (1.95–3.50 µm)', 'G235M/F170LP (3.30–5.25 µm)']):
        fs_csv  = f'{OUTPUT_DIR}/coeffs_fs_v1_{grating}.csv'
        ifu_csv = f'{IFU_CAL_DIR}/coeffs_ifu_v1_{grating}.csv'

        if os.path.exists(fs_csv):
            d = np.loadtxt(fs_csv,  delimiter=',', skiprows=1)
            ax.plot(d[:, 0], d[:, 1], color='k',    lw=2.0, label='FS v1 k(λ)')
            ax.fill_between(d[:, 0], d[:, 1] - 0.05, d[:, 1] + 0.05,
                            color='k', alpha=0.08)
        if os.path.exists(ifu_csv):
            d = np.loadtxt(ifu_csv, delimiter=',', skiprows=1)
            ax.plot(d[:, 0], d[:, 1], color='navy', lw=2.0, ls='--',
                    label='IFU v1 k(λ)')

        ax.axhline(1.0, color='0.5', ls=':', lw=0.8)
        ax.set_xlabel('Wavelength (µm)', fontsize=11)
        ax.set_ylabel('k(λ)' if ax == axes[0] else '', fontsize=11)
        ax.set_title(label, fontsize=10)
        ax.set_ylim(0, 2.0)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.2)

    plt.tight_layout()
    outpath = f'{OUTPUT_DIR}/fs_vs_ifu_k_comparison.png'
    plt.savefig(outpath, dpi=200, bbox_inches='tight')
    plt.close()
    print(f'Summary comparison saved: {outpath}')


# ── Diagnostics ────────────────────────────────────────────────────────────────
def print_file_status():
    """Print a table of which data files exist."""
    print('\n=== Data file status ===')
    for pid, src in SOURCES.items():
        name = src['name']
        for key in ['nrs2_g140m', 'nrs2_g235m', 'l3_g140m', 'l3_g235m', 'l3_g395m']:
            path = src[key]
            exists = '✅' if os.path.exists(path) else '❌'
            size   = f'{os.path.getsize(path) / 1e6:.1f} MB' if os.path.exists(path) else '--'
            print(f'  {name} {key:15s}: {exists}  {size:8s}  {os.path.basename(path)}')


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    print('=== FS Calibration v1 Solver (Parlanti model) ===')
    print(f'Output: {OUTPUT_DIR}')
    print_file_status()

    results = {}

    for grating, wav_grid, nrs2_lo in [
        ('G140M', GRID_G140M, 1.93),
        ('G235M', GRID_G235M, 3.25),
    ]:
        ks, alphas, betas, all_data = solve_fs(grating, wav_grid)

        if ks is not None:
            n_src = len(all_data)
            k_med = np.nanmedian(ks)
            k_min = np.nanmin(ks)
            k_max = np.nanmax(ks)
            a_max = np.nanmax(alphas)

            print(f'\n  {grating} results ({n_src} sources):')
            print(f'    k(λ): median={k_med:.3f}, range=[{k_min:.3f}, {k_max:.3f}]')
            print(f'    α(λ): max={a_max:.4f}')
            print(f'    β(λ): max={np.nanmax(betas):.4f}')

            save_csv(grating, wav_grid, ks, alphas, betas)
            plot_results(grating, wav_grid, ks, alphas, betas, all_data, nrs2_lo)
            results[grating] = (ks, alphas, betas, all_data, wav_grid, nrs2_lo)
        else:
            print(f'  {grating}: SKIPPED (insufficient data)')

    if results:
        summary_comparison()

    print('\n=== Done ===')


if __name__ == '__main__':
    main()
