"""
Parlanti et al. (2025) k/α/β solver — IFU version 1 (ifu_v1)

Uses the correct data sources:
  - S_obs(λ):  IFU stage3_ext NRS2 portion (Jy, from cube-extracted x1d)
  - S_truth(λ): IFU stage3 nominal of the NEXT grating (same calibration system)
    * G140M NRS2 truth: G235M stage3 nominal  (overlap 1.87–3.17 µm)
    * G235M NRS2 truth: G395M FS nominal       (overlap 3.15–5.14 µm)

Three standard stars: P330E (1538), G191-B2B (1537), J1743045 (1536)

Model (Parlanti eq. 1):
    S_obs(λ) = k(λ)·f(λ) + α(λ)·f(λ/2) + β(λ)·f(λ/3)

where f(λ) is the intrinsic reference (CALSPEC).

Solver: NNLS bootstrap, point-wise, with 40-channel box smoothing.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import interp1d
from scipy.optimize import nnls

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE        = '/Users/dcoe/NIRSpec/wavext'
IFU_DIR     = f'{BASE}/data/IFU'
FS_DIR      = f'{BASE}/data'
CALSPEC_DIR = f'{BASE}/data/CALSPEC'
OUTPUT_DIR  = f'{BASE}/nirspec_wavext_work/plots/Parlanti/cal/ifu_v1'
os.makedirs(OUTPUT_DIR, exist_ok=True)

C_ANG_S = 2.99792458e18  # Å/s

# ── Source table ──────────────────────────────────────────────────────────────
SOURCES = {
    '1538': {
        'name':     'P330E',
        'calspec':  'p330e_mod_008.fits',
        'ifu_dir':  f'{IFU_DIR}/PID1538_P330E',
        'g395m_fs': f'{FS_DIR}/PID1538_P330E/jw01538-o160_t002-s000000001_nirspec_f290lp-g395m-s1600a1-sub2048_x1d.fits',
    },
    '1537': {
        'name':     'G191-B2B',
        'calspec':  'g191b2b_mod_012.fits',
        'ifu_dir':  f'{IFU_DIR}/PID1537_G191-B2B',
        'g395m_fs': f'{FS_DIR}/PID1537_G191-B2B/jw01537-o007_t001-s000000001_nirspec_f290lp-g395m-s1600a1-sub2048_x1d.fits',
    },
    '1536': {
        'name':     'J1743045',
        'calspec':  '1743045_mod_007.fits',
        'ifu_dir':  f'{IFU_DIR}/PID1536_J1743045',
        'g395m_fs': f'{FS_DIR}/PID1536_J1743045/jw01536-o002_t004-s000000001_nirspec_f290lp-g395m-s1600a1-sub2048_x1d.fits',
    },
}

# Wavelength grids for solving (NRS2 region only, in µm)
GRID_G140M = np.linspace(1.87, 3.55, 200)   # G140M NRS2 range
GRID_G235M = np.linspace(3.15, 5.40, 200)   # G235M NRS2 range

# ── Helpers ────────────────────────────────────────────────────────────────────
def load_spec(path, label=''):
    if not os.path.exists(path):
        print(f'  WARNING: missing {os.path.basename(path)}' + (f' [{label}]' if label else ''))
        return np.array([]), np.array([])
    with fits.open(path) as h:
        d  = h[1].data
        wl = d['WAVELENGTH'].astype(float)
        fl = d['FLUX'].astype(float)
        dq = d['DQ'].astype(int) if 'DQ' in d.dtype.names else np.zeros(len(wl), int)
    good = (wl > 0) & np.isfinite(fl) & (fl > 0) & (dq == 0)
    return wl[good], fl[good]

def load_calspec_jy(fname):
    with fits.open(f'{CALSPEC_DIR}/{fname}') as h:
        d    = h[1].data
        wl_a = d['WAVELENGTH'].astype(float)
        flam = d['FLUX'].astype(float)
    fnu_jy = flam * (wl_a**2) / C_ANG_S * 1e23
    order  = np.argsort(wl_a)
    return wl_a[order] / 1e4, fnu_jy[order]

def box_smooth(y, box=40):
    s = np.convolve(y, np.ones(box) / box, mode='same')
    s[:box//2]  = y[:box//2]
    s[-box//2:] = y[-box//2:]
    return s


# ── Core solver ────────────────────────────────────────────────────────────────
def solve_ifu(grating, wav_grid):
    """
    Solve for k(λ), α(λ), β(λ) using IFU stage3_ext NRS2 vs next-grating truth.

    For G140M:
        S_obs  = G140M stage3_ext NRS2 portion
        S_truth = G235M stage3 nominal
        f(λ)   = CALSPEC model

    For G235M:
        S_obs  = G235M stage3_ext NRS2 portion
        S_truth = G395M FS nominal
        f(λ)   = CALSPEC model

    At each wavelength λ, NNLS solves:
        S_obs(λ) ≈ k·f(λ) + α·f(λ/2) + β·f(λ/3)
    across three stars simultaneously.
    """
    print(f'\n--- Solving {grating} (IFU v1, 3 sources) ---')

    all_data = []

    for pid, src in SOURCES.items():
        name    = src['name']
        ifu_dir = src['ifu_dir']

        # Load CALSPEC as intrinsic reference
        wl_cs, fj_cs = load_calspec_jy(src['calspec'])
        f_ref = interp1d(wl_cs, fj_cs, bounds_error=False, fill_value=np.nan)

        if grating == 'G140M':
            # Observed: IFU G140M stage3_ext (full NRS1+NRS2 in Jy)
            x1d_ext  = f'{ifu_dir}/stage3_ext/f100lp_g140m-f100lp_x1d.fits'
            # Truth: IFU G235M stage3 nominal (Jy)
            x1d_truth = f'{ifu_dir}/stage3/f170lp_g235m-f170lp_x1d.fits'
            # NRS2 boundary for G140M
            nrs2_lo = 1.87

        else:  # G235M
            x1d_ext   = f'{ifu_dir}/stage3_ext/f170lp_g235m-f170lp_x1d.fits'
            x1d_truth = src['g395m_fs']
            nrs2_lo   = 3.15

        w_full, f_full = load_spec(x1d_ext, f'{name} stage3_ext')
        w_truth, f_truth = load_spec(x1d_truth, f'{name} truth')

        if len(w_full) < 10 or len(w_truth) < 10:
            print(f'  Skipping {name}: insufficient data')
            continue

        # Extract only the NRS2 portion of stage3_ext
        nrs2_mask = w_full >= nrs2_lo
        w_e = w_full[nrs2_mask]
        f_e = f_full[nrs2_mask]

        if len(w_e) < 10:
            print(f'  Skipping {name}: no NRS2 pixels above {nrs2_lo} µm')
            continue

        interp_ext   = interp1d(w_e,      f_e,      bounds_error=False, fill_value=np.nan)
        interp_truth = interp1d(w_truth,  f_truth,  bounds_error=False, fill_value=np.nan)

        print(f'  {name}: NRS2 {w_e.min():.3f}–{w_e.max():.3f} µm ({len(w_e)} px), '
              f'truth {w_truth.min():.3f}–{w_truth.max():.3f} µm')

        all_data.append({
            'name':        name,
            'ext':         interp_ext,
            'truth':       interp_truth,
            'calspec':     f_ref,
        })

    if len(all_data) == 0:
        raise RuntimeError(f'No valid source data for {grating}')

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

            S_obs = d['ext'](L)
            if not np.isfinite(S_obs) or S_obs <= 0:
                continue

            # Use truth to verify the calibration scale at this wavelength
            # (truth is used only to assess quality, not in the NNLS solve itself)
            # The model is: S_obs = k*f_cs + alpha*f_cs(lambda/2) + beta*f_cs(lambda/3)
            # Normalise by CALSPEC to make it dimensionless
            norm = f_cs
            A_rows.append([1.0, f_cs2 / norm, f_cs3 / norm])
            b_vec.append(S_obs / norm)

        if len(A_rows) == 0:
            ks.append(1.0); alphas.append(0.0); betas.append(0.0)
            continue

        # Regularized NNLS: anchor k≈1.0, allow ghost coefficients
        k_anchor = 2.0   # softer anchor — let data speak (IFU NRS2 is ~2x under)
        g_reg    = 0.01  # weak ghost regularization

        A_ext = np.vstack([
            np.array(A_rows),
            [k_anchor, 0.0, 0.0],   # prior: k ≈ 1.0 (with weight k_anchor)
            [0.0, g_reg, 0.0],
            [0.0, 0.0, g_reg],
        ])
        b_ext = np.hstack([
            np.array(b_vec),
            [k_anchor * 1.0],        # k prior = 1.0
            [0.0],
            [0.0],
        ])

        sol, _ = nnls(A_ext, b_ext)
        ks.append(sol[0])
        alphas.append(sol[1])
        betas.append(sol[2])

    # ── Smooth coefficients (Parlanti: 40-channel box) ─────────────────────────
    ks_s     = box_smooth(np.array(ks),     40)
    alphas_s = box_smooth(np.array(alphas), 40)
    betas_s  = box_smooth(np.array(betas),  40)

    return ks_s, alphas_s, betas_s, all_data


# ── Diagnostic plot ────────────────────────────────────────────────────────────
def plot_results(grating, wav_grid, ks, alphas, betas, all_data, nrs2_lo):
    fig, axes = plt.subplots(2, 1, figsize=(12, 9), sharex=True)
    fig.suptitle(f'IFU Calibration v1 — {grating} NRS2\nk/α/β coefficients + truth comparison',
                 fontsize=13)

    colors = {'P330E': '#f39c12', 'G191-B2B': '#2980b9', 'J1743045': '#27ae60'}

    # ── Top: k(λ) ──────────────────────────────────────────────────────────────
    ax = axes[0]
    ax.axhline(1.0, color='0.6', lw=0.8, ls='--')
    ax.plot(wav_grid, ks, color='k', lw=2.0, label='k(λ)  [smoothed, IFU v1]')
    ax.set_ylabel('k(λ)', fontsize=11)
    ax.set_ylim(0, 3.5)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.2)
    ax.set_title('First-order throughput factor', fontsize=10)

    # Also plot alpha / beta
    ax2 = ax.twinx()
    ax2.plot(wav_grid, alphas, color='darkorange', lw=1.2, ls=':', label='α(λ)')
    ax2.plot(wav_grid, betas,  color='purple',     lw=1.2, ls=':', label='β(λ)')
    ax2.set_ylabel('α / β', fontsize=10, color='darkorange')
    ax2.set_ylim(-0.1, 1.0)
    ax2.legend(fontsize=8, loc='upper right')

    # ── Bottom: S_obs vs S_truth after applying k correction ──────────────────
    ax = axes[1]
    k_interp = interp1d(wav_grid, ks, bounds_error=False, fill_value=np.nan)
    a_interp = interp1d(wav_grid, alphas, bounds_error=False, fill_value=np.nan)
    b_interp = interp1d(wav_grid, betas,  bounds_error=False, fill_value=np.nan)

    for d in all_data:
        name  = d['name']
        color = colors.get(name, 'gray')
        lo_wl = nrs2_lo + 0.05
        hi_wl = wav_grid.max() - 0.05
        w_eval = np.linspace(lo_wl, hi_wl, 400)

        f_obs   = d['ext'](w_eval)
        f_truth = d['truth'](w_eval)
        f_cs    = d['calspec'](w_eval)
        f_cs2   = d['calspec'](w_eval / 2)
        f_cs3   = d['calspec'](w_eval / 3)

        ki = k_interp(w_eval)
        ai = a_interp(w_eval)
        bi = b_interp(w_eval)

        ok = np.isfinite(f_obs) & np.isfinite(f_truth) & (f_obs > 0) & (f_truth > 0)
        if ok.sum() < 5:
            continue

        # Raw ratio (before correction)
        ratio_raw = f_obs[ok] / f_truth[ok]
        # Corrected: apply full Parlanti correction, then compare to truth
        f_corr = (f_obs - ai * f_cs2 - bi * f_cs3) / ki
        ratio_corr = f_corr[ok] / f_truth[ok]

        from scipy.ndimage import median_filter
        def sm(y): return median_filter(y, size=21)

        ax.plot(w_eval[ok], sm(ratio_raw),  color=color, lw=1.0, alpha=0.4, ls='--', label=f'{name} raw')
        ax.plot(w_eval[ok], sm(ratio_corr), color=color, lw=1.5, alpha=0.9, label=f'{name} corrected')

    ax.axhline(1.0, color='k', lw=1.0)
    ax.axhspan(0.90, 1.10, color='limegreen', alpha=0.07)
    ax.set_ylabel('S_corrected / S_truth', fontsize=11)
    ax.set_xlabel('Wavelength (µm)', fontsize=11)
    ax.set_ylim(0.3, 2.0)
    ax.legend(fontsize=8, ncol=2)
    ax.grid(True, alpha=0.2)
    ax.set_title('Ratio to truth reference (1.0 = perfect calibration)', fontsize=10)

    plt.tight_layout()
    outpath = f'{OUTPUT_DIR}/ifu_v1_{grating.lower()}_coeffs.png'
    plt.savefig(outpath, dpi=200, bbox_inches='tight')
    plt.close()
    print(f'  Plot saved: {outpath}')


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    print('=== IFU Calibration v1 Solver ===')
    print(f'Output: {OUTPUT_DIR}')

    results = {}

    for grating, wav_grid, nrs2_lo in [
        ('G140M', GRID_G140M, 1.87),
        ('G235M', GRID_G235M, 3.15),
    ]:
        ks, alphas, betas, all_data = solve_ifu(grating, wav_grid)

        # Save coefficients
        data = np.column_stack([wav_grid, ks, alphas, betas])
        csv_path = f'{OUTPUT_DIR}/coeffs_ifu_v1_{grating}.csv'
        np.savetxt(csv_path, data, delimiter=',', header='wav,k,alpha,beta', comments='')
        print(f'  Saved: {csv_path}')

        # Print summary statistics
        ok = np.isfinite(ks)
        print(f'  k({grating}): median={np.median(ks[ok]):.3f}  '
              f'range=[{ks[ok].min():.3f}, {ks[ok].max():.3f}]')
        print(f'  α({grating}): median={np.median(alphas[ok]):.4f}')
        print(f'  β({grating}): median={np.median(betas[ok]):.4f}')

        results[grating] = (wav_grid, ks, alphas, betas)

        plot_results(grating, wav_grid, ks, alphas, betas, all_data, nrs2_lo)

    # Combined coefficients summary plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for ax, (grating, (wav_grid, ks, alphas, betas)) in zip(axes, results.items()):
        ax.plot(wav_grid, ks,     lw=2.0, label='k(λ)')
        ax.plot(wav_grid, alphas, lw=1.5, ls='--', label='α(λ)')
        ax.plot(wav_grid, betas,  lw=1.5, ls=':',  label='β(λ)')
        ax.axhline(1.0, color='0.7', lw=0.8)
        ax.axhline(0.0, color='0.7', lw=0.5, ls=':')
        ax.set_title(f'{grating} NRS2', fontsize=12)
        ax.set_xlabel('Wavelength (µm)', fontsize=11)
        ax.set_ylabel('Coefficient value', fontsize=11)
        ax.set_ylim(-0.1, 3.5)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.25)
    fig.suptitle('IFU v1 — Parlanti k/α/β coefficients (3 standard stars)', fontsize=13)
    plt.tight_layout()
    summary_path = f'{OUTPUT_DIR}/ifu_v1_coefficients_summary.png'
    plt.savefig(summary_path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f'\n  Summary plot: {summary_path}')
    print('\nDone.')


if __name__ == '__main__':
    main()
