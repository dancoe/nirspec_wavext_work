"""
Derive Parlanti et al. (2025) k(λ), α(λ), β(λ) calibration coefficients from
JWST NIRSpec IFU observations of spectrophotometric standards.

Physical model (per source i, per wavelength λ):
    S_obs_i(λ) = k(λ)·S_cal_i(λ) + α(λ)·S_cal_i(λ/2) + β(λ)·S_cal_i(λ/3)

where S_cal_i(λ) is the CALSPEC model SED interpolated to the pipeline wavelength
grid, and S_obs_i(λ) is the extracted x1d flux.  For each λ this is a linear
system in (k, α, β), solved by least squares across the 4 available sources.

Sources:
    G191-B2B   → PID 1537  (hot WD)
    P330E      → PID 1538  (G-star)
    J1743045   → PID 1536  (solar analog)
    P330E      → PID 6645  (independent observation)

Outputs:
    results/calibration/calib_g140m_f100lp.fits
    results/calibration/calib_g235m_f170lp.fits

The output FITS tables match the Parlanti calibration_functions_*.fits format:
columns [WAVELENGTH (µm), K, ALPHA, BETA].
"""
import os
import numpy as np
import astropy.io.fits as fits
from scipy.ndimage import uniform_filter1d
from scipy.interpolate import interp1d

# ── Paths ───────────────────────────────────────────────────────────────────
BASE          = '/Users/dcoe/NIRSpec/wavext'
IFU_DIR       = f'{BASE}/data/IFU'
CALSPEC_DIR   = f'{BASE}/data/CALSPEC'
PARLANTI_DIR  = f'{BASE}/data/parlanti_repo/calibration_files'
OUTDIR        = f'{BASE}/results/calibration'
os.makedirs(OUTDIR, exist_ok=True)

# ── Source configuration ─────────────────────────────────────────────────────
# Maps PID → (target name, CALSPEC filename)
# NOTE: PID 6645 targets "GSC 02581-02323" (a faint field star near P330E),
#       NOT P330E itself. Excluded from calibration solver.
SOURCES = {
    '1537': ('G191-B2B', 'g191b2b_mod_012.fits'),   # hot WD      V=11.8
    '1538': ('P330E',    'p330e_mod_008.fits'),       # G2V solar   V=12.9
    '1536': ('J1743045', '1743045_mod_007.fits'),     # A8III       V=13.5
}

# ── Grating configurations ────────────────────────────────────────────────────
GRATING_CONFIGS = {
    'g140m_f100lp': {
        'label':     'G140M/F100LP',
        'x1d_stem':  'f100lp_g140m-f100lp_x1d.fits',
        'stage':     'stage3_ext',
        'parlanti':  'calibration_functions_g140m_f100lp.fits',
        'pid_dirs': {
            '1537': 'PID1537_G191-B2B',
            '1538': 'PID1538_P330E',
            '1536': 'PID1536_J1743045',
        },
    },
    'g235m_f170lp': {
        'label':     'G235M/F170LP',
        'x1d_stem':  'f170lp_g235m-f170lp_x1d.fits',
        'stage':     'stage3_ext',
        'parlanti':  'calibration_functions_g235m_f170lp.fits',
        'pid_dirs': {
            '1537': 'PID1537_G191-B2B',
            '1538': 'PID1538_P330E',
            '1536': 'PID1536_J1743045',
        },
    },
}

# Smoothing half-width in wavelength channels (Parlanti use 40 channels)
SMOOTH_WIDTH = 40

# ── Unit conversion ───────────────────────────────────────────────────────────
C_ANGSTROM_PER_S = 2.99792458e18   # speed of light in Å/s

def flam_to_jy(wavelength_angstrom, flux_flam):
    """Convert FLAM [erg/s/cm²/Å] → Fν [Jy], wavelength in Å."""
    # F_ν [erg/s/cm²/Hz] = F_λ × λ² / c
    # 1 Jy = 1e-23 erg/s/cm²/Hz
    fnu_cgs = flux_flam * wavelength_angstrom**2 / C_ANGSTROM_PER_S
    return fnu_cgs * 1e23   # Jy


def load_calspec(calspec_file):
    """Load a CALSPEC FITS file and return (wavelength_um, flux_jy) arrays."""
    hdul = fits.open(calspec_file)
    d = hdul[1].data
    wl_ang = d['WAVELENGTH'].astype(np.float64)   # Å
    f_flam  = d['FLUX'].astype(np.float64)         # erg/s/cm²/Å
    hdul.close()
    # Convert to µm and Jy
    wl_um  = wl_ang / 1e4
    f_jy   = flam_to_jy(wl_ang, f_flam)
    # Sort (should already be sorted, but ensure)
    order  = np.argsort(wl_um)
    return wl_um[order], f_jy[order]


def load_x1d(x1d_file):
    """Load an x1d FITS file. Return (wavelength_um, flux_jy, dq) arrays."""
    hdul = fits.open(x1d_file)
    d = hdul[1].data
    wl   = d['WAVELENGTH'].astype(np.float64)    # µm
    flux = d['FLUX'].astype(np.float64)           # Jy
    dq   = d['DQ'].astype(np.int32)
    hdul.close()
    return wl, flux, dq


def interp_calspec(calspec_wl, calspec_flux, target_wl):
    """Interpolate CALSPEC onto target wavelength grid (all in µm)."""
    f = interp1d(calspec_wl, calspec_flux,
                 kind='linear', bounds_error=False, fill_value=np.nan)
    return f(target_wl)


def solve_calib(wl_grid, obs_matrix, cal_matrix):
    """
    At each wavelength, solve:
        obs_matrix[:, i] = k[i]*cal_matrix[:, 3*i] + alpha[i]*cal_matrix[:, 3*i+1] + beta[i]*cal_matrix[:, 3*i+2]

    Parameters
    ----------
    wl_grid : (N,) wavelength array in µm
    obs_matrix  : (N_src, N) observed flux at each wavelength [Jy]
    cal_matrix  : (N_src, N, 3) CALSPEC at λ, λ/2, λ/3 [Jy]

    Returns
    -------
    k, alpha, beta : (N,) arrays of calibration coefficients
    valid          : (N,) bool mask of pixels with valid solutions
    """
    N   = len(wl_grid)
    k   = np.full(N, np.nan)
    alpha  = np.full(N, np.nan)
    beta   = np.full(N, np.nan)
    valid  = np.zeros(N, dtype=bool)

    for i in range(N):
        y  = obs_matrix[:, i]        # (N_src,) observed fluxes
        A  = cal_matrix[:, i, :]     # (N_src, 3) design matrix

        # Require at least MIN_SRC valid sources per wavelength
        finite = np.isfinite(y) & np.all(np.isfinite(A), axis=1) & (y > 0) & (np.all(A > 0, axis=1))
        if finite.sum() < 2:
            continue

        y_fit = y[finite]
        A_fit = A[finite]

        # Non-negative least squares (physically: k, alpha, beta >= 0)
        try:
            from scipy.optimize import lsq_linear
            bounds = ([0, 0, 0], [np.inf, np.inf, np.inf])
            res = lsq_linear(A_fit, y_fit, bounds=bounds, method='bvls')
            result = res.x
        except Exception:
            try:
                result, _, _, _ = np.linalg.lstsq(A_fit, y_fit, rcond=None)
            except np.linalg.LinAlgError:
                continue

        k_i, alpha_i, beta_i = result
        # Basic sanity: k must be positive
        if k_i > 0:
            k[i]     = k_i
            alpha[i] = alpha_i
            beta[i]  = beta_i
            valid[i] = True

    return k, alpha, beta, valid


def smooth_coeffs(arr, width, valid):
    """Apply uniform (box) smoothing of width `width` channels, ignoring NaNs."""
    arr_copy = arr.copy()
    arr_copy[~valid] = np.nan
    # Fill NaN with interpolated values for smooth filter, then restore
    from scipy.interpolate import interp1d as itp
    x = np.arange(len(arr_copy))
    ok = valid & np.isfinite(arr_copy)
    if ok.sum() < 2:
        return arr_copy
    interper = itp(x[ok], arr_copy[ok], kind='linear', bounds_error=False,
                   fill_value=(arr_copy[ok][0], arr_copy[ok][-1]))
    filled = interper(x)
    smoothed = uniform_filter1d(filled, size=width, mode='nearest')
    smoothed[~valid] = np.nan
    return smoothed


def save_calib_fits(outpath, wl_um, k, alpha, beta, grating_label):
    """Save calibration coefficients to a FITS binary table."""
    col_wl    = fits.Column(name='WAVELENGTH', format='D', unit='um',    array=wl_um)
    col_k     = fits.Column(name='K',          format='D', unit='',      array=k)
    col_alpha = fits.Column(name='ALPHA',      format='D', unit='',      array=alpha)
    col_beta  = fits.Column(name='BETA',       format='D', unit='',      array=beta)

    hdr = fits.Header()
    hdr['GRATING']  = grating_label
    hdr['ORIGIN']   = 'solve_parlanti_calib.py'
    hdr['DESCRIP']  = 'NIRSpec IFU higher-order contamination correction'
    hdr['SMOOTH']   = (SMOOTH_WIDTH, 'channel smoothing full-width')
    hdr['NSOURCES'] = (len(SOURCES), 'number of calibration sources')

    table_hdu = fits.BinTableHDU.from_columns([col_wl, col_k, col_alpha, col_beta],
                                               header=hdr)
    primary = fits.PrimaryHDU()
    hdul = fits.HDUList([primary, table_hdu])
    hdul.writeto(outpath, overwrite=True)
    print(f'  → written {outpath}')


def run_grating(key, cfg):
    """Run the full solve pipeline for one grating/filter configuration."""
    label = cfg['label']
    print(f'\n{"="*60}')
    print(f'{label}')
    print('='*60)

    # ── Load x1d spectra ────────────────────────────────────────────────────
    x1d_data = {}   # pid → (wl, flux, dq)
    wl_ref = None   # reference wavelength grid (from first loaded x1d)
    for pid, pid_dir in cfg['pid_dirs'].items():
        x1d_path = f'{IFU_DIR}/{pid_dir}/{cfg["stage"]}/{cfg["x1d_stem"]}'
        if not os.path.exists(x1d_path):
            print(f'  WARNING: missing x1d for PID {pid}: {x1d_path}')
            continue
        wl, flux, dq = load_x1d(x1d_path)
        x1d_data[pid] = (wl, flux, dq)
        print(f'  PID {pid} ({SOURCES[pid][0]}): {len(wl)} channels, '
              f'{wl.min():.4f}-{wl.max():.4f} µm')
        if wl_ref is None:
            wl_ref = wl

    if wl_ref is None or len(x1d_data) == 0:
        print('  ERROR: No x1d data loaded. Skipping.')
        return

    N      = len(wl_ref)
    n_src  = len(x1d_data)
    pids = sorted(x1d_data.keys())

    # ── Load CALSPEC models ──────────────────────────────────────────────────
    calspec_cache = {}  # filename → (wl_um, flux_jy)
    for pid in pids:
        _, calspec_fname = SOURCES[pid]
        if calspec_fname not in calspec_cache:
            cpath = f'{CALSPEC_DIR}/{calspec_fname}'
            calspec_cache[calspec_fname] = load_calspec(cpath)
            print(f'  Loaded CALSPEC {calspec_fname}: '
                  f'{calspec_cache[calspec_fname][0].min():.3f}-'
                  f'{calspec_cache[calspec_fname][0].max():.3f} µm')

    # ── Build obs_matrix and cal_matrix ─────────────────────────────────────
    # obs_matrix : (n_src, N)    observed Jy at wl_ref
    # cal_matrix : (n_src, N, 3) CALSPEC at λ, λ/2, λ/3 [Jy]
    obs_matrix = np.full((n_src, N), np.nan)
    cal_matrix = np.full((n_src, N, 3), np.nan)

    for src_idx, pid in enumerate(pids):
        wl_i, flux_i, dq_i  = x1d_data[pid]
        _, calspec_fname = SOURCES[pid]
        cs_wl, cs_flux = calspec_cache[calspec_fname]

        # Interpolate observed x1d onto reference grid (usually identical)
        if np.allclose(wl_i, wl_ref, atol=1e-7):
            obs = flux_i.copy()
        else:
            obs = interp_calspec(wl_i, flux_i, wl_ref)

        # Mask bad DQ pixels and non-positive flux
        bad = (dq_i > 0) | ~np.isfinite(obs) | (obs <= 0)
        obs[bad] = np.nan
        obs_matrix[src_idx] = obs

        # Interpolate CALSPEC at λ, λ/2, λ/3
        for order_idx, divisor in enumerate([1, 2, 3]):
            wl_query = wl_ref / divisor
            cal_matrix[src_idx, :, order_idx] = interp_calspec(cs_wl, cs_flux, wl_query)

    # Report coverage
    valid_obs = np.sum(np.isfinite(obs_matrix), axis=0)
    print(f'  Wavelength bins with ≥3 valid sources: '
          f'{np.sum(valid_obs >= 3)}/{N}')

    # ── Solve per wavelength ─────────────────────────────────────────────────
    print('  Solving for k, α, β ...')
    k_raw, alpha_raw, beta_raw, valid_mask = solve_calib(wl_ref, obs_matrix, cal_matrix)
    n_valid = valid_mask.sum()
    print(f'  Valid solutions: {n_valid}/{N} ({100*n_valid/N:.1f}%)')

    if n_valid < 10:
        print('  ERROR: too few valid solutions, check inputs.')
        return

    # ── Smooth ───────────────────────────────────────────────────────────────
    print(f'  Smoothing with width={SMOOTH_WIDTH} channels ...')
    k_smooth     = smooth_coeffs(k_raw,     SMOOTH_WIDTH, valid_mask)
    alpha_smooth = smooth_coeffs(alpha_raw, SMOOTH_WIDTH, valid_mask)
    beta_smooth  = smooth_coeffs(beta_raw,  SMOOTH_WIDTH, valid_mask)

    # ── Print statistics ─────────────────────────────────────────────────────
    ok = valid_mask
    print(f'  k:     [{np.nanmin(k_smooth[ok]):.3f}, {np.nanmax(k_smooth[ok]):.3f}]  '
          f'median={np.nanmedian(k_smooth[ok]):.3f}')
    print(f'  alpha: [{np.nanmin(alpha_smooth[ok]):.4f}, {np.nanmax(alpha_smooth[ok]):.4f}]  '
          f'median={np.nanmedian(alpha_smooth[ok]):.4f}')
    print(f'  beta:  [{np.nanmin(beta_smooth[ok]):.5f}, {np.nanmax(beta_smooth[ok]):.5f}]  '
          f'median={np.nanmedian(beta_smooth[ok]):.5f}')

    # ── Save ─────────────────────────────────────────────────────────────────
    outpath = f'{OUTDIR}/calib_{key}.fits'
    save_calib_fits(outpath, wl_ref, k_smooth, alpha_smooth, beta_smooth, label)

    # Also save unsmoothed for diagnostics
    outpath_raw = f'{OUTDIR}/calib_{key}_raw.fits'
    save_calib_fits(outpath_raw, wl_ref, k_raw, alpha_raw, beta_raw,
                    label + ' (unsmoothed)')

    return {
        'wl':     wl_ref,
        'k':      k_smooth,
        'alpha':  alpha_smooth,
        'beta':   beta_smooth,
        'k_raw':  k_raw,
        'alpha_raw': alpha_raw,
        'beta_raw':  beta_raw,
        'valid':  valid_mask,
    }


def main():
    print('Parlanti k/α/β solver')
    print(f'Sources: {list(SOURCES.keys())}')
    print(f'Output: {OUTDIR}')

    results = {}
    for key, cfg in GRATING_CONFIGS.items():
        results[key] = run_grating(key, cfg)

    print('\nDone.')
    return results


if __name__ == '__main__':
    main()
