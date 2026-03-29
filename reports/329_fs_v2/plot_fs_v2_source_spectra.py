"""
FS v2 source spectra validation plots.

For each standard star (P330E, G191-B2B, J1743045), compares:
  - Observed FS NRS2 per-exposure x1d (gray)
  - Truth: MAST L3 G235M (for G140M NRS2) or G395M (for G235M NRS2) (blue dashed)
  - Recalibrated: (S_obs - α·f(λ/2) - β·f(λ/3)) / k(λ)  [FS v2 coefficients] (red)
  - CALSPEC model (grey dotted)

Saves one PNG per source.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import interp1d

BASE        = '/Users/dcoe/NIRSpec/wavext'
FS_DIR      = f'{BASE}/data'
CALSPEC_DIR = f'{BASE}/data/CALSPEC'
COEFFS_DIR  = f'{BASE}/nirspec_wavext_work/plots/Parlanti/cal/fs_v2'
OUTPUT_DIR  = f'{BASE}/nirspec_wavext_work/reports/329_fs_v2'
os.makedirs(OUTPUT_DIR, exist_ok=True)

C_ANG_S = 2.99792458e18  # Å/s

SOURCES = {
    '1538': {
        'name':       'P330E',
        'calspec':    'p330e_mod_008.fits',
        'nrs2_g140m': f'{FS_DIR}/PID1538_P330E/nrs2_spec2_cal/jw01538160001_06101_00001_nrs2_x1d.fits',
        'nrs2_g235m': f'{FS_DIR}/PID1538_P330E/nrs2_spec2_cal/jw01538160001_08101_00005_nrs2_x1d.fits',
        'l3_g235m':   f'{FS_DIR}/PID1538_P330E/jw01538-o160_t002-s000000001_nirspec_f170lp-g235m-s1600a1-sub2048_x1d.fits',
        'l3_g395m':   f'{FS_DIR}/PID1538_P330E/jw01538-o160_t002-s000000001_nirspec_f290lp-g395m-s1600a1-sub2048_x1d.fits',
    },
    '1537': {
        'name':       'G191-B2B',
        'calspec':    'g191b2b_mod_012.fits',
        'nrs2_g140m': f'{FS_DIR}/PID1537_G191-B2B/nrs2_spec2_cal/jw01537007001_07101_00003_nrs2_x1d.fits',
        'nrs2_g235m': f'{FS_DIR}/PID1537_G191-B2B/nrs2_spec2_cal/jw01537007001_09101_00003_nrs2_x1d.fits',
        'l3_g235m':   f'{FS_DIR}/PID1537_G191-B2B/jw01537-o007_t001-s000000001_nirspec_f170lp-g235m-s1600a1-sub2048_x1d.fits',
        'l3_g395m':   f'{FS_DIR}/PID1537_G191-B2B/jw01537-o007_t001-s000000001_nirspec_f290lp-g395m-s1600a1-sub2048_x1d.fits',
    },
    '1536': {
        'name':       'J1743045',
        'calspec':    '1743045_mod_007.fits',
        'nrs2_g140m': f'{FS_DIR}/PID1536_J1743045/nrs2_spec2_cal/jw01536002001_03104_00003_nrs2_x1d.fits',
        'nrs2_g235m': f'{FS_DIR}/PID1536_J1743045/nrs2_spec2_cal/jw01536002001_03106_00001_nrs2_x1d.fits',
        'l3_g235m':   f'{FS_DIR}/PID1536_J1743045/jw01536-o002_t004-s000000001_nirspec_f170lp-g235m-s1600a1-sub2048_x1d.fits',
        'l3_g395m':   f'{FS_DIR}/PID1536_J1743045/jw01536-o002_t004-s000000001_nirspec_f290lp-g395m-s1600a1-sub2048_x1d.fits',
    },
}


def load_spec(path):
    if not os.path.exists(path):
        return np.array([]), np.array([])
    with fits.open(path) as h:
        ext = h[1].data
        wl  = ext['WAVELENGTH'].astype(float)
        fl  = ext['FLUX'].astype(float)
        dq  = ext['DQ'].astype(int) if 'DQ' in ext.dtype.names else np.zeros(len(wl), int)
    good = (wl > 0) & np.isfinite(fl) & (fl > 0) & (dq == 0)
    return wl[good], fl[good]


def load_calspec_jy(fname):
    with fits.open(f'{CALSPEC_DIR}/{fname}') as h:
        d    = h[1].data
        wl_a = d['WAVELENGTH'].astype(float)
        flam = d['FLUX'].astype(float)
    fnu   = flam * wl_a**2 / C_ANG_S * 1e23
    order = np.argsort(wl_a)
    return interp1d(wl_a[order] / 1e4, fnu[order],
                    bounds_error=False, fill_value=np.nan)


def load_coeffs(grating):
    csv_path = f'{COEFFS_DIR}/coeffs_fs_v2_{grating}.csv'
    data = np.loadtxt(csv_path, delimiter=',', skiprows=1)
    wl, k, a, b = data[:, 0], data[:, 1], data[:, 2], data[:, 3]
    ki = interp1d(wl, k, bounds_error=False, fill_value=np.nan)
    ai = interp1d(wl, np.nan_to_num(a, nan=0.0), bounds_error=False, fill_value=0.0)
    bi = interp1d(wl, np.nan_to_num(b, nan=0.0), bounds_error=False, fill_value=0.0)
    return ki, ai, bi


def plot_source(sid, src):
    name = src['name']
    print(f'  Plotting {name}...')

    f_ref = load_calspec_jy(src['calspec'])

    fig, axes = plt.subplots(2, 1, figsize=(14, 10), sharex=False)
    fig.suptitle(f'FS v2 Calibration (329) — {name} Spectra Comparison', fontsize=14)

    for i, grating in enumerate(['G140M', 'G235M']):
        ax = axes[i]
        k_interp, a_interp, b_interp = load_coeffs(grating)

        if grating == 'G140M':
            nrs2_path    = src['nrs2_g140m']
            truth_path   = src['l3_g235m']
            nrs2_lo      = 1.93
            truth_label  = 'Truth (Nominal G235M)'
            truth_lo, truth_hi = 1.95, 3.50
        else:
            nrs2_path    = src['nrs2_g235m']
            truth_path   = src['l3_g395m']
            nrs2_lo      = 3.25
            truth_label  = 'Truth (Nominal G395M)'
            truth_lo, truth_hi = 3.30, 5.25

        w_all, f_all = load_spec(nrs2_path)
        if len(w_all) < 5:
            ax.text(0.5, 0.5, f'{grating} data not available',
                    transform=ax.transAxes, ha='center')
            continue

        mask = w_all >= nrs2_lo
        w_obs, f_obs = w_all[mask], f_all[mask]
        if len(w_obs) < 5:
            ax.text(0.5, 0.5, 'No NRS2 data', transform=ax.transAxes, ha='center')
            continue

        # Truth
        w_tr, f_tr = load_spec(truth_path)
        if len(w_tr) > 5:
            tr_mask = (w_tr > truth_lo) & (w_tr < truth_hi)
            w_tr, f_tr = w_tr[tr_mask], f_tr[tr_mask]

        # Apply FS v2 correction
        ki    = k_interp(w_obs)
        ai    = a_interp(w_obs)
        bi    = b_interp(w_obs)
        f_cs2 = f_ref(w_obs / 2)
        f_cs3 = f_ref(w_obs / 3)
        f_cs2 = np.where(np.isfinite(f_cs2) & (f_cs2 > 0), f_cs2, 0.0)
        f_cs3 = np.where(np.isfinite(f_cs3) & (f_cs3 > 0), f_cs3, 0.0)

        ki_safe = np.where(np.isfinite(ki) & (ki > 0.01), ki, np.nan)
        f_corr  = (f_obs - ai * f_cs2 - bi * f_cs3) / ki_safe

        # Y-limits from truth + CALSPEC only
        w_cal = np.linspace(w_obs.min(), w_obs.max(), 500)
        f_cal = f_ref(w_cal)

        ref_vals = []
        if len(w_tr) > 5:
            ref_vals.append(f_tr)
        ref_vals.append(f_cal[np.isfinite(f_cal) & (f_cal > 0)])
        all_ref = np.concatenate(ref_vals)
        ymin = np.nanmin(all_ref) * 0.4
        ymax = np.nanmax(all_ref) * 3.0

        ax.plot(w_obs, f_obs,  color='0.6', lw=0.8, alpha=0.6,
                label=f'Observed {grating} (Extended NRS2)')
        if len(w_tr) > 5:
            ax.plot(w_tr, f_tr, color='navy', lw=1.8, ls='--', alpha=0.8,
                    label=truth_label)
        ax.plot(w_obs, f_corr, color='red', lw=1.5,
                label='Recalibrated (FS v2)')
        ax.plot(w_cal, f_cal, '0.5', lw=0.8, ls=':', label='CALSPEC model')

        ax.set_yscale('log')
        ax.set_ylim(ymin, ymax)
        ax.set_ylabel('Flux [Jy]', fontsize=11)
        ax.set_xlabel('Wavelength [µm]', fontsize=11)
        ax.set_title(f'{grating} NRS2 region', fontsize=11)
        ax.grid(True, which='both', alpha=0.1)
        ax.legend(fontsize=9, ncol=2)

    plt.tight_layout()
    outpath = f'{OUTPUT_DIR}/fs_v2_spectra_{name}.png'
    plt.savefig(outpath, dpi=200, bbox_inches='tight')
    plt.close()
    print(f'    Saved: {outpath}')


if __name__ == '__main__':
    print('Generating FS v2 source spectra plots...')
    for sid, src in SOURCES.items():
        plot_source(sid, src)
    print('Done.')
