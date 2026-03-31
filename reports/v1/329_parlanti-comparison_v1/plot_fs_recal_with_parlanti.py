import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import interp1d

# Paths
BASE        = '/Users/dcoe/NIRSpec/wavext'
DATA_DIR    = f'{BASE}/data'
CALSPEC_DIR = f'{BASE}/data/CALSPEC'
PARLANTI_DIR = f'{BASE}/data/parlanti_repo/calibration_files'
OUTPUT_DIR   = f'{BASE}/nirspec_wavext_work/reports/329_parlanti-comparison_v1'
os.makedirs(OUTPUT_DIR, exist_ok=True)

C_ANG_S = 2.99792458e18

# Source data - Use the EXTENDED NRS2 files produced by our internal pipeline
SOURCES = {
    '1538': {
        'name':     'P330E',
        'calspec':  'p330e_mod_008.fits',
        'nrs2_g140m': f'{DATA_DIR}/PID1538_P330E/nrs2_spec2_cal/jw01538160001_06101_00001_nrs2_x1d.fits',
        'nrs2_g235m': f'{DATA_DIR}/PID1538_P330E/nrs2_spec2_cal/jw01538160001_08101_00005_nrs2_x1d.fits',
        'truth_g235m': f'{DATA_DIR}/PID1538_P330E/jw01538-o160_t002-s000000001_nirspec_f170lp-g235m-s1600a1-sub2048_x1d.fits',
        'truth_g395m': f'{DATA_DIR}/PID1538_P330E/jw01538-o160_t002-s000000001_nirspec_f290lp-g395m-s1600a1-sub2048_x1d.fits',
    },
    '1537': {
        'name':     'G191-B2B',
        'calspec':  'g191b2b_mod_012.fits',
        'nrs2_g140m': f'{DATA_DIR}/PID1537_G191-B2B/nrs2_spec2_cal/jw01537007001_07101_00003_nrs2_x1d.fits',
        'nrs2_g235m': f'{DATA_DIR}/PID1537_G191-B2B/nrs2_spec2_cal/jw01537007001_09101_00003_nrs2_x1d.fits',
        'truth_g235m': f'{DATA_DIR}/PID1537_G191-B2B/jw01537-o007_t001-s000000001_nirspec_f170lp-g235m-s1600a1-sub2048_x1d.fits',
        'truth_g395m': f'{DATA_DIR}/PID1537_G191-B2B/jw01537-o007_t001-s000000001_nirspec_f290lp-g395m-s1600a1-sub2048_x1d.fits',
    },
    '1536': {
        'name':     'J1743045',
        'calspec':  '1743045_mod_007.fits',
        'nrs2_g140m': f'{DATA_DIR}/PID1536_J1743045/nrs2_spec2_cal/jw01536002001_03104_00003_nrs2_x1d.fits',
        'nrs2_g235m': f'{DATA_DIR}/PID1536_J1743045/nrs2_spec2_cal/jw01536002001_03106_00001_nrs2_x1d.fits',
        'truth_g235m': f'{DATA_DIR}/PID1536_J1743045/jw01536-o002_t004-s000000001_nirspec_f170lp-g235m-s1600a1-sub2048_x1d.fits',
        'truth_g395m': f'{DATA_DIR}/PID1536_J1743045/jw01536-o002_t004-s000000001_nirspec_f290lp-g395m-s1600a1-sub2048_x1d.fits',
    },
}

def load_spec(path):
    if not os.path.exists(path):
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

def plot_recal_with_parlanti(sid, src):
    name = src['name']
    print(f'Recalibrating {name} spectra with Parlanti coeffs...')
    
    fig, axes = plt.subplots(2, 1, figsize=(14, 10), sharex=False)
    fig.suptitle(f'FS Recalibration — Parlanti Original Coeffs applied to {name}', fontsize=14)

    wl_cs, fj_cs = load_calspec_jy(src['calspec'])
    f_ref = interp1d(wl_cs, fj_cs, bounds_error=False, fill_value=np.nan)

    gratings = ['G140M', 'G235M']
    fnames_parlanti = ['calibration_functions_g140m_f100lp.fits', 'calibration_functions_g235m_f170lp.fits']
    
    for i, (grating, p_fname) in enumerate(zip(gratings, fnames_parlanti)):
        ax = axes[i]
        
        with fits.open(f'{PARLANTI_DIR}/{p_fname}') as h:
            d_p = h[1].data
            k_p_interp = interp1d(d_p['wavelength'], d_p['k'],     bounds_error=False, fill_value=np.nan)
            a_p_interp = interp1d(d_p['wavelength'], d_p['alpha'], bounds_error=False, fill_value=np.nan)
            b_p_interp = interp1d(d_p['wavelength'], d_p['beta'],  bounds_error=False, fill_value=np.nan)

        obs_path = src[f'nrs2_{grating.lower()}']
        truth_path = src[f'truth_{"g235m" if grating == "G140M" else "g395m"}']
        
        w_obs, f_obs = load_spec(obs_path)
        if len(w_obs) == 0:
            print(f"No extended data for {grating}")
            continue
            
        w_truth, f_truth = load_spec(truth_path)
        if grating == 'G140M': mask = (w_truth > 1.95) & (w_truth < 3.50)
        else: mask = (w_truth > 3.30) & (w_truth < 5.25)
        w_truth_overlap, f_truth_overlap = w_truth[mask], f_truth[mask]

        ki = k_p_interp(w_obs)
        ai = a_p_interp(w_obs)
        bi = b_p_interp(w_obs)
        f_cs2 = f_ref(w_obs / 2)
        f_cs3 = f_ref(w_obs / 3)
        
        f_corr = (f_obs - ai * f_cs2 - bi * f_cs3) / ki

        ax.plot(w_obs,  f_obs,  'k', lw=1.0, alpha=0.3, label=f'Observed {grating} (NRS2 Extended)')
        if len(w_truth_overlap)>0:
            ax.plot(w_truth_overlap, f_truth_overlap, 'navy', lw=1.5, ls='--', alpha=0.7, label='Truth (Next Grating)')
        ax.plot(w_obs,  f_corr,  '#e67e22', lw=1.8, label='Recalibrated (Parlanti Original Coeffs)')
        
        w_cs_plot = np.linspace(w_obs.min(), w_obs.max(), 500)
        ax.plot(w_cs_plot, f_ref(w_cs_plot), '0.5', lw=0.8, ls=':', label='CALSPEC model')

        f_cs_overlap = f_ref(w_obs)
        vals = []
        if len(w_truth_overlap)>0: vals.append(f_truth_overlap)
        vals.append(f_cs_overlap[np.isfinite(f_cs_overlap) & (f_cs_overlap > 0)])
        if vals:
            all_vals = np.concatenate(vals)
            ax.set_ylim(np.nanmin(all_vals) * 0.5, np.nanmax(all_vals) * 2.5)

        ax.set_yscale('log')
        ax.set_ylabel('Flux [Jy]', fontsize=12)
        ax.set_xlabel('Wavelength [µm]', fontsize=12)
        ax.set_title(f'FS {grating} (Extended)', fontsize=12)
        ax.grid(True, which='both', alpha=0.1)
        ax.legend(fontsize=9, ncol=2)

    plt.tight_layout()
    outpath = f'{OUTPUT_DIR}/parlanti_recal_spectra_{name}.png'
    plt.savefig(outpath, dpi=200, bbox_inches='tight')
    plt.close()
    print(f'Saved: {outpath}')

if __name__ == '__main__':
    for sid, src in SOURCES.items():
        plot_recal_with_parlanti(sid, src)
