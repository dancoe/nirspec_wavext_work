import os, sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import interp1d
import pandas as pd

# Paths
DATA_DIR = '/Users/dcoe/NIRSpec/wavext/data'
COEFF_DIR = '/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/plots/Parlanti/cal/153678_v3'
OUTPUT_DIR = COEFF_DIR

def load_spec(path):
    if not os.path.exists(path): return np.array([]), np.array([])
    try:
        with fits.open(path) as hdul:
            d = hdul[1].data
            w, f = d['wavelength'], d['flux']
            msk = (w > 0) & np.isfinite(f)
            return w[msk], f[msk]
    except: return np.array([]), np.array([])

def make_summary_v3(target_name):
    print(f"Generating summary plot V3 for {target_name}...")
    
    # Load coeffs
    c140 = pd.read_csv(f'{COEFF_DIR}/coeffs_G140M.csv')
    c235 = pd.read_csv(f'{COEFF_DIR}/coeffs_G235M.csv')
    
    # Target specific paths
    if target_name == 'P330E':
        dir = f'{DATA_DIR}/PID1538_P330E'
        files = {
            'G140M_ext': f'{dir}/jw01538160001_06101_00001_nrs2_extract_1d.fits',
            'G140M_nom': f'{dir}/jw01538-o160_t002-s000000001_nirspec_f100lp-g140m-s1600a1-sub2048_x1d.fits',
            'G235M_ext': f'{dir}/jw01538160001_08101_00005_nrs2_extract_1d.fits',
            'G235M_nom': f'{dir}/jw01538-o160_t002-s000000001_nirspec_f170lp-g235m-s1600a1-sub2048_x1d.fits',
            'G395M_nom': f'{dir}/jw01538-o160_t002-s000000001_nirspec_f290lp-g395m-s1600a1-sub2048_x1d.fits',
            'PRISM': f'{dir}/jw01538-o160_t002-s000000001_nirspec_clear-prism-s1600a1-sub2048_x1d.fits'
        }
    else: # G191-B2B
        dir = f'{DATA_DIR}/PID1537_G191-B2B'
        files = {
            'G140M_ext': f'{dir}/jw01537007001_07101_00003_nrs2_extract_1d.fits',
            'G140M_nom': f'{dir}/jw01537-o007_t001-s000000001_nirspec_f100lp-g140m-s1600a1-sub2048_x1d.fits',
            'G235M_ext': f'{dir}/jw01537007001_09101_00003_nrs2_extract_1d.fits',
            'G235M_nom': f'{dir}/jw01537-o007_t001-s000000001_nirspec_f170lp-g235m-s1600a1-sub2048_x1d.fits',
            'G395M_nom': f'{dir}/jw01537-o007_t001-s000000001_nirspec_f290lp-g395m-s1600a1-sub2048_x1d.fits',
            'PRISM': f'{dir}/jw01537-o007_t001-s000000001_nirspec_clear-prism-s1600a1-sub2048_x1d.fits'
        }

    plt.figure(figsize=(14, 12))
    plt.suptitle(f'NIRSpec Multi-Source Bootstrap Calibration (V3 - Kappa Fixed) — {target_name}', fontsize=16)

    # G140M Subplot
    plt.subplot(2,1,1)
    w_e, f_e = load_spec(files['G140M_ext'])
    w_n, f_n = load_spec(files['G140M_nom'])
    w_t, f_t = load_spec(files['G235M_nom']) # Primary truth
    w_t2, f_t2 = load_spec(files['G395M_nom']) # Secondary truth
    w_p, f_p = load_spec(files['PRISM']) # Reference

    # Scaling
    # Use a broader range or the start of the extension to anchor the scale
    m_e = (w_e > 1.8) & (w_e < 3.0) # Any overlap
    m_t = (w_t > 1.8) & (w_t < 3.0)
    
    # We want a narrow window near the start of the extension
    w_start = np.max([np.nanmin(w_e), np.nanmin(w_t)])
    window = 0.05
    mask_e = (w_e >= w_start) & (w_e < w_start + window)
    mask_t = (w_t >= w_start) & (w_t < w_start + window)
    
    if np.any(mask_e) and np.any(mask_t):
        scale = np.nanmedian(f_t[mask_t]) / np.nanmedian(f_e[mask_e])
    else:
        # Fallback: broad median
        scale = np.nanmedian(f_t) / np.nanmedian(f_e)
    
    f_e *= scale
    
    interp_e = interp1d(w_e, f_e, bounds_error=False, fill_value=0)
    interp_n = interp1d(w_n, f_n, bounds_error=False, fill_value=0)
    
    f_cal = (interp_e(c140['wav']) - c140['alpha']*interp_n(c140['wav']/2) - c140['beta']*interp_n(c140['wav']/3)) / c140['k']
    
    plt.plot(w_n, f_n, 'k-', lw=1.5, label='Nominal G140M')
    plt.plot(w_t, f_t, 'k-', alpha=0.3, label='Observed G235M (Truth)')
    plt.plot(w_t2, f_t2, 'k-', alpha=0.1, label='Observed G395M')
    plt.plot(w_e, f_e, 'g--', alpha=0.4, label='Extended G140M (Raw)')
    plt.plot(c140['wav'], f_cal, 'r-', lw=2, label='Extended G140M (Recalibrated)')
    
    # Arrows indicating extension
    plt.annotate('', xy=(3.5, np.nanmedian(f_n)), xytext=(1.8, np.nanmedian(f_n)), arrowprops=dict(arrowstyle='<->', color='blue', lw=2, alpha=0.3))
    plt.text(2.6, np.nanmedian(f_n)*1.2, 'Extended G140M', color='blue', weight='bold', alpha=0.5)

    plt.yscale('log'); plt.title('G140M Extension: Overlap with Observed G235M'); plt.legend(loc='upper right'); plt.grid(alpha=0.1)
    plt.xlim(0.8, 3.8); plt.ylabel('Flux [Jy]')

    # G235M Subplot
    plt.subplot(2,1,2)
    w_e, f_e = load_spec(files['G235M_ext'])
    w_n, f_n = load_spec(files['G235M_nom'])
    w_t, f_t = load_spec(files['G395M_nom']) # Primary truth
    
    # Scaling for G235M
    w_start = np.max([np.nanmin(w_e), np.nanmin(w_t)])
    window = 0.1
    mask_e = (w_e >= w_start) & (w_e <= w_start + window)
    mask_t = (w_t >= w_start) & (w_t <= w_start + window)
    
    if np.any(mask_e) and np.any(mask_t):
        scale = np.nanmedian(f_t[mask_t]) / np.nanmedian(f_e[mask_e])
    else:
        scale = np.nanmedian(f_t) / np.nanmedian(f_e)
    f_e *= scale
    
    interp_e = interp1d(w_e, f_e, bounds_error=False, fill_value=0)
    interp_n = interp1d(w_n, f_n, bounds_error=False, fill_value=0)
    f_cal = (interp_e(c235['wav']) - c235['alpha']*interp_n(c235['wav']/2) - c235['beta']*interp_n(c235['wav']/3)) / c235['k']
    
    plt.plot(w_n, f_n, 'k-', lw=1.5, label='Nominal G235M')
    plt.plot(w_t, f_t, 'k-', alpha=0.3, label='Observed G395M (Truth)')
    plt.plot(w_e, f_e, 'g--', alpha=0.4, label='Extended G235M (Raw)')
    plt.plot(c235['wav'], f_cal, 'r-', lw=2, label='Extended G235M (Recalibrated)')
    
    plt.annotate('', xy=(5.3, np.nanmedian(f_n)), xytext=(3.1, np.nanmedian(f_n)), arrowprops=dict(arrowstyle='<->', color='#B8860B', lw=2, alpha=0.3))
    plt.text(4.2, np.nanmedian(f_n)*1.2, 'Extended G235M', color='#B8860B', weight='bold', alpha=0.5)

    plt.yscale('log'); plt.title('G235M Extension: Overlap with Observed G395M'); plt.legend(loc='upper right'); plt.grid(alpha=0.1)
    plt.xlim(1.5, 5.5); plt.xlabel('Wavelength [µm]'); plt.ylabel('Flux [Jy]')

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(f'{OUTPUT_DIR}/SUMMARY_V3_{target_name}.png')
    plt.close()

make_summary_v3('P330E')
make_summary_v3('G191-B2B')
