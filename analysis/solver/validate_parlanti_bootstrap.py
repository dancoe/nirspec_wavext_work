import os, sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import interp1d
import pandas as pd

# Paths
DATA_DIR = '/Users/dcoe/NIRSpec/wavext/data'
COEFF_DIR = '/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/plots/Parlanti/cal/153678'
OUTPUT_DIR = '/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/plots/Parlanti/cal/153678'
os.makedirs(OUTPUT_DIR, exist_ok=True)

def load_spec(path):
    if not os.path.exists(path):
        return np.array([]), np.array([])
    with fits.open(path) as hdul:
        data = hdul[1].data
        w = data['wavelength']
        f = data['flux']
        msk = (w > 0) & np.isfinite(f) & (f > 0)
        return w[msk], f[msk]

# Load coefficients
c140 = pd.read_csv(f'{COEFF_DIR}/coeffs_G140M.csv')
c235 = pd.read_csv(f'{COEFF_DIR}/coeffs_G235M.csv')

interp_k140 = interp1d(c140['wav'], c140['k'], bounds_error=False, fill_value=0)
interp_a140 = interp1d(c140['wav'], c140['alpha'], bounds_error=False, fill_value=0)
interp_b140 = interp1d(c140['wav'], c140['beta'], bounds_error=False, fill_value=0)

interp_k235 = interp1d(c235['wav'], c235['k'], bounds_error=False, fill_value=0)
interp_a235 = interp1d(c235['wav'], c235['alpha'], bounds_error=False, fill_value=0)
interp_b235 = interp1d(c235['wav'], c235['beta'], bounds_error=False, fill_value=0)

targets = ['G191-B2B', 'P330E']
pids = {'G191-B2B': 'PID1537_G191-B2B', 'P330E': 'PID1538_P330E'}

# File mappings (manually corrected)
specs = {
    'G191-B2B': {
        'G140M_ext': f'{DATA_DIR}/PID1537_G191-B2B/jw01537007001_07101_00003_nrs2_extract_1d.fits',
        'G140M_nom': f'{DATA_DIR}/PID1537_G191-B2B/jw01537-o007_t001-s000000001_nirspec_f100lp-g140m-s1600a1-sub2048_x1d.fits',
        'G235M_nom': f'{DATA_DIR}/PID1537_G191-B2B/jw01537-o007_t001-s000000001_nirspec_f170lp-g235m-s1600a1-sub2048_x1d.fits',
        'G395M_nom': f'{DATA_DIR}/PID1537_G191-B2B/jw01537-o007_t001-s000000001_nirspec_f290lp-g395m-s1600a1-sub2048_x1d.fits'
    },
    'P330E': {
        'G140M_ext': f'{DATA_DIR}/PID1538_P330E/jw01538160001_06101_00001_nrs2_extract_1d.fits',
        'G140M_nom': f'{DATA_DIR}/PID1538_P330E/jw01538-o160_t002-s000000001_nirspec_f100lp-g140m-s1600a1-sub2048_x1d.fits',
        'G235M_ext': f'{DATA_DIR}/PID1538_P330E/jw01538160001_08101_00005_nrs2_extract_1d.fits',
        'G235M_nom': f'{DATA_DIR}/PID1538_P330E/jw01538-o160_t002-s000000001_nirspec_f170lp-g235m-s1600a1-sub2048_x1d.fits',
        'G395M_nom': f'{DATA_DIR}/PID1538_P330E/jw01538-o160_t002-s000000001_nirspec_f290lp-g395m-s1600a1-sub2048_x1d.fits'
    }
}

for tgt in targets:
    if 'G140M_ext' not in specs[tgt]: continue
    
    # 1. Calibrate G140M
    w_ext, f_ext = load_spec(specs[tgt]['G140M_ext'])
    w_nom, f_nom = load_spec(specs[tgt]['G140M_nom'])
    w_truth, f_truth = load_spec(specs[tgt]['G235M_nom'])
    
    if len(w_ext) == 0 or len(w_truth) == 0: continue
    
    # Scale DN/s to Jy
    scale_wav = 2.0
    mask_truth = (w_truth > scale_wav-0.05) & (w_truth < scale_wav+0.05)
    mask_ext = (w_ext > scale_wav-0.05) & (w_ext < scale_wav+0.05)
    scale = np.median(f_truth[mask_truth]) / np.median(f_ext[mask_ext])
    f_ext_jy = f_ext * scale
    
    # Apply coefficients
    # S_obs = k f_truth + a f_clean(L/2) + b f_clean(L/3)
    # => f_cal = (S_obs - a f_clean(L/2) - b f_clean(L/3)) / k
    interp_clean = interp1d(w_nom, f_nom, bounds_error=False, fill_value=0)
    
    w_vec = w_ext
    k = interp_k140(w_vec)
    a = interp_a140(w_vec)
    b = interp_b140(w_vec)
    
    f_cal = (f_ext_jy - a * interp_clean(w_vec/2) - b * interp_clean(w_vec/3)) / np.maximum(k, 0.05)
    
    plt.figure(figsize=(10,6))
    plt.plot(w_truth, f_truth, 'k-', alpha=0.5, label='G235M (Truth)')
    plt.plot(w_nom, f_nom, 'b-', label='G140M (Nominal)')
    plt.plot(w_ext, f_ext_jy, 'g--', alpha=0.3, label='G140M (Extended Raw Scaled)')
    plt.plot(w_ext, f_cal, 'r-', label='G140M (Calibrated Bootstrap)')
    plt.xlim(1.5, 3.5); plt.ylim(0, np.percentile(f_truth, 99) * 1.5)
    plt.title(f'{tgt} - G140M Calibration Validation')
    plt.legend(); plt.grid()
    plt.savefig(f'{OUTPUT_DIR}/{tgt}_G140M_cal_val.png')
    plt.close()

print(f"Validation plots saved to {OUTPUT_DIR}")
