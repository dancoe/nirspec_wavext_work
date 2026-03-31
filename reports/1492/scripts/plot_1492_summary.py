import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import interp1d

# --- Paths ---
BASE = '/Users/dcoe/NIRSpec/wavext'
DATA_DIR = f'{BASE}/data/PID1492'
RESULTS_DIR = f'{BASE}/results/v7'
PLOT_DIR = f'{BASE}/nirspec_wavext_work/reports/1492/plots'
os.makedirs(PLOT_DIR, exist_ok=True)

C_ANG_S = 2.99792458e18

# For PID 1492, Truth comes from adjacent grating overlap (stitching)
PID1492_MAP = {
    'g140m': f'{DATA_DIR}/spec3_nrs1_nom/jw01492-o001_t001-s000000001_nirspec_f170lp-g235m-s200a1-subs200a1_x1d.fits',
    'g235m': f'{DATA_DIR}/spec3_nrs1_nom/jw01492-o001_t001-s000000001_nirspec_f290lp-g395m-s200a1-subs200a1_x1d.fits',
}

def load_spec(path):
    if not os.path.exists(path): return None, None
    with fits.open(path) as h:
        d = h[1].data
        wl, fl = d['WAVELENGTH'].astype(float), d['FLUX'].astype(float)
        dq = d['DQ'] if 'DQ' in d.dtype.names else np.zeros_like(wl)
    mask = (wl > 0) & np.isfinite(fl) & (dq == 0)
    return wl[mask], fl[mask]

def plot_1492_validation(grating):
    # Load v7 coefficients
    cal_file = os.path.join(RESULTS_DIR, f'calib_v7_fs_{grating}_all.fits')
    if not os.path.exists(cal_file):
        print(f"v7 calibration file not found: {cal_file}")
        return
    with fits.open(cal_file) as h:
        c = h[1].data
        wl_grid, k, alpha, beta = c['WAVELENGTH'], c['K'], c['ALPHA'], c['BETA']
    
    # Load Observed NRS2 Jy data
    obs_path = f'{DATA_DIR}/spec3_nrs1_nom/*_x1d.fits' # This might be named differently
    # Let's find any nrs2 file in the data dir
    hits = glob.glob(os.path.join(DATA_DIR, '**', '*_x1d.fits'), recursive=True)
    nrs2_obs = [h for h in hits if 'nrs2' in h.lower() and grating.lower() in h.lower()]
    
    if not nrs2_obs:
        print(f"No NRS2 observed files found for {grating}")
        return

    wl_o, fl_o = load_spec(nrs2_obs[0])
    if wl_o is None: return

    # Load "Truth" (adjacent grating)
    truth_path = PID1492_MAP.get(grating.lower())
    wl_t, fl_t = load_spec(truth_path)
    if wl_t is None:
        print(f"Truth file not found for {grating}: {truth_path}")
        return

    f_t_interp = interp1d(wl_t, fl_t, bounds_error=False, fill_value=0)
    
    # Model: S_model = k*f1 + a*f2 + b*f3
    k_i = np.interp(wl_o, wl_grid, k, left=0, right=0)
    a_i = np.interp(wl_o, wl_grid, alpha, left=0, right=0)
    b_i = np.interp(wl_o, wl_grid, beta, left=0, right=0)
    
    model = k_i * f_t_interp(wl_o) + a_i * f_t_interp(wl_o/2) + b_i * f_t_interp(wl_o/3)
    
    plt.figure(figsize=(10,6))
    plt.plot(wl_o, fl_o, 'k-', alpha=0.3, label=f'Observed (NRS2 Jy - v7)')
    plt.plot(wl_o, model, 'r--', label='v7 Model (k*f1 + α*f2 + β*f3)')
    plt.plot(wl_o, k_i * f_t_interp(wl_o), 'b-', alpha=0.5, label='Predicted 1st Order (k*f1)')
    
    # Truth itself for context
    # plt.plot(wl_t, fl_t, 'g:', label='Truth (Adjacent Grating)')

    plt.title(f'PID 1492 (IRAS-05248) - {grating.upper()} NRS2 Validation (v7)')
    plt.xlabel('Wavelength (µm)')
    plt.ylabel('Flux (Jy)')
    plt.legend()
    plt.grid(alpha=0.2)
    plt.savefig(os.path.join(PLOT_DIR, f'1492_nrs2_validation_{grating}.png'), dpi=150)
    plt.close()

if __name__ == '__main__':
    for g in ['g140m', 'g235m']:
        plot_1492_validation(g)
