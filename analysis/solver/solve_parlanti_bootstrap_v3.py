import os, sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import interp1d
from scipy.optimize import nnls
import pandas as pd

# Paths
DATA_DIR = '/Users/dcoe/NIRSpec/wavext/data'
OUTPUT_DIR = '/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/plots/Parlanti/cal/153678_v3'
os.makedirs(OUTPUT_DIR, exist_ok=True)

def load_spec(path):
    if not os.path.exists(path): return np.array([]), np.array([])
    try:
        with fits.open(path) as hdul:
            data = hdul[1].data
            w, f = data['wavelength'], data['flux']
            msk = (w > 0) & np.isfinite(f) & (f > 0)
            return w[msk], f[msk]
    except: return np.array([]), np.array([])

def solve_bootstrap_v6(grating, wav_grid):
    print(f"\n--- Solving for {grating} (V6 - Multi-Truth North Star) ---")
    
    samples = ['G191-B2B', 'P330E', 'IRAS-05248']
    all_data = []

    for name in samples:
        if name == 'G191-B2B':
            d = f'{DATA_DIR}/PID1537_G191-B2B'
            ext = f'{d}/jw01537007001_{"07101_00003" if grating == "G140M" else "09101_00003"}_nrs2_extract_1d.fits'
            nom = f'{d}/jw01537-o007_t001-s000000001_nirspec_f100lp-g140m-s1600a1-sub2048_x1d.fits' if grating == 'G140M' else f'{d}/jw01537-o007_t001-s000000001_nirspec_f170lp-g235m-s1600a1-sub2048_x1d.fits'
            t1 = f'{d}/jw01537-o007_t001-s000000001_nirspec_f170lp-g235m-s1600a1-sub2048_x1d.fits'
            t2 = f'{d}/jw01537-o007_t001-s000000001_nirspec_f290lp-g395m-s1600a1-sub2048_x1d.fits'
            t3 = f'{d}/jw01537-o007_t001-s000000001_nirspec_clear-prism-s1600a1-sub2048_x1d.fits'
        elif name == 'P330E':
            d = f'{DATA_DIR}/PID1538_P330E'
            ext = f'{d}/jw01538160001_{"06101_00001" if grating == "G140M" else "08101_00005"}_nrs2_extract_1d.fits'
            nom = f'{d}/jw01538-o160_t002-s000000001_nirspec_f100lp-g140m-s1600a1-sub2048_x1d.fits' if grating == 'G140M' else f'{d}/jw01538-o160_t002-s000000001_nirspec_f170lp-g235m-s1600a1-sub2048_x1d.fits'
            t1 = f'{d}/jw01538-o160_t002-s000000001_nirspec_f170lp-g235m-s1600a1-sub2048_x1d.fits'
            t2 = f'{d}/jw01538-o160_t002-s000000001_nirspec_f290lp-g395m-s1600a1-sub2048_x1d.fits'
            t3 = f'{d}/jw01538-o160_t002-s000000001_nirspec_clear-prism-s1600a1-sub2048_x1d.fits'
        else: # IRAS
            d = f'{DATA_DIR}/PID1492'
            ext = f'{d}/jw01492003001_{"03102_00005_nrs2" if grating == "G140M" else "03104_00004_nrs2_g235m"}_extract_1d.fits'
            nom = f'{d}/jw01492-o001_t001-s000000001_nirspec_f100lp-g140m-s200a1-subs200a1_x1d.fits' if grating == 'G140M' else f'{d}/jw01492-o001_t001-s000000001_nirspec_f170lp-g235m-s200a1-subs200a1_x1d.fits'
            t1 = f'{d}/jw01492-o001_t001-s000000001_nirspec_f170lp-g235m-s200a1-subs200a1_x1d.fits'
            t2 = f'{d}/jw01492-o001_t001-s000000001_nirspec_f290lp-g395m-s200a1-subs200a1_x1d.fits'
            t3 = f'{d}/jw01492-o001_t001-s000000001_nirspec_clear-prism-s200a1-subs200a1_x1d.fits'

        w_e, f_e = load_spec(ext)
        w_n, f_n = load_spec(nom)
        w_t1, f_t1 = load_spec(t1)
        w_t2, f_t2 = load_spec(t2)
        w_t3, f_t3 = load_spec(t3)
        
        if len(w_e) == 0: continue
        
        # Merge truths (t1, t2, t3) into a single master truth for this star
        # Preference: G235M > G395M > PRISM for G140M extension
        # G395M > PRISM for G235M extension
        master_w = np.unique(np.concatenate([w_t1, w_t2, w_t3]))
        master_f = np.zeros_like(master_w)
        
        i1 = interp1d(w_t1, f_t1, bounds_error=False, fill_value=0) if len(w_t1) > 0 else (lambda x: 0)
        i2 = interp1d(w_t2, f_t2, bounds_error=False, fill_value=0) if len(w_t2) > 0 else (lambda x: 0)
        i3 = interp1d(w_t3, f_t3, bounds_error=False, fill_value=0) if len(w_t3) > 0 else (lambda x: 0)
        
        for i, w_val in enumerate(master_w):
            f1, f2, f3 = (i1(w_val) if len(w_t1)>0 else 0), (i2(w_val) if len(w_t2)>0 else 0), (i3(w_val) if len(w_t3)>0 else 0)
            if f1 > 0: master_f[i] = f1
            elif f2 > 0: master_f[i] = f2
            else: master_f[i] = f3
            
        interp_master_truth = interp1d(master_w, master_f, bounds_error=False, fill_value=0)
        interp_nom = interp1d(w_n, f_n, bounds_error=False, fill_value=0)
        
        # Scaling at boundary (MUST match nominal range at crossover)
        # Handle gaps by finding the first available overlap window
        w_all = np.concatenate([w_e, master_w])
        w_start = np.max([np.nanmin(w_e), np.nanmin(master_w)])
        window = 0.05
        mask_e = (w_e >= w_start) & (w_e < w_start + window)
        mask_t = (master_w >= w_start) & (master_w < w_start + window)
        
        if np.any(mask_e) and np.any(mask_t):
            scale = np.nanmedian(master_f[mask_t]) / np.nanmedian(f_e[mask_e])
        else:
            scale = np.nanmedian(master_f) / np.nanmedian(f_e)
        
        all_data.append({
            'ext': interp1d(w_e, f_e * scale, bounds_error=False, fill_value=0),
            'nom': interp_nom,
            'truth': interp_master_truth
        })

    ks, alphas, betas = [], [], []
    for L in wav_grid:
        A, C = [], []
        for d in all_data:
            f_t = d['truth'](L)
            if f_t > 1e-12:
                r2 = d['nom'](L/2) / f_t
                r3 = d['nom'](L/3) / f_t
                S_rel = d['ext'](L) / f_t
                if np.isfinite([r2, r3, S_rel]).all():
                    A.append([1.0, r2, r3])
                    C.append(S_rel)
        
        if len(A) == 0:
            ks.append(1.0); alphas.append(0); betas.append(0)
            continue
            
        # Regularized solve - Strength of priors
        k_wt = 5.0    # Anchor k=1 strongly
        g_wt = 0.01   # Allow ghosts to grow if data forces it
        
        A_ext = np.vstack([np.array(A), [[k_wt, 0, 0]], [[0, g_wt, 0]], [[0, 0, g_wt]]])
        C_ext = np.hstack([np.array(C), [k_wt * 1.0], [0], [0]])
        sol, _ = nnls(A_ext, C_ext)
        ks.append(sol[0]); alphas.append(sol[1]); betas.append(sol[2])

    def smooth(y, box): return np.convolve(y, np.ones(box)/box, mode='same')
    ks_s = smooth(ks, 40) # Parlanti uses 40
    # Ensure k doesn't hit zero at ends
    ks_s[:20] = ks[:20]; ks_s[-20:] = ks[-20:]
    return ks_s, smooth(alphas, 40), smooth(betas, 40)

grid140 = np.linspace(1.8, 3.4, 200)
grid235 = np.linspace(3.0, 5.3, 200)
k1, a1, b1 = solve_bootstrap_v6('G140M', grid140)
k2, a2, b2 = solve_bootstrap_v6('G235M', grid235)

pd.DataFrame({'wav': grid140, 'k': k1, 'alpha': a1, 'beta': b1}).to_csv(f'{OUTPUT_DIR}/coeffs_G140M.csv', index=False)
pd.DataFrame({'wav': grid235, 'k': k2, 'alpha': a2, 'beta': b2}).to_csv(f'{OUTPUT_DIR}/coeffs_G235M.csv', index=False)
