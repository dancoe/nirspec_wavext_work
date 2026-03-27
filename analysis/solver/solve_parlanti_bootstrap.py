import os, sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import interp1d
from scipy.optimize import nnls
import pandas as pd

# Paths
DATA_DIR = '/Users/dcoe/NIRSpec/wavext/data'
OUTPUT_DIR = '/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/plots/Parlanti/cal/153678_v2'
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Formatting from PARLANTI_PLOTS.md
# G140M: Blue / G235M: Dark Goldenrod
COLORS = {
    'G140M': {'main': 'blue', '2nd': 'cornflowerblue', '3rd': 'lightblue'},
    'G235M': {'main': '#B8860B', '2nd': '#DAA520', '3rd': '#FFD700'} # Dark Goldenrod, Goldenrod, Gold
}

def load_spec(path):
    if not os.path.exists(path):
        return np.array([]), np.array([])
    try:
        with fits.open(path) as hdul:
            data = hdul[1].data
            w = data['wavelength']
            f = data['flux']
            msk = (w > 0) & np.isfinite(f) & (f > 0)
            return w[msk], f[msk]
    except:
        return np.array([]), np.array([])

def solve_bootstrap_v4(grating, wav_grid):
    print(f"\n--- Solving for {grating} (Bootstrap V4 - Priority 1st Order) ---")
    
    # Sources
    samples = {
        'G191-B2B': {
            'ext': f'{DATA_DIR}/PID1537_G191-B2B/jw01537007001_{"07101_00003" if grating == "G140M" else "09101_00003"}_nrs2_extract_1d.fits',
            'clean_nom': f'{DATA_DIR}/PID1537_G191-B2B/jw01537-o007_t001-s000000001_nirspec_f100lp-g140m-s1600a1-sub2048_x1d.fits' if grating == 'G140M' else f'{DATA_DIR}/PID1537_G191-B2B/jw01537-o007_t001-s000000001_nirspec_f170lp-g235m-s1600a1-sub2048_x1d.fits',
            'truth': f'{DATA_DIR}/PID1537_G191-B2B/jw01537-o007_t001-s000000001_nirspec_f170lp-g235m-s1600a1-sub2048_x1d.fits' if grating == 'G140M' else f'{DATA_DIR}/PID1537_G191-B2B/jw01537-o007_t001-s000000001_nirspec_f290lp-g395m-s1600a1-sub2048_x1d.fits'
        },
        'P330E': {
            'ext': f'{DATA_DIR}/PID1538_P330E/jw01538160001_{"06101_00001" if grating == "G140M" else "08101_00005"}_nrs2_extract_1d.fits',
            'clean_nom': f'{DATA_DIR}/PID1538_P330E/jw01538-o160_t002-s000000001_nirspec_f100lp-g140m-s1600a1-sub2048_x1d.fits' if grating == 'G140M' else f'{DATA_DIR}/PID1538_P330E/jw01538-o160_t002-s000000001_nirspec_f170lp-g235m-s1600a1-sub2048_x1d.fits',
            'truth': f'{DATA_DIR}/PID1538_P330E/jw01538-o160_t002-s000000001_nirspec_f170lp-g235m-s1600a1-sub2048_x1d.fits' if grating == 'G140M' else f'{DATA_DIR}/PID1538_P330E/jw01538-o160_t002-s000000001_nirspec_f290lp-g395m-s1600a1-sub2048_x1d.fits'
        },
        'IRAS-05248': {
            'ext': f'{DATA_DIR}/PID1492/jw01492003001_{"03102_00005_nrs2" if grating == "G140M" else "03104_00004_nrs2_g235m"}_extract_1d.fits',
            'clean_nom': f'{DATA_DIR}/PID1492/jw01492-o001_t001-s000000001_nirspec_f100lp-g140m-s200a1-subs200a1_x1d.fits' if grating == 'G140M' else f'{DATA_DIR}/PID1492/jw01492-o001_t001-s000000001_nirspec_f170lp-g235m-s200a1-subs200a1_x1d.fits',
            'truth': f'{DATA_DIR}/PID1492/jw01492-o001_t001-s000000001_nirspec_f170lp-g235m-s200a1-subs200a1_x1d.fits' if grating == 'G140M' else f'{DATA_DIR}/PID1492/jw01492-o001_t001-s000000001_nirspec_f290lp-g395m-s200a1-subs200a1_x1d.fits'
        }
    }

    all_As, all_Cs = [], []
    valid_sources = []

    for name, paths in samples.items():
        w_ext, f_ext = load_spec(paths['ext'])
        w_clean, f_clean = load_spec(paths['clean_nom']) 
        w_truth, f_truth = load_spec(paths['truth']) 
        
        if len(w_ext) == 0 or len(w_truth) == 0:
            print(f"Skipping {name} due to missing data.")
            continue
            
        valid_sources.append(name)
        
        # Scaling ext (DN/s) to truth (Jy) at crossover
        # For G140M, nominal ends at 1.88. G235M truth starts at 1.7.
        # Let's crossover at 1.85
        scale_wav = 1.85 if grating == 'G140M' else 3.1
        mask_truth = (w_truth > scale_wav - 0.05) & (w_truth < scale_wav + 0.05)
        mask_ext = (w_ext > scale_wav - 0.05) & (w_ext < scale_wav + 0.05)
        
        if len(f_ext[mask_ext]) == 0 or len(f_truth[mask_truth]) == 0:
             # Try a different range if 1.85/3.1 is not populated
             scale_wav = np.max(w_truth[w_truth < np.max(w_ext)]) - 0.1
             mask_truth = (w_truth > scale_wav - 0.05) & (w_truth < scale_wav + 0.05)
             mask_ext = (w_ext > scale_wav - 0.05) & (w_ext < scale_wav + 0.05)
        
        if len(f_ext[mask_ext]) == 0 or len(f_truth[mask_truth]) == 0:
            scale = 1.0; print(f"Source {name} scale fallback")
        else:
            scale = np.nanmedian(f_truth[mask_truth]) / np.nanmedian(f_ext[mask_ext])
        
        f_ext_scaled = f_ext * scale
        interp_truth = interp1d(w_truth, f_truth, bounds_error=False, fill_value=0)
        interp_clean = interp1d(w_clean, f_clean, bounds_error=False, fill_value=0)
        interp_ext_scaled = interp1d(w_ext, f_ext_scaled, bounds_error=False, fill_value=0)
        
        source_A, source_C = [], []
        for L in wav_grid:
            tr_val = interp_truth(L)
            if tr_val > 0:
                source_A.append([1.0, interp_clean(L/2)/tr_val, interp_clean(L/3)/tr_val])
                source_C.append(interp_ext_scaled(L)/tr_val)
            else:
                source_A.append([0, 0, 0])
                source_C.append(0)
        all_As.append(source_A)
        all_Cs.append(source_C)

    ks, alphas, betas = [], [], []
    # Priors
    k_prior_val = 0.8
    k_prior_weight = 1.0 # Stronger weight to keep k near 0.8 during cutoff
    ghost_penalty = 2.0 # Stronger penalty on ghosts
    
    for i in range(len(wav_grid)):
        A_rows = []
        C_vals = []
        for j in range(len(all_As)):
            if all_As[j][i][0] > 0:
                A_rows.append(all_As[j][i])
                C_vals.append(all_Cs[j][i])
        
        if len(A_rows) < 1:
            ks.append(k_prior_val); alphas.append(0); betas.append(0)
            continue
            
        A_mat = np.array(A_rows)
        C_vec = np.array(C_vals)
        
        A_ext = np.vstack([
            A_mat,
            np.array([[k_prior_weight, 0, 0]]),
            np.array([[0, ghost_penalty, 0]]),
            np.array([[0, 0, ghost_penalty]])
        ])
        C_ext = np.hstack([
            C_vec,
            [k_prior_weight * k_prior_val],
            [0], [0]
        ])
        
        try:
            sol, _ = nnls(A_ext, C_ext)
            ks.append(sol[0]); alphas.append(sol[1]); betas.append(sol[2])
        except:
            ks.append(k_prior_val); alphas.append(0); betas.append(0)

    def smooth(y, box_pts):
        box = np.ones(box_pts)/box_pts
        return np.convolve(y, box, mode='same')
    
    ks = smooth(np.array(ks), 15)
    alphas = smooth(np.array(alphas), 15)
    betas = smooth(np.array(betas), 15)
    
    return ks, alphas, betas

grid_g140m = np.linspace(1.8, 3.4, 150)
grid_g235m = np.linspace(3.1, 5.3, 150)

k1, a1, b1 = solve_bootstrap_v4('G140M', grid_g140m)
k2, a2, b2 = solve_bootstrap_v4('G235M', grid_g235m)

# Save
pd.DataFrame({'wav': grid_g140m, 'k': k1, 'alpha': a1, 'beta': b1}).to_csv(f'{OUTPUT_DIR}/coeffs_G140M.csv', index=False)
pd.DataFrame({'wav': grid_g235m, 'k': k2, 'alpha': a2, 'beta': b2}).to_csv(f'{OUTPUT_DIR}/coeffs_G235M.csv', index=False)

# Plotting
plt.figure(figsize=(12,12))

# G140M
plt.subplot(2,1,1)
c140 = COLORS['G140M']
plt.plot(grid_g140m, k1, color=c140['main'], label='k(λ) (1st)', linewidth=2.5)
plt.plot(grid_g140m, a1, color=c140['2nd'], label='α̃(λ) (2nd)', alpha=0.8)
plt.plot(grid_g140m, b1, color=c140['3rd'], label='β̃(λ) (3rd)', alpha=0.8)
plt.yscale('log')
plt.ylim(1e-4, 5)
plt.title('G140M Coefficients (Bootstrap V4 - Priority 1st)')
plt.ylabel('Value'); plt.legend(); plt.grid(True, which="both", ls="-", alpha=0.2)

# G235M
plt.subplot(2,1,2)
c235 = COLORS['G235M']
plt.plot(grid_g235m, k2, color=c235['main'], label='k(λ) (1st)', linewidth=2.5)
plt.plot(grid_g235m, a2, color=c235['2nd'], label='α̃(λ) (2nd)', alpha=0.8)
plt.plot(grid_g235m, b2, color=c235['3rd'], label='β̃(λ) (3rd)', alpha=0.8)
plt.yscale('log')
plt.ylim(1e-4, 5)
plt.title('G235M Coefficients (Bootstrap V4 - Priority 1st)')
plt.xlabel('Wavelength [µm]'); plt.ylabel('Value'); plt.legend(); plt.grid(True, which="both", ls="-", alpha=0.2)

plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/coeffs_bootstrap_v2.png')
print(f"Results saved to {OUTPUT_DIR}/coeffs_bootstrap_v2.png")
