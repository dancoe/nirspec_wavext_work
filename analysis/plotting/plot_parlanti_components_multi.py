import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import interp1d

# Paths
DATA_DIR = '/Users/dcoe/NIRSpec/wavext/data'
COEFF_DIR = '/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/plots/Parlanti/cal/153678_v3'
OUTPUT_DIR = '/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/plots/Parlanti/cal/153678_v3'
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Target: P330E
TARGET = 'P330E'
TARGET_DIR = f'{DATA_DIR}/PID1538_P330E'

# File Mapping
FILES = {
    'PRISM': f'{TARGET_DIR}/jw01538-o107_t002-s000000001_nirspec_clear-prism-s1600a1-sub512_x1d.fits',
    'G140M_NOM': f'{TARGET_DIR}/jw01538-o160_t002-s000000001_nirspec_f100lp-g140m-s1600a1-sub2048_x1d.fits',
    'G235M_NOM': f'{TARGET_DIR}/jw01538-o160_t002-s000000001_nirspec_f170lp-g235m-s1600a1-sub2048_x1d.fits',
    'G395M_NOM': f'{TARGET_DIR}/jw01538-o160_t002-s000000001_nirspec_f290lp-g395m-s1600a1-sub2048_x1d.fits',
    'G140M_EXT': f'{TARGET_DIR}/jw01538160001_06101_00001_nrs2_extract_1d.fits',
    'G235M_EXT': f'{TARGET_DIR}/jw01538160001_08101_00005_nrs2_extract_1d.fits',
}

def load_spec(path):
    if not os.path.exists(path):
        print(f"Warning: Missing {path}")
        return np.array([]), np.array([])
    with fits.open(path) as hdul:
        data = hdul[1].data
        w, f = data['wavelength'], data['flux']
        msk = (w > 0) & np.isfinite(f) & (f > 0)
        return w[msk], f[msk]

# Load All
specs = {}
for key, path in FILES.items():
    w, f = load_spec(path)
    if len(w) > 0:
        specs[key] = (w, f)

# Interpolate Reference f(lambda) from PRISM (best global reference for P330E)
if 'PRISM' in specs:
    ref_w, ref_f = specs['PRISM']
    f_interp = interp1d(ref_w, ref_f, bounds_error=False, fill_value=0)
else:
    print("Error: PRISM reference data missing.")
    sys.exit(1)

# Load Coefficients
def load_coeffs(path):
    if not os.path.exists(path): return None
    # Using numpy to load CSV (skipping header)
    data = np.loadtxt(path, delimiter=',', skiprows=1)
    return {'wav': data[:,0], 'k': data[:,1], 'alpha': data[:,2], 'beta': data[:,3]}

c140 = load_coeffs(f'{COEFF_DIR}/coeffs_G140M.csv')
c235 = load_coeffs(f'{COEFF_DIR}/coeffs_G235M.csv')

# --- Plotting ---
plt.style.use('default') # White background
fig, ax = plt.subplots(figsize=(12, 8))

# 1. Plot Reference
ax.plot(ref_w, ref_f, color='black', alpha=0.3, lw=1.5, label=f'Intrinsic Reference f(λ) [{TARGET}]', zorder=1)

# Shifted References (Visualizing where contamination ORIGINATES)
grid_full = np.linspace(0.6, 5.6, 2000)
ax.plot(grid_full, f_interp(grid_full/2), color='#2ecc71', alpha=0.4, lw=1.2, label='Potential Ghost Source f(λ/2)', linestyle='--')
ax.plot(grid_full, f_interp(grid_full/3), color='#e056fd', alpha=0.4, lw=1.2, label='Potential Ghost Source f(λ/3)', linestyle='--')

# 2. Plot Extended Obs and Calibrated Components
for grating, coeff_data, color, nom_key, ext_key in [
    ('G140M', c140, '#3498db', 'G140M_NOM', 'G140M_EXT'),
    ('G235M', c235, '#f1c40f', 'G235M_NOM', 'G235M_EXT')
]:
    # a. Plot Nominal Range
    if nom_key in specs:
        w_n, f_n = specs[nom_key]
        msk_n = (w_n > 1.2) & (w_n < 1.7) if grating == 'G140M' else (w_n > 2.0) & (w_n < 3.0)
        sc_n = np.nanmedian(f_interp(w_n[msk_n])) / np.nanmedian(f_n[msk_n])
        ax.plot(w_n, f_n * sc_n, color=color, lw=1.5, alpha=0.7, label=f'{grating} Nominal')

    # b. Plot Extended Range
    if ext_key in specs and coeff_data is not None:
        w_obs, f_obs = specs[ext_key]
        w_grid, k, alpha, beta = coeff_data['wav'], coeff_data['k'], coeff_data['alpha'], coeff_data['beta']
        
        # Scaling: Use the NOMINAL range as the anchor to preserve ghosts in the extension
        if nom_key in specs:
            w_n, f_n = specs[nom_key]
            # Use overlap between extended and nominal (clean region)
            msk_ovl = (w_obs > np.nanmin(w_n) + 0.1) & (w_obs < np.nanmax(w_n) - 0.1)
            f_nom_interp = interp1d(w_n, f_n, bounds_error=False)(w_obs[msk_ovl])
            scale = np.nanmedian(f_nom_interp) / np.nanmedian(f_obs[msk_ovl])
        else:
            # Fallback to scaling at the boundary if nominal is missing
            msk_ovl = (w_obs > w_grid[0]) & (w_obs < w_grid[0] + 0.1)
            scale = np.nanmedian(f_interp(w_obs[msk_ovl])) / np.nanmedian(f_obs[msk_ovl])
            
        S_obs = f_obs * scale
        
        # 1. Plot Unscaled Raw Data (NOW prominent at high flux units)
        ax.plot(w_obs, f_obs, color=color, alpha=0.5, lw=0.6, label=f'{grating} Raw (Unscaled)')
        
        # 2. Plot Scaled Observation (Anchored to Physics)
        ax.plot(w_obs, S_obs, color=color, alpha=0.7, lw=0.8, label=f'{grating} Calibrated (Jy)')
        
        # Calculate Components
        f_L, f_L2, f_L3 = f_interp(w_grid), f_interp(w_grid/2), f_interp(w_grid/3)
        comp_2nd = alpha * f_L2
        comp_3rd = beta * f_L3
        model_sum = (k * f_L) + comp_2nd + comp_3rd
        
        # Model Overlay
        ax.plot(w_grid, model_sum, color='red', ls='--', lw=1.5, label=f'Model Sum [{grating}]' if grating == 'G140M' else "")
        
        # Corrected Data
        f_obs_interp = interp1d(w_obs, S_obs, bounds_error=False, fill_value=np.nan)(w_grid)
        corrected = f_obs_interp - comp_2nd - comp_3rd
        ax.plot(w_grid, corrected, color='black', ls=(0, (1, 1)), lw=1.8, alpha=0.9, label=f'Corrected Data [{grating}]' if grating == 'G140M' else "")

# Formatting
ax.set_yscale('log')
ax.set_ylim(1e-3, 10000) # Catch artifacts/spikes without excessive scale 1e6
ax.set_xlim(0.6, 5.6)
ax.set_xlabel('Wavelength [µm]', fontsize=14)
ax.set_ylabel('Flux Units (Jy for Calibrated, ADU/s? for Raw)', fontsize=12)
ax.set_title(f'Spectral Decomposition: Ghost Contribution & Unit Check ({TARGET})', fontsize=16, fontweight='bold', pad=15)
ax.legend(loc='lower left', fontsize=8, ncol=2, frameon=True, framealpha=1.0)
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/CAL_COMPONENTS_{TARGET}.png', dpi=300)
print(f"Plot saved to {OUTPUT_DIR}/CAL_COMPONENTS_{TARGET}.png")
