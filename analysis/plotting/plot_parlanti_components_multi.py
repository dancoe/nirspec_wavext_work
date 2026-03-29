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
plt.style.use('dark_background')
fig, ax = plt.subplots(figsize=(12, 8))

# 1. Plot Reference and Component Potentials
ax.plot(ref_w, ref_f, color='white', alpha=0.5, lw=1.5, label=f'Baseline f(λ) [{TARGET}]', zorder=1)

# Shifted References (Potential Contaminants)
grid_full = np.linspace(0.6, 5.6, 1000)
ax.plot(grid_full, f_interp(grid_full/2), color='#2ecc71', alpha=0.3, lw=1, label='Reference f(λ/2) [Ghost Source 2nd]', linestyle='--')
ax.plot(grid_full, f_interp(grid_full/3), color='#e056fd', alpha=0.3, lw=1, label='Reference f(λ/3) [Ghost Source 3rd]', linestyle='--')

# 2. Plot Extended Obs and Calibrated Components
for grating, coeff_data, color, ext_key in [
    ('G140M', c140, '#3498db', 'G140M_EXT'),
    ('G235M', c235, '#f1c40f', 'G235M_EXT')
]:
    if coeff_data is None: continue
    w_grid = coeff_data['wav']
    k = coeff_data['k']
    alpha = coeff_data['alpha']
    beta = coeff_data['beta']
    
    f_L = f_interp(w_grid)
    f_L2 = f_interp(w_grid/2)
    f_L3 = f_interp(w_grid/3)
    
    if ext_key in specs:
        w_obs, f_obs = specs[ext_key]
        scale_msk = (w_obs > w_grid[0]) & (w_obs < w_grid[0] + 0.1)
        scale = np.nanmedian(f_interp(w_obs[scale_msk])) / np.nanmedian(f_obs[scale_msk])
        f_obs_scaled = f_obs * scale
        
        ax.plot(w_obs, f_obs_scaled, color=color, alpha=0.2, lw=0.8, label=f'{grating} Raw Extended')
        
        comp_1st = k * f_L
        comp_2nd = alpha * f_L2
        comp_3rd = beta * f_L3
        comp_sum = comp_1st + comp_2nd + comp_3rd
        
        ax.plot(w_grid, comp_2nd, color='#2ecc71', lw=1.5, label=f'α(λ)·f(λ/2) [{grating}]' if grating == 'G140M' else "")
        ax.plot(w_grid, comp_3rd, color='#e056fd', lw=1.5, label=f'β(λ)·f(λ/3) [{grating}]' if grating == 'G140M' else "")
        
        ax.plot(w_grid, comp_sum, color='red', lw=1, label=f'Model Sum (k f + α f/2 + β f/3)' if grating == 'G140M' else "")

# Formatting
ax.set_yscale('log')
ax.set_ylim(1e-3, 5)
ax.set_xlim(0.6, 5.6)
ax.set_xlabel('Wavelength [µm]', fontsize=12)
ax.set_ylabel('Flux [Jy]', fontsize=12)
ax.set_title(f'NIRSpec Wavelength Extension: Ghost Component Decomposition ({TARGET})', fontsize=14, fontweight='bold', pad=20)
ax.legend(loc='lower right', fontsize=9, ncol=2)
ax.grid(alpha=0.1)

ax.annotate('Contamination comes from light at λ/2 and λ/3', 
            xy=(4.0, 1e-1), xytext=(4.2, 5e-1),
            arrowprops=dict(facecolor='white', shrink=0.05, width=1, headwidth=5),
            fontsize=10, color='white')

plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/CAL_COMPONENTS_{TARGET}.png', dpi=300)
print(f"Plot saved to {OUTPUT_DIR}/CAL_COMPONENTS_{TARGET}.png")
