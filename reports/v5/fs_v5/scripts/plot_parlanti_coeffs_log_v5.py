import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits

# Paths
BASE        = '/Users/dcoe/NIRSpec/wavext'
COEFFS_DIR  = f'{BASE}/results/v5'
OUTPUT_DIR  = f'{BASE}/nirspec_wavext_work/reports/v5/fs_v5/plots'
os.makedirs(OUTPUT_DIR, exist_ok=True)

CONFIGS = [
    ('G140M/F100LP', 'g140m_f100lp'),
    ('G235M/F170LP', 'g235m_f170lp'),
]

# Style matches v2/Parlanti-style log plots
colors = {'k': 'black', 'alpha': 'darkorange', 'beta': 'royalblue'}

fig, axes = plt.subplots(1, 2, figsize=(16, 7))
fig.suptitle('NIRSpec FS Calibration v5 — Joint 4-Source Solve (Log Scale)\n'
             'k(λ) ~ 0.96 | α(λ) ~ 0.003–0.012 (Degeneracy Broken)', 
             fontsize=16, fontweight='bold')

for i, (label, key) in enumerate(CONFIGS):
    ax = axes[i]
    fits_path = f'{COEFFS_DIR}/calib_v5_{key}.fits'
    
    if not os.path.exists(fits_path):
        ax.text(0.5, 0.5, f'Missing {key}', transform=ax.transAxes, ha='center')
        continue
        
    with fits.open(fits_path) as h:
        d = h[1].data
        wl = d['WAVELENGTH']
        k  = d['K']
        a  = d['ALPHA']
        b  = d['BETA']
        
        # Log floor for zero/negative values (though v5 NNLS should be positive)
        a_plot = np.where(a > 1e-6, a, 1e-6)
        b_plot = np.where(b > 1e-6, b, 1e-6)
        
        ax.plot(wl, k,      color=colors['k'],     lw=2.5, label=r'$k(\lambda)$ [v5]')
        ax.plot(wl, a_plot, color=colors['alpha'], lw=2.0, label=r'$\tilde{\alpha}(\lambda)$ [v5]')
        ax.plot(wl, b_plot, color=colors['beta'],  lw=2.0, label=r'$\tilde{\beta}(\lambda)$ [v5]')
        
        ax.axhline(1.0, color='red', ls='--', alpha=0.3, label='k=1.0')
        
        # Parlanti Comparison (if reference available)
        # For now, just the v5 solution
        
        ax.set_yscale('log')
        ax.set_ylim(1e-4, 5.0)
        ax.set_xlabel('Wavelength [µm]', fontsize=13)
        ax.set_ylabel('Coefficient Value', fontsize=13)
        ax.set_title(label, fontsize=14, fontweight='bold')
        ax.grid(True, which='both', alpha=0.2)
        ax.legend(loc='lower left', fontsize=10)

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
outpath = f'{OUTPUT_DIR}/fs_v5_coeffs_log_2panel.png'
plt.savefig(outpath, dpi=200, bbox_inches='tight')
plt.close()
print(f'Saved: {outpath}')
