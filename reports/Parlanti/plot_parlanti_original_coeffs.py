import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

# Paths
BASE        = '/Users/dcoe/NIRSpec/wavext'
PARLANTI_DIR = f'{BASE}/data/parlanti_repo/calibration_files'
OUTPUT_DIR  = f'{BASE}/nirspec_wavext_work/reports/Parlanti'
os.makedirs(OUTPUT_DIR, exist_ok=True)

def plot_parlanti_coeffs():
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
    fig.suptitle('Original Parlanti et al. (2025) Calibration — k/α/β coefficients', fontsize=14)

    colors = {'k': 'violet', 'alpha': 'darkorange', 'beta': '#3498db'}

    for ax, grating, fname in zip(axes, ['G140M', 'G235M'],
                                  ['calibration_functions_g140m_f100lp.fits',
                                   'calibration_functions_g235m_f170lp.fits']):
        path = f'{PARLANTI_DIR}/{fname}'
        if not os.path.exists(path):
            print(f'Warning: {path} not found')
            continue
            
        with fits.open(path) as h:
            data = h[1].data
            wl = data['wavelength']
            k = data['k']
            a = data['alpha']
            b = data['beta']

        ax.plot(wl, k, color=colors['k'], lw=2, label=r'$k(\lambda)$')
        ax.plot(wl, a, color=colors['alpha'], lw=1.5, label=r'$\tilde{\alpha}(\lambda)$')
        ax.plot(wl, b, color=colors['beta'], lw=1.5, label=r'$\tilde{\beta}(\lambda)$')

        ax.set_yscale('log')
        ax.set_ylim(1e-4, 2.0)
        ax.set_xlabel('Wavelength [µm]', fontsize=12)
        if ax == axes[0]:
            ax.set_ylabel('Value', fontsize=12)
        
        ax.set_title(f'{grating}', fontsize=12)
        ax.grid(True, which='both', alpha=0.2)
        ax.legend(fontsize=10)

    plt.tight_layout()
    outpath = f'{OUTPUT_DIR}/parlanti_original_coeffs_log.png'
    plt.savefig(outpath, dpi=200, bbox_inches='tight')
    plt.close()
    print(f'Saved: {outpath}')

if __name__ == '__main__':
    plot_parlanti_coeffs()
