"""
Compare v7 and v6 Parlanti coefficients for NRS2 extended regions.
"""
import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

BASE = '/Users/dcoe/NIRSpec/wavext'
V6_DIR = f'{BASE}/results/v6'
V7_DIR = f'{BASE}/results/v7'
PLOT_DIR = f'{BASE}/nirspec_wavext_work/reports/v7/fs_v7/plots'
os.makedirs(PLOT_DIR, exist_ok=True)

def compare_coeffs(grating, filt):
    key6 = f'{grating}_{filt}'
    path6 = os.path.join(V6_DIR, f'calib_v6_fs_{key6}.fits')
    path7 = os.path.join(V7_DIR, f'calib_v7_fs_{grating}_all.fits')
    
    if not (os.path.exists(path6) and os.path.exists(path7)):
        print(f'Missing {grating}/{filt} results for comparison')
        return

    with fits.open(path6) as h6, fits.open(path7) as h7:
        d6, d7 = h6[1].data, h7[1].data
        
    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    
    # K (throughput)
    axes[0].plot(d6['WAVELENGTH'], d6['K'], 'k--', label='v6 (DN/s scaled)')
    axes[0].plot(d7['WAVELENGTH'], d7['K'], 'r-',  label='v7 (NRS2 Jy)')
    axes[0].set_ylabel('k (Throughput)')
    axes[0].legend()
    axes[0].set_title(f'v7 vs v6 Comparison: FS {grating.upper()}/{filt.upper()}')
    
    # Alpha (2nd order)
    axes[1].plot(d6['WAVELENGTH'], d6['ALPHA'], 'k--')
    axes[1].plot(d7['WAVELENGTH'], d7['ALPHA'], 'b-')
    axes[1].set_ylabel(r'$\alpha$ (2nd order)')
    
    # Beta (3rd order)
    axes[2].plot(d6['WAVELENGTH'], d6['BETA'], 'k--')
    axes[2].plot(d7['WAVELENGTH'], d7['BETA'], 'g-')
    axes[2].set_ylabel(r'$\beta$ (3rd order)')
    axes[2].set_xlabel('Wavelength (µm)')
    
    plt.tight_layout()
    out = os.path.join(PLOT_DIR, f'fs_v7_v6_comparison_{grating}.png')
    plt.savefig(out)
    print(f'Saved: {out}')

if __name__ == '__main__':
    plt.switch_backend('Agg')
    compare_coeffs('g140m', 'f100lp')
    compare_coeffs('g235m', 'f170lp')
