import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

# Paths
BASE        = '/Users/dcoe/NIRSpec/wavext'
PARLANTI_DIR = f'{BASE}/data/parlanti_repo/calibration_files'
FS_DIR      = f'{BASE}/nirspec_wavext_work/plots/Parlanti/cal/fs_v1'
IFU_DIR      = f'{BASE}/nirspec_wavext_work/plots/Parlanti/cal/ifu_v1'
OUTPUT_DIR   = f'{BASE}/nirspec_wavext_work/reports/329_parlanti-comparison_v1'
os.makedirs(OUTPUT_DIR, exist_ok=True)

def plot_coeff_comparison():
    # 3x2 grid: (k, alpha, beta) rows, (G140M, G235M) columns
    fig, axes = plt.subplots(3, 2, figsize=(15, 14), sharex='col')
    fig.suptitle('Instrumental Coefficient Comparison: Parlanti vs. FS v1 vs. IFU v1', fontsize=18)

    gratings = ['G140M', 'G235M']
    fnames_parlanti = ['calibration_functions_g140m_f100lp.fits', 'calibration_functions_g235m_f170lp.fits']
    
    # Define colors for SOURCES
    colors = {'Parlanti': '0.3', 'FS_v1': 'red', 'IFU_v1': 'c'}
    lstyles = {'Parlanti': '-', 'FS_v1': '--', 'IFU_v1': '-.'}

    for i, (grating, p_fname) in enumerate(zip(gratings, fnames_parlanti)):
        # -- Load Data --
        # Parlanti
        with fits.open(f'{PARLANTI_DIR}/{p_fname}') as h:
            d_p = h[1].data
            wl_p, k_p, a_p, b_p = d_p['wavelength'], d_p['k'], d_p['alpha'], d_p['beta']
        
        # FS v1
        d_fs = np.loadtxt(f'{FS_DIR}/coeffs_fs_v1_{grating}.csv', delimiter=',', skiprows=1)
        wl_fs, k_fs, a_fs, b_fs = d_fs[:, 0], d_fs[:, 1], d_fs[:, 2], d_fs[:, 3]
        
        # IFU v1
        d_ifu = np.loadtxt(f'{IFU_DIR}/coeffs_ifu_v1_{grating}.csv', delimiter=',', skiprows=1)
        wl_ifu, k_ifu, a_ifu, b_ifu = d_ifu[:, 0], d_ifu[:, 1], d_ifu[:, 2], d_ifu[:, 3]
        
        # -- Row 0: k(lambda) --
        ax = axes[0, i]
        ax.plot(wl_p,   k_p,   color=colors['Parlanti'], lw=2.5, label='Parlanti Original')
        ax.plot(wl_fs,  k_fs,  color=colors['FS_v1'],   lw=2.0, ls='-', label='FS v1 (329)')
        ax.plot(wl_ifu, k_ifu, color=colors['IFU_v1'],  lw=2.0, ls='-', label='IFU v1 (329)')
        ax.set_ylabel(r'$k(\lambda)$ (1st order)', fontsize=12)
        ax.set_title(f'{grating}', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.2)
        ax.legend(fontsize=10)
        ax.set_ylim(0, 2.0)
        
        # -- Row 1: alpha(lambda) --
        ax = axes[1, i]
        ax.plot(wl_p,   a_p,   color=colors['Parlanti'], lw=2.0, label='Parlanti Original')
        ax.plot(wl_fs,  a_fs,  color=colors['FS_v1'],   lw=1.5, label='FS v1 (329)')
        ax.plot(wl_ifu, a_ifu, color=colors['IFU_v1'],  lw=1.5, label='IFU v1 (329)')
        ax.set_yscale('log')
        ax.set_ylim(1e-4, 1.0)
        ax.set_ylabel(r'$\tilde{\alpha}(\lambda)$ (2nd order)', fontsize=12)
        ax.grid(True, which='both', alpha=0.1)
        ax.legend(fontsize=10)

        # -- Row 2: beta(lambda) --
        ax = axes[2, i]
        ax.plot(wl_p,   b_p,   color=colors['Parlanti'], lw=2.0, label='Parlanti Original')
        ax.plot(wl_fs,  b_fs,  color=colors['FS_v1'],   lw=1.5, label='FS v1 (329)')
        ax.plot(wl_ifu, b_ifu, color=colors['IFU_v1'],  lw=1.5, label='IFU v1 (329)')
        ax.set_yscale('log')
        ax.set_ylim(1e-4, 1.0)
        ax.set_ylabel(r'$\tilde{\beta}(\lambda)$ (3rd order)', fontsize=12)
        ax.set_xlabel('Wavelength [µm]', fontsize=12)
        ax.grid(True, which='both', alpha=0.1)
        ax.legend(fontsize=10)

    plt.tight_layout()
    outpath = f'{OUTPUT_DIR}/coeff_comparison_by_type.png'
    plt.savefig(outpath, dpi=200, bbox_inches='tight')
    plt.close()
    print(f'Saved: {outpath}')

if __name__ == '__main__':
    plot_coeff_comparison()
