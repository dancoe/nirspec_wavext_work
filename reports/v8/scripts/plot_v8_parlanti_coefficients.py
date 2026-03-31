import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from astropy.io import fits

# --- Paths ---
BASE = '/Users/dcoe/NIRSpec/wavext'
RESULTS_DIR = f'{BASE}/results/v7'
PAR_DIR = f'{BASE}/data/parlanti_repo/calibration_files'
OUT_DIR = f'{BASE}/nirspec_wavext_work/reports/v8/parlanti'

os.makedirs(OUT_DIR, exist_ok=True)

# Colors and Mapping
COLORS = {
    'k': '#9400D3',      # Violet
    'alpha': '#FF8C00',  # Orange
    'beta': '#1E90FF'    # Blue
}

LABELS = {
    'k': r'$k(\lambda)$',
    'alpha': r'$\tilde{\alpha}(\lambda)$',
    'beta': r'$\tilde{\beta}(\lambda)$'
}

# Requested X-Ranges
XRANGES = {
    'g140m': (0.85, 3.7),
    'g235m': (1.5, 5.5)
}

def load_par_fits(grating):
    """Load original Parlanti coefficients from the repository FITS files."""
    fname = (f'calibration_functions_{grating.lower()}_'
             f'{"f100lp" if grating.upper() == "G140M" else "f170lp"}.fits')
    path = os.path.join(PAR_DIR, fname)
    if not os.path.exists(path):
        return None, None, None, None
    with fits.open(path) as h:
        d = h[1].data
        return d['wavelength'], d['k'], d['alpha'], d['beta']

def plot_v8_coefficients(mode, out_filename, add_faint_ref=False):
    """Plot our v8 coefficients (2 panels: G140M, G235M) with all 3 vars."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
    plt.subplots_adjust(wspace=0.1)

    gratings = ['g140m', 'g235m']
    titles = ['G140M/F100LP', 'G235M/F170LP']

    for i, (g, title) in enumerate(zip(gratings, titles)):
        ax = axes[i]
        ax.axhline(1.0, color='black', lw=0.6, alpha=0.3)

        if add_faint_ref:
            wl_ref, k_ref, a_ref, b_ref = load_par_fits(g)
            if wl_ref is not None:
                ax.plot(wl_ref, k_ref, color=COLORS['k'], lw=1.0, alpha=0.3, ls='--')
                ax.plot(wl_ref, a_ref, color=COLORS['alpha'], lw=1.0, alpha=0.3, ls='--')
                ax.plot(wl_ref, b_ref, color=COLORS['beta'], lw=1.0, alpha=0.3, ls='--')

        fname = f'calib_v7_{mode}_{g}_all.fits'
        path = os.path.join(RESULTS_DIR, fname)
        if os.path.exists(path):
            with fits.open(path) as h:
                data = h[1].data
                wl, k, alpha, beta = data['WAVELENGTH'], data['K'], data['ALPHA'], data['BETA']
            
            ax.plot(wl, k, color=COLORS['k'], lw=1.5)
            ax.plot(wl, alpha, color=COLORS['alpha'], lw=1.5)
            ax.plot(wl, beta, color=COLORS['beta'], lw=1.5)
        else:
            ax.text(0.5, 0.5, f"Data not found:\n{fname}", ha='center', va='center', transform=ax.transAxes)

        ax.set_yscale('log')
        ax.set_ylim(1e-4, 10.0)
        ax.set_xlim(XRANGES[g])
        ax.set_xlabel(r'Wavelength [$\mu$m]', fontsize=12)
        ax.set_title(f'{mode.upper()} {title} coeffs (This work v8)', fontsize=14)
        ax.grid(False) 
        
        if i == 0:
            ax.set_ylabel(f'Value ({mode.upper()})', fontsize=12)
            legend_elements = [
                Line2D([0], [0], color=COLORS['k'], lw=1.5, label=LABELS['k']),
                Line2D([0], [0], color=COLORS['alpha'], lw=1.5, label=LABELS['alpha']),
                Line2D([0], [0], color=COLORS['beta'], lw=1.5, label=LABELS['beta']),
            ]
            if add_faint_ref:
                legend_elements.extend([
                    Line2D([0], [0], color='black', lw=1.5, ls='-', label='This work (v8)'),
                    Line2D([0], [0], color='black', lw=1.0, ls='--', alpha=0.3, label='Parlanti+25'),
                ])
            ax.legend(handles=legend_elements, loc='lower left', fontsize=10)

    out_path = os.path.join(OUT_DIR, out_filename)
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out_path}")

def plot_single_coefficient_comparison(mode, var_name):
    """Plot a single coefficient comparison (This work v8 vs Parlanti+25) across 2 panels."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
    plt.subplots_adjust(wspace=0.1)

    gratings = ['g140m', 'g235m']
    titles = ['G140M/F100LP', 'G235M/F170LP']
    
    color = COLORS[var_name]
    label = LABELS[var_name]

    for i, (g, title) in enumerate(zip(gratings, titles)):
        ax = axes[i]
        ax.axhline(1.0, color='black', lw=0.6, alpha=0.3)

        # Parlanti+25 Reference - Increased alpha and fill beneath
        wl_ref, k_ref, a_ref, b_ref = load_par_fits(g)
        ref_data = {'k': k_ref, 'alpha': a_ref, 'beta': b_ref}
        val_ref = ref_data[var_name]
        
        if wl_ref is not None:
            # Cranked up alpha slightly (0.4 vs previous 0.3) and filled beneath
            ax.plot(wl_ref, val_ref, color=color, lw=1.0, alpha=0.45, ls='--')
            ax.fill_between(wl_ref, 1e-6, val_ref, color=color, alpha=0.08)

        # Our v8 derivation
        fname = f'calib_v7_{mode}_{g}_all.fits'
        path = os.path.join(RESULTS_DIR, fname)
        if os.path.exists(path):
            with fits.open(path) as h:
                data = h[1].data
                wl = data['WAVELENGTH']
                v8_data = {'k': data['K'], 'alpha': data['ALPHA'], 'beta': data['BETA']}
                val_v8 = v8_data[var_name]
            
            ax.plot(wl, val_v8, color=color, lw=2.0) # slightly thicker for contrast

        ax.set_yscale('log')
        ax.set_ylim(1e-4, 10.0)
        ax.set_xlim(XRANGES[g])
        ax.set_xlabel(r'Wavelength [$\mu$m]', fontsize=12)
        ax.set_title(f'{mode.upper()} {title} {label} (This work v8)', fontsize=14)
        ax.grid(False) 
        
        if i == 0:
            ax.set_ylabel(f'Value ({label})', fontsize=12)
            legend_elements = [
                Line2D([0], [0], color=color, lw=2.0, ls='-', label=f'This work (v8)'),
                Line2D([0], [0], color=color, lw=1.0, ls='--', alpha=0.45, label='Parlanti+25'),
            ]
            ax.legend(handles=legend_elements, loc='lower left', fontsize=10)

    out_filename = f'parlanti_coefficients_v8_{mode}_{var_name}_comparison.png'
    out_path = os.path.join(OUT_DIR, out_filename)
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out_path}")

def plot_literature_coefficients():
    """Plot original Parlanti (2025) coefficients alone."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
    plt.subplots_adjust(wspace=0.1)

    gratings = ['g140m', 'g235m']
    titles = ['G140M/F100LP', 'G235M/F170LP']

    for i, (g, title) in enumerate(zip(gratings, titles)):
        ax = axes[i]
        ax.axhline(1.0, color='black', lw=0.6, alpha=0.3)
        wl, k, alpha, beta = load_par_fits(g)
        if wl is not None:
            ax.plot(wl, k, color=COLORS['k'], label=LABELS['k'], lw=1.5)
            ax.plot(wl, alpha, color=COLORS['alpha'], label=LABELS['alpha'], lw=1.5)
            ax.plot(wl, beta, color=COLORS['beta'], label=LABELS['beta'], lw=1.5)

        ax.set_yscale('log')
        ax.set_ylim(1e-4, 10.0)
        ax.set_xlim(XRANGES[g])
        ax.set_xlabel(r'Wavelength [$\mu$m]', fontsize=12)
        ax.set_title(f'{title} coeffs (Parlanti+25)', fontsize=14)
        ax.grid(False) 
        if i == 0:
            ax.set_ylabel('Value (Ref)', fontsize=12)
            ax.legend(loc='lower left', fontsize=12)

    out_path = os.path.join(OUT_DIR, 'parlanti_coefficients_literature_ref.png')
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out_path}")

if __name__ == '__main__':
    # Full multi-variable plots
    plot_v8_coefficients('fs', 'parlanti_coefficients_v8_fs.png', add_faint_ref=False)
    plot_v8_coefficients('ifu', 'parlanti_coefficients_v8_ifu.png', add_faint_ref=False)
    plot_v8_coefficients('fs', 'parlanti_coefficients_v8_fs_with_ref.png', add_faint_ref=True)
    plot_v8_coefficients('ifu', 'parlanti_coefficients_v8_ifu_with_ref.png', add_faint_ref=True)
    plot_literature_coefficients()
    for mode in ['fs', 'ifu']:
        for var in ['k', 'alpha', 'beta']:
            plot_single_coefficient_comparison(mode, var)
