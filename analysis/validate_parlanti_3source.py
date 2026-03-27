import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import interp1d

# Paths
DATA_DIR = '/Users/dcoe/NIRSpec/wavext/data'
WAVEXT_DIR = '/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work'

def load_spec(path):
    with fits.open(path) as hdul:
        if 'EXTRACT1D' in hdul:
            data = hdul['EXTRACT1D'].data
            w, f = data['wavelength'], data['flux']
            msk = (w > 0) & np.isfinite(f) & (f > 0)
            return w[msk], f[msk]
    return np.array([]), np.array([])

def validate_grating(grating, boundary_wav):
    print(f"\n--- Validating {grating} ---")
    
    # Load coefficients
    coeff_file = os.path.join(WAVEXT_DIR, f'plots/Parlanti/cal/direct_coefficients_{grating}.txt')
    data = np.loadtxt(coeff_file)
    wav_grid = data[:, 0]
    ks = data[:, 1]
    alphas = data[:, 2]
    betas = data[:, 3]

    # Sources configuration (same as solver)
    sources = {
        'IRAS-05248': {
            'prism': f'{DATA_DIR}/PID1492/jw01492-o001_t001-s000000001_nirspec_clear-prism-s200a1-subs200a1_x1d.fits',
            'ext_n2': f'{DATA_DIR}/PID1492/jw01492003001_03102_00005_nrs2_g140m_extract_1d.fits' if grating == 'G140M' else f'{DATA_DIR}/PID1492/jw01492003001_03104_00004_nrs2_g235m_extract_1d.fits'
        },
        'G191-B2B': {
            'prism': f'{DATA_DIR}/PID1537_G191-B2B/jw01537-o009_t001-s000000001_nirspec_clear-prism-s1600a1-sub512_x1d.fits',
            'ext_n2': f'{DATA_DIR}/PID1537_G191-B2B/jw01537007001_07101_00003_nrs2_extract_1d.fits' if grating == 'G140M' else f'{DATA_DIR}/PID1537_G191-B2B/jw01537007001_09101_00003_nrs2_extract_1d.fits'
        },
        'P330E': {
            'prism': f'{DATA_DIR}/PID1538_P330E/jw01538-o107_t002-s000000001_nirspec_clear-prism-s1600a1-sub512_x1d.fits',
            'ext_n2': f'{DATA_DIR}/PID1538_P330E/jw01538160001_06101_00001_nrs2_extract_1d.fits' if grating == 'G140M' else f'{DATA_DIR}/PID1538_P330E/jw01538160001_08101_00005_nrs2_extract_1d.fits'
        }
    }

    fig, axes = plt.subplots(3, 1, figsize=(12, 18), sharex=True)
    
    for ax, (name, paths) in zip(axes, sources.items()):
        w_prism, f_prism = load_spec(paths['prism'])
        w_ext, c_ext = load_spec(paths['ext_n2'])
        
        # Scaling (consistent with solver)
        mask_prism = (w_prism > boundary_wav - 0.1) & (w_prism < boundary_wav + 0.1)
        mask_ext = (w_ext > boundary_wav) & (w_ext < boundary_wav + 0.2)
        scale = np.median(f_prism[mask_prism]) / np.median(c_ext[mask_ext])
        c_ext_scaled = c_ext * scale
        
        interp_prism = interp1d(w_prism, f_prism, bounds_error=False, fill_value=0)
        
        # Calculate components
        f1 = interp_prism(wav_grid)
        f2 = interp_prism(wav_grid/2)
        f3 = interp_prism(wav_grid/3)
        
        comp1 = ks * f1
        comp2 = alphas * f2
        comp3 = betas * f3
        total_pred = comp1 + comp2 + comp3
        
        # Plotting
        ax.step(w_ext, c_ext_scaled, where='mid', color='gray', alpha=0.3, label='Observed (NRS2 Scaled)')
        ax.plot(wav_grid, total_pred, 'k-', lw=1.5, label='Total Predicted Model')
        ax.plot(wav_grid, comp1, 'r--', lw=1, label=r'$k \cdot f(\lambda)$ (1st order)')
        ax.plot(wav_grid, comp2, 'g--', lw=1, label=r'$\tilde{\alpha} \cdot f(\lambda/2)$ (2nd order)')
        ax.plot(wav_grid, comp3, 'b--', lw=1, label=r'$\tilde{\beta} \cdot f(\lambda/3)$ (3rd order)')
        
        # Baseline ground truth
        ax.step(w_prism, f_prism, where='mid', color='cyan', alpha=0.5, label='Ground Truth (PRISM)')
        
        ax.set_title(f'Source-Specific Validation: {name} ({grating})')
        ax.set_ylabel('Flux (scaled)')
        ax.legend(loc='upper right', fontsize=8)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(wav_grid[0]-0.1, wav_grid[-1]+0.1)

    axes[-1].set_xlabel('Wavelength (um)')
    plt.tight_layout()
    
    out = os.path.join(WAVEXT_DIR, f'plots/Parlanti/cal/validation_results_{grating}.png')
    plt.savefig(out, dpi=200)
    plt.close()
    print(f"SUCCESS: Validation plot saved to {out}")

if __name__ == "__main__":
    validate_grating('G235M', 3.3)
    validate_grating('G140M', 1.9)
