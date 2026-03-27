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

def plot_summary():
    # Load coefficients for G140M and G235M
    coeffs = {}
    for g in ['G140M', 'G235M']:
        path = os.path.join(WAVEXT_DIR, f'plots/Parlanti/cal/3source/direct_coefficients_{g}.txt')
        if os.path.exists(path):
            coeffs[g] = np.loadtxt(path)

    # Pick a representative target for the "Calibrated" view - IRAS-05248 or G191-B2B
    target = 'P330E'
    sources = {
        'G140M': {
            'prism': f'{DATA_DIR}/PID1538_P330E/jw01538-o107_t002-s000000001_nirspec_clear-prism-s1600a1-sub512_x1d.fits',
            'ext_n2': f'{DATA_DIR}/PID1538_P330E/jw01538160001_06101_00001_nrs2_extract_1d.fits',
            'boundary': 1.9,
            'ylim': [1e-4, 0.5]
        },
        'G235M': {
            'prism': f'{DATA_DIR}/PID1538_P330E/jw01538-o107_t002-s000000001_nirspec_clear-prism-s1600a1-sub512_x1d.fits',
            'ext_n2': f'{DATA_DIR}/PID1538_P330E/jw01538160001_08101_00005_nrs2_extract_1d.fits',
            'boundary': 3.3,
            'ylim': [1e-4, 0.5]
        }
    }

    fig, axes = plt.subplots(2, 1, figsize=(10, 10), sharex=False)

    # --- TOP PANEL: COEFFICIENTS (Global) ---
    ax = axes[0]
    for g, c in coeffs.items():
        w, k, a, b = c[:,0], c[:,1], c[:,2], c[:,3]
        ls = '-' if g == 'G235M' else '--'
        ax.plot(w, k, color='C0', ls=ls, label=f'k ({g})', lw=2)
        ax.plot(w, a, color='C1', ls=ls, label=f'alpha ({g})', lw=1.5)
        ax.plot(w, b, color='C2', ls=ls, label=f'beta ({g})', lw=1.5)
    
    ax.set_ylim(-0.05, 1.1)
    ax.set_title('Global Direct Solver Coefficients ($k, \\tilde{\\alpha}, \\tilde{\\beta}$)')
    ax.set_xlabel('Wavelength (µm)')
    ax.set_ylabel('Throughput')
    ax.legend(ncol=2, fontsize=9)
    ax.grid(True, alpha=0.3)

    # --- BOTTOM PANEL: CALIBRATED SPECTRA (Representative target P330E) ---
    ax = axes[1]
    all_cal = []
    for grating in ['G140M', 'G235M']:
        if grating not in sources: continue
        p = sources[grating]
        w_p, f_p = load_spec(p['prism'])
        w_ext, c_ext = load_spec(p['ext_n2'])
        
        # Scaling
        mask_prism = (w_p > p['boundary'] - 0.1) & (w_p < p['boundary'] + 0.1)
        mask_ext = (w_ext > p['boundary']) & (w_ext < p['boundary'] + 0.2)
        scale = np.median(f_p[mask_prism]) / np.median(c_ext[mask_ext])
        c_ext_scaled = c_ext * scale
        
        # Apply the universal calibration
        c = coeffs[grating]
        w_c, ks, alphas, betas = c[:,0], c[:,1], c[:,2], c[:,3]
        interp_prism = interp1d(w_p, f_p, bounds_error=False, fill_value=0)
        
        # Predicted components
        comp2 = np.interp(w_c, w_c, alphas) * interp_prism(w_c/2)
        comp3 = np.interp(w_c, w_c, betas) * interp_prism(w_c/3)
        
        # S_cal = (S_obs - alpha*f/2 - beta*f/3) / k
        # Interpolate observed spectrum onto the coefficient grid
        s_obs_interp = np.interp(w_c, w_ext, c_ext_scaled)
        # Avoid division by zero
        k_safe = np.maximum(ks, 0.01)
        s_cal = (s_obs_interp - comp2 - comp3) / k_safe
        
        # Plot
        ax.plot(w_ext, c_ext_scaled, color='lightgreen', alpha=0.5, label='_nolegend_') # Uncalibrated
        ax.plot(w_c, s_cal, color='red', lw=1.5, label='_nolegend_') # Calibrated
        ax.plot(w_p, f_p, 'k--', alpha=0.5, label='_nolegend_') # Ground truth
        
    ax.plot([], [], 'lightgreen', label='Uncalibrated (NRS2 Scaled)')
    ax.plot([], [], 'red', lw=1.5, label='Calibrated (Inverse Solve)')
    ax.plot([], [], 'k--', label='Ground Truth (PRISM)')

    ax.set_yscale('log')
    ax.set_ylim(2e-4, 0.5)
    ax.set_xlim(1.5, 5.6)
    ax.set_title(f'Calibration Application Example: {target}')
    ax.set_xlabel('Wavelength (µm)')
    ax.set_ylabel('Flux (scaled)')
    ax.legend(fontsize=9)
    ax.grid(True, which='both', alpha=0.3)

    plt.tight_layout()
    out = os.path.join(WAVEXT_DIR, 'plots/Parlanti/cal/3source/FS_3source_cal.png')
    plt.savefig(out, dpi=200)
    plt.close()
    print(f"SUCCESS: Plot saved to {out}")

if __name__ == "__main__":
    plot_summary()
