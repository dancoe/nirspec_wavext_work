import os
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import interp1d

# Paths
BASE        = '/Users/dcoe/NIRSpec/wavext'
DATA_DIR    = f'{BASE}/data'
COEFFS_DIR  = f'{BASE}/results/v5'
OUTPUT_DIR  = f'{BASE}/nirspec_wavext_work/reports/v5/ifu_v5/plots'
os.makedirs(OUTPUT_DIR, exist_ok=True)

SOURCES = {
    '2186': {
        'name': 'UGC-5101',
        'path': f'{DATA_DIR}/PID2186_UGC-5101',
        'z': 0.039,
        'gratings': ['G235M', 'G395M']
    },
    '2654': {
        'name': 'SDSSJ0841',
        'path': f'{DATA_DIR}/PID2654_SDSSJ0841',
        'z': 2.96,
        'gratings': ['G140M']
    },
}

GRATING_COLORS = {
    'G140M': 'tab:blue',
    'G235M': 'tab:green',
    'G395M': 'tab:red',
    'PRISM': 'tab:orange',
}

def load_spec(path):
    if not os.path.exists(path): return None, None, None, None
    try:
        with fits.open(path) as h:
            hdr = h[0].header
            d = h[1].data
            wl = d['WAVELENGTH'].astype(float)
            fl = d['FLUX'].astype(float)
            dq = d['DQ'].astype(int) if 'DQ' in d.dtype.names else np.zeros(len(wl))
            good = (wl > 0) & np.isfinite(fl) & (dq == 0)
            # Find grating/detector in headers
            grating = hdr.get('GRATING') or h[1].header.get('GRATING')
            detector = hdr.get('DETECTOR') or h[1].header.get('DETECTOR')
            return wl[good], fl[good], grating, detector
    except: return None, None, None, None

def load_coeffs(grating):
    filt = 'f100lp' if grating == 'G140M' else 'f170lp'
    path = f'{COEFFS_DIR}/calib_v5_{grating.lower()}_{filt}.fits'
    if not os.path.exists(path): return None, None, None
    with fits.open(path) as h:
        d = h[1].data
        wl, k, a, b = d['WAVELENGTH'], d['K'], d['ALPHA'], d['BETA']
        return interp1d(wl, k, bounds_error=False), interp1d(wl, a, bounds_error=False), interp1d(wl, b, bounds_error=False)

def plot_full_spectra_ifu(pid, src):
    name = src['name']
    z = src['z']
    print(f'  Plotting full spectrum for {name} (z={z})...')
    
    fig, ax = plt.subplots(figsize=(16, 8))
    fig.suptitle(f'IFU v5 Validation — {name} Full Merged Spectrum (z={z})', fontsize=16)

    plotted_labels = set()
    
    # Discovery logic: Stage 3 extended and Stage 3 nominal
    search_paths = [
        f"{src['path']}/stage3_ext/**/*x1d.fits",
        f"{src['path']}/stage3_nom/**/*x1d.fits",
        f"{src['path']}/**/*x1d.fits" # Fallback
    ]
    
    files = []
    for p in search_paths:
        files.extend(glob.glob(p, recursive=True))
    files = list(set(files)) # Unique
    
    all_wl = []
    all_fl = []

    for f in sorted(files):
        fname = os.path.basename(f)
        if '-o00' in fname: continue 
            
        w, f_val, grating, detector = load_spec(f)
        if w is None or len(w) < 10: continue

        base_color = GRATING_COLORS.get(grating, '0.5')
        # For IFU cubes, NRS1/NRS2 are often combined, but x1d might still have detector info or be 'BOTH'
        alpha = 0.4
        lw = 1.0
        
        ax.plot(w, f_val, color=base_color, lw=lw, alpha=alpha, label='_nolegend_')
        
        # Add Label
        label_key = f"{grating}"
        if label_key not in plotted_labels:
            if len(w) > 50:
                mid_idx = len(w) // 2
                ax.text(w[mid_idx], np.nanmedian(f_val) * 1.3, f"{grating}", color=base_color, fontsize=10, ha='center', weight='bold')
                plotted_labels.add(label_key)

        # Recalibration (v5)
        # We apply correction to the grating being validated (G140M or G235M)
        if grating in ['G140M', 'G235M']:
            ki, ai, bi = load_coeffs(grating)
            if ki:
                # We need a "truth" estimate for the ghosts. 
                # For science targets, we use the observed spectrum itself (self-correction).
                lo = 1.93 if grating == 'G140M' else 3.25
                mask = w >= lo
                if mask.any():
                    w_e, f_e = w[mask], f_val[mask]
                    
                    # Interpolator for self-correction
                    interp_obs = interp1d(w, f_val, bounds_error=False, fill_value=0.0)
                    
                    f_gh2 = interp_obs(w_e/2)
                    f_gh3 = interp_obs(w_e/3)
                    
                    k_vals = ki(w_e)
                    k_safe = np.where(k_vals > 0.05, k_vals, np.nan)
                    f_corr = (f_e - ai(w_e)*f_gh2 - bi(w_e)*f_gh3) / k_safe
                    
                    ax.plot(w_e, f_corr, color='magenta', lw=1.5, alpha=0.9, label='_nolegend_')
                    if f"Corr {grating}" not in plotted_labels:
                        mid_e = len(w_e) // 2
                        ax.text(w_e[mid_e], np.nanmedian(f_corr) * 1.8, f"{grating} v5 Corr", color='magenta', fontsize=11, ha='center', weight='bold')
                        plotted_labels.add(f"Corr {grating}")

    # Mark key lines for context
    lines = {
        'H-alpha': 0.6563,
        '[OIII]': 0.5007,
        'H-beta': 0.4861,
        'MgII': 0.2799,
        'Pa-alpha': 1.875,
        '3.3 PAH': 3.3
    }
    
    for lname, lwl in lines.items():
        obs_wl = lwl * (1 + z)
        if 0.6 < obs_wl < 5.6:
            ax.axvline(obs_wl, color='gray', ls=':', alpha=0.5)
            ax.text(obs_wl, ax.get_ylim()[1]*0.1 if ax.get_yscale() == 'linear' else 1e-5, lname, rotation=90, verticalalignment='bottom', fontsize=8, color='gray')

    ax.set_yscale('log')
    ax.set_xlim(0.6, 5.6)
    ax.set_xlabel('Wavelength [µm]', fontsize=14)
    ax.set_ylabel('Flux [Jy]', fontsize=14)
    ax.grid(True, which='both', alpha=0.1)

    outpath = f'{OUTPUT_DIR}/ifu_v5_full_spectrum_{name}.png'
    plt.savefig(outpath, dpi=200, bbox_inches='tight')
    plt.close()
    print(f'    Saved: {outpath}')

if __name__ == '__main__':
    for pid, src in SOURCES.items():
        plot_full_spectra_ifu(pid, src)
