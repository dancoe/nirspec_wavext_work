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
IFU_DIR     = f'{BASE}/data/IFU'
CALSPEC_DIR = f'{BASE}/data/CALSPEC'
COEFFS_DIR  = f'{BASE}/results/v7'
OUTPUT_DIR  = f'{BASE}/nirspec_wavext_work/reports/v7/ifu_v7/plots'
os.makedirs(OUTPUT_DIR, exist_ok=True)

C_ANG_S = 2.99792458e18

SOURCES = {
    '1537': {'name': 'G191-B2B',  'calspec': 'g191b2b_mod_012.fits', 'ifu_dir': 'PID1537_G191-B2B'},
    '1538': {'name': 'P330E',     'calspec': 'p330e_mod_008.fits',    'ifu_dir': 'PID1538_P330E'},
    '1536': {'name': 'J1743045',  'calspec': '1743045_mod_007.fits',  'ifu_dir': 'PID1536_J1743045'},
}

GRATING_COLORS = {
    'G140M': 'tab:blue',
    'G235M': 'goldenrod',
    'G395M': 'tab:red',
    'G395H': 'tab:red',
    'PRISM': 'm',
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
            return wl[good], fl[good], hdr.get('GRATING'), hdr.get('DETECTOR')
    except: return None, None, None, None

def load_calspec_jy(fname):
    with fits.open(f'{CALSPEC_DIR}/{fname}') as h:
        d = h[1].data
        wl_a, flam = d['WAVELENGTH'].astype(float), d['FLUX'].astype(float)
    fnu = flam * wl_a**2 / C_ANG_S * 1e23
    return interp1d(wl_a / 1e4, fnu, bounds_error=False, fill_value=np.nan)

def load_coeffs(grating):
    fits_path = f'{COEFFS_DIR}/calib_v7_ifu_{grating.lower()}_all.fits'
    if not os.path.exists(fits_path): return None, None, None
    with fits.open(fits_path) as h:
        d = h[1].data
        wl, k, a, b = d['WAVELENGTH'], d['K'], d['ALPHA'], d['BETA']
        return interp1d(wl, k, bounds_error=False, fill_value=1.0), \
               interp1d(wl, a, bounds_error=False, fill_value=0.0), \
               interp1d(wl, b, bounds_error=False, fill_value=0.0)

def plot_full_spectra(pid, src):
    name = src['name']
    print(f'  Plotting full spectrum for {name}...')
    f_ref = load_calspec_jy(src['calspec']) if src['calspec'] else None if src['calspec'] else None
    
    # Discovery logic
    files = glob.glob(f"{IFU_DIR}/{src['ifu_dir']}/**/*x1d.fits", recursive=True)
    
    fig, ax = plt.subplots(figsize=(18, 9))
    fig.suptitle(f'IFU v7 Calibration — {name} Full Merged Spectrum', fontsize=16)

    plotted_labels = set()

    # Pre-load coefficients
    coeffs = {
        'G140M': load_coeffs('G140M'),
        'G235M': load_coeffs('G235M'),
    }

    for f in files:
        fname = os.path.basename(f)
        w, f_val, grating, detector = load_spec(f)
        if w is None or len(w) < 10: continue

        base_color = GRATING_COLORS.get(grating, '0.5')
        alpha = 0.4 if detector == 'NRS2' else 0.7
        lw = 0.8
        
        # New Style for PRISM vs others
        if grating == 'PRISM':
            ax.plot(w, f_val, color='green', lw=5, alpha=0.2, zorder=2)
        else:
            # Standard Style: Dots at alpha=0.6, Line at alpha=0.2
            ax.plot(w, f_val, '.', markersize=1, color=base_color, alpha=0.6, zorder=2)
            ax.plot(w, f_val, color=base_color, lw=0.5, alpha=0.2, zorder=2)
        
        # Add Label at Mid-point
        label_key = f"{grating}"
        if label_key not in plotted_labels:
            if len(w) > 50:
                mid_idx = len(w) // 2
                l_color = 'green' if grating == 'PRISM' else base_color
                ax.text(w[mid_idx], f_val[mid_idx] * 1.5, f"{grating}", color=l_color, fontsize=14, ha='center', weight='bold')
                plotted_labels.add(label_key)

        # Recalibration (Extension logic) - Only for IFU NRS2
        if detector == 'NRS2' and grating in ['G140M', 'G235M']:
            ki, ai, bi = coeffs[grating]
            if ki is not None:
                lo = 1.87 if grating == 'G140M' else 3.15 
                mask = w >= lo
                if mask.any():
                    w_e, f_e = w[mask], f_val[mask]
                    f_c2, f_c3 = f_ref(w_e/2), f_ref(w_e/3)
                    f_c2 = np.where(np.isfinite(f_c2), f_c2, 0)
                    f_c3 = np.where(np.isfinite(f_c3), f_c3, 0)
                    
                    k_vals = ki(w_e)
                    a_vals = ai(w_e)
                    b_vals = bi(w_e)
                    
                    # Correction: f_corr = (S_obs - a*f2 - b*f3) / k
                    f_corr = (f_e - a_vals*f_c2 - b_vals*f_c3) / k_vals
                    
                    ax.plot(w_e, f_corr, color='magenta', lw=1.5, alpha=0.9, zorder=10)
                    if f"Rec {grating}" not in plotted_labels:
                        mid_e = len(w_e) // 2
                        ax.text(w_e[mid_e], f_corr[mid_e] * 2.5, f"{grating} v7 Corr", color='magenta', fontsize=14, ha='center', weight='bold')
                        plotted_labels.add(f"Rec {grating}")

    # Truth baseline (NRS1)
    for f in files:
        w, f_val, grating, detector = load_spec(f)
        if w is not None and detector == 'NRS1':
            ax.plot(w, f_val, color='black', lw=1.0, alpha=0.5, ls='--', zorder=4)

    # CALSPEC
    w_range = np.linspace(0.6, 5.6, 3000)
    f_truth = f_ref(w_range)
    ax.plot(w_range, f_truth, color='black', lw=2.0, ls='-', alpha=0.5, zorder=1)
    ax.text(1.0, f_ref(1.0) * 0.7, 'CALSPEC', color='black', weight='bold', fontsize=14, ha='center', alpha=0.8)
    
    # Ghost curves with labels
    # Ghost curves with labels
    # Use the new requested style: pixels (alpha 0.5/0.4), line at alpha=0.2
    # 2nd Order
    ax.plot(w_range, f_gh2, ',', markersize=1, color='black', alpha=0.5, zorder=1)
    ax.plot(w_range, f_gh2, color='black', lw=0.5, ls='--', alpha=0.2, zorder=1)
    ax.text(2.2, f_ref(2.2/2) * 1.1, '2nd order (λ/2)', color='black', fontsize=14, 
            alpha=0.8, weight='bold', zorder=20)
    
    # 3rd Order
    ax.plot(w_range, f_gh3, ',', markersize=1, color='black', alpha=0.4, zorder=1)
    ax.plot(w_range, f_gh3, color='black', lw=0.5, ls=':',  alpha=0.2, zorder=1)
    ax.text(3.7, f_ref(3.7/3) * 1.1, '3rd order (λ/3)', color='black', fontsize=14, 
            alpha=0.7, weight='bold', zorder=20)
    
    ax.set_yscale('log')
    ax.set_xlim(0.6, 5.6)
    
    # Auto-scale Y based on truth
    f_mask = f_truth[np.isfinite(f_truth) & (f_truth > 0)]
    if len(f_mask) > 0:
        ax.set_ylim(np.min(f_mask)*0.5, np.max(f_mask)*3.0)

    ax.set_xlabel('Wavelength [µm]', fontsize=14)
    ax.set_ylabel('Flux [Jy]', fontsize=14)
    ax.grid(True, which='both', alpha=0.1)

    outpath = f'{OUTPUT_DIR}/full_spectrum_merged_v7_{name}.png'
    plt.savefig(outpath, dpi=200, bbox_inches='tight')
    plt.close()
    print(f'    Saved: {outpath}')

if __name__ == '__main__':
    for pid, src in SOURCES.items():
        plot_full_spectra(pid, src)
