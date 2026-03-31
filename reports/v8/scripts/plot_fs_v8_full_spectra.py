
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
DATA_DIR    = f'{BASE}/data/FS'
CALSPEC_DIR = f'{BASE}/data/CALSPEC'
COEFFS_DIR  = f'{BASE}/results/v7'
OUTPUT_DIR  = f'{BASE}/nirspec_wavext_work/reports/v8/fs_v8/plots'
os.makedirs(OUTPUT_DIR, exist_ok=True)

C_ANG_S = 2.99792458e18

SOURCES = {
    '1538': {'name': 'P330E',       'calspec': 'p330e_mod_008.fits',    'dirs': ['PID1538_P330E']},
    '1537': {'name': 'G191-B2B',   'calspec': 'g191b2b_mod_012.fits',  'dirs': ['PID1537_G191-B2B']},
    '1536': {'name': 'J1743045',    'calspec': '1743045_mod_007.fits',  'dirs': ['PID1536_J1743045']},
    '6644': {'name': 'NGC2506-G31', 'calspec': 'ngc2506g31_mod_003.fits', 'dirs': ['PID6644_NGC2506G31', 'PID6644_NGC2506-G31']},
    '1492': {'name': 'IRAS-05248',  'calspec': None,                    'dirs': ['PID1492']},
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
    fits_path = f'{COEFFS_DIR}/calib_v7_fs_{grating.lower()}_all.fits'
    if not os.path.exists(fits_path): return None, None, None
    with fits.open(fits_path) as h:
        d = h[1].data
        wl, k, a, b = d['WAVELENGTH'], d['K'], d['ALPHA'], d['BETA']
        k = np.where(np.isfinite(k), k, 1.0)
        a = np.where(np.isfinite(a), a, 0.0)
        b = np.where(np.isfinite(b), b, 0.0)
        return interp1d(wl, k, bounds_error=False, fill_value=1.0), \
               interp1d(wl, a, bounds_error=False, fill_value=0.0), \
               interp1d(wl, b, bounds_error=False, fill_value=0.0)

def plot_full_spectra(pid, src):
    name = src['name']
    print(f'  Plotting full spectrum for {name}...')
    f_ref = load_calspec_jy(src['calspec']) if src['calspec'] else None
    
    files = []
    for d in src['dirs']:
        files.extend(glob.glob(f"{DATA_DIR}/{d}/**/*x1d.fits", recursive=True))
    
    fig, ax = plt.subplots(figsize=(18, 9))
    fig.suptitle(f'{name} — PID {pid} — FS wavext v8', fontsize=16)

    plotted_labels = set()

    coeffs = {
        'G140M': load_coeffs('G140M'),
        'G235M': load_coeffs('G235M'),
    }

    for f in files:
        fname = os.path.basename(f)
        if '_s2d' in fname: continue
            
        w, f_val, grating, detector = load_spec(f)
        if w is None or len(w) < 10: continue

        base_color = GRATING_COLORS.get(grating, '0.5')
        
        # PRISM special styling
        if grating == 'PRISM':
            ax.plot(w, f_val, color='green', lw=7, alpha=0.3, zorder=2)
        else:
            # Dots at alpha=0.6, Line at alpha=0.2
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

        # Recalibration
        path_lower = f.lower()
        if detector == 'NRS2' and grating in ['G140M', 'G235M'] and ('nrs2_spec2_cal' in path_lower or 'v7' in path_lower or 'ext' in path_lower) and f_ref:
            ki, ai, bi = coeffs[grating]
            if ki is not None:
                lo = 1.90 if grating == 'G140M' else 3.15 
                mask = w >= lo
                if mask.any():
                    w_e, f_e = w[mask], f_val[mask]
                    f_c2, f_c3 = f_ref(w_e/2), f_ref(w_e/3)
                    f_c2 = np.where(np.isfinite(f_c2), f_c2, 0)
                    f_c3 = np.where(np.isfinite(f_c3), f_c3, 0)
                    
                    k_vals = ki(w_e)
                    a_vals = ai(w_e)
                    b_vals = bi(w_e)
                    
                    f_corr = (f_e - a_vals*f_c2 - b_vals*f_c3) / k_vals
                    
                    ax.plot(w_e, f_corr, color='magenta', lw=0.5, alpha=0.9, zorder=10)
                    ax.fill_between(w_e, f_e, f_corr, color='magenta', alpha=0.1, zorder=5)
                    if f"Rec {grating}" not in plotted_labels:
                        mid_e = len(w_e) // 2
                        ax.text(w_e[mid_e], f_corr[mid_e] * 2.5, f"{grating} Corr", color='magenta', fontsize=14, ha='center', weight='bold')
                        plotted_labels.add(f"Rec {grating}")

    # Truth L3 (MAST) if available
    truth_patterns = ["**/*-s00001*_x1d.fits", "**/*_s3d_x1d.fits"]
    for pat in truth_patterns:
        for d in src['dirs']:
            truth_files = glob.glob(f"{DATA_DIR}/{d}/{pat}", recursive=True)
            for tf in truth_files:
                wt, ft, g, det = load_spec(tf)
                if wt is not None and det == 'NRS1':
                    ax.plot(wt, ft, color='black', lw=1, alpha=0.6, ls='--', zorder=4)

    # CALSPEC Model & Ghosts
    w_range = np.linspace(0.6, 5.6, 3000)
    if f_ref:
        f_truth = f_ref(w_range)
        f_gh2   = f_ref(w_range/2)
        f_gh3   = f_ref(w_range/3)
    
        ax.plot(w_range, f_truth, color='black', lw=1, ls='-', alpha=0.5, zorder=1)
        ax.text(1.0, f_ref(1.0) * 0.7, 'CALSPEC', color='black', weight='bold', fontsize=14, ha='center', alpha=0.8)
    
        ax.plot(w_range, f_gh2, ',', markersize=1, color='black', alpha=0.5, zorder=1)
        ax.text(2.2, f_ref(2.2/2) * 1.1, '2nd order (λ/2)', color='0.50', fontsize=14, 
                alpha=0.8, weight='bold', zorder=20)
    
        ax.plot(w_range, f_gh3, ',', markersize=1, color='black', alpha=0.4, zorder=1)
        ax.text(3.7, f_ref(3.7/3) * 1.1, '3rd order (λ/3)', color='0.70', fontsize=14, 
                alpha=0.7, weight='bold', zorder=20)

    ax.set_yscale('log')
    ax.set_xlim(0.6, 5.6)
    
    if f_ref:
        f_mask = f_truth[np.isfinite(f_truth) & (f_truth > 0)]
        if len(f_mask) > 0:
            ax.set_ylim(np.min(f_mask)*0.5, np.max(f_mask)*1.4)
    else:
        ax.autoscale(axis='y', tight=False)

    ax.set_xlabel('Wavelength [µm]', fontsize=14)
    ax.set_ylabel('Flux [Jy]', fontsize=14)
    ax.grid(True, which='both', alpha=0.1)

    outpath = f'{OUTPUT_DIR}/full_spectrum_merged_v8_{name}.png'
    plt.savefig(outpath, dpi=200, bbox_inches='tight')
    plt.close()
    print(f'    Saved: {outpath}')

if __name__ == '__main__':
    for pid, src in SOURCES.items():
        try:
            plot_full_spectra(pid, src)
        except Exception as e:
            print(f"Failed {pid}: {e}")
