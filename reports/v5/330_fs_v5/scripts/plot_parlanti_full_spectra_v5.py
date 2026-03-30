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
FS_DIR      = f'{BASE}/data'
CALSPEC_DIR = f'{BASE}/data/CALSPEC'
COEFFS_DIR  = f'{BASE}/results/v5'
OUTPUT_DIR  = f'{BASE}/nirspec_wavext_work/reports/v5/330_fs_v5/plots'
os.makedirs(OUTPUT_DIR, exist_ok=True)

C_ANG_S = 2.99792458e18

SOURCES = {
    '1538': {'name': 'P330E',       'calspec': 'p330e_mod_008.fits',    'path': f'{FS_DIR}/PID1538_P330E'},
    '1537': {'name': 'G191-B2B',   'calspec': 'g191b2b_mod_012.fits',  'path': f'{FS_DIR}/PID1537_G191-B2B'},
    '1536': {'name': 'J1743045',    'calspec': '1743045_mod_007.fits',  'path': f'{FS_DIR}/PID1536_J1743045'},
    '6644': {'name': 'NGC2506-G31', 'calspec': 'ngc2506g31_mod_003.fits', 'path': f'{FS_DIR}/PID6644_NGC2506G31'},
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
            good = (wl > 0) & np.isfinite(fl) & (fl > 0) & (dq == 0)
            return wl[good], fl[good], hdr.get('GRATING'), hdr.get('DETECTOR')
    except: return None, None, None, None

def load_calspec_jy(fname):
    with fits.open(f'{CALSPEC_DIR}/{fname}') as h:
        d = h[1].data
        wl_a, flam = d['WAVELENGTH'].astype(float), d['FLUX'].astype(float)
    fnu = flam * wl_a**2 / C_ANG_S * 1e23
    return interp1d(wl_a / 1e4, fnu, bounds_error=False, fill_value=np.nan)

def load_coeffs(grating):
    # Map grating to v5 config key
    key = 'g140m_f100lp' if grating == 'G140M' else 'g235m_f170lp'
    fits_path = f'{COEFFS_DIR}/calib_v5_{key}.fits'
    if not os.path.exists(fits_path): return None, None, None
    with fits.open(fits_path) as h:
        d = h[1].data
        wl, k, a, b = d['WAVELENGTH'], d['K'], d['ALPHA'], d['BETA']
        return interp1d(wl, k, bounds_error=False), interp1d(wl, a, bounds_error=False), interp1d(wl, b, bounds_error=False)

def plot_full_spectra(pid, src):
    name = src['name']
    print(f'  Plotting full spectrum for {name}...')
    f_ref = load_calspec_jy(src['calspec'])
    
    # Discovery logic similar to v5 solver
    files = glob.glob(f"{src['path']}/**/*x1d.fits", recursive=True)
    # Filter for actually relevant files (avoid those with -o prefix if needed, but recursive usually finds the right ones)
    
    fig, ax = plt.subplots(figsize=(18, 9))
    fig.suptitle(f'FS v5 Calibration (4-Source Break) — {name} Full Merged Spectrum', fontsize=16)

    plotted_labels = set()

    for f in files:
        # Avoid some types of duplicates if necessary
        fname = os.path.basename(f)
        if '-o00' in fname: continue # Skip some intermediate products
            
        w, f_val, grating, detector = load_spec(f)
        if w is None or len(w) < 10: continue

        base_color = GRATING_COLORS.get(grating, '0.5')
        alpha = 0.2 if detector == 'NRS2' else 0.6
        lw = 0.7
        
        ax.plot(w, f_val, color=base_color, lw=lw, alpha=alpha)
        
        # Add Label at Mid-point
        label_key = f"{grating}/{detector}"
        if label_key not in plotted_labels:
            if len(w) > 50:
                mid_idx = len(w) // 2
                ax.text(w[mid_idx], f_val[mid_idx] * 1.5, f"{grating} {detector}", color=base_color, fontsize=9, ha='center', weight='bold')
                plotted_labels.add(label_key)

        # Recalibration (The v5 Core)
        if detector == 'NRS2' and grating in ['G140M', 'G235M']:
            ki, ai, bi = load_coeffs(grating)
            if ki:
                lo = 1.93 if grating == 'G140M' else 3.25 # Start where ghosting hits
                mask = w >= lo
                if mask.any():
                    w_e, f_e = w[mask], f_val[mask]
                    f_c2, f_c3 = f_ref(w_e/2), f_ref(w_e/3)
                    f_c2 = np.where(np.isfinite(f_c2), f_c2, 0)
                    f_c3 = np.where(np.isfinite(f_c3), f_c3, 0)
                    k_vals = ki(w_e)
                    k_safe = np.where(k_vals > 0.05, k_vals, 1.0)
                    f_corr = (f_e - ai(w_e)*f_c2 - bi(w_e)*f_c3) / k_safe
                    
                    ax.plot(w_e, f_corr, color='magenta', lw=1.2, alpha=0.8)
                    if f"Rec {grating}" not in plotted_labels:
                        mid_e = len(w_e) // 2
                        ax.text(w_e[mid_e], f_corr[mid_e] * 2.5, f"{grating} v5 Corr", color='magenta', fontsize=10, ha='center', weight='bold')
                        plotted_labels.add(f"Rec {grating}")

    # Truth L3 (MAST)
    truth_files = glob.glob(f"{src['path']}/**/*-s000000001*_x1d.fits", recursive=True)
    for tf in truth_files:
        wt, ft, g, d = load_spec(tf)
        if wt is not None:
            ax.plot(wt, ft, color='black', lw=1.0, alpha=0.9, ls='--')
            if f"Truth {g}" not in plotted_labels:
                 mid_t = len(wt) // 2
                 ax.text(wt[mid_t], ft[mid_t] * 0.4, f"{g} MAST L3", color='black', fontsize=9, ha='center')
                 plotted_labels.add(f"Truth {g}")

    # CALSPEC Model & Ghosts (Unscaled) - Same Color (Goldenrod), Layered Transparency
    w_range = np.linspace(0.6, 5.6, 2000)
    f_truth = f_ref(w_range)
    f_gh2   = f_ref(w_range/2)
    f_gh3   = f_ref(w_range/3)
    
    # 1st order Truth (higher zorder, higher alpha)
    ax.plot(w_range, f_truth, color='goldenrod', lw=1.5, ls='-', alpha=0.6, zorder=5, label='_nolegend_')
    ax.text(1.0, f_ref(1.0) * 1.2, 'CALSPEC', color='goldenrod', weight='bold', fontsize=12, ha='center', zorder=10)
    
    # 2nd order Truth (lower zorder, lower alpha, dashed)
    ax.plot(w_range, f_gh2, color='goldenrod', lw=1.0, ls='--', alpha=0.4, zorder=2, label='_nolegend_')
    ax.text(2.0, f_ref(2.0/2) * 1.2, r'2nd order ($\lambda/2$)', color='goldenrod', weight='bold', fontsize=11, ha='center', alpha=0.8, zorder=10)
    
    # 3rd order Truth (lowest zorder, lowest alpha, dotted)
    ax.plot(w_range, f_gh3, color='goldenrod', lw=1.0, ls=':',  alpha=0.3, zorder=1, label='_nolegend_')
    ax.text(3.0, f_ref(3.0/3) * 1.2, r'3rd order ($\lambda/3$)', color='goldenrod', weight='bold', fontsize=10, ha='center', alpha=0.7, zorder=10)

    ax.set_yscale('log')
    ax.set_xlim(0.6, 5.6)
    
    f_cal = f_ref(w_range)
    f_cal = f_cal[np.isfinite(f_cal) & (f_cal > 0)]
    if len(f_cal) > 0:
        ax.set_ylim(np.min(f_cal)*0.2, np.max(f_cal)*10.0) # Increased top margin for labels

    ax.set_xlabel('Wavelength [µm]', fontsize=14)
    ax.set_ylabel('Flux [Jy]', fontsize=14)
    ax.grid(True, which='both', alpha=0.1)

    outpath = f'{OUTPUT_DIR}/full_spectrum_merged_v5_{name}.png'
    plt.savefig(outpath, dpi=200, bbox_inches='tight')
    plt.close()
    print(f'    Saved: {outpath}')

if __name__ == '__main__':
    for pid, src in SOURCES.items():
        plot_full_spectra(pid, src)
