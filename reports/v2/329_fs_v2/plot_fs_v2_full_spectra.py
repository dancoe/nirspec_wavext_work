import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import interp1d

# Paths
BASE        = '/Users/dcoe/NIRSpec/wavext'
FS_DIR      = f'{BASE}/data'
CALSPEC_DIR = f'{BASE}/data/CALSPEC'
COEFFS_DIR  = f'{BASE}/nirspec_wavext_work/plots/Parlanti/cal/fs_v2'
OUTPUT_DIR  = f'{BASE}/nirspec_wavext_work/reports/329_fs_v2'
os.makedirs(OUTPUT_DIR, exist_ok=True)

C_ANG_S = 2.99792458e18

SOURCES = {
    '1538': {'name': 'P330E',    'calspec': 'p330e_mod_008.fits', 'path': f'{FS_DIR}/PID1538_P330E'},
    '1537': {'name': 'G191-B2B', 'calspec': 'g191b2b_mod_012.fits', 'path': f'{FS_DIR}/PID1537_G191-B2B'},
    '1536': {'name': 'J1743045', 'calspec': '1743045_mod_007.fits', 'path': f'{FS_DIR}/PID1536_J1743045'},
}

GRATING_COLORS = {
    'G140M': 'tab:blue',
    'G235M': 'tab:green',
    'G395M': 'tab:red',
    'PRISM': 'tab:orange',
}

def load_spec(path):
    if not os.path.exists(path): return None, None, None, None
    with fits.open(path) as h:
        hdr = h[0].header
        try:
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
    csv = f'{COEFFS_DIR}/coeffs_fs_v2_{grating}.csv'
    if not os.path.exists(csv): return None, None, None
    data = np.loadtxt(csv, delimiter=',', skiprows=1)
    wl, k, a, b = data[:, 0], data[:, 1], data[:, 2], data[:, 3]
    return interp1d(wl, k, bounds_error=False), interp1d(wl, a, bounds_error=False), interp1d(wl, b, bounds_error=False)

def add_label(ax, x, y, text, color):
    # Position slightly above the line
    ax.text(x, y * 1.5, text, color=color, fontsize=10, ha='center', va='bottom', weight='bold')

def plot_full_spectra(sid, src):
    name = src['name']
    print(f'  Plotting full spectrum for {name}...')
    f_ref = load_calspec_jy(src['calspec'])
    
    files = glob.glob(f"{src['path']}/**/*_x1d.fits", recursive=True)
    fig, ax = plt.subplots(figsize=(18, 9))
    fig.suptitle(f'FS v2 Calibration (329) — {name} Full Merged Spectrum', fontsize=16)

    plotted_labels = set()

    for f in files:
        if 'jw0' + sid + '-' in os.path.basename(f) or 'jw01538-o' in os.path.basename(f) or 'jw01537-o' in os.path.basename(f) or 'jw01536-o' in os.path.basename(f):
            continue
            
        w, f_val, grating, detector = load_spec(f)
        if w is None or len(w) < 10: continue

        base_color = GRATING_COLORS.get(grating, '0.5')
        color = base_color if detector == 'NRS1' else base_color # Maybe same?
        alpha = 0.4 if detector == 'NRS2' else 0.7
        lw = 0.8
        
        ax.plot(w, f_val, color=color, lw=lw, alpha=alpha)
        
        # Add Label at Mid-point
        if f"{grating}/{detector}" not in plotted_labels:
            if len(w) > 20:
                mid_idx = len(w) // 2
                mid_x   = w[mid_idx]
                mid_y   = np.median(f_val[max(0, mid_idx-5):min(len(w), mid_idx+5)])
                add_label(ax, mid_x, mid_y, f"{grating} {detector}", color)
                plotted_labels.add(f"{grating}/{detector}")

        # Recalibration
        if detector == 'NRS2' and grating in ['G140M', 'G235M']:
            ki, ai, bi = load_coeffs(grating)
            if ki:
                lo = 1.93 if grating == 'G140M' else 3.25
                mask = w >= lo
                if mask.any():
                    w_e, f_e = w[mask], f_val[mask]
                    f_c2, f_c3 = f_ref(w_e/2), f_ref(w_e/3)
                    f_c2 = np.where(np.isfinite(f_c2) & (f_c2 > 0), f_c2, 0.0)
                    f_c3 = np.where(np.isfinite(f_c3) & (f_c3 > 0), f_c3, 0.0)
                    k_vals = ki(w_e)
                    k_safe = np.where(np.isfinite(k_vals) & (k_vals > 0.05), k_vals, np.nan)
                    f_corr = (f_e - ai(w_e)*f_c2 - bi(w_e)*f_c3) / k_safe
                    
                    ax.plot(w_e, f_corr, color='magenta', lw=1.5, alpha=0.9)
                    if f"Rec {grating}" not in plotted_labels:
                        mid_idx_e = len(w_e) // 2
                        add_label(ax, w_e[mid_idx_e], np.median(f_corr[max(0, mid_idx_e-5):min(len(w_e), mid_idx_e+5)]), f"{grating} Recal", 'magenta')
                        plotted_labels.add(f"Rec {grating}")

    # Truth L3
    truth_files = glob.glob(f"{src['path']}/*-s000000001*_x1d.fits")
    for tf in truth_files:
        wt, ft, g, d = load_spec(tf)
        if wt is not None:
            ax.plot(wt, ft, color='black', lw=1.2, alpha=0.9, ls='--')
            if f"Truth {g}" not in plotted_labels:
                 mid_idx_t = len(wt) // 2
                 # Move Truth higher than nominal labels
                 add_label(ax, wt[mid_idx_t], np.median(ft[max(0, mid_idx_t-5):min(len(wt), mid_idx_t+5)]) * 1.8, f"{g} MAST Level 3", 'black')
                 plotted_labels.add(f"Truth {g}")

    # CALSPEC
    w_range = np.linspace(0.6, 5.6, 1000)
    ax.plot(w_range, f_ref(w_range), color='goldenrod', lw=0.8, ls=':', alpha=0.9, label='CALSPEC model')
    # Add CALSPEC label near the shorter wavelengths (e.g. 1.0 um)
    add_label(ax, 1.0, f_ref(1.0) * 1.8, 'CALSPEC Model', 'goldenrod')

    ax.set_yscale('log')
    ax.set_xlim(0.6, 5.6)
    
    f_cal_range = f_ref(w_range)
    f_cal_range = f_cal_range[np.isfinite(f_cal_range) & (f_cal_range > 0)]
    if len(f_cal_range) > 0:
        # Tighter y-limits based on CALSPEC (0.4x to 3x)
        ax.set_ylim(np.min(f_cal_range)*0.4, np.max(f_cal_range)*3.0)

    ax.set_xlabel('Wavelength [µm]', fontsize=12)
    ax.set_ylabel('Flux [Jy]', fontsize=12)
    ax.grid(True, which='both', alpha=0.15)

    outpath = f'{OUTPUT_DIR}/fs_v2_full_spectrum_{name}.png'
    plt.savefig(outpath, dpi=200, bbox_inches='tight')
    plt.close()
    print(f'    Saved: {outpath}')

if __name__ == '__main__':
    for sid, src in SOURCES.items():
        plot_full_spectra(sid, src)
