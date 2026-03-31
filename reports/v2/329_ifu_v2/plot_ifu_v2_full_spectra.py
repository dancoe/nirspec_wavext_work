import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import interp1d

# Paths
BASE        = '/Users/dcoe/NIRSpec/wavext'
IFU_DIR     = f'{BASE}/data/IFU'
FS_DIR      = f'{BASE}/data'
CALSPEC_DIR = f'{BASE}/data/CALSPEC'
COEFFS_DIR  = f'{BASE}/nirspec_wavext_work/plots/Parlanti/cal/ifu_v2'
OUTPUT_DIR  = f'{BASE}/nirspec_wavext_work/reports/329_ifu_v2'
os.makedirs(OUTPUT_DIR, exist_ok=True)

C_ANG_S = 2.99792458e18

SOURCES = {
    '1538': {'name': 'P330E',    'calspec': 'p330e_mod_008.fits', 'ifu_dir': f'{IFU_DIR}/PID1538_P330E', 'fs_truth': f'{FS_DIR}/PID1538_P330E/jw01538-o160_t002-s000000001_nirspec_f290lp-g395m-s1600a1-sub2048_x1d.fits'},
    '1537': {'name': 'G191-B2B', 'calspec': 'g191b2b_mod_012.fits', 'ifu_dir': f'{IFU_DIR}/PID1537_G191-B2B', 'fs_truth': f'{FS_DIR}/PID1537_G191-B2B/jw01537-o007_t001-s000000001_nirspec_f290lp-g395m-s1600a1-sub2048_x1d.fits'},
    '1536': {'name': 'J1743045', 'calspec': '1743045_mod_007.fits', 'ifu_dir': f'{IFU_DIR}/PID1536_J1743045', 'fs_truth': f'{FS_DIR}/PID1536_J1743045/jw01536-o002_t004-s000000001_nirspec_f290lp-g395m-s1600a1-sub2048_x1d.fits'},
}

GRATING_COLORS = {
    'G140M': 'tab:blue',
    'G235M': 'tab:green',
    'G395M': 'tab:red',
}

def load_spec(path):
    if not os.path.exists(path): return None, None, 'Unknown'
    with fits.open(path) as h:
        hdr = h[0].header
        try:
            d = h[1].data
            wl = d['WAVELENGTH'].astype(float)
            fl = d['FLUX'].astype(float)
            dq = d['DQ'].astype(int) if 'DQ' in d.dtype.names else np.zeros(len(wl))
            good = (wl > 0) & np.isfinite(fl) & (fl > 0) & (dq == 0)
            g_hdr = hdr.get('GRATING', 'Unknown')
            if g_hdr == 'Unknown':
                if 'g140m' in path.lower(): g_hdr = 'G140M'
                elif 'g235m' in path.lower(): g_hdr = 'G235M'
            return wl[good], fl[good], g_hdr
        except: return None, None, 'Unknown'

def load_calspec_jy(fname):
    with fits.open(f'{CALSPEC_DIR}/{fname}') as h:
        d = h[1].data
        wl_a, flam = d['WAVELENGTH'].astype(float), d['FLUX'].astype(float)
    fnu = flam * wl_a**2 / C_ANG_S * 1e23
    return interp1d(wl_a / 1e4, fnu, bounds_error=False, fill_value=np.nan)

def load_coeffs(grating):
    csv = f'{COEFFS_DIR}/coeffs_ifu_v2_{grating}.csv'
    if not os.path.exists(csv): return None, None, None
    data = np.loadtxt(csv, delimiter=',', skiprows=1)
    wl, k, a, b = data[:, 0], data[:, 1], data[:, 2], data[:, 3]
    return interp1d(wl, k, bounds_error=False), interp1d(wl, a, bounds_error=False), interp1d(wl, b, bounds_error=False)

def add_label(ax, x, y, text, color):
    ax.text(x, y * 1.5, text, color=color, fontsize=10, ha='center', va='bottom', weight='bold')

def plot_full_spectra_ifu(sid, src):
    name = src['name']
    print(f'  Plotting full spectrum for {name}...')
    f_ref = load_calspec_jy(src['calspec'])
    
    fig, ax = plt.subplots(figsize=(18, 9))
    fig.suptitle(f'IFU v2 Calibration (329) — {name} Full Merged Spectrum', fontsize=16)

    gratings = ['G140M', 'G235M']
    plotted_labels = set()

    for g in gratings:
        color = GRATING_COLORS.get(g, '0.5')
        path_ext = glob.glob(os.path.join(src['ifu_dir'], 'stage3_ext', f'*_{g.lower()}*x1d.fits'))
        for f in path_ext:
            w, f_val, _ = load_spec(f)
            if w is None or len(w) < 10: continue

            ax.plot(w, f_val, color=color, lw=0.8, alpha=0.4)
            if f"Obs {g}" not in plotted_labels:
                 mid_idx = len(w) // 2
                 add_label(ax, w[mid_idx], np.median(f_val[max(0, mid_idx-5):min(len(w), mid_idx+5)]), f"Obs {g}", color)
                 plotted_labels.add(f"Obs {g}")

            # Recalibrate
            ki, ai, bi = load_coeffs(g)
            if ki:
                lo_ext = 1.87 if g == 'G140M' else 3.15
                mask = w >= lo_ext
                if mask.any():
                    w_e, f_e = w[mask], f_val[mask]
                    f_c2, f_c3 = f_ref(w_e/2), f_ref(w_e/3)
                    f_c2 = np.where(np.isfinite(f_c2) & (f_c2 > 0), f_c2, 0.0)
                    f_c3 = np.where(np.isfinite(f_c3) & (f_c3 > 0), f_c3, 0.0)
                    k_vals = ki(w_e)
                    k_safe = np.where(np.isfinite(k_vals) & (k_vals > 0.05), k_vals, np.nan)
                    f_corr = (f_e - ai(w_e)*f_c2 - bi(w_e)*f_c3) / k_safe
                    
                    ax.plot(w_e, f_corr, color='magenta', lw=1.5, alpha=0.9)
                    if f"Rec {g}" not in plotted_labels:
                        mid_idx_e = len(w_e) // 2
                        add_label(ax, w_e[mid_idx_e], np.median(f_corr[max(0, mid_idx_e-5):min(len(w_e), mid_idx_e+5)]), f"{g} Recal", 'magenta')
                        plotted_labels.add(f"Rec {g}")

    # Truth L3 IFU
    truth_matches = glob.glob(os.path.join(src['ifu_dir'], 'stage3', '*_x1d.fits'))
    for tf in truth_matches:
        wt, ft, tg = load_spec(tf)
        if wt is not None:
             ax.plot(wt, ft, color='black', lw=1.2, alpha=0.9, ls='--')
             if f"Truth IFU {tg}" not in plotted_labels:
                  mid_idx_t = len(wt) // 2
                  # Move Truth higher than nominal labels
                  add_label(ax, wt[mid_idx_t], np.median(ft[max(0, mid_idx_t-5):min(len(wt), mid_idx_t+5)]) * 1.8, f"{tg} MAST Level 3", 'black')
                  plotted_labels.add(f"Truth IFU {tg}")
             
    # FS Truth (G395M)
    wt, ft, tg = load_spec(src['fs_truth'])
    if wt is not None:
        ax.plot(wt, ft, color='darkgreen', lw=1.2, alpha=0.8, ls='-.')
        if 'Truth FS G395M' not in plotted_labels:
             mid_idx_t = len(wt) // 2
             # Move Truth higher
             add_label(ax, wt[mid_idx_t], np.median(ft[max(0, mid_idx_t-5):min(len(wt), mid_idx_t+5)]) * 1.8, "G395M FS Truth", 'darkgreen')
             plotted_labels.add('Truth FS G395M')

    # CALSPEC
    w_range = np.linspace(0.6, 5.6, 1000)
    ax.plot(w_range, f_ref(w_range), color='goldenrod', lw=0.8, ls=':', alpha=0.9, label='CALSPEC model')
    # Label near the shorter end
    add_label(ax, 1.0, f_ref(1.0) * 1.8, 'CALSPEC Model', 'goldenrod')

    ax.set_yscale('log')
    ax.set_xlim(0.6, 5.6)
    
    f_cal_range = f_ref(w_range)
    f_cal_range = f_cal_range[np.isfinite(f_cal_range) & (f_cal_range > 0)]
    if len(f_cal_range) > 0:
        # Tighter y-limits (0.4x to 3x)
        ax.set_ylim(np.min(f_cal_range)*0.4, np.max(f_cal_range)*3.0)

    ax.set_xlabel('Wavelength [µm]', fontsize=12)
    ax.set_ylabel('Flux [Jy]', fontsize=12)
    ax.grid(True, which='both', alpha=0.15)

    outpath = f'{OUTPUT_DIR}/ifu_v2_full_spectrum_{name}.png'
    plt.savefig(outpath, dpi=200, bbox_inches='tight')
    plt.close()
    print(f'    Saved: {outpath}')

if __name__ == '__main__':
    for sid, src in SOURCES.items():
        plot_full_spectra_ifu(sid, src)
