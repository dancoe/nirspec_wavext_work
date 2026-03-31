
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
FS_DIR      = f'{BASE}/data/FS'
CALSPEC_DIR = f'{BASE}/data/CALSPEC'
COEFFS_DIR  = f'{BASE}/results/v7'
OUTPUT_DIR  = f'{BASE}/nirspec_wavext_work/reports/v9/ifu_v9/plots'
os.makedirs(OUTPUT_DIR, exist_ok=True)

C_ANG_S = 2.99792458e18

# v9 IFU Targets (Uniform 0.5" Aperture)
# SDSSJ0841 (2654) EXCLUDED
TARGETS = [
    {'pid': '1536', 'name': 'J1743045',  'calspec': '1743045_mod_007.fits',  'dir': f'{IFU_DIR}/PID1536_J1743045/stage3_ext'},
    {'pid': '1537', 'name': 'G191-B2B',  'calspec': 'g191b2b_mod_012.fits', 'dir': f'{IFU_DIR}/PID1537_G191-B2B/stage3_ext'},
    {'pid': '1538', 'name': 'P330E',     'calspec': 'p330e_mod_008.fits',    'dir': f'{IFU_DIR}/PID1538_P330E/stage3_ext'},
    {'pid': '2186', 'name': 'UGC-5101',  'calspec': None,                    'dir': f'{FS_DIR}/PID2186_UGC-5101', 'science': True},
    {'pid': '6645', 'name': 'P330E-C3',  'calspec': 'p330e_mod_008.fits',    'dir': f'{IFU_DIR}/PID6645_P330E-C3/stage3_ext'},
]

GRATING_COLORS = {
    'G140M': 'tab:blue',
    'G235M': 'goldenrod',
    'G395M': 'tab:red',
}

def load_calspec_jy(fname):
    if not fname: return None
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
        k = np.where(np.isfinite(k), k, 1.0)
        a = np.where(np.isfinite(a), a, 0.0)
        b = np.where(np.isfinite(b), b, 0.0)
        return interp1d(wl, k, bounds_error=False, fill_value=1.0), \
               interp1d(wl, a, bounds_error=False, fill_value=0.0), \
               interp1d(wl, b, bounds_error=False, fill_value=0.0)

def extract_05_radius(s3d_file):
    try:
        with fits.open(s3d_file) as hdul:
            sci = hdul['SCI'].data
            hdr = hdul['SCI'].header
            wl = np.arange(hdr['NAXIS3']) * hdr['CDELT3'] + hdr['CRVAL3']
            pixar_sr = hdr.get('PIXAR_SR', 1.0)
            conv = 1e6 * pixar_sr if hdr.get('BUNIT') == 'MJy/sr' else 1.0
            
            white_light = np.nanmean(sci, axis=0)
            peak_y, peak_x = np.unravel_index(np.nanargmax(white_light), white_light.shape)
            yy, xx = np.indices(white_light.shape)
            dist_arcsec = np.sqrt(((xx - peak_x)*hdr['CDELT1'])**2 + ((yy - peak_y)*hdr['CDELT2'])**2) * 3600
            mask = dist_arcsec <= 0.5
            
            spec = np.nansum(sci[:, mask], axis=1) * conv
            return wl, spec
    except Exception as e:
        print(f"Error extracting {s3d_file}: {e}")
        return None, None

def plot_full_spectra(target):
    name = target['name']
    pid = target['pid']
    print(f'  Plotting v9 full spectrum for {name} (PID {pid})...')
    f_ref = load_calspec_jy(target['calspec'])
    
    fig, ax = plt.subplots(figsize=(18, 9))
    fig.suptitle(f'{name} — PID {pid} — IFU wavext v9 — Standard 0.5" Extraction', fontsize=16)

    ax.set_yscale('log')

    y_min_target, y_max_target = 1e-3, 10.0
    if f_ref:
        w_range_pre = np.linspace(0.6, 5.6, 1000)
        f_truth_pre = f_ref(w_range_pre)
        f_mask_pre = f_truth_pre[np.isfinite(f_truth_pre) & (f_truth_pre > 0)]
        if len(f_mask_pre) > 0:
            y_min_target = np.nanmin(f_mask_pre) * 0.5
            y_max_target = np.nanmax(f_mask_pre) * 1.5
            ax.set_ylim(y_min_target, y_max_target)

    plotted_labels = set()
    gratings = ['G140M', 'G235M', 'G395M']
    
    coeffs = {
        'G140M': load_coeffs('G140M'),
        'G235M': load_coeffs('G235M'),
    }

    for g_name in gratings:
        s3d_files = glob.glob(f"{target['dir']}/**/*{g_name.lower()}*s3d.fits", recursive=True)
        if not s3d_files: continue
        
        # Unified for v9: Always extract 0.5" radius
        w, f_val = extract_05_radius(s3d_files[0])
        
        if w is None: continue

        base_color = GRATING_COLORS.get(g_name, '0.5')
        
        # Gap aware plotting
        if g_name == 'G140M':
            m1 = w < 2.17
            m2 = w > 2.28
            ax.plot(w[m1], f_val[m1], '.', markersize=1, color=base_color, alpha=0.6, zorder=2)
            ax.plot(w[m1], f_val[m1], color=base_color, lw=0.5, alpha=0.2, zorder=2)
            ax.plot(w[m2], f_val[m2], '.', markersize=1, color=base_color, alpha=0.6, zorder=2)
            ax.plot(w[m2], f_val[m2], color=base_color, lw=0.5, alpha=0.2, zorder=2)
        elif g_name == 'G235M':
            m1 = w < 3.66
            m2 = w > 3.82
            ax.plot(w[m1], f_val[m1], '.', markersize=1, color=base_color, alpha=0.6, zorder=2)
            ax.plot(w[m1], f_val[m1], color=base_color, lw=0.5, alpha=0.2, zorder=2)
            ax.plot(w[m2], f_val[m2], '.', markersize=1, color=base_color, alpha=0.6, zorder=2)
            ax.plot(w[m2], f_val[m2], color=base_color, lw=0.5, alpha=0.2, zorder=2)
        else:
            ax.plot(w, f_val, '.', markersize=1, color=base_color, alpha=0.6, zorder=2)
            ax.plot(w, f_val, color=base_color, lw=0.5, alpha=0.2, zorder=2)
        
        # Correction logic
        if g_name in ['G140M', 'G235M'] and f_ref:
            ki, ai, bi = coeffs[g_name]
            if ki is not None:
                lo = 1.90 if g_name == 'G140M' else 3.15 
                mask = w >= lo
                if mask.any():
                    w_e, f_e = w[mask], f_val[mask]
                    f_c2, f_c3 = f_ref(w_e/2), f_ref(w_e/3)
                    f_c2 = np.where(np.isfinite(f_c2), f_c2, 0)
                    f_c3 = np.where(np.isfinite(f_c3), f_c3, 0)
                    k_vals, a_vals, b_vals = ki(w_e), ai(w_e), bi(w_e)
                    f_corr = (f_e - a_vals*f_c2 - b_vals*f_c3) / k_vals
                    
                    f_e_fill    = np.clip(f_e,    y_min_target, y_max_target * 10)
                    f_corr_fill = np.clip(f_corr, y_min_target, y_max_target * 10)

                    ax.plot(w_e, f_corr, color='magenta', lw=0.5, alpha=0.9, zorder=10)
                    ax.fill_between(w_e, f_e_fill, f_corr_fill, color='magenta', alpha=0.1, zorder=5)
                    if f"Rec {g_name}" not in plotted_labels:
                        mid_e = len(w_e) // 2
                        ax.text(w_e[mid_e], f_corr[mid_e] * 2.5, f"{g_name} Corr", color='magenta', fontsize=14, ha='center', weight='bold')
                        plotted_labels.add(f"Rec {g_name}")
        
        # Add Label at Mid-point
        if g_name not in plotted_labels and len(w) > 50:
            mid_idx = len(w) // 2
            ax.text(w[mid_idx], f_val[mid_idx] * 1.5, f"{g_name}", color=base_color, fontsize=14, ha='center', weight='bold')
            plotted_labels.add(g_name)

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

    ax.set_xlim(0.6, 5.6)
    if not f_ref:
        ax.autoscale(axis='y', tight=False)

    ax.set_xlabel('Wavelength [µm]', fontsize=14)
    ax.set_ylabel('Flux [Jy]', fontsize=14)
    ax.grid(True, which='both', alpha=0.1)
    ax.text(0.02, 0.02, "Uniform v9 Extraction: r=0.5\" fixed circular", transform=ax.transAxes, fontsize=12, color='blue', alpha=0.7)

    outpath = f'{OUTPUT_DIR}/full_spectrum_merged_v9_{name}.png'
    plt.savefig(outpath, dpi=200, bbox_inches='tight')
    plt.close()
    print(f'    Saved: {outpath}')

if __name__ == '__main__':
    for target in TARGETS:
        try:
            plot_full_spectra(target)
        except Exception as e:
            print(f"Failed {target['pid']}: {e}")
