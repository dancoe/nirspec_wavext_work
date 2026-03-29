import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import interp1d
from scipy.ndimage import median_filter

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE        = '/Users/dcoe/NIRSpec/wavext'
FS_DIR      = f'{BASE}/data'
CALSPEC_DIR = f'{BASE}/data/CALSPEC'
COEFFS_DIR  = f'{BASE}/nirspec_wavext_work/plots/Parlanti/cal/fs_v1'
OUTPUT_DIR  = f'{BASE}/nirspec_wavext_work/reports/329_fs_v1'
os.makedirs(OUTPUT_DIR, exist_ok=True)

C_ANG_S = 2.99792458e18  # Å/s

# ── Source Table ──────────────────────────────────────────────────────────────
SOURCES = {
    '1538': {
        'name':     'P330E',
        'calspec':  'p330e_mod_008.fits',
        'nrs2_g140m': f'{FS_DIR}/PID1538_P330E/nrs2_spec2_cal/jw01538160001_06101_00001_nrs2_x1d.fits',
        'nrs2_g235m': f'{FS_DIR}/PID1538_P330E/nrs2_spec2_cal/jw01538160001_08101_00005_nrs2_x1d.fits',
        'l3_g235m': f'{FS_DIR}/PID1538_P330E/jw01538-o160_t002-s000000001_nirspec_f170lp-g235m-s1600a1-sub2048_x1d.fits',
        'l3_g395m': f'{FS_DIR}/PID1538_P330E/jw01538-o160_t002-s000000001_nirspec_f290lp-g395m-s1600a1-sub2048_x1d.fits',
    },
    '1537': {
        'name':     'G191-B2B',
        'calspec':  'g191b2b_mod_012.fits',
        'nrs2_g140m': f'{FS_DIR}/PID1537_G191-B2B/nrs2_spec2_cal/jw01537007001_07101_00003_nrs2_x1d.fits',
        'nrs2_g235m': f'{FS_DIR}/PID1537_G191-B2B/nrs2_spec2_cal/jw01537007001_09101_00003_nrs2_x1d.fits',
        'l3_g235m': f'{FS_DIR}/PID1537_G191-B2B/jw01537-o007_t001-s000000001_nirspec_f170lp-g235m-s1600a1-sub2048_x1d.fits',
        'l3_g395m': f'{FS_DIR}/PID1537_G191-B2B/jw01537-o007_t001-s000000001_nirspec_f290lp-g395m-s1600a1-sub2048_x1d.fits',
    },
    '1536': {
        'name':     'J1743045',
        'calspec':  '1743045_mod_007.fits',
        'nrs2_g140m': f'{FS_DIR}/PID1536_J1743045/nrs2_spec2_cal/jw01536002001_03104_00003_nrs2_x1d.fits',
        'nrs2_g235m': f'{FS_DIR}/PID1536_J1743045/nrs2_spec2_cal/jw01536002001_03106_00001_nrs2_x1d.fits',
        'l3_g235m': f'{FS_DIR}/PID1536_J1743045/jw01536-o002_t004-s000000001_nirspec_f170lp-g235m-s1600a1-sub2048_x1d.fits',
        'l3_g395m': f'{FS_DIR}/PID1536_J1743045/jw01536-o002_t004-s000000001_nirspec_f290lp-g395m-s1600a1-sub2048_x1d.fits',
    },
}

# ── Helpers ────────────────────────────────────────────────────────────────────
def load_spec(path):
    if not os.path.exists(path):
        return np.array([]), np.array([])
    with fits.open(path) as h:
        ext = None
        for name in ['EXTRACT1D', 'EXTRACT1D']:
            if name in [hh.name for hh in h]:
                ext = h[name].data
                break
        if ext is None:
            ext = h[1].data
        wl = ext['WAVELENGTH'].astype(float)
        fl = ext['FLUX'].astype(float)
        dq = ext['DQ'].astype(int) if 'DQ' in ext.dtype.names else np.zeros(len(wl), int)
    good = (wl > 0) & np.isfinite(fl) & (fl > 0) & (dq == 0)
    return wl[good], fl[good]

def load_calspec_jy(fname):
    with fits.open(f'{CALSPEC_DIR}/{fname}') as h:
        d    = h[1].data
        wl_a = d['WAVELENGTH'].astype(float)
        flam = d['FLUX'].astype(float)
    fnu_jy = flam * (wl_a**2) / C_ANG_S * 1e23
    order  = np.argsort(wl_a)
    return wl_a[order] / 1e4, fnu_jy[order]   # µm, Jy

def plot_source_spectra(sid, src):
    name = src['name']
    print(f'Plotting {name}...')
    
    fig, axes = plt.subplots(2, 1, figsize=(14, 10), sharex=False)
    fig.suptitle(f'FS v1 Calibration (329) — {name} Spectra Comparison', fontsize=14)

    # CALSPEC Reference
    wl_cs, fj_cs = load_calspec_jy(src['calspec'])
    f_ref = interp1d(wl_cs, fj_cs, bounds_error=False, fill_value=np.nan)

    for i, grating in enumerate(['G140M', 'G235M']):
        ax = axes[i]
        
        # Load Coefficients
        csv_path = f'{COEFFS_DIR}/coeffs_fs_v1_{grating}.csv'
        data = np.loadtxt(csv_path, delimiter=',', skiprows=1)
        wl_grid = data[:, 0]
        ks      = data[:, 1]
        alphas  = data[:, 2]
        betas   = data[:, 3]
        
        k_interp = interp1d(wl_grid, ks,     bounds_error=False, fill_value=np.nan)
        a_interp = interp1d(wl_grid, alphas, bounds_error=False, fill_value=np.nan)
        b_interp = interp1d(wl_grid, betas,  bounds_error=False, fill_value=np.nan)

        # Load Data
        nrs2_path  = src['nrs2_g140m' if grating == 'G140M' else 'nrs2_g235m']
        truth_path = src['l3_g235m' if grating == 'G140M' else 'l3_g395m']
        
        w_nrs2,  f_nrs2  = load_spec(nrs2_path)
        w_truth, f_truth = load_spec(truth_path)

        # Truncate truth to typical overlap area
        if grating == 'G140M':
            mask = (w_truth > 1.95) & (w_truth < 3.50)
        else:
            mask = (w_truth > 3.30) & (w_truth < 5.25)
        w_truth_overlap, f_truth_overlap = w_truth[mask], f_truth[mask]

        # Calculate Corrected Spec
        ki = k_interp(w_nrs2)
        ai = a_interp(w_nrs2)
        bi = b_interp(w_nrs2)
        f_cs2 = f_ref(w_nrs2 / 2)
        f_cs3 = f_ref(w_nrs2 / 3)
        
        f_corr = (f_nrs2 - ai * f_cs2 - bi * f_cs3) / ki

        # Plotting
        ax.plot(w_nrs2,  f_nrs2,  'k', lw=1.0, alpha=0.3, label=f'Observed {grating} (Extended NRS2)')
        if len(w_truth_overlap) > 0:
            truth_label = f'Truth (Nominal {"G235M" if grating == "G140M" else "G395M"})'
            ax.plot(w_truth_overlap, f_truth_overlap, 'navy', lw=1.5, ls='--', alpha=0.7, label=truth_label)
        ax.plot(w_nrs2,  f_corr,  'r', lw=1.5, label=f'Recalibrated {grating} (Parlanti Corr)')
        
        # Plot CALSPEC for general reference
        w_cs_plot = np.linspace(w_nrs2.min(), w_nrs2.max(), 500)
        ax.plot(w_cs_plot, f_ref(w_cs_plot), '0.5', lw=0.8, ls=':', label='CALSPEC model')

        # Y-range scaling (CALSPEC and Truth only, ignore recalibrated outliers)
        f_cs_overlap = f_ref(w_nrs2)
        vals = []
        if len(w_truth_overlap) > 0:
            vals.append(f_truth_overlap)
        vals.append(f_cs_overlap[np.isfinite(f_cs_overlap) & (f_cs_overlap > 0)])
        
        if vals:
            all_vals = np.concatenate(vals)
            ymin = np.nanmin(all_vals) * 0.5
            ymax = np.nanmax(all_vals) * 2.5
            ax.set_ylim(ymin, ymax)

        ax.set_yscale('log')
        ax.set_ylabel('Flux [Jy]', fontsize=12)
        ax.set_xlabel('Wavelength [µm]', fontsize=12)
        ax.set_title(f'{grating} region', fontsize=12)
        ax.grid(True, which='both', alpha=0.1)
        ax.legend(fontsize=9, ncol=2)

    plt.tight_layout()
    outpath = f'{OUTPUT_DIR}/fs_v1_spectra_{name}.png'
    plt.savefig(outpath, dpi=200, bbox_inches='tight')
    plt.close()
    print(f'Saved: {outpath}')

if __name__ == '__main__':
    for sid, src in SOURCES.items():
        plot_source_spectra(sid, src)
