"""
Coefficient comparison and spectral validation plots for v5 calibration.
Features: 
  - 4-source simultaneous solve results (v5)
  - Comparison against Parlanti (2025) and v4 (3-star) results
  - Spectral residuals for all 4 standards
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import glob
from scipy.interpolate import interp1d
from scipy.ndimage import median_filter

# ── Paths ────────────────────────────────────────────────────────────────────
BASE         = '/Users/dcoe/NIRSpec/wavext'
DATA_DIR     = f'{BASE}/data'
CALSPEC_DIR  = f'{BASE}/data/CALSPEC'
PARLANTI_DIR = f'{BASE}/data/parlanti_repo/calibration_files'
V4_DIR       = f'{BASE}/results/calibration'
V5_DIR       = f'{BASE}/results/v5'
PLOT_DIR     = f'{BASE}/nirspec_wavext_work/reports/v5/fs_v5/plots'
os.makedirs(PLOT_DIR, exist_ok=True)

C_ANG_S = 2.99792458e18  # Å/s

# ── Source configuration ─────────────────────────────────────────────────────
SOURCES = {
    '1537': ('G191-B2B',   'g191b2b_mod_012.fits', 'PID1537_G191-B2B', 'blue'),
    '1538': ('P330E',      'p330e_mod_008.fits',    'PID1538_P330E',    'goldenrod'),
    '1536': ('J1743045',   '1743045_mod_007.fits',  'PID1536_J1743045', 'forestgreen'),
    '6644': ('NGC2506-G31','ngc2506g31_mod_003.fits',  'PID6644_NGC2506G31', 'red'),
}

CONFIGS = {
    'g140m_f100lp': {
        'label': 'G140M/F100LP',
        'short': 'G140M',
        'v5': f'{V5_DIR}/calib_v5_g140m_f100lp.fits',
        'v4': f'{V4_DIR}/calib_g140m_f100lp.fits',
        'par': f'{PARLANTI_DIR}/calibration_functions_g140m_f100lp.fits',
        'wl_range': (0.97, 3.30),
    },
    'g235m_f170lp': {
        'label': 'G235M/F170LP',
        'short': 'G235M',
        'v5': f'{V5_DIR}/calib_v5_g235m_f170lp.fits',
        'v4': f'{V4_DIR}/calib_g235m_f170lp.fits',
        'par': f'{PARLANTI_DIR}/calibration_functions_g235m_f170lp.fits',
        'wl_range': (1.66, 5.50),
    }
}

# ── Utilities ────────────────────────────────────────────────────────────────
def load_calib(path):
    if not os.path.exists(path): return None
    with fits.open(path) as h:
        d = h[1].data
        wl = d['WAVELENGTH'].astype(float)
        k = d['K'].astype(float)
        try:
            a = d['ALPHA'].astype(float)
            b = d['BETA'].astype(float)
        except:
            a = d['alpha'].astype(float)
            b = d['beta'].astype(float)
        return wl, k, a, b

def plot_coefficients():
    print("Plotting coefficients...")
    for key, cfg in CONFIGS.items():
        res_v5 = load_calib(cfg['v5'])
        res_v4 = load_calib(cfg['v4'])
        res_par = load_calib(cfg['par'])
        
        fig, axes = plt.subplots(3, 1, figsize=(10, 12), sharex=True)
        fig.suptitle(f'Coefficient Comparison: {cfg["label"]}\nv5 (4-source) vs v4 (3-source) vs Parlanti', fontsize=14)
        
        labels = [r'$k(\lambda)$', r'$\alpha(\lambda)$', r'$\beta(\lambda)$']
        ylims = [(0.4, 1.2), (-0.01, 0.15), (-0.005, 0.03)]
        
        for i, (label, ylim) in enumerate(zip(labels, ylims)):
            ax = axes[i]
            if res_par: ax.plot(res_par[0], res_par[i+1], color='gray', alpha=0.5, label='Parlanti 2025', lw=2)
            if res_v4:  ax.plot(res_v4[0], res_v4[i+1], color='red', alpha=0.4, label='v4 (3-star)', ls='--')
            if res_v5:  ax.plot(res_v5[0], res_v5[i+1], color='blue', alpha=0.9, label='v5 (4-star)', lw=1.5)
            
            ax.set_ylabel(label, fontsize=12)
            ax.set_ylim(*ylim)
            ax.grid(True, alpha=0.2)
            if i == 0: ax.legend()
        
        axes[2].set_xlabel('Wavelength (µm)', fontsize=12)
        axes[2].set_xlim(*cfg['wl_range'])
        
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        out = f"{PLOT_DIR}/coeffs_{key}.png"
        plt.savefig(out, dpi=150)
        plt.close()
        print(f"  Saved {out}")

# ── Spectral Loading (Recycled from solver) ──────────────────────────────────
def load_l3_x1d(pid, grating, detector):
    potential_roots = [pid]
    if pid == '1538': potential_roots.append('PID1538_P330E')
    if pid == '1537': potential_roots.append('PID1537_G191-B2B')
    if pid == '1536': potential_roots.append('PID1536_J1743045')
    if pid == '6644': potential_roots.append('PID6644_NGC2506-G31')
    if pid == '6644': potential_roots.append('PID6644_NGC2506G31')
    for root in potential_roots:
        root_path = os.path.join(DATA_DIR, root)
        if not os.path.exists(root_path): continue
        all_x1ds = glob.glob(os.path.join(root_path, '**', '*x1d.fits'), recursive=True)
        files = [f for f in all_x1ds if grating.lower() in os.path.basename(f).lower()]
        if detector == 'nrs2':
            files = [f for f in files if 'nrs2' in os.path.basename(f).lower()]
        else:
            files = [f for f in files if 'nrs2' not in os.path.basename(f).lower()]
            files = sorted(files, key=lambda x: 'spec3' in x.lower(), reverse=True)
        if files:
            with fits.open(files[0]) as hdul:
                d = hdul[1].data
                return d['WAVELENGTH'].copy(), d['FLUX'].copy()
    return None, None

def stitch_nrs1_nrs2(pid, grating):
    wl1, fl1 = load_l3_x1d(pid, grating, 'nrs1')
    wl2, fl2 = load_l3_x1d(pid, grating, 'nrs2')
    if wl1 is None or wl2 is None: return None, None
    merged_wl = np.concatenate([wl1, wl2])
    merged_fl = np.concatenate([fl1, fl2])
    order = np.argsort(merged_wl)
    return merged_wl[order], merged_fl[order]

def load_calspec(fname):
    with fits.open(f'{CALSPEC_DIR}/{fname}') as h:
        d = h[1].data
        wl = d['WAVELENGTH'].astype(float) / 1e4
        fl = d['FLUX'].astype(float) * (wl*1e4)**2 / 2.99792458e18 * 1e23 # Jy
        return wl, fl

def apply_v5(wl, obs, k_terp, a_terp, b_terp):
    f_obs_terp = interp1d(wl, obs, bounds_error=False, fill_value=0)
    obs_half = f_obs_terp(wl/2)
    obs_third = f_obs_terp(wl/3)
    return (obs - a_terp * obs_half - b_terp * obs_third ) / k_terp

def plot_residuals():
    print("Plotting v5 residual summary...")
    import glob
    for key, cfg in CONFIGS.items():
        res_v5 = load_calib(cfg['v5'])
        if not res_v5: continue
        wl_c, k_c, a_c, b_c = res_v5
        
        fig, axes = plt.subplots(2, 1, figsize=(12, 10), sharex=True)
        fig.suptitle(f'v5 Spectral Validation — {cfg["label"]}', fontsize=15)
        for ax, title in zip(axes, ['Before ghost correction (S_obs/S_cal)', 'After v5 ghost correction (S_corr/S_cal)']):
            ax.axhline(1.0, color='black', alpha=0.5)
            ax.set_ylabel('Ratio', fontsize=12)
            ax.set_title(title, fontsize=12)
            ax.grid(True, alpha=0.15)
            ax.set_ylim(0.5, 1.5)

        for pid, (star_name, cs_file, dir_name, color) in SOURCES.items():
            wl, obs = stitch_nrs1_nrs2(pid, key.split('_')[0])
            if wl is None: continue
            cs_wl, cs_fl = load_calspec(cs_file)
            cs_terp = interp1d(cs_wl, cs_fl, bounds_error=False, fill_value=np.nan)(wl)
            
            # Calibration interpolants
            k_t = interp1d(wl_c, k_c, bounds_error=False, fill_value=1.0)(wl)
            a_t = interp1d(wl_c, a_c, bounds_error=False, fill_value=0.0)(wl)
            b_t = interp1d(wl_c, b_c, bounds_error=False, fill_value=0.0)(wl)
            
            # Correction
            corr = apply_v5(wl, obs, k_t, a_t, b_t)
            
            # Smooth for plotting
            ratio_raw = median_filter(obs/cs_terp, size=51)
            ratio_corr = median_filter(corr/cs_terp, size=51)
            
            axes[0].plot(wl, ratio_raw, color=color, alpha=0.8, label=f'{star_name}')
            axes[1].plot(wl, ratio_corr, color=color, alpha=0.9, label=f'{star_name}')

        axes[1].set_xlabel('Wavelength (µm)', fontsize=12)
        axes[1].set_xlim(*cfg['wl_range'])
        axes[1].legend(ncol=4, loc='upper center', bbox_to_anchor=(0.5, -0.1))
        
        plt.tight_layout(rect=[0, 0.05, 1, 0.95])
        out = f"{PLOT_DIR}/residuals_{key}.png"
        plt.savefig(out, dpi=150)
        plt.close()
        print(f"  Saved {out}")

def main():
    plot_coefficients()
    plot_residuals()
    print("Done.")

if __name__ == '__main__':
    main()
