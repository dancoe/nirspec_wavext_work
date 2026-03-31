import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import interp1d

# --- Paths ---
BASE = '/Users/dcoe/NIRSpec/wavext'
DATA_DIR = f'{BASE}/data'
CALSPEC_DIR = f'{BASE}/data/CALSPEC'
RESULTS_DIR = f'{BASE}/results/v7'
PLOT_DIR = f'{BASE}/nirspec_wavext_work/reports/v7/fs_v7/plots'
os.makedirs(PLOT_DIR, exist_ok=True)

C_ANG_S = 2.99792458e18

# --- Sources ---
SOURCES = {
    '1537': {'name': 'G191-B2B',    'calspec': 'g191b2b_mod_012.fits',   'dirs': ['PID1537_G191-B2B']},
    '1538': {'name': 'P330E',       'calspec': 'p330e_mod_008.fits',     'dirs': ['PID1538_P330E']},
    '1536': {'name': 'J1743045',    'calspec': '1743045_mod_007.fits',   'dirs': ['PID1536_J1743045']},
    '6644': {'name': 'NGC2506-G31', 'calspec': 'ngc2506g31_mod_003.fits', 'dirs': ['PID6644_NGC2506G31', 'PID6644_NGC2506-G31']},
}

PID1492_MAP = {
    'g140m': {'truth': f'{DATA_DIR}/PID1492/jw01492-o001_t001-s000000001_nirspec_f170lp-g235m-s200a1-subs200a1_x1d.fits'},
    'g235m': {'truth': f'{DATA_DIR}/PID1492/jw01492-o001_t001-s000000001_nirspec_f290lp-g395m-s200a1-subs200a1_x1d.fits'},
}

def load_spec(path):
    if not os.path.exists(path): return None, None
    with fits.open(path) as h:
        d = h[1].data
        wl, fl = d['WAVELENGTH'].astype(float), d['FLUX'].astype(float)
        dq = d['DQ'] if 'DQ' in d.dtype.names else np.zeros_like(wl)
    mask = (wl > 0) & np.isfinite(fl) & (dq == 0)
    return wl[mask], fl[mask]

def get_truth(pid, grating):
    if pid == '1492':
        return load_spec(PID1492_MAP[grating]['truth'])
    s = SOURCES[pid]
    with fits.open(os.path.join(CALSPEC_DIR, s['calspec'])) as h:
        d = h[1].data
        wl_a, flam = d['WAVELENGTH'].astype(float), d['FLUX'].astype(float)
    fnu = flam * wl_a**2 / C_ANG_S * 1e23
    return wl_a/1e4, fnu

def plot_full_spectra(grating):
    cal_file = os.path.join(RESULTS_DIR, f'calib_v7_fs_{grating}_all.fits')
    if not os.path.exists(cal_file): return
    with fits.open(cal_file) as h:
        c = h[1].data
        wl_grid, k, alpha, beta = c['WAVELENGTH'], c['K'], c['ALPHA'], c['BETA']
    
    # Simple correction function
    def apply_corr(wl, fl):
        f_corr = np.zeros_like(fl)
        # Interpolate coeffs to observed wavelength
        k_i = np.interp(wl, wl_grid, k, left=np.nan, right=np.nan)
        a_i = np.interp(wl, wl_grid, alpha, left=np.nan, right=np.nan)
        b_i = np.interp(wl, wl_grid, beta, left=np.nan, right=np.nan)
        
        # We need f(lam/2) and f(lam/3) to subtract ghosts.
        # For simplicity in this validation plot, we'll just plot corrected vs truth
        # S_obs = k*f1 + a*f2 + b*f3  => f1 = (S_obs - a*f2 - b*f3) / k
        # Since we are plotting "Full Spectrum", we usually plot the Corrected Flux.
        # But wait, to plot the full spectrum correctly, we need the truth model to do the subtraction.
        # A better way for these plots is: plot Observed vs Model(k*f1 + a*f2 + b*f3).
        return k_i, a_i, b_i

    # Let's iterate over sources
    pids = list(SOURCES.keys()) + ['1492']
    for pid in pids:
        name = SOURCES[pid]['name'] if pid in SOURCES else 'PID 1492'
        
        # Find NRS2 Jy files
        roots = SOURCES[pid]['dirs'] if pid in SOURCES else ['PID1492']
        hits = []
        for r in roots:
            files = glob.glob(os.path.join(DATA_DIR, r, 'nrs2_spec2_cal', f'*_x1d.fits'))
            for f in files:
                with fits.open(f) as h:
                    if h[0].header.get('GRATING','').lower() == grating.lower():
                        hits.append(f)
        
        if not hits: continue
        
        wl_obs, fl_obs = [], []
        for h in hits:
            w, f = load_spec(h)
            if w is not None:
                wl_obs.append(w); fl_obs.append(f)
        
        if not wl_obs: continue
        wl_o = np.concatenate(wl_obs); fl_o = np.concatenate(fl_obs)
        idx = np.argsort(wl_o); wl_o = wl_o[idx]; fl_o = fl_o[idx]
        
        # Truth
        wl_t, fl_t = get_truth(pid, grating)
        f_t_interp = interp1d(wl_t, fl_t, bounds_error=False, fill_value=0)
        
        # Model: k*f(l) + a*f(l/2) + b*f(l/3)
        k_i = np.interp(wl_o, wl_grid, k, left=0, right=0)
        a_i = np.interp(wl_o, wl_grid, alpha, left=0, right=0)
        b_i = np.interp(wl_o, wl_grid, beta, left=0, right=0)
        
        model = k_i * f_t_interp(wl_o) + a_i * f_t_interp(wl_o/2) + b_i * f_t_interp(wl_o/3)
        
        plt.figure(figsize=(10,6))
        # New Style: Small dots '.' at alpha=0.6, but the line at only alpha=0.2.
        plt.plot(wl_o, fl_o, '.', markersize=2, color='black', alpha=0.6, label='Observed (NRS2 Jy)')
        plt.plot(wl_o, fl_o, 'k-', alpha=0.2, lw=0.5)
        
        plt.plot(wl_o, model, 'r--', lw=1.5, label='Model (k*f1 + α*f2 + β*f3)')
        plt.plot(wl_o, k_i * f_t_interp(wl_o), 'b-', label='Predicted 1st Order (k*f1)')
        
        # Restore 2nd and 3rd order components
        f2_comp = a_i * f_t_interp(wl_o/2)
        f3_comp = b_i * f_t_interp(wl_o/3)
        plt.plot(wl_o, f2_comp, ',', markersize=1, color='black', alpha=0.5, label='2nd Order (α*f2) λ/2')
        plt.plot(wl_o, f2_comp, color='black', lw=0.5, alpha=0.2)
        plt.plot(wl_o, f3_comp, ',', markersize=1, color='black',   alpha=0.4, label='3rd Order (β*f3) λ/3')
        plt.plot(wl_o, f3_comp, color='black',   lw=0.5, alpha=0.2)
        
        # Scaling logic: Keep it closer to CALSPEC/Model max
        y_max = np.nanmax(model) * 1.5
        y_min = np.nanmin(k_i * f_t_interp(wl_o)) * 0.5
        if np.isfinite(y_min) and np.isfinite(y_max):
            plt.ylim(y_min, y_max)
        
        # Also plot NRS1 if available for context
        nrs1_path = os.path.join(DATA_DIR, roots[0], 'nrs1_spec3_cal', f'*_x1d.fits') if pid != '1492' else PID1492_MAP[grating].get('nrs1','')
        # This is a bit complex to find the exact file, skipping for now to keep it short.
        
        plt.title(f'{name} - {grating.upper()} Full Spectrum Validation (v7)', fontsize=16)
        plt.xlabel('Wavelength (µm)', fontsize=14)
        plt.ylabel('Flux (Jy)', fontsize=14)
        plt.legend(fontsize=12)
        plt.grid(alpha=0.2)
        plt.savefig(os.path.join(PLOT_DIR, f'fs_v7_full_spec_{grating}_{pid}.png'), dpi=150)
        plt.close()

if __name__ == '__main__':
    for g in ['g140m', 'g235m']:
        plot_full_spectra(g)
