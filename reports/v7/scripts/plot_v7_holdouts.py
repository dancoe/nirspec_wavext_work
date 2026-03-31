"""
NIRSpec FS v7 — Hold-out Validation Plotting
Compare hold-out source (observed) vs ghost model prediction (from other sources).
"""
import os
import glob
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

BASE        = '/Users/dcoe/NIRSpec/wavext'
DATA_DIR    = f'{BASE}/data'
CALSPEC_DIR = f'{BASE}/data/CALSPEC'
V7_DIR      = f'{BASE}/results/v7'
PLOT_DIR    = f'{BASE}/nirspec_wavext_work/reports/v7/fs_v7/plots'
os.makedirs(PLOT_DIR, exist_ok=True)

C_ANG_S = 2.99792458e18

SOURCES = {
    '1537': {'name': 'G191-B2B',    'calspec': 'g191b2b_mod_012.fits',   'dirs': ['PID1537_G191-B2B']},
    '1538': {'name': 'P330E',       'calspec': 'p330e_mod_008.fits',     'dirs': ['PID1538_P330E']},
    '1536': {'name': 'J1743045',    'calspec': '1743045_mod_007.fits',   'dirs': ['PID1536_J1743045']},
    '6644': {'name': 'NGC2506-G31', 'calspec': 'ngc2506g31_mod_003.fits', 'dirs': ['PID6644_NGC2506-G31', 'PID6644_NGC2506G31']},
}

PID1492_MAP = {
    'g140m': {
        'truth':     f'{DATA_DIR}/PID1492/jw01492-o001_t001-s000000001_nirspec_f170lp-g235m-s200a1-subs200a1_x1d.fits',
        'nrs1':      f'{DATA_DIR}/PID1492/jw01492-o001_t001-s000000001_nirspec_f100lp-g140m-s200a1-subs200a1_x1d.fits',
    },
    'g235m': {
        'truth':     f'{DATA_DIR}/PID1492/jw01492-o001_t001-s000000001_nirspec_f290lp-g395m-s200a1-subs200a1_x1d.fits',
        'nrs1':      f'{DATA_DIR}/PID1492/jw01492-o001_t001-s000000001_nirspec_f170lp-g235m-s200a1-subs200a1_x1d.fits',
    },
}

def load_spec(path):
    if not os.path.exists(path): return np.array([]), np.array([])
    with fits.open(path) as h:
        d = h[1].data
        wl, fl = d['WAVELENGTH'].astype(float), d['FLUX'].astype(float)
    return wl, fl

def calspec_jy_interp(fname):
    with fits.open(os.path.join(CALSPEC_DIR, fname)) as h:
        d = h[1].data
        wl_a, flam = d['WAVELENGTH'].astype(float), d['FLUX'].astype(float)
    fnu = flam * wl_a**2 / C_ANG_S * 1e23
    return interp1d(wl_a / 1e4, fnu, bounds_error=False, fill_value=np.nan)

def find_v7_jy_x1ds(pid, grating):
    roots = SOURCES.get(pid, {}).get('dirs', [f'PID{pid}'])
    hits = []
    for root in roots:
        path = os.path.join(DATA_DIR, root, 'nrs2_spec2_cal', '*_x1d.fits')
        for f in glob.glob(path):
            with fits.open(f) as h:
                if h[0].header.get('GRATING', '').upper() == grating.upper():
                    hits.append(f)
    return hits

def plot_holdout(grating):
    print(f'\nPlotting Hold-outs for {grating.upper()}')
    sources = list(SOURCES.keys()) + ['1492']
    
    fig, axes = plt.subplots(len(sources), 1, figsize=(10, 3*len(sources)), sharex=True)
    if len(sources) == 1: axes = [axes]

    for ax, pid in zip(axes, sources):
        sname = SOURCES.get(pid, {}).get('name', 'PID1492')
        cal_path = os.path.join(V7_DIR, f'calib_v7_fs_{grating}_holdout_{pid}.fits')
        if not os.path.exists(cal_path):
            print(f'  SKIP {sname}: hold-out cal file missing')
            continue

        with fits.open(cal_path) as ch:
            cwl = ch[1].data['WAVELENGTH']
            ck  = ch[1].data['K']
            ca  = ch[1].data['ALPHA']
            cb  = ch[1].data['BETA']
        
        # Reference flux for prediction
        if pid == '1492':
            m = PID1492_MAP[grating]
            w1, f1 = load_spec(m['nrs1'])
            wt, ft = load_spec(m['truth'])
            stitch = wt.min()
            comb_w = np.concatenate([w1[w1 < stitch], wt[wt >= stitch]])
            comb_f = np.concatenate([f1[w1 < stitch], ft[wt >= stitch]])
            f_ref = interp1d(comb_w, comb_f, bounds_error=False, fill_value=np.nan)
        else:
            f_ref = calspec_jy_interp(SOURCES[pid]['calspec'])

        # Predict observed: f_pred = k*f(l) + a*f(l/2) + b*f(l/3)
        f_pred = ck*f_ref(cwl) + ca*f_ref(cwl/2.0) + cb*f_ref(cwl/3.0)
        
        # Load OBSERVED (v7 Jy)
        hits = find_v7_jy_x1ds(pid, grating)
        if not hits: continue
        obs_wl, obs_fl = [], []
        for path in hits:
            w, f = load_spec(path)
            obs_wl.append(w); obs_fl.append(f)
        owl = np.concatenate(obs_wl)
        ofl = np.concatenate(obs_fl)
        
        # Plot
        ax.plot(owl, ofl, 'k.', alpha=0.1, label='Observed (NRS2 Jy)')
        ax.plot(cwl, f_pred, 'r-', lw=2, label='Predicted (Hold-out Ghost Model)')
        ax.plot(cwl, ck*f_ref(cwl), 'b--', alpha=0.5, label='Predicted (Ghost-free)')
        
        ax.set_ylabel('Flux (Jy)')
        ax.set_title(f'Hold-out Validation: {sname} ({grating.upper()})')
        if pid == sources[0]: ax.legend()
        
    axes[-1].set_xlabel('Wavelength (µm)')
    plt.tight_layout()
    out = os.path.join(PLOT_DIR, f'fs_v7_holdout_validation_{grating}.png')
    plt.savefig(out)
    print(f'Saved: {out}')

if __name__ == '__main__':
    plt.switch_backend('Agg')
    for g in ('g140m', 'g235m'):
        plot_holdout(g)
