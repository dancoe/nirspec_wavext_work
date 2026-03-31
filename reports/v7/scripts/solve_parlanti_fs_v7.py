"""
NIRSpec Wavelength Extension — FS v7 Solver
Transitioning to Full Flux Calibration (Jy)

Key upgrades over FS v6:
-----------------------
1.  **Full Flux Calibration (Jy)**: Reads natively flux-calibrated Jy spectra from 
    the Spec2 nrs2_spec2_cal/ directories, produced with the custom extended-photom
    reference file. Manual DN/s→Jy scaling is removed.
2.  **Leave-One-Out Cross-Validation**: Systematically holds out one source at a 
    time to assess coefficient stability and per-source performance.

Inputs:
-------
  - *_x1d.fits from data/PIDXXXX/nrs2_spec2_cal/ (Jy units)

Outputs:
--------
  - results/v7/calib_v7_fs_g140m_f100lp_all.fits
  - results/v7/calib_v7_fs_g140m_f100lp_holdout_{SOURCE}.fits
  - results/v7/calib_v7_fs_g235m_f170lp_all.fits
  - results/v7/calib_v7_fs_g235m_f170lp_holdout_{SOURCE}.fits
"""
import os
import glob
import numpy as np
from astropy.io import fits
from scipy.interpolate import interp1d
from scipy.optimize import nnls

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE        = '/Users/dcoe/NIRSpec/wavext'
DATA_DIR    = f'{BASE}/data'
CALSPEC_DIR = f'{BASE}/data/CALSPEC'
OUTDIR      = f'{BASE}/results/v7'
os.makedirs(OUTDIR, exist_ok=True)

C_ANG_S = 2.99792458e18  # Å/s

# ── Sources ───────────────────────────────────────────────────────────────────
SOURCES = {
    '1537': {'name': 'G191-B2B',    'sptype': 'WD-DA.8', 'calspec': 'g191b2b_mod_012.fits',   'dirs': ['PID1537_G191-B2B']},
    '1538': {'name': 'P330E',       'sptype': 'G2V',     'calspec': 'p330e_mod_008.fits',     'dirs': ['PID1538_P330E']},
    '1536': {'name': 'J1743045',    'sptype': 'A8III',   'calspec': '1743045_mod_007.fits',   'dirs': ['PID1536_J1743045']},
    '6644': {'name': 'NGC2506-G31', 'sptype': 'G1V',     'calspec': 'ngc2506g31_mod_003.fits', 'dirs': ['PID6644_NGC2506G31', 'PID6644_NGC2506-G31']},
}

# PID 1492 (Science Target) — Truth comes from cross-grating STITCHING
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

# NRS2 extended lower boundary (µm)
NRS2_LO = {'g140m': 1.87, 'g235m': 3.15}

# Wavelength grid for solver
GRID = {
    'g140m': np.linspace(1.87, 3.20, 400),
    'g235m': np.linspace(3.15, 5.15, 400),
}

SMOOTH_BOX = 40


# ── Helpers ────────────────────────────────────────────────────────────────────
def load_spec(path, label=''):
    """Load x1d FITS → (wl, flux); bad-pixel mask applied."""
    if not os.path.exists(path):
        return np.array([]), np.array([])
    with fits.open(path) as h:
        d  = h[1].data
        wl = d['WAVELENGTH'].astype(float)
        fl = d['FLUX'].astype(float)
        dq = (d['DQ'].astype(int) if 'DQ' in d.dtype.names else np.zeros(len(wl), int))
    good = (wl > 0.3) & np.isfinite(fl) & (fl > 0) & (dq == 0)
    return wl[good], fl[good]


def calspec_jy_interp(fname):
    """Return µm→Jy interpolator for a CALSPEC model FITS."""
    with fits.open(os.path.join(CALSPEC_DIR, fname)) as h:
        d    = h[1].data
        wl_a = d['WAVELENGTH'].astype(float)
        flam = d['FLUX'].astype(float)
    fnu = flam * wl_a**2 / C_ANG_S * 1e23  # Jy
    idx = np.argsort(wl_a)
    return interp1d(wl_a[idx] / 1e4, fnu[idx], bounds_error=False, fill_value=np.nan)


def make_interp(wl, fl, extrapolate='nearest'):
    if len(wl) < 2: return lambda x: np.full_like(x, np.nan)
    idx = np.argsort(wl)
    w, f = wl[idx], fl[idx]
    fv = (float(f[0]), float(f[-1])) if extrapolate == 'nearest' else np.nan
    return interp1d(w, f, bounds_error=False, fill_value=fv)


def box_smooth(y, box=40):
    s = np.convolve(y, np.ones(box) / box, mode='same')
    s[:box // 2]  = y[:box // 2]
    s[-box // 2:] = y[-box // 2:]
    return s


def find_v7_jy_x1ds(pid, grating):
    """Find all Level-3 NRS2 Jy-calibrated x1d produced by v7-Spec2 for specific grating."""
    roots = SOURCES.get(pid, {}).get('dirs', [f'PID{pid}'])
    grating = grating.upper()
    hits = []
    for root in roots:
        path = os.path.join(DATA_DIR, root, 'nrs2_spec2_cal', '*_x1d.fits')
        files = glob.glob(path)
        for f in files:
            with fits.open(f) as h:
                if h[0].header.get('GRATING', '').upper() == grating:
                    hits.append(f)
    return hits


# ── NNLS Core ─────────────────────────────────────────────────────────────────
def run_nnls(grid, sources, hold_out=None):
    """Solve the 3-column system across the grid using provided sources."""
    active = [s for s in sources if s['id'] != hold_out]
    n = len(grid)
    k_raw, a_raw, b_raw = np.full(n, np.nan), np.full(n, np.nan), np.full(n, np.nan)
    
    for i, lam in enumerate(grid):
        A, b = [], []
        for s in active:
            s_obs = float(s['f_obs'](lam))
            f1    = float(s['f_ref'](lam))
            f2    = float(s['f_ref'](lam/2.0))
            f3    = float(s['f_ref'](lam/3.0))
            if not (np.isfinite(s_obs) and s_obs > 0 and 
                    np.isfinite(f1) and f1 > 0 and 
                    np.isfinite(f2) and f2 > 0 and 
                    np.isfinite(f3) and f3 > 0):
                continue
            A.append([f1, f2, f3])
            b.append(s_obs)
        
        if len(b) < 2: continue
        x, _ = nnls(np.array(A), np.array(b))
        k_raw[i], a_raw[i], b_raw[i] = x

    # Interp and Smooth
    for arr in (k_raw, a_raw, b_raw):
        bad = ~np.isfinite(arr)
        if bad.any() and (~bad).sum() > 1:
            arr[bad] = np.interp(np.where(bad)[0], np.where(~bad)[0], arr[~bad])
    
    k = np.clip(box_smooth(np.clip(k_raw, 0, None), SMOOTH_BOX), 0.01, None)
    a = np.clip(box_smooth(np.clip(a_raw, 0, None), SMOOTH_BOX), 0, None)
    b = np.clip(box_smooth(np.clip(b_raw, 0, None), SMOOTH_BOX), 0, None)
    
    return k, a, b


def solve_grating(grating):
    print(f'\n=== FS v7 Solve: {grating.upper()} ===')
    grid = GRID[grating]
    lo   = NRS2_LO[grating]
    
    source_list = []
    
    # 1. Standard Stars (Jy)
    for pid, s in SOURCES.items():
        hits = find_v7_jy_x1ds(pid, grating)
        if not hits:
            print(f'  SKIP {s["name"]}: no v7 Jy file found')
            continue
        
        # Load and stack all matching exposures for this source
        all_wl, all_fl = [], []
        for path in hits:
            wl, fl = load_spec(path)
            mask = wl >= lo
            if mask.sum() < 20: continue
            all_wl.append(wl[mask])
            all_fl.append(fl[mask])
        
        if not all_wl: continue
        
        # Median combine (or concatenate)
        wl_cat = np.concatenate(all_wl)
        fl_cat = np.concatenate(all_fl)
        
        source_list.append({
            'id': pid, 'name': s['name'],
            'f_obs': make_interp(wl_cat, fl_cat),
            'f_ref': calspec_jy_interp(s['calspec'])
        })
        print(f'  Added: {s["name"]} ({len(hits)} exposures stack)')

    # 2. PID 1492 (Stitched Jy)
    hits_1492 = find_v7_jy_x1ds('1492', grating)
    if hits_1492:
        all_wl, all_fl = [], []
        for path in hits_1492:
            wl, fl = load_spec(path)
            mask = wl >= lo
            if mask.sum() > 50:
                all_wl.append(wl[mask])
                all_fl.append(fl[mask])
        
        if all_wl:
            wl_cat, fl_cat = np.concatenate(all_wl), np.concatenate(all_fl)
            map_1492 = PID1492_MAP[grating]
            wl_tr, fl_tr = load_spec(map_1492['truth'])     # next grating nominal
            wl_n1, fl_n1 = load_spec(map_1492['nrs1'])      # current grating nominal (nrs1)
            
            stitch = wl_tr.min()
            m1, mt = wl_n1 < stitch, wl_tr >= stitch
            wl_comb = np.concatenate([wl_n1[m1], wl_tr[mt]])
            fl_comb = np.concatenate([fl_n1[m1], fl_tr[mt]])
            
            source_list.append({
                'id': '1492', 'name': 'PID1492',
                'f_obs': make_interp(wl_cat, fl_cat),
                'f_ref': make_interp(wl_comb, fl_comb)
            })
            print(f'  Added: PID1492 (Stitched Jy, {len(hits_1492)} stacked)')

    if len(source_list) < 2:
        print(f'  ERROR: insufficient sources for {grating}')
        return

    # Solve ALL
    print(f'  Solving with ALL {len(source_list)} sources...')
    k_all, a_all, b_all = run_nnls(grid, source_list)
    save_result(grating, grid, k_all, a_all, b_all, 'all', [s['name'] for s in source_list])
    
    # Solve Hold-outs
    for s in source_list:
        print(f'  Solving hold-out: {s["name"]}...')
        k, a, b = run_nnls(grid, source_list, hold_out=s['id'])
        save_result(grating, grid, k, a, b, f'holdout_{s["id"]}', [x['name'] for x in source_list if x['id'] != s['id']])

def save_result(grating, grid, k, a, b, tag, sources):
    fname = os.path.join(OUTDIR, f'calib_v7_fs_{grating}_{tag}.fits')
    cols = [
        fits.Column('WAVELENGTH', 'D', array=grid),
        fits.Column('K', 'D', array=k),
        fits.Column('ALPHA', 'D', array=a),
        fits.Column('BETA', 'D', array=b)
    ]
    hdr = fits.Header()
    hdr['VERSION'] = 'v7'
    hdr['GRATING'] = grating.upper()
    hdr['TAG'] = tag
    hdr['SOURCES'] = ','.join(sources)
    fits.HDUList([fits.PrimaryHDU(header=hdr), fits.BinTableHDU.from_columns(cols)]).writeto(fname, overwrite=True)

if __name__ == '__main__':
    for g in ('g140m', 'g235m'):
        solve_grating(g)
    print('\nFS v7 solving complete.')
