"""
Derive NIRSpec Wavelength Extension k/α/β coefficients (v5 Analysis).

This version:
  1. Targets Fixed Slit (FS) datasets (S1600A1).
  2. Includes 4 sources to break the k/α degeneracy:
     - G191-B2B (PID 1537)   — hot WD
     - P330E (PID 1538)      — G2V standard
     - J1743045 (PID 1536)   — A8III solar analog
     - NGC2506-G31 (PID 6644) — G1V cool star (The Degeneracy Breaker)
  3. Uses our Level-3 reduction for both NRS1 (nominal) and NRS2 (extended).
  4. Solves the Parlanti et al. (2025) 3nd-order throughput model:
     S_obs(λ) = k(λ)·S_cal(λ) + α(λ)·S_cal(λ/2) + β(λ)·S_cal(λ/3)

Output:
  results/v5/calib_v5_g140m_f100lp.fits
  results/v5/calib_v5_g235m_f170lp.fits
"""
import os
import glob
import numpy as np
import astropy.io.fits as fits
from scipy.ndimage import uniform_filter1d
from scipy import optimize
from scipy.interpolate import interp1d

# ── Configuration ─────────────────────────────────────────────────────────────
BASE          = '/Users/dcoe/NIRSpec/wavext'
DATA_DIR      = f'{BASE}/data'
CALSPEC_DIR   = f'{BASE}/data/CALSPEC'
PARLANTI_DIR  = f'{BASE}/data/parlanti_repo/calibration_files'
OUTDIR        = f'{BASE}/results/v5'
os.makedirs(OUTDIR, exist_ok=True)

# Map PID → (target display name, CALSPEC filename)
SOURCES = {
    '1537': ('G191-B2B',   'g191b2b_mod_012.fits'),
    '1538': ('P330E',      'p330e_mod_008.fits'),
    '1536': ('J1743045',   '1743045_mod_007.fits'),
    '6644': ('NGC2506-G31', 'ngc2506g31_mod_003.fits'),
}

# Grating configurations for v5 analysis
CONFIGS = {
    'g140m_f100lp': {
        'label': 'G140M/F100LP',
        'slits': ['s1600a1'],
    },
    'g235m_f170lp': {
        'label': 'G235M/F170LP',
        'slits': ['s1600a1'],
    },
}

# Smoothing full-width (channels)
SMOOTH_WIDTH = 40

# ── Data Loading ──────────────────────────────────────────────────────────────

def flam_to_jy(wl_ang, flux_flam):
    """Convert erg/s/cm²/Å → Jy."""
    c_ang_s = 2.99792458e18
    fnu = flux_flam * wl_ang**2 / c_ang_s
    return fnu * 1e23

def load_calspec(pid):
    """Load CALSPEC for spectral type reference."""
    name, fname = SOURCES[pid]
    cpath = os.path.join(CALSPEC_DIR, fname)
    with fits.open(cpath) as hdul:
        d = hdul[1].data
        wl_um = d['WAVELENGTH'].astype(np.float64) / 1e4
        flux_jy = flam_to_jy(d['WAVELENGTH'], d['FLUX'].astype(np.float64))
    return wl_um, flux_jy

def load_l3_x1d(pid, grating, detector):
    """
    Find and load Level-3 x1d.fits from nrs1/nrs2 reduction subdirs.
    detector should be 'nrs1' or 'nrs2'.
    """
    subdir = 'nrs1_spec3_nom' if detector == 'nrs1' else 'nrs2_spec3_ext'
    potential_roots = [pid]
    if pid == '1538': potential_roots.append('PID1538_P330E')
    if pid == '1537': potential_roots.append('PID1537_G191-B2B')
    if pid == '1536': potential_roots.append('PID1536_J1743045')
    if pid == '6644': potential_roots.append('PID6644_NGC2506-G31')
    if pid == '6644': potential_roots.append('PID6644_NGC2506G31')

    for root in potential_roots:
        root_path = os.path.join(DATA_DIR, root)
        if not os.path.exists(root_path): continue
        
        # Get all x1d files in this PID's directory
        all_x1ds = glob.glob(os.path.join(root_path, '**', '*x1d.fits'), recursive=True)
        
        # Filter by grating
        files = [f for f in all_x1ds if grating.lower() in os.path.basename(f).lower()]
        
        # Filter by detector
        if detector == 'nrs2':
            # Looking for nrs2 extension
            files = [f for f in files if 'nrs2' in os.path.basename(f).lower()]
        else:
            # Looking for nrs1 nominal
            # It might be our nrs1_... style or MAST root style (no nrs2)
            files = [f for f in files if 'nrs2' not in os.path.basename(f).lower()]
            # Prefer our spec3 style if both exist
            files = sorted(files, key=lambda x: 'spec3' in x.lower(), reverse=True)
            
        if files:
            with fits.open(files[0]) as hdul:
                d = hdul[1].data
                return d['WAVELENGTH'].copy(), d['FLUX'].copy(), d['DQ'].copy()
    
    return None

def stitch_nrs1_nrs2(pid, grating):
    """Stitch NRS1 and NRS2 detections for a given PID/grating."""
    res1 = load_l3_x1d(pid, grating, 'nrs1')
    res2 = load_l3_x1d(pid, grating, 'nrs2')
    
    if res1 is None or res2 is None:
        return None, None, None
        
    wl1, fl1, dq1 = res1
    wl2, fl2, dq2 = res2
    
    # Simple stitch: NRS1 below split, NRS2 above split
    split = 2.9 if 'g140' in grating.lower() else 3.5
    
    mask1 = wl1 < split
    mask2 = wl2 >= split
    
    wl_all = np.concatenate([wl1[mask1], wl2[mask2]])
    fl_all = np.concatenate([fl1[mask1], fl2[mask2]])
    dq_all = np.concatenate([dq1[mask1], dq2[mask2]])
    
    return wl_all, fl_all, dq_all

# ── Solver ────────────────────────────────────────────────────────────────────

def solve_nnls_coeffs(wl_ref, obs_matrix, cal_matrix):
    """
    Solve S_obs = k*S_λ + α*S_λ/2 + β*S_λ/3 using NNLS per wavelength.
    """
    n_wl = len(wl_ref)
    n_src = obs_matrix.shape[0]
    
    k = np.full(n_wl, np.nan)
    alpha = np.full(n_wl, np.nan)
    beta = np.full(n_wl, np.nan)
    valid = np.zeros(n_wl, dtype=bool)
    
    for i in range(n_wl):
        y = obs_matrix[:, i]
        A = cal_matrix[:, i, :]
        
        # Require all 4 sources ideally, at least 3 for k,alpha,beta
        finite = np.isfinite(y) & np.all(np.isfinite(A), axis=1) & (y > 0)
        if finite.sum() < 3: continue
        
        y_fit = y[finite]
        A_fit = A[finite]
        
        # NNLS ensures non-negative coefficients
        # A_fit is (n_src_finite, 3)
        result, _ = optimize.nnls(A_fit, y_fit)
        
        k[i], alpha[i], beta[i] = result
        valid[i] = True
        
    return k, alpha, beta, valid

def smooth(arr, width, valid):
    """Interpolate across gaps and apply boxcar smoothing."""
    arr_copy = arr.copy()
    x = np.arange(len(arr))
    ok = valid & np.isfinite(arr)
    if ok.sum() < 2: return arr
    f = interp1d(x[ok], arr[ok], kind='linear', bounds_error=False,
                 fill_value=(arr[ok][0], arr[ok][-1]))
    filled = f(x)
    smoothed = uniform_filter1d(filled, size=width, mode='nearest')
    smoothed[~valid] = np.nan
    return smoothed

# ── Main ──────────────────────────────────────────────────────────────────────

def run_v5_solve(key, cfg):
    print(f'\n=== V5 Solve: {cfg["label"]} ===')
    
    # 1. Load data
    wl_ref = None
    all_obs = {}
    pids = list(SOURCES.keys())
    
    for pid in pids:
        wl, fl, dq = stitch_nrs1_nrs2(pid, key.split('_')[0])
        if wl is None:
            print(f'  [!] Missing data for PID {pid}')
            continue
        
        all_obs[pid] = (wl, fl, dq)
        if wl_ref is None: wl_ref = wl # Use first source as grid
        print(f'  PID {pid} Loaded. Range: {wl.min():.2f}-{wl.max():.2f} µm')

    if not all_obs: return
    
    n_wl = len(wl_ref)
    n_pids = len(pids)
    obs_matrix = np.full((n_pids, n_wl), np.nan)
    cal_matrix = np.full((n_pids, n_wl, 3), np.nan)
    
    # 2. Preparation: Interpolate onto reference grid
    for idx, pid in enumerate(pids):
        if pid not in all_obs: continue
        wl_i, fl_i, dq_i = all_obs[pid]
        
        # Interpolate observed
        f_obs = interp1d(wl_i, fl_i, bounds_error=False, fill_value=np.nan)
        f_dq = interp1d(wl_i, dq_i, kind='nearest', bounds_error=False, fill_value=1)
        
        obs_interp = f_obs(wl_ref)
        dq_interp = f_dq(wl_ref)
        obs_interp[dq_interp > 0] = np.nan
        obs_matrix[idx] = obs_interp
        
        # Interpolate CALSPEC truth
        cs_wl, cs_fl = load_calspec(pid)
        f_cs = interp1d(cs_wl, cs_fl, bounds_error=False, fill_value=np.nan)
        
        cal_matrix[idx, :, 0] = f_cs(wl_ref)      # order 1
        cal_matrix[idx, :, 1] = f_cs(wl_ref / 2.0) # order 2
        cal_matrix[idx, :, 2] = f_cs(wl_ref / 3.0) # order 3

    # 3. Solve
    k_raw, a_raw, b_raw, valid = solve_nnls_coeffs(wl_ref, obs_matrix, cal_matrix)
    print(f'  Valid raw solutions: {valid.sum()}/{n_wl}')
    
    # 4. Smooth
    k_sm = smooth(k_raw, SMOOTH_WIDTH, valid)
    a_sm = smooth(a_raw, SMOOTH_WIDTH, valid)
    b_sm = smooth(b_raw, SMOOTH_WIDTH, valid)
    
    # 5. Save
    outname = f'{OUTDIR}/calib_v5_{key}.fits'
    col_wl = fits.Column(name='WAVELENGTH', format='D', array=wl_ref)
    col_k = fits.Column(name='K', format='D', array=k_sm)
    col_a = fits.Column(name='ALPHA', format='D', array=a_sm)
    col_b = fits.Column(name='BETA', format='D', array=b_sm)
    col_kraw = fits.Column(name='K_RAW', format='D', array=k_raw)
    col_araw = fits.Column(name='ALPHA_RAW', format='D', array=a_raw)
    col_braw = fits.Column(name='BETA_RAW', format='D', array=b_raw)
    
    primary = fits.PrimaryHDU()
    table = fits.BinTableHDU.from_columns([col_wl, col_k, col_a, col_b, col_kraw, col_araw, col_braw])
    table.header['GRATING'] = cfg['label']
    table.header['V5_CAL'] = True
    table.header['SOURCES'] = ','.join(pids)
    
    hdul = fits.HDUList([primary, table])
    hdul.writeto(outname, overwrite=True)
    print(f'  → Results saved to {outname}')

if __name__ == '__main__':
    for key, cfg in CONFIGS.items():
        run_v5_solve(key, cfg)
