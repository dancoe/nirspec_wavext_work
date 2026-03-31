"""
NIRSpec Wavelength Extension — IFU v7 Solver
Apply the Parlanti correction model to IFU x1d spectra.

Key differences vs FS v7:
--------------------------
1. **Input files**: IFU Stage 3 x1d from data/IFU/PID{pid}_{target}/stage3_ext/
   (units: Jy, produced by Spec3Pipeline cube_build + extract_1d)
2. **No photom override needed**: IFU photom_0016 has relresponse=1.0 everywhere;
   flux calibration is fully handled by the pipeline through cube_build.
3. **Truth for calibrators**: CALSPEC model spectra (interpolated to Jy).
4. **Hold-out CV**: Same leave-one-out approach as FS v7.
5. **PID 2186 (UGC5101)**: Science validation target; NRS1 truth = adjacent grating.

Inputs:
-------
  data/IFU/PID{pid}_{target}/stage3_ext/*_{grating}_*x1d.fits  (Jy)
  data/CALSPEC/{model}.fits

Outputs:
--------
  results/v7/calib_v7_ifu_g140m_all.fits
  results/v7/calib_v7_ifu_g140m_holdout_{SOURCE}.fits
  results/v7/calib_v7_ifu_g235m_all.fits
  results/v7/calib_v7_ifu_g235m_holdout_{SOURCE}.fits
"""
import os
import glob
import numpy as np
from astropy.io import fits
from scipy.interpolate import interp1d
from scipy.optimize import nnls

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE        = '/Users/dcoe/NIRSpec/wavext'
IFU_DIR     = f'{BASE}/data/IFU'
DATA_DIR    = f'{BASE}/data'
CALSPEC_DIR = f'{BASE}/data/CALSPEC'
OUTDIR      = f'{BASE}/results/v7'
os.makedirs(OUTDIR, exist_ok=True)

C_ANG_S = 2.99792458e18  # Å/s

# ── IFU Calibration Sources (CALSPEC standards with IFU observations) ─────────
SOURCES = {
    '1537': {'name': 'G191-B2B',  'sptype': 'WD-DA.8', 'calspec': 'g191b2b_mod_012.fits',
             'ifu_dir': 'PID1537_G191-B2B'},
    '1538': {'name': 'P330E',     'sptype': 'G2V',     'calspec': 'p330e_mod_008.fits',
             'ifu_dir': 'PID1538_P330E'},
    '1536': {'name': 'J1743045',  'sptype': 'A8III',   'calspec': '1743045_mod_007.fits',
             'ifu_dir': 'PID1536_J1743045'},
}

# NRS2 extended lower boundary (µm) — where the extended range begins
NRS2_LO = {'g140m': 1.87, 'g235m': 3.15}

# Wavelength grid for solver (extended NRS2 range only)
GRID = {
    'g140m': np.linspace(1.87, 3.20, 400),
    'g235m': np.linspace(3.15, 5.15, 400),
}

# x1d filename pattern within stage3_ext/
X1D_PATTERNS = {
    'g140m': ['*g140m*x1d.fits', '*G140M*x1d.fits', '*f100lp*x1d.fits'],
    'g235m': ['*g235m*x1d.fits', '*G235M*x1d.fits', '*f170lp*x1d.fits'],
}

SMOOTH_BOX = 40


# ── Helpers ────────────────────────────────────────────────────────────────────
def load_spec(path, lo=None, hi=None):
    """Load IFU x1d FITS → (wl_um, flux_jy); applies DQ mask and optional wl clip."""
    if not path or not os.path.exists(path):
        return np.array([]), np.array([])
    with fits.open(path) as h:
        # IFU x1d uses EXTRACT1D extension
        ext = 'EXTRACT1D' if 'EXTRACT1D' in [e.name for e in h] else 1
        d  = h[ext].data
        wl = d['WAVELENGTH'].astype(float)
        fl = d['FLUX'].astype(float)
        dq = (d['DQ'].astype(int) if 'DQ' in d.dtype.names else np.zeros(len(wl), int))
    good = (wl > 0.3) & np.isfinite(fl) & (fl > 0) & (dq == 0)
    if lo is not None:
        good &= (wl >= lo)
    if hi is not None:
        good &= (wl <= hi)
    return wl[good], fl[good]


def calspec_jy_interp(fname):
    """Return µm→Jy interpolator for a CALSPEC model FITS."""
    path = os.path.join(CALSPEC_DIR, fname)
    if not os.path.exists(path):
        print(f'  WARNING: CALSPEC model not found: {fname}')
        return lambda x: np.full_like(x, np.nan)
    with fits.open(path) as h:
        d    = h[1].data
        wl_a = d['WAVELENGTH'].astype(float)
        flam = d['FLUX'].astype(float)
    fnu = flam * wl_a**2 / C_ANG_S * 1e23  # erg/s/cm²/Å → Jy
    idx = np.argsort(wl_a)
    return interp1d(wl_a[idx] / 1e4, fnu[idx], bounds_error=False, fill_value=np.nan)


def make_interp(wl, fl, extrapolate='nearest'):
    """Linear interpolator with nearest-edge fill for extrapolation."""
    if len(wl) < 2:
        return lambda x: np.full_like(np.asarray(x, float), np.nan)
    idx = np.argsort(wl)
    w, f = wl[idx], fl[idx]
    fv = (float(f[0]), float(f[-1])) if extrapolate == 'nearest' else np.nan
    return interp1d(w, f, bounds_error=False, fill_value=fv)


def box_smooth(y, box=40):
    """Box-car smooth; preserve edges."""
    s = np.convolve(y, np.ones(box) / box, mode='same')
    s[:box // 2]  = y[:box // 2]
    s[-box // 2:] = y[-box // 2:]
    return s


def find_ifu_x1d(pid, grating):
    """Find IFU stage3_ext x1d file for a given PID and grating."""
    info = SOURCES.get(pid)
    if info is None:
        return None
    ifu_base = os.path.join(IFU_DIR, info['ifu_dir'], 'stage3_ext')
    for pat in X1D_PATTERNS[grating]:
        hits = sorted(glob.glob(os.path.join(ifu_base, pat)))
        if hits:
            return hits[0]
    return None


# ── NNLS Core ─────────────────────────────────────────────────────────────────
def run_nnls(grid, sources, hold_out=None):
    """
    Solve the 3-parameter Parlanti system at each wavelength grid point.
    Model: S_obs(λ) = k(λ)*f(λ) + α(λ)*f(λ/2) + β(λ)*f(λ/3)
    """
    active = [s for s in sources if s['id'] != hold_out]
    n = len(grid)
    k_raw = np.full(n, np.nan)
    a_raw = np.full(n, np.nan)
    b_raw = np.full(n, np.nan)

    for i, lam in enumerate(grid):
        A_rows, b_vec = [], []
        for s in active:
            s_obs = float(s['f_obs'](lam))
            f1    = float(s['f_ref'](lam))
            f2    = float(s['f_ref'](lam / 2.0))
            f3    = float(s['f_ref'](lam / 3.0))
            if not (np.isfinite(s_obs) and s_obs > 0 and
                    np.isfinite(f1) and f1 > 0 and
                    np.isfinite(f2) and f2 > 0 and
                    np.isfinite(f3) and f3 > 0):
                continue
            A_rows.append([f1, f2, f3])
            b_vec.append(s_obs)

        if len(b_vec) < 2:
            continue
        x, _ = nnls(np.array(A_rows), np.array(b_vec))
        k_raw[i], a_raw[i], b_raw[i] = x

    # Fill NaN gaps by linear interpolation, then smooth
    for arr in (k_raw, a_raw, b_raw):
        bad = ~np.isfinite(arr)
        if bad.any() and (~bad).sum() > 1:
            arr[bad] = np.interp(np.where(bad)[0], np.where(~bad)[0], arr[~bad])

    k = np.clip(box_smooth(np.clip(k_raw, 0, None), SMOOTH_BOX), 0.01, None)
    a = np.clip(box_smooth(np.clip(a_raw, 0, None), SMOOTH_BOX), 0, None)
    b = np.clip(box_smooth(np.clip(b_raw, 0, None), SMOOTH_BOX), 0, None)

    return k, a, b


# ── Per-grating Solve ─────────────────────────────────────────────────────────
def solve_grating(grating):
    print(f'\n{"="*60}')
    print(f'IFU v7 Solve: {grating.upper()}')
    print(f'{"="*60}')
    grid = GRID[grating]
    lo   = NRS2_LO[grating]

    source_list = []

    for pid, s in SOURCES.items():
        x1d_path = find_ifu_x1d(pid, grating)
        if x1d_path is None:
            print(f'  SKIP {s["name"]} (PID {pid}): no x1d found')
            continue

        wl, fl = load_spec(x1d_path, lo=lo)
        if len(wl) < 20:
            print(f'  SKIP {s["name"]} (PID {pid}): insufficient extended-range points')
            continue

        calspec_interp = calspec_jy_interp(s['calspec'])
        # Quick sanity check — does the CALSPEC model cover the extended range?
        test_wl = np.array([lo, lo + 0.5, grid[-1]])
        test_fl = calspec_interp(test_wl)
        if not np.all(np.isfinite(test_fl) & (test_fl > 0)):
            print(f'  SKIP {s["name"]}: CALSPEC model has bad values in extended range')
            continue

        source_list.append({
            'id':    pid,
            'name':  s['name'],
            'f_obs': make_interp(wl, fl),
            'f_ref': calspec_interp,
        })
        print(f'  Added: {s["name"]} (PID {pid}) — {len(wl)} pts, '
              f'wl=[{wl.min():.3f},{wl.max():.3f}]µm, '
              f'med={np.nanmedian(fl):.3e}Jy')
        print(f'    x1d: {os.path.relpath(x1d_path, BASE)}')

    if len(source_list) < 2:
        print(f'  ERROR: only {len(source_list)} source(s) — need ≥2 for NNLS. Skipping {grating}.')
        return

    # Solve with ALL sources
    print(f'\n  Solving with ALL {len(source_list)} sources ...')
    k_all, a_all, b_all = run_nnls(grid, source_list)
    save_result(grating, grid, k_all, a_all, b_all, 'all',
                [s['name'] for s in source_list])

    # Leave-one-out cross-validation
    for s in source_list:
        remaining = [x['name'] for x in source_list if x['id'] != s['id']]
        if len(remaining) < 1:
            continue
        print(f'  Hold-out: {s["name"]} ...')
        k, a, b = run_nnls(grid, source_list, hold_out=s['id'])
        save_result(grating, grid, k, a, b, f'holdout_{s["id"]}', remaining)


def save_result(grating, grid, k, a, b, tag, sources):
    fname = os.path.join(OUTDIR, f'calib_v7_ifu_{grating}_{tag}.fits')
    cols = [
        fits.Column('WAVELENGTH', 'D', array=grid),
        fits.Column('K',          'D', array=k),
        fits.Column('ALPHA',      'D', array=a),
        fits.Column('BETA',       'D', array=b),
    ]
    hdr = fits.Header()
    hdr['VERSION'] = 'v7-ifu'
    hdr['GRATING'] = grating.upper()
    hdr['MODE']    = 'IFU'
    hdr['TAG']     = tag
    hdr['SOURCES'] = ','.join(sources)
    fits.HDUList([fits.PrimaryHDU(header=hdr),
                  fits.BinTableHDU.from_columns(cols)]).writeto(fname, overwrite=True)
    print(f'  Saved: {os.path.relpath(fname, BASE)}')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='IFU v7 Parlanti NNLS solver')
    parser.add_argument('--grating', default='all',
                        help='Grating to solve: g140m, g235m, or all')
    args = parser.parse_args()

    gratings = ['g140m', 'g235m'] if args.grating == 'all' else [args.grating.lower()]
    for g in gratings:
        solve_grating(g)

    print('\nIFU v7 solving complete.')
