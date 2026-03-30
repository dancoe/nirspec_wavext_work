"""
NIRSpec Wavelength Extension — FS v6 Solver

Key upgrade over FS v5
-----------------------
v5 used 4 CALSPEC standard stars (P330E, G191-B2B, J1743045, NGC2506-G31).
v6 adds PID 1492 as a 5th source via cross-grating truth: the MAST-pipeline
(nominal) x1d from the ADJACENT grating serves as the truth spectrum for
the NRS2 extended region.

G140M:  4 CALSPEC stds + PID 1492 (cross-grating: G235M MAST as truth)
G235M:  4 CALSPEC stds + PID 1492 (cross-grating: G395M MAST as truth)

PID 1492 data details
----------------------
  G140M NRS2 ext  :  data/PID1492/jw01492003001_03102_00005_nrs2_g140m_extract_1d.fits
                     Units: DN/s  → normalised to pseudo-Jy via overlap scale
  G235M MAST truth:  data/PID1492/jw01492-o001_t001-..._nirspec_f170lp-g235m-..._x1d.fits
                     Units: Jy   (nominal pipeline, covers 1.70–3.07 µm)

  G235M NRS2 ext  :  data/PID1492/jw01492003001_03104_00004_nrs2_g235m_extract_1d.fits
                     Units: DN/s  → normalised to pseudo-Jy via overlap scale
  G395M MAST truth:  data/PID1492/jw01492-o001_t001-..._nirspec_f290lp-g395m-..._x1d.fits
                     Units: Jy   (nominal pipeline, covers 2.88–5.09 µm)

Normalisation:  scale = median(MAST_Jy / NRS2_DNs)  over the wavelength
overlap between the NRS2 extended spectrum and the MAST truth. This converts
the DN/s spectrum to the same Jy scale, making it compatible with the
simultaneous NNLS solve.

Outputs
--------
  results/v6/calib_v6_fs_g140m_f100lp.fits
  results/v6/calib_v6_fs_g235m_f170lp.fits
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
OUTDIR      = f'{BASE}/results/v6'
os.makedirs(OUTDIR, exist_ok=True)

C_ANG_S = 2.99792458e18  # Å/s

# ── CALSPEC standard star sources (same as v5) ────────────────────────────────
CALSPEC_SOURCES = {
    '1537': ('G191-B2B',    'g191b2b_mod_012.fits'),
    '1538': ('P330E',       'p330e_mod_008.fits'),
    '1536': ('J1743045',    '1743045_mod_007.fits'),
    '6644': ('NGC2506-G31', 'ngc2506g31_mod_003.fits'),
}

# Root directories to search per PID
PID_DIRS = {
    '1537': ['PID1537_G191-B2B'],
    '1538': ['PID1538_P330E'],
    '1536': ['PID1536_J1743045'],
    '6644': ['PID6644_NGC2506G31', 'PID6644_NGC2506-G31'],
}

# NRS2 extended lower boundary (µm)
NRS2_LO = {'g140m': 1.87, 'g235m': 3.15}

# Wavelength grid for solver
GRID = {
    'g140m': np.linspace(1.87, 3.20, 400),   # FS NRS2 extent
    'g235m': np.linspace(3.15, 5.15, 400),
}

SMOOTH_BOX = 40

# ── PID 1492 file map ─────────────────────────────────────────────────────────
PID1492_DIR = f'{DATA_DIR}/PID1492'

PID1492_FILES = {
    'g140m': {
        'nrs2_ext':  f'{PID1492_DIR}/jw01492003001_03102_00005_nrs2_g140m_extract_1d.fits',
        'truth':     f'{PID1492_DIR}/jw01492-o001_t001-s000000001_nirspec_f170lp-g235m-s200a1-subs200a1_x1d.fits',
        'nrs1_truth': f'{PID1492_DIR}/jw01492-o001_t001-s000000001_nirspec_f100lp-g140m-s200a1-subs200a1_x1d.fits',
    },
    'g235m': {
        'nrs2_ext':  f'{PID1492_DIR}/jw01492003001_03104_00004_nrs2_g235m_extract_1d.fits',
        'truth':     f'{PID1492_DIR}/jw01492-o001_t001-s000000001_nirspec_f290lp-g395m-s200a1-subs200a1_x1d.fits',
        'nrs1_truth': f'{PID1492_DIR}/jw01492-o001_t001-s000000001_nirspec_f170lp-g235m-s200a1-subs200a1_x1d.fits',
    },
}


# ── Helpers ────────────────────────────────────────────────────────────────────
def load_spec(path, label=''):
    """Load x1d FITS → (wl, flux); bad-pixel mask applied."""
    if not os.path.exists(path):
        print(f'  MISSING: {os.path.basename(path)}  [{label}]')
        return np.array([]), np.array([])
    with fits.open(path) as h:
        d  = h[1].data
        wl = d['WAVELENGTH'].astype(float)
        fl = d['FLUX'].astype(float)
        dq = (d['DQ'].astype(int) if 'DQ' in d.dtype.names
              else np.zeros(len(wl), int))
    good = (wl > 0.3) & np.isfinite(fl) & (fl > 0) & (dq == 0)
    return wl[good], fl[good]


def calspec_jy_interp(fname):
    """Return µm→Jy interpolator for a CALSPEC model FITS."""
    with fits.open(os.path.join(CALSPEC_DIR, fname)) as h:
        d    = h[1].data
        wl_a = d['WAVELENGTH'].astype(float)
        flam = d['FLUX'].astype(float)
    fnu = flam * wl_a**2 / C_ANG_S * 1e23
    idx = np.argsort(wl_a)
    return interp1d(wl_a[idx] / 1e4, fnu[idx],
                    bounds_error=False, fill_value=np.nan)


def make_interp(wl, fl, extrapolate='nearest'):
    idx = np.argsort(wl)
    w, f = wl[idx], fl[idx]
    fv = (float(f[0]), float(f[-1])) if extrapolate == 'nearest' else np.nan
    return interp1d(w, f, bounds_error=False, fill_value=fv)


def box_smooth(y, box=40):
    s = np.convolve(y, np.ones(box) / box, mode='same')
    s[:box // 2]  = y[:box // 2]
    s[-box // 2:] = y[-box // 2:]
    return s


def find_nrs2_ext(pid, grating_name):
    """Find Level-3 NRS2 extended x1d for a standard-star PID."""
    for root in PID_DIRS.get(pid, [pid]):
        root_path = os.path.join(DATA_DIR, root)
        if not os.path.isdir(root_path):
            continue
        pattern = os.path.join(root_path, '**', f'*{grating_name}*x1d.fits')
        hits = [f for f in glob.glob(pattern, recursive=True)
                if 'nrs2' in os.path.basename(f).lower()]
        if hits:
            return hits[0]
    return None


def find_nrs1_nom(pid, grating_name):
    """Find Level-3 NRS1 nominal x1d for a standard-star PID."""
    for root in PID_DIRS.get(pid, [pid]):
        root_path = os.path.join(DATA_DIR, root)
        if not os.path.isdir(root_path):
            continue
        # First try high-level MAST x1d
        mast = glob.glob(os.path.join(root_path, f'*{grating_name}*x1d.fits'))
        mast = [f for f in mast if 'nrs2' not in os.path.basename(f).lower()]
        if mast:
            return mast[0]
        # Then nrs1_spec3_nom style
        pattern = os.path.join(root_path, '**', f'*{grating_name}*x1d.fits')
        hits = [f for f in glob.glob(pattern, recursive=True)
                if 'nrs2' not in os.path.basename(f).lower()]
        if hits:
            return sorted(hits, key=lambda x: 'spec3' in x)[-1]
    return None


# ── Build PID 1492 cross-grating source ───────────────────────────────────────
def build_pid1492_source(grating):
    """
    Build an active-source dict for PID 1492 for the given grating.

    Returns dict with 'name', 'f_obs', 'f_ref' or None on failure.

    Strategy
    ---------
    1. Load NRS2 extended (DN/s) and MAST truth (Jy) spectra.
    2. Compute scale = median(MAST_Jy / NRS2_DNs) over wavelength overlap,
       which normalises the NRS2 spectrum to pseudo-Jy.
    3. Build a stitched truth spectrum:
         λ  < stitch : NRS1 MAST (lower grating nominal, Jy)
         λ >= stitch : MAST truth (next grating nominal, Jy)
       for use in direct and ghost (λ/2, λ/3) columns of A.
    """
    fmap = PID1492_FILES[grating]
    nrs2_lo = NRS2_LO[grating]

    # Load NRS2 extended (DN/s)
    wl_ext, fl_ext = load_spec(fmap['nrs2_ext'], f'PID1492 {grating.upper()} NRS2 ext')
    if len(wl_ext) < 10:
        return None

    # Restrict to extended region only
    mask = wl_ext >= nrs2_lo
    wl_nrs, fl_nrs = wl_ext[mask], fl_ext[mask]
    if len(wl_nrs) < 10:
        return None

    # Load next-grating MAST truth (Jy)
    wl_truth, fl_truth = load_spec(fmap['truth'], f'PID1492 {grating.upper()} MAST truth')
    if len(wl_truth) < 10:
        return None

    # Overlap between NRS2 ext and MAST truth
    lo = max(wl_nrs.min(), wl_truth.min())
    hi = min(wl_nrs.max(), wl_truth.max())
    if hi <= lo:
        print(f'  PID1492 {grating.upper()}: no wavelength overlap → skip')
        return None

    # Compute normalisation scale (DN/s → pseudo-Jy)
    wl_grid = np.linspace(lo, hi, 50)
    truth_interp = make_interp(wl_truth, fl_truth)
    nrs2_interp  = make_interp(wl_nrs, fl_nrs)
    ratio = truth_interp(wl_grid) / nrs2_interp(wl_grid)
    good  = np.isfinite(ratio) & (ratio > 0)
    if good.sum() < 5:
        print(f'  PID1492 {grating.upper()}: insufficient overlap for scale → skip')
        return None
    scale = float(np.nanmedian(ratio[good]))
    print(f'  PID1492 {grating.upper()}: overlap {lo:.2f}–{hi:.2f} µm, '
          f'DN/s→Jy scale = {scale:.4e}')

    # Normalised NRS2 observed spectrum
    fl_nrs_jy = fl_nrs * scale
    f_obs = make_interp(wl_nrs, fl_nrs_jy)

    # Build stitched truth spectrum for ghost terms (λ/2, λ/3)
    # Use NRS1 MAST (lower grating nominal) for short wavelengths
    wl_nrs1, fl_nrs1 = load_spec(fmap['nrs1_truth'], f'PID1492 {grating.upper()} NRS1')
    stitch = wl_truth.min()  # stitch at start of next-grating coverage
    if len(wl_nrs1) > 10:
        # Combine NRS1 (short) + next-grating MAST (long) into single truth
        lo_part = wl_nrs1 < stitch
        hi_part = wl_truth >= stitch
        wl_comb = np.concatenate([wl_nrs1[lo_part], wl_truth[hi_part]])
        fl_comb = np.concatenate([fl_nrs1[lo_part], fl_truth[hi_part]])
    else:
        wl_comb, fl_comb = wl_truth, fl_truth

    idx = np.argsort(wl_comb)
    f_ref = make_interp(wl_comb[idx], fl_comb[idx], extrapolate='nearest')

    return {
        'name':   f'PID1492 [cross-grating]',
        'sptype': 'IRAS-05248',
        'f_obs':  f_obs,
        'f_ref':  f_ref,
    }


# ── Core NNLS solver ──────────────────────────────────────────────────────────
def solve_grating(grating, extra_sources=None):
    print(f'\n{"="*60}')
    print(f'FS v6 NNLS — {grating.upper()}')
    print(f'{"="*60}')

    wav_grid = GRID[grating]
    nrs2_lo  = NRS2_LO[grating]
    active   = []

    # 1. CALSPEC standard stars
    for pid, (name, calspec_fname) in CALSPEC_SOURCES.items():
        path_nrs2 = find_nrs2_ext(pid, grating)
        if path_nrs2 is None:
            print(f'  SKIP {name}: no NRS2 ext x1d found')
            continue
        wl2, fl2 = load_spec(path_nrs2, f'{name} NRS2 ext')
        mask = (wl2 >= nrs2_lo)
        wl_nrs, fl_nrs = wl2[mask], fl2[mask]
        if len(wl_nrs) < 10:
            print(f'  SKIP {name}: <10 NRS2 pixels')
            continue

        f_cs  = calspec_jy_interp(calspec_fname)
        r_mid = (float(np.interp(np.median(wl_nrs), wl_nrs, fl_nrs))
                 / float(f_cs(np.median(wl_nrs))))
        if not (r_mid > 0.01):
            print(f'  SKIP {name}: R_mid={r_mid:.4f} implausible')
            continue

        active.append({
            'name':   name,
            'sptype': '—',
            'f_obs':  make_interp(wl_nrs, fl_nrs),
            'f_ref':  f_cs,
        })
        print(f'  {name} (PID {pid}): NRS2 '
              f'{wl_nrs.min():.3f}–{wl_nrs.max():.3f} µm  R_mid={r_mid:.3f}')

    # 2. Cross-grating science source(s)
    if extra_sources:
        for es in extra_sources:
            if es is not None:
                active.append(es)
                print(f'  {es["name"]} [cross-grating]: added')

    if len(active) < 2:
        raise RuntimeError(f'Only {len(active)} source(s) for {grating}; need ≥2')
    print(f'  Total sources: {len(active)}')

    # NNLS
    n = len(wav_grid)
    k_raw = np.full(n, np.nan)
    a_raw = np.full(n, np.nan)
    b_raw = np.full(n, np.nan)
    cond  = np.full(n, np.nan)

    for j, lam in enumerate(wav_grid):
        rows_A, rows_b = [], []
        for d in active:
            s_obs = float(d['f_obs'](lam))
            f1    = float(d['f_ref'](lam))
            f2    = float(d['f_ref'](lam / 2.0))
            f3    = float(d['f_ref'](lam / 3.0))
            if not (np.isfinite(s_obs) and s_obs > 0
                    and np.isfinite(f1) and f1 > 0
                    and np.isfinite(f2) and f2 > 0
                    and np.isfinite(f3) and f3 > 0):
                continue
            rows_A.append([f1, f2, f3])
            rows_b.append(s_obs)

        if len(rows_b) < 2:
            continue

        A  = np.array(rows_A, dtype=float)
        bv = np.array(rows_b, dtype=float)
        try:
            sv = np.linalg.svd(A, compute_uv=False)
            cond[j] = sv[0] / sv[-1] if sv[-1] > 0 else np.inf
        except np.linalg.LinAlgError:
            pass
        x, _ = nnls(A, bv)
        k_raw[j], a_raw[j], b_raw[j] = x

    for arr in (k_raw, a_raw, b_raw):
        bad = ~np.isfinite(arr)
        if bad.any():
            xi = np.where(~bad)[0]
            if len(xi) > 1:
                arr[bad] = np.interp(np.where(bad)[0], xi, arr[xi])

    k_raw = np.clip(k_raw, 0.01, None)
    a_raw = np.clip(a_raw, 0.0,  None)
    b_raw = np.clip(b_raw, 0.0,  None)

    ks  = np.clip(box_smooth(k_raw, SMOOTH_BOX), 0.01, None)
    als = np.clip(box_smooth(a_raw, SMOOTH_BOX), 0.0,  None)
    bts = np.clip(box_smooth(b_raw, SMOOTH_BOX), 0.0,  None)

    med_cond = np.nanmedian(cond)
    print(f'\n  Condition number: median={med_cond:.1f}, max={np.nanmax(cond):.1f}')
    print(f'  k:  median={np.nanmedian(ks):.3f}, range={ks.min():.3f}–{ks.max():.3f}')
    print(f'  α:  median={np.nanmedian(als):.4f}, max={als.max():.4f}')
    print(f'  β:  median={np.nanmedian(bts):.4f}, max={bts.max():.4f}')

    return wav_grid, ks, als, bts, k_raw, a_raw, b_raw, [d['name'] for d in active]


def save_fits(wav, k, a, b, k_raw, a_raw, b_raw, key, sources_used):
    outpath = os.path.join(OUTDIR, f'calib_v6_fs_{key}.fits')
    cols = [
        fits.Column('WAVELENGTH', 'D', array=wav),
        fits.Column('K',          'D', array=k),
        fits.Column('ALPHA',      'D', array=a),
        fits.Column('BETA',       'D', array=b),
        fits.Column('K_RAW',      'D', array=k_raw),
        fits.Column('ALPHA_RAW',  'D', array=a_raw),
        fits.Column('BETA_RAW',   'D', array=b_raw),
    ]
    primary = fits.PrimaryHDU()
    tbl = fits.BinTableHDU.from_columns(cols)
    tbl.header['VERSION'] = 'v6'
    tbl.header['MODE']    = 'FS'
    tbl.header['GRATING'] = key.upper()
    tbl.header['SOURCES'] = ','.join(sources_used[:5])   # truncate for FITS
    tbl.header['HISTORY'] = 'FS v6: 4 CALSPEC stds + PID1492 cross-grating'
    fits.HDUList([primary, tbl]).writeto(outpath, overwrite=True)
    print(f'  Saved: {outpath}')


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    for grating in ('g140m', 'g235m'):
        pid1492_src = build_pid1492_source(grating)
        if pid1492_src is None:
            print(f'  PID 1492 not usable for {grating.upper()} — running with 4 stds only')

        wav, k, a, b, kr, ar, br, src_names = solve_grating(
            grating, extra_sources=[pid1492_src])
        key = f'{grating}_{"f100lp" if grating == "g140m" else "f170lp"}'
        save_fits(wav, k, a, b, kr, ar, br, key, src_names)

    print('\nFS v6 solving complete.')


if __name__ == '__main__':
    main()
