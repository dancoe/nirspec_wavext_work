"""
NIRSpec Wavelength Extension — IFU v6 Solver

Key upgrade over IFU v4/v5
---------------------------
v4/v5 used only the 3 CALSPEC standard stars (P330E, G191-B2B, J1743045).
v6 adds UGC-5101 (PID 2186) as a 4th constraint source for G235M using a
cross-grating truth: the G395M/F290LP NOMINAL stage3 x1d acts as the truth
spectrum for the G235M NRS2 extended region (3.15–5.27 µm overlap). This
science target data contributes to the coefficient determination regardless
of its non-standard status.

G140M:  3 CALSPEC sources  (P330E + G191-B2B + J1743045)
G235M:  3 CALSPEC sources + UGC-5101 cross-grating  (4 sources total)

Cross-grating truth for UGC-5101
---------------------------------
  observed  = G235M stage3_ext x1d, NRS2 part (λ > 3.15 µm), Jy
  truth(λ)  = G395M stage3_nom x1d interpolated at λ        [Jy]
  truth(λ/2)= G235M stage3_ext x1d (NRS1 part, λ/2 ~ 1.6–2.6) [Jy]
  truth(λ/3)= G235M stage3_ext x1d (NRS1 part, λ/3 ~ 1.0–1.8) [Jy]

All spectra already Jy — no unit conversion needed.

Outputs
--------
  results/v6/calib_v6_ifu_g140m_f100lp.fits
  results/v6/calib_v6_ifu_g235m_f170lp.fits
"""
import os
import numpy as np
from astropy.io import fits
from scipy.interpolate import interp1d
from scipy.optimize import nnls

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE        = '/Users/dcoe/NIRSpec/wavext'
IFU_DIR     = f'{BASE}/data/IFU'
DATA_DIR    = f'{BASE}/data'
CALSPEC_DIR = f'{BASE}/data/CALSPEC'
OUTDIR      = f'{BASE}/results/v6'
os.makedirs(OUTDIR, exist_ok=True)

C_ANG_S = 2.99792458e18  # speed of light in Å/s

# ── CALSPEC standard star sources ─────────────────────────────────────────────
CALSPEC_SOURCES = [
    {
        'name':    'P330E',
        'sptype':  'G2V',
        'calspec': 'p330e_mod_008.fits',
        'ifu_dir': f'{IFU_DIR}/PID1538_P330E',
    },
    {
        'name':    'G191-B2B',
        'sptype':  'WD-DA.8',
        'calspec': 'g191b2b_mod_012.fits',
        'ifu_dir': f'{IFU_DIR}/PID1537_G191-B2B',
    },
    {
        'name':    'J1743045',
        'sptype':  'A8III',
        'calspec': '1743045_mod_007.fits',
        'ifu_dir': f'{IFU_DIR}/PID1536_J1743045',
    },
]

# Stage3_ext x1d filename tag per grating
EXT_TAG = {
    'G140M': 'f100lp_g140m-f100lp_x1d.fits',
    'G235M': 'f170lp_g235m-f170lp_x1d.fits',
}

# NRS2 extended lower wavelength boundary (µm)
NRS2_LO = {'G140M': 1.87, 'G235M': 3.15}

# Wavelength grid for solver
GRID = {
    'G140M': np.linspace(1.87, 3.55, 400),
    'G235M': np.linspace(3.15, 5.27, 400),
}

SMOOTH_BOX = 40  # channels, identical to Parlanti 2025

# ── Helpers ────────────────────────────────────────────────────────────────────
def load_spec(path, label=''):
    """Load x1d FITS, return (wl, flux) with bad pixels masked."""
    if not os.path.exists(path):
        print(f'  MISSING: {os.path.basename(path)} [{label}]')
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
    """Return interpolator (µm → Jy) for a CALSPEC model."""
    path = os.path.join(CALSPEC_DIR, fname)
    with fits.open(path) as h:
        d    = h[1].data
        wl_a = d['WAVELENGTH'].astype(float)
        flam = d['FLUX'].astype(float)
    fnu_jy = flam * wl_a**2 / C_ANG_S * 1e23
    idx = np.argsort(wl_a)
    return interp1d(wl_a[idx] / 1e4, fnu_jy[idx],
                    bounds_error=False, fill_value=np.nan)


def make_interp(wl, fl, extrapolate='nan'):
    """Build linear interp1d; optionally extrapolate as nearest or nan."""
    idx = np.argsort(wl)
    w, f = wl[idx], fl[idx]
    if extrapolate == 'nearest':
        fv = (float(f[0]), float(f[-1]))
    else:
        fv = np.nan
    return interp1d(w, f, bounds_error=False, fill_value=fv)


def box_smooth(y, box=40):
    """Simple boxcar smooth, preserving edge values."""
    s = np.convolve(y, np.ones(box) / box, mode='same')
    s[:box // 2]  = y[:box // 2]
    s[-box // 2:] = y[-box // 2:]
    return s


# ── Build science-target (cross-grating) truth spectrum ───────────────────────
def build_ugc5101_truth_g235m():
    """
    Build a stitched truth interpolator for UGC-5101 G235M correction.

    Stitch:
      λ  < 3.0 µm : G235M stage3_ext NRS1 part  (reliable nominal calibration)
      λ >= 3.0 µm : G395M stage3_nom              (nominal G395M coverage)

    Returns interp1d covering 0.9–5.5 µm for use as truth(λ), truth(λ/2),
    and truth(λ/3).
    """
    path_g235m = f'{DATA_DIR}/PID2186_UGC-5101/stage3_ext/g235m_f170lp_g235m-f170lp_x1d.fits'
    path_g395m = f'{DATA_DIR}/PID2186_UGC5101/stage3_nom/g395m_f290lp_g395m-f290lp_x1d.fits'

    wl235, fl235 = load_spec(path_g235m, 'UGC5101 G235M ext')
    wl395, fl395 = load_spec(path_g395m, 'UGC5101 G395M nom')

    if len(wl235) < 10 or len(wl395) < 10:
        return None

    # Use G235M for λ < 3.0 µm and G395M for λ >= 3.0 µm
    # Scale G395M to match G235M in the 3.0–3.15 µm overlap region
    mask_235_lo = wl235 < 3.15
    mask_235_xg = (wl235 >= 2.87) & (wl235 < 3.15)   # overlap zone
    mask_395_xg = (wl395 >= 2.87) & (wl395 < 3.15)   # overlap zone

    if mask_235_xg.sum() < 5 or mask_395_xg.sum() < 5:
        # No overlap to derive scale: use both as-is
        scale = 1.0
        print('  UGC5101: no G235M/G395M overlap for scale — using scale=1')
    else:
        f235_med = np.nanmedian(np.interp(wl395[mask_395_xg],
                                          wl235[mask_235_xg.nonzero()[0]],
                                          fl235[mask_235_xg]))
        f395_med = np.nanmedian(fl395[mask_395_xg])
        if f395_med > 0 and f235_med > 0:
            scale = f235_med / f395_med
        else:
            scale = 1.0
        print(f'  UGC5101: G235M/G395M overlap scale = {scale:.4f}')

    # Combine: G235M up to 3.0 µm, scaled G395M from 3.0 µm
    STITCH = 3.0
    wl_lo = wl235[wl235 < STITCH]
    fl_lo = fl235[wl235 < STITCH]
    wl_hi = wl395[wl395 >= STITCH]
    fl_hi = fl395[wl395 >= STITCH] * scale

    wl_all = np.concatenate([wl_lo, wl_hi])
    fl_all = np.concatenate([fl_lo, fl_hi])
    idx = np.argsort(wl_all)
    wl_all, fl_all = wl_all[idx], fl_all[idx]

    return make_interp(wl_all, fl_all, extrapolate='nearest')


# ── Core NNLS solver ──────────────────────────────────────────────────────────
def solve_grating(grating, extra_sources=None):
    """
    Solve k/α/β by NNLS for a given grating.

    extra_sources : list of dicts, each with keys:
        'name'   : str label
        'f_obs'  : interp1d of observed NRS2 extended spectrum (Jy)
        'f_ref'  : interp1d of truth spectrum (Jy) for direct + ghost terms
    """
    print(f'\n{"="*60}')
    print(f'IFU v6 NNLS — {grating}')
    print(f'{"="*60}')

    wav_grid = GRID[grating]
    nrs2_lo  = NRS2_LO[grating]

    active = []

    # 1. Load CALSPEC sources
    for src in CALSPEC_SOURCES:
        path = os.path.join(src['ifu_dir'], 'stage3_ext', EXT_TAG[grating])
        wl_all, fl_all = load_spec(path, f"{src['name']} stage3_ext {grating}")
        if len(wl_all) < 10:
            continue
        # Use only NRS2 extended part as observed
        mask = wl_all >= nrs2_lo
        wl_nrs, fl_nrs = wl_all[mask], fl_all[mask]
        if len(wl_nrs) < 10:
            print(f'  SKIP {src["name"]}: <10 NRS2 pixels')
            continue

        f_cs = calspec_jy_interp(src['calspec'])
        r_mid = (float(np.interp(np.median(wl_nrs), wl_nrs, fl_nrs))
                 / float(f_cs(np.median(wl_nrs))))
        if not (r_mid > 0.01):
            print(f'  SKIP {src["name"]}: implausible ratio R_mid={r_mid:.4f}')
            continue

        active.append({
            'name':   src['name'],
            'sptype': src['sptype'],
            'f_obs':  make_interp(wl_nrs, fl_nrs),
            'f_ref':  f_cs,
        })
        print(f'  {src["name"]} ({src["sptype"]}): NRS2 '
              f'{wl_nrs.min():.3f}–{wl_nrs.max():.3f} µm  R_mid={r_mid:.3f}')

    # 2. Add cross-grating science sources
    if extra_sources:
        for es in extra_sources:
            active.append(es)
            print(f'  {es["name"]} [cross-grating]: added as {grating} source')

    if len(active) < 2:
        raise RuntimeError(f'Only {len(active)} source(s) for {grating}; need ≥2')
    print(f'  Total sources: {len(active)}')

    # NNLS at each grid point
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

    # Fill isolated gaps by linear interpolation
    for arr in (k_raw, a_raw, b_raw):
        bad = ~np.isfinite(arr)
        if bad.any():
            xi = np.where(~bad)[0]
            if len(xi) > 1:
                arr[bad] = np.interp(np.where(bad)[0], xi, arr[xi])

    k_raw = np.clip(k_raw, 0.01, None)
    a_raw = np.clip(a_raw, 0.0,  None)
    b_raw = np.clip(b_raw, 0.0,  None)

    # 40-channel boxcar smooth
    ks  = np.clip(box_smooth(k_raw, SMOOTH_BOX), 0.01, None)
    als = np.clip(box_smooth(a_raw, SMOOTH_BOX), 0.0,  None)
    bts = np.clip(box_smooth(b_raw, SMOOTH_BOX), 0.0,  None)

    med_cond = np.nanmedian(cond)
    print(f'\n  Condition number: median={med_cond:.1f}, max={np.nanmax(cond):.1f}')
    print(f'  k:  median={np.nanmedian(ks):.3f}, range={ks.min():.3f}–{ks.max():.3f}')
    print(f'  α:  median={np.nanmedian(als):.4f}, max={als.max():.4f}')
    print(f'  β:  median={np.nanmedian(bts):.4f}, max={bts.max():.4f}')

    return wav_grid, ks, als, bts, k_raw, a_raw, b_raw


def save_fits(wav, k, a, b, k_raw, a_raw, b_raw, key, sources_used):
    """Save coefficient table to FITS."""
    outpath = os.path.join(OUTDIR, f'calib_v6_ifu_{key}.fits')
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
    tbl.header['MODE']    = 'IFU'
    tbl.header['GRATING'] = key.upper()
    tbl.header['SOURCES'] = ','.join(sources_used)
    tbl.header['HISTORY'] = 'IFU v6: CALSPEC standards + cross-grating science targets'
    fits.HDUList([primary, tbl]).writeto(outpath, overwrite=True)
    print(f'  Saved: {outpath}')


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    # ── G140M: 3 CALSPEC sources only (no cross-grating partner for G140M IFU) ─
    wav, k, a, b, kr, ar, br = solve_grating('G140M', extra_sources=None)
    save_fits(wav, k, a, b, kr, ar, br,
              'g140m_f100lp',
              [s['name'] for s in CALSPEC_SOURCES])

    # ── G235M: 3 CALSPEC sources + UGC-5101 cross-grating ──────────────────────
    print('\nBuilding UGC-5101 cross-grating truth for G235M...')
    f_truth_ugc = build_ugc5101_truth_g235m()

    if f_truth_ugc is not None:
        # Load UGC-5101 G235M extended observed spectrum (NRS2 part)
        path_obs = (f'{DATA_DIR}/PID2186_UGC-5101/stage3_ext/'
                    'g235m_f170lp_g235m-f170lp_x1d.fits')
        wl_ugc, fl_ugc = load_spec(path_obs, 'UGC5101 G235M ext observed')
        nrs2_lo = NRS2_LO['G235M']
        mask = wl_ugc >= nrs2_lo
        wl_nrs, fl_nrs = wl_ugc[mask], fl_ugc[mask]

        if len(wl_nrs) >= 10:
            # Sanity: median ratio obs/truth in extended region
            truth_med = float(f_truth_ugc(np.median(wl_nrs)))
            obs_med   = float(np.nanmedian(fl_nrs))
            r_mid     = obs_med / truth_med if truth_med > 0 else 0
            print(f'  UGC5101 NRS2 obs/truth ratio at median λ: {r_mid:.3f}')

            ugc_source = {
                'name':   'UGC5101 [cross-grating G395M truth]',
                'sptype': 'ULIRG-z0.039',
                'f_obs':  make_interp(wl_nrs, fl_nrs),
                'f_ref':  f_truth_ugc,
            }
            extra = [ugc_source]
        else:
            print('  UGC5101 G235M NRS2: too few pixels, skipping')
            extra = None
    else:
        print('  UGC5101 truth not available — G235M solved with 3 CALSPEC only')
        extra = None

    wav, k, a, b, kr, ar, br = solve_grating('G235M', extra_sources=extra)
    sources_g235 = [s['name'] for s in CALSPEC_SOURCES]
    if extra:
        sources_g235.append('UGC5101-xgrating')
    save_fits(wav, k, a, b, kr, ar, br, 'g235m_f170lp', sources_g235)

    print('\nIFU v6 solving complete.')


if __name__ == '__main__':
    main()
