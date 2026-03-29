"""
Parlanti et al. (2025) k/α/β solver — Fixed Slit version 2 (fs_v2)

Mirrors the IFU v2 algorithm (solve_parlanti_ifu_v2.py) but uses per-exposure
FS NRS2 x1d extractions in place of IFU stage3_ext cubes.

v2 algorithm
------------
The FS v1 solver used a regularised NNLS with a k-prior anchored to ~1.0,
which biased k too high (FS v1 medians: 0.846 / 0.976).

v2 uses the same two-step approach proven in IFU v2:

  Step 1 – k(λ)
    For each k-source, compute the observed ratio
        R_i(λ) = S_obs_i(λ) / f_CALSPEC_i(λ)
    then subtract the Parlanti second/third-order contributions:
        k_i(λ) = R_i(λ) − α_P(λ)·r₂_i(λ) − β_P(λ)·r₃_i(λ)
    k(λ) = median of k_i, then 40-channel box-smoothed.

  k-sources: P330E (PID 1538), J1743045 (PID 1536).
  G191-B2B excluded: carries large source-dependent photometric excess at
  λ > 2.5 µm that contaminates the k estimate.

  Step 2 – α̃(λ), β̃(λ)
    Loaded directly from the Parlanti published FITS calibration files:
      data/parlanti_repo/calibration_files/calibration_functions_g140m_f100lp.fits
      data/parlanti_repo/calibration_files/calibration_functions_g235m_f170lp.fits
    These are pipeline-independent (detector-optics based).

Data files used:
  S_obs:   data/PID{pid}/nrs2_spec2_cal/*_nrs2_x1d.fits  (FS per-exposure NRS2)
  f(λ):    CALSPEC model SEDs (Jy)
  α,β:     data/parlanti_repo/calibration_files/calibration_functions_*.fits

Outputs:
  CSV  plots/Parlanti/cal/fs_v2/coeffs_fs_v2_{G140M|G235M}.csv
  PNG  plots/Parlanti/cal/fs_v2/fs_v2_{G140M|G235M}_coeffs.png
"""
import os

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import interp1d

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE         = '/Users/dcoe/NIRSpec/wavext'
FS_DIR       = f'{BASE}/data'
CALSPEC_DIR  = f'{BASE}/data/CALSPEC'
PARLANTI_CAL = f'{BASE}/data/parlanti_repo/calibration_files'
OUTPUT_DIR   = f'{BASE}/nirspec_wavext_work/plots/Parlanti/cal/fs_v2'
os.makedirs(OUTPUT_DIR, exist_ok=True)

C_ANG_S = 2.99792458e18   # Å/s

# ── Source tables ──────────────────────────────────────────────────────────────
# Sources used for k derivation per grating.
# G191-B2B excluded from k due to known hot-WD pipeline photometric excess.
K_SOURCES = {
    'G140M': [
        {'name': 'P330E',    'calspec': 'p330e_mod_008.fits',
         'nrs2_path': f'{FS_DIR}/PID1538_P330E/nrs2_spec2_cal/jw01538160001_06101_00001_nrs2_x1d.fits'},
        {'name': 'J1743045', 'calspec': '1743045_mod_007.fits',
         'nrs2_path': f'{FS_DIR}/PID1536_J1743045/nrs2_spec2_cal/jw01536002001_03104_00003_nrs2_x1d.fits'},
    ],
    'G235M': [
        {'name': 'P330E',    'calspec': 'p330e_mod_008.fits',
         'nrs2_path': f'{FS_DIR}/PID1538_P330E/nrs2_spec2_cal/jw01538160001_08101_00005_nrs2_x1d.fits'},
        {'name': 'J1743045', 'calspec': '1743045_mod_007.fits',
         'nrs2_path': f'{FS_DIR}/PID1536_J1743045/nrs2_spec2_cal/jw01536002001_03106_00001_nrs2_x1d.fits'},
    ],
}

# All sources (for validation plot; includes G191-B2B)
ALL_SOURCES = [
    {'name': 'P330E',    'calspec': 'p330e_mod_008.fits',
     'nrs2_g140m': f'{FS_DIR}/PID1538_P330E/nrs2_spec2_cal/jw01538160001_06101_00001_nrs2_x1d.fits',
     'nrs2_g235m': f'{FS_DIR}/PID1538_P330E/nrs2_spec2_cal/jw01538160001_08101_00005_nrs2_x1d.fits',
     'l3_g235m': f'{FS_DIR}/PID1538_P330E/jw01538-o160_t002-s000000001_nirspec_f170lp-g235m-s1600a1-sub2048_x1d.fits',
     'l3_g395m': f'{FS_DIR}/PID1538_P330E/jw01538-o160_t002-s000000001_nirspec_f290lp-g395m-s1600a1-sub2048_x1d.fits',
    },
    {'name': 'G191-B2B', 'calspec': 'g191b2b_mod_012.fits',
     'nrs2_g140m': f'{FS_DIR}/PID1537_G191-B2B/nrs2_spec2_cal/jw01537007001_07101_00003_nrs2_x1d.fits',
     'nrs2_g235m': f'{FS_DIR}/PID1537_G191-B2B/nrs2_spec2_cal/jw01537007001_09101_00003_nrs2_x1d.fits',
     'l3_g235m': f'{FS_DIR}/PID1537_G191-B2B/jw01537-o007_t001-s000000001_nirspec_f170lp-g235m-s1600a1-sub2048_x1d.fits',
     'l3_g395m': f'{FS_DIR}/PID1537_G191-B2B/jw01537-o007_t001-s000000001_nirspec_f290lp-g395m-s1600a1-sub2048_x1d.fits',
    },
    {'name': 'J1743045', 'calspec': '1743045_mod_007.fits',
     'nrs2_g140m': f'{FS_DIR}/PID1536_J1743045/nrs2_spec2_cal/jw01536002001_03104_00003_nrs2_x1d.fits',
     'nrs2_g235m': f'{FS_DIR}/PID1536_J1743045/nrs2_spec2_cal/jw01536002001_03106_00001_nrs2_x1d.fits',
     'l3_g235m': f'{FS_DIR}/PID1536_J1743045/jw01536-o002_t004-s000000001_nirspec_f170lp-g235m-s1600a1-sub2048_x1d.fits',
     'l3_g395m': f'{FS_DIR}/PID1536_J1743045/jw01536-o002_t004-s000000001_nirspec_f290lp-g395m-s1600a1-sub2048_x1d.fits',
    },
]

NRS2_LO = {'G140M': 1.93, 'G235M': 3.25}
GRID = {
    'G140M': np.linspace(1.95, 3.50, 300),
    'G235M': np.linspace(3.30, 5.25, 300),
}

# ── Helpers ────────────────────────────────────────────────────────────────────
def load_spec(path, label=''):
    if not os.path.exists(path):
        if label:
            print(f'  WARNING missing: {os.path.basename(path)} [{label}]')
        return np.array([]), np.array([])
    with fits.open(path) as h:
        ext = h[1].data
        wl  = ext['WAVELENGTH'].astype(float)
        fl  = ext['FLUX'].astype(float)
        dq  = (ext['DQ'].astype(int)
               if 'DQ' in ext.dtype.names
               else np.zeros(len(wl), int))
    good = (wl > 0) & np.isfinite(fl) & (fl > 0) & (dq == 0)
    return wl[good], fl[good]


def calspec_jy(fname):
    """CALSPEC in Jy vs µm."""
    with fits.open(f'{CALSPEC_DIR}/{fname}') as h:
        d    = h[1].data
        wl_a = d['WAVELENGTH'].astype(float)
        flam = d['FLUX'].astype(float)
    fnu   = flam * wl_a**2 / C_ANG_S * 1e23
    order = np.argsort(wl_a)
    return interp1d(wl_a[order] / 1e4, fnu[order],
                    bounds_error=False, fill_value=np.nan)


def box_smooth(y, box=40):
    """40-channel box smooth (Parlanti §2.2)."""
    s = np.convolve(y, np.ones(box) / box, mode='same')
    s[:box // 2]  = y[:box // 2]
    s[-box // 2:] = y[-box // 2:]
    return s


def load_parlanti_coeffs(grating):
    """Return (k_interp, alpha_interp, beta_interp) from Parlanti FITS."""
    fname = {
        'G140M': 'calibration_functions_g140m_f100lp.fits',
        'G235M': 'calibration_functions_g235m_f170lp.fits',
    }[grating]
    path = os.path.join(PARLANTI_CAL, fname)
    if not os.path.exists(path):
        print(f'  WARNING: Parlanti cal file not found: {path}')
        return (lambda x: np.ones_like(x, float),
                lambda x: np.zeros_like(x, float),
                lambda x: np.zeros_like(x, float))
    with fits.open(path) as h:
        d   = h[1].data
        wlp = d['wavelength']
        kp  = d['k']
        ap  = d['alpha']
        bp  = d['beta']
    ki = interp1d(wlp, kp, bounds_error=False, fill_value=np.nan)
    ai = interp1d(wlp, np.clip(ap, 0, None), bounds_error=False, fill_value=0.0)
    bi = interp1d(wlp, np.clip(bp, 0, None), bounds_error=False, fill_value=0.0)
    return ki, ai, bi


# ── Core solver ────────────────────────────────────────────────────────────────
def solve_grating(grating):
    """
    Derive k(λ), α(λ), β(λ) for one grating using the FS v2 algorithm:

    1.  For each k-source, compute the observed ratio
            R_i(λ) = S_obs_i(λ) / f_CALSPEC_i(λ)
        and the Parlanti-corrected ratio
            k_i(λ) = R_i(λ) − α_P(λ)·r2_i(λ) − β_P(λ)·r3_i(λ)

    2.  k(λ) = median of valid k_i(λ), then 40-channel box-smoothed.

    3.  α(λ) = Parlanti α_P(λ),  β(λ) = Parlanti β_P(λ).
    """
    print(f'\n{"="*60}')
    print(f'Solving {grating} (FS v2: k from median R; α,β from Parlanti)')
    print(f'{"="*60}')

    wav_grid  = GRID[grating]
    nrs2_lo   = NRS2_LO[grating]
    srcs      = K_SOURCES[grating]
    k_parlanti, a_parlanti, b_parlanti = load_parlanti_coeffs(grating)

    active = []
    for src in srcs:
        w_all, f_all = load_spec(src['nrs2_path'], f"{src['name']} FS NRS2")
        if len(w_all) < 10:
            continue
        mask = w_all >= nrs2_lo
        w_nrs, f_nrs = w_all[mask], f_all[mask]
        if len(w_nrs) < 10:
            print(f'  WARNING: {src["name"]} — no NRS2 data above {nrs2_lo} µm')
            continue

        f_cs  = calspec_jy(src['calspec'])
        r_mid = float(np.interp(np.median(w_nrs), w_nrs, f_nrs)) / float(f_cs(np.median(w_nrs)))
        if not (r_mid > 0.05):
            print(f'  WARNING: {src["name"]} — R at median λ = {r_mid:.4f} (too low; skipping)')
            continue

        f_obs_i = interp1d(w_nrs, f_nrs, bounds_error=False, fill_value=np.nan)
        print(f'  {src["name"]}: NRS2 {w_nrs.min():.3f}–{w_nrs.max():.3f} µm '
              f'({len(w_nrs)} px)  R_mid={r_mid:.3f}')
        active.append({'name': src['name'],
                       'f_obs': f_obs_i,
                       'f_ref': calspec_jy(src['calspec'])})

    if len(active) < 1:
        raise RuntimeError(f'No valid sources for {grating}')

    # ── Compute k_i(λ) = R_i − α_P·r2_i − β_P·r3_i at each grid point ───────
    k_per_source = []
    for d in active:
        k_i = []
        for lam in wav_grid:
            f_ref_val = float(d['f_ref'](lam))
            s_obs_val = float(d['f_obs'](lam))
            if not (np.isfinite(s_obs_val) and s_obs_val > 0 and
                    np.isfinite(f_ref_val) and f_ref_val > 0):
                k_i.append(np.nan)
                continue
            R_i   = s_obs_val / f_ref_val
            f_sec = float(d['f_ref'](lam / 2))
            f_thi = float(d['f_ref'](lam / 3))
            r2_i  = (f_sec / f_ref_val
                     if (np.isfinite(f_sec) and f_sec > 0) else 0.0)
            r3_i  = (f_thi / f_ref_val
                     if (np.isfinite(f_thi) and f_thi > 0) else 0.0)
            alpha_p = float(a_parlanti(lam))
            beta_p  = float(b_parlanti(lam))
            k_i.append(R_i - alpha_p * r2_i - beta_p * r3_i)
        k_per_source.append(np.array(k_i))

    # ── Median k across sources (NaN-safe) ────────────────────────────────────
    stack    = np.vstack(k_per_source)
    k_raw    = np.nanmedian(stack, axis=0)
    mask_nan = ~np.isfinite(k_raw)
    if mask_nan.any():
        xi = np.where(~mask_nan)[0]
        k_raw[mask_nan] = np.interp(np.where(mask_nan)[0], xi, k_raw[xi])

    k_raw = np.clip(k_raw, 0.01, None)
    ks    = box_smooth(k_raw, 40)
    ks    = np.clip(ks, 0.01, None)

    als = np.nan_to_num(np.clip(a_parlanti(wav_grid).astype(float), 0, None), nan=0.0)
    bts = np.nan_to_num(np.clip(b_parlanti(wav_grid).astype(float), 0, None), nan=0.0)

    # ── Diagnostics ───────────────────────────────────────────────────────────
    print(f'\n  {grating} k(λ): median={np.median(ks):.3f}, '
          f'range={np.min(ks):.3f}–{np.max(ks):.3f}')
    print(f'  {grating} α(λ) [Parlanti]: max={np.nanmax(als):.4f}')
    print(f'  {grating} β(λ) [Parlanti]: max={np.nanmax(bts):.4f}')
    for d, ki_arr in zip(active, k_per_source):
        print(f'    k_{d["name"]}: median={np.nanmedian(ki_arr):.3f}  '
              f'range={np.nanmin(ki_arr):.3f}–{np.nanmax(ki_arr):.3f}')

    # ── Save CSV ───────────────────────────────────────────────────────────────
    csv_path = os.path.join(OUTPUT_DIR, f'coeffs_fs_v2_{grating}.csv')
    np.savetxt(csv_path,
               np.column_stack([wav_grid, ks, als, bts]),
               delimiter=',', header='wavelength_um,k,alpha,beta', comments='')
    print(f'  CSV saved: {csv_path}')

    return wav_grid, ks, als, bts, active


# ── Quick validation plot ──────────────────────────────────────────────────────
def plot_coeffs(grating, wav_grid, ks, als, bts, active):
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(wav_grid, ks,  color='k',          lw=2.0, label=r'$k(\lambda)$ FS v2')
    ax.plot(wav_grid, als, color='darkorange',  lw=1.5, label=r'$\tilde{\alpha}(\lambda)$ [Parlanti]')
    ax.plot(wav_grid, bts, color='steelblue',   lw=1.5, label=r'$\tilde{\beta}(\lambda)$ [Parlanti]')
    ax.axhline(1.0, color='r', lw=0.8, ls='--', alpha=0.4, label='k=1')
    ax.set_yscale('log')
    ax.set_ylim(1e-4, 3.0)
    ax.set_xlabel('Wavelength [µm]', fontsize=11)
    ax.set_ylabel('Value', fontsize=11)
    ax.set_title(f'FS v2 — {grating} calibration coefficients', fontsize=12)
    ax.legend(fontsize=9)
    ax.grid(True, which='both', alpha=0.2)
    plt.tight_layout()
    outpath = os.path.join(OUTPUT_DIR, f'fs_v2_{grating.lower()}_coeffs.png')
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Quick-check plot saved: {outpath}')


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    print('=== FS Calibration v2 Solver (Parlanti model, v2 approach) ===')
    print(f'Output: {OUTPUT_DIR}')
    for grating in ['G140M', 'G235M']:
        wav_grid, ks, als, bts, active = solve_grating(grating)
        plot_coeffs(grating, wav_grid, ks, als, bts, active)
    print('\n=== Done ===')


if __name__ == '__main__':
    main()
