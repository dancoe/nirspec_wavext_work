"""
Parlanti et al. (2025) k/α/β solver — IFU version 3 (ifu_v3)

Upgrade over IFU v2
--------------------
IFU v2 already used Level-3 stage3_ext products (Spec3Pipeline output combining
NRS1+NRS2 with Parlanti cubepar). IFU v3 makes this explicit, adds MAST Level-3
comparison plots (stage3 NRS1-only vs stage3_ext NRS1+NRS2), and is labelled
consistently with the FS v3 analysis.

Custom pipeline tweaks applied (WAVEXT.md):
  1. jwst/assign_wcs/nirspec.py: NRS2 M-grating hard-coded error removed
  2. wavelengthrange_0008.asdf: extended range for G140M/G235M
  3. Parlanti sflat/fflat/cubepar overrides (NRS2 extended coverage)
  PYTHONPATH not needed — fork installed via pip install -e .

Data sources
-------------
  IFU extended (our pipeline):  stage3_ext/{grating}-x1d.fits
  MAST L3 nominal (NRS1 only):  stage3/{grating}-x1d.fits

Algorithm (identical to v2)
-----------------------------
  Step 1 – k(λ) from median unbiased ratio:
    k_i(λ) = R_i(λ) − α_P(λ)·r₂_i(λ) − β_P(λ)·r₃_i(λ)
    k(λ) = median_i[k_i(λ)], 40-channel box-smoothed.
    k-sources: P330E (PID 1538), J1743045 (PID 1536).
    G191-B2B excluded: large source-dependent photometric excess at λ > 2.5 µm.
    P330E-C3 (PID 6645): no valid G140M or G235M NRS2 data; excluded.

  Step 2 – α̃(λ), β̃(λ): Parlanti published calibration FITS.

Outputs
--------
  CSV  plots/Parlanti/cal/ifu_v3/coeffs_ifu_v3_{G140M|G235M}.csv
  PNG  plots/Parlanti/cal/ifu_v3/  (coeffs, spectra, MAST vs L3 ext, full)
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import interp1d

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE         = '/Users/dcoe/NIRSpec/wavext'
IFU_DIR      = f'{BASE}/data/IFU'
FS_DIR       = f'{BASE}/data'
CALSPEC_DIR  = f'{BASE}/data/CALSPEC'
PARLANTI_CAL = f'{BASE}/data/parlanti_repo/calibration_files'
OUTPUT_DIR   = f'{BASE}/nirspec_wavext_work/plots/Parlanti/cal/ifu_v3'
os.makedirs(OUTPUT_DIR, exist_ok=True)

C_ANG_S = 2.99792458e18  # Å/s

# ── Source tables ──────────────────────────────────────────────────────────────
K_SOURCES = {
    'G140M': [
        {'name': 'P330E',    'pid': '1538', 'calspec': 'p330e_mod_008.fits',
         'ifu_dir': f'{IFU_DIR}/PID1538_P330E'},
        {'name': 'J1743045', 'pid': '1536', 'calspec': '1743045_mod_007.fits',
         'ifu_dir': f'{IFU_DIR}/PID1536_J1743045'},
    ],
    'G235M': [
        {'name': 'P330E',    'pid': '1538', 'calspec': 'p330e_mod_008.fits',
         'ifu_dir': f'{IFU_DIR}/PID1538_P330E'},
        {'name': 'J1743045', 'pid': '1536', 'calspec': '1743045_mod_007.fits',
         'ifu_dir': f'{IFU_DIR}/PID1536_J1743045'},
    ],
}

ALL_SOURCES = [
    {
        'name':      'P330E',
        'pid':       '1538',
        'calspec':   'p330e_mod_008.fits',
        'ifu_dir':   f'{IFU_DIR}/PID1538_P330E',
        'mast_g395m_fs': f'{FS_DIR}/PID1538_P330E/jw01538-o160_t002-s000000001_nirspec_f290lp-g395m-s1600a1-sub2048_x1d.fits',
    },
    {
        'name':      'G191-B2B',
        'pid':       '1537',
        'calspec':   'g191b2b_mod_012.fits',
        'ifu_dir':   f'{IFU_DIR}/PID1537_G191-B2B',
        'mast_g395m_fs': f'{FS_DIR}/PID1537_G191-B2B/jw01537-o007_t001-s000000001_nirspec_f290lp-g395m-s1600a1-sub2048_x1d.fits',
    },
    {
        'name':      'J1743045',
        'pid':       '1536',
        'calspec':   '1743045_mod_007.fits',
        'ifu_dir':   f'{IFU_DIR}/PID1536_J1743045',
        'mast_g395m_fs': f'{FS_DIR}/PID1536_J1743045/jw01536-o002_t004-s000000001_nirspec_f290lp-g395m-s1600a1-sub2048_x1d.fits',
    },
]

NRS2_LO = {'G140M': 1.87, 'G235M': 3.15}
GRID = {
    'G140M': np.linspace(1.87, 3.55, 400),
    'G235M': np.linspace(3.15, 5.27, 400),
}

# Mapping from grating to stage3_ext x1d filename tag
EXT_TAG = {
    'G140M': 'f100lp_g140m-f100lp_x1d.fits',
    'G235M': 'f170lp_g235m-f170lp_x1d.fits',
}
# Mapping from grating to stage3 (MAST L3 nominal, NRS1 only) filename tag
NOM_TAG = {
    'G140M': 'f100lp_g140m-f100lp_x1d.fits',
    'G235M': 'f170lp_g235m-f170lp_x1d.fits',
}

# ── Helpers ────────────────────────────────────────────────────────────────────
def load_spec(path, label=''):
    if not os.path.exists(path):
        if label:
            print(f'  WARNING missing: {os.path.basename(path)} [{label}]')
        return np.array([]), np.array([])
    with fits.open(path) as h:
        d  = h[1].data
        wl = d['WAVELENGTH'].astype(float)
        fl = d['FLUX'].astype(float)
        dq = (d['DQ'].astype(int) if 'DQ' in d.dtype.names
              else np.zeros(len(wl), int))
    good = (wl > 0) & np.isfinite(fl) & (fl > 0) & (dq == 0)
    return wl[good], fl[good]


def calspec_jy(fname):
    with fits.open(f'{CALSPEC_DIR}/{fname}') as h:
        d    = h[1].data
        wl_a = d['WAVELENGTH'].astype(float)
        flam = d['FLUX'].astype(float)
    fnu   = flam * wl_a**2 / C_ANG_S * 1e23
    order = np.argsort(wl_a)
    return interp1d(wl_a[order] / 1e4, fnu[order],
                    bounds_error=False, fill_value=np.nan)


def box_smooth(y, box=40):
    s = np.convolve(y, np.ones(box) / box, mode='same')
    s[:box // 2]  = y[:box // 2]
    s[-box // 2:] = y[-box // 2:]
    return s


def load_parlanti_coeffs(grating):
    fname = {
        'G140M': 'calibration_functions_g140m_f100lp.fits',
        'G235M': 'calibration_functions_g235m_f170lp.fits',
    }[grating]
    path = os.path.join(PARLANTI_CAL, fname)
    if not os.path.exists(path):
        print(f'  WARNING: {path} not found')
        return (lambda x: np.ones_like(x, float),
                lambda x: np.zeros_like(x, float),
                lambda x: np.zeros_like(x, float))
    with fits.open(path) as h:
        d = h[1].data
        wlp, kp, ap, bp = d['wavelength'], d['k'], d['alpha'], d['beta']
    return (interp1d(wlp, kp,               bounds_error=False, fill_value=np.nan),
            interp1d(wlp, np.clip(ap, 0, None), bounds_error=False, fill_value=0.0),
            interp1d(wlp, np.clip(bp, 0, None), bounds_error=False, fill_value=0.0))


# ── Core solver ────────────────────────────────────────────────────────────────
def solve_grating(grating):
    print(f'\n{"="*60}')
    print(f'IFU v3 Solving {grating} (Level-3 stage3_ext x1d; same algorithm as v2)')
    print(f'{"="*60}')

    wav_grid = GRID[grating]
    nrs2_lo  = NRS2_LO[grating]
    _, a_parlanti, b_parlanti = load_parlanti_coeffs(grating)

    active = []
    for src in K_SOURCES[grating]:
        path = os.path.join(src['ifu_dir'], 'stage3_ext', EXT_TAG[grating])
        w_all, f_all = load_spec(path, f"{src['name']} stage3_ext")
        if len(w_all) < 10:
            continue
        mask = w_all >= nrs2_lo
        w_nrs, f_nrs = w_all[mask], f_all[mask]
        if len(w_nrs) < 10:
            print(f'  WARNING: {src["name"]} — <10 NRS2 pixels above {nrs2_lo} µm')
            continue

        f_cs  = calspec_jy(src['calspec'])
        r_mid = float(np.interp(np.median(w_nrs), w_nrs, f_nrs)) / float(f_cs(np.median(w_nrs)))
        if not (r_mid > 0.05):
            print(f'  WARNING: {src["name"]} — R_mid={r_mid:.4f} too low; skipping')
            continue

        f_obs_i = interp1d(w_nrs, f_nrs, bounds_error=False, fill_value=np.nan)
        print(f'  {src["name"]}: NRS2 {w_nrs.min():.3f}–{w_nrs.max():.3f} µm '
              f'({len(w_nrs)} px)  R_mid={r_mid:.3f}')
        active.append({'name': src['name'],
                       'f_obs': f_obs_i,
                       'f_ref': calspec_jy(src['calspec'])})

    if not active:
        raise RuntimeError(f'No valid k-sources for {grating}')

    k_per_source = []
    for d in active:
        k_i = []
        for lam in wav_grid:
            f_ref_val = float(d['f_ref'](lam))
            s_obs_val = float(d['f_obs'](lam))
            if not (np.isfinite(s_obs_val) and s_obs_val > 0 and
                    np.isfinite(f_ref_val) and f_ref_val > 0):
                k_i.append(np.nan); continue
            R_i   = s_obs_val / f_ref_val
            f_sec = float(d['f_ref'](lam / 2))
            f_thi = float(d['f_ref'](lam / 3))
            r2_i  = (f_sec / f_ref_val if (np.isfinite(f_sec) and f_sec > 0) else 0.0)
            r3_i  = (f_thi / f_ref_val if (np.isfinite(f_thi) and f_thi > 0) else 0.0)
            alpha_p = float(a_parlanti(lam))
            beta_p  = float(b_parlanti(lam))
            k_i.append(R_i - alpha_p * r2_i - beta_p * r3_i)
        k_per_source.append(np.array(k_i))

    stack = np.vstack(k_per_source)
    k_raw = np.nanmedian(stack, axis=0)
    mask_nan = ~np.isfinite(k_raw)
    if mask_nan.any():
        xi = np.where(~mask_nan)[0]
        k_raw[mask_nan] = np.interp(np.where(mask_nan)[0], xi, k_raw[xi])
    k_raw = np.clip(k_raw, 0.01, None)
    ks    = box_smooth(k_raw, 40)
    ks    = np.clip(ks, 0.01, None)

    als = np.nan_to_num(np.clip(a_parlanti(wav_grid).astype(float), 0, None), nan=0.0)
    bts = np.nan_to_num(np.clip(b_parlanti(wav_grid).astype(float), 0, None), nan=0.0)

    print(f'\n  {grating} k(λ): median={np.median(ks):.3f}, '
          f'range={np.min(ks):.3f}–{np.max(ks):.3f}')
    print(f'  {grating} α(λ) [Parlanti]: max={np.nanmax(als):.4f}')
    print(f'  {grating} β(λ) [Parlanti]: max={np.nanmax(bts):.4f}')
    for d, ki_arr in zip(active, k_per_source):
        print(f'    k_{d["name"]}: median={np.nanmedian(ki_arr):.3f}  '
              f'range={np.nanmin(ki_arr):.3f}–{np.nanmax(ki_arr):.3f}')

    csv_path = os.path.join(OUTPUT_DIR, f'coeffs_ifu_v3_{grating}.csv')
    np.savetxt(csv_path,
               np.column_stack([wav_grid, ks, als, bts]),
               delimiter=',', header='wav,k,alpha,beta', comments='')
    print(f'  CSV saved: {csv_path}')

    return wav_grid, ks, als, bts, k_raw, active


# ── Plot 1: coefficients (3-panel) ────────────────────────────────────────────
def plot_coeffs(grating, wav_grid, k_raw, ks, als, bts, active):
    fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
    fig.suptitle(
        f'IFU v3 — {grating} NRS2 calibration coefficients\n'
        f'Level-3 stage3_ext x1d; k = median(R − Parlanti contributions); 40-ch smoothed',
        fontsize=12)

    ax = axes[0]
    ax.plot(wav_grid, k_raw, color='0.75', lw=0.8, label='k raw (median)')
    ax.plot(wav_grid, ks,    color='k',    lw=2.0, label='k smoothed')
    ax.axhline(1.0, color='r', lw=0.8, ls='--', alpha=0.5, label='k=1')
    ax.set_ylabel('k(λ)'); ax.set_ylim(-0.05, 2.2)
    ax.legend(fontsize=8); ax.grid(True, alpha=0.2)
    ax.set_title(f'{grating} first-order efficiency (IFU v3, Level-3)', fontsize=10)

    axes[1].plot(wav_grid, als, color='darkorange', lw=2.0, label='α̃ [Parlanti]')
    axes[1].axhline(0, color='0.5', lw=0.6)
    axes[1].set_ylabel('α̃(λ)'); axes[1].legend(fontsize=8); axes[1].grid(True, alpha=0.2)

    axes[2].plot(wav_grid, bts, color='steelblue', lw=2.0, label='β̃ [Parlanti]')
    axes[2].axhline(0, color='0.5', lw=0.6)
    axes[2].set_ylabel('β̃(λ)'); axes[2].set_xlabel('Wavelength [µm]')
    axes[2].legend(fontsize=8); axes[2].grid(True, alpha=0.2)

    plt.tight_layout()
    out = os.path.join(OUTPUT_DIR, f'ifu_v3_{grating.lower()}_coeffs.png')
    plt.savefig(out, dpi=150, bbox_inches='tight'); plt.close()
    print(f'  Saved: {out}')


# ── Plot 2: per-source NRS2 validation ────────────────────────────────────────
def plot_source_spectra(grating, wav_grid, ks, als, bts):
    k_itp  = interp1d(wav_grid, ks,  bounds_error=False, fill_value=np.nan)
    al_itp = interp1d(wav_grid, als, bounds_error=False, fill_value=0.0)
    bt_itp = interp1d(wav_grid, bts, bounds_error=False, fill_value=0.0)

    for src in ALL_SOURCES:
        name  = src['name']
        f_cs  = calspec_jy(src['calspec'])
        nrs2_lo = NRS2_LO[grating]

        ext_path   = os.path.join(src['ifu_dir'], 'stage3_ext', EXT_TAG[grating])
        truth_path = src['mast_g395m_fs'] if grating == 'G235M' else \
                     os.path.join(src['ifu_dir'], 'stage3_ext',
                                  'f170lp_g235m-f170lp_x1d.fits')

        w_all, f_all = load_spec(ext_path, f'{name} stage3_ext')
        w_tru, f_tru = load_spec(truth_path, f'{name} truth')

        fig, ax = plt.subplots(figsize=(10, 5))
        ax.set_title(f'IFU v3 — {name} / {grating}  [Level-3 stage3_ext]', fontsize=12)

        wav_fine = np.linspace(min(wav_grid), max(wav_grid), 800)
        ax.plot(wav_fine, f_cs(wav_fine), 'k:', lw=1.0, alpha=0.6, label='CALSPEC model')

        if len(w_tru) > 5:
            truth_label = ('MAST L3 G395M FS (truth)' if grating == 'G235M'
                           else 'IFU stage3_ext G235M (truth)')
            ax.plot(w_tru, f_tru, 'b--', lw=1.5, alpha=0.8, label=truth_label)

        if len(w_all) > 5:
            m = w_all >= nrs2_lo
            w_e, f_e = w_all[m], f_all[m]
            ax.plot(w_e, f_e, color='0.6', lw=1.0, alpha=0.7,
                    label='IFU stage3_ext NRS2 (raw)')
            if len(w_e) > 10:
                k_v  = k_itp(w_e)
                al_v = al_itp(w_e)
                bt_v = bt_itp(w_e)
                f_sc = f_cs(w_e / 2)
                f_tc = f_cs(w_e / 3)
                f_c  = (f_e - al_v * np.where(np.isfinite(f_sc), f_sc, 0)
                            - bt_v * np.where(np.isfinite(f_tc), f_tc, 0)) / k_v
                good = np.isfinite(f_c) & (k_v > 0.01)
                ax.plot(w_e[good], f_c[good], 'r-', lw=2.0,
                        label='IFU stage3_ext NRS2 corrected (v3)')

        ax.set_xlabel('Wavelength [µm]', fontsize=11)
        ax.set_ylabel('Flux [Jy]', fontsize=11)
        ax.legend(fontsize=9); ax.grid(True, alpha=0.2); ax.set_ylim(bottom=0)
        plt.tight_layout()
        out = os.path.join(OUTPUT_DIR, f'ifu_v3_spectra_{name}_{grating}.png')
        plt.savefig(out, dpi=150, bbox_inches='tight'); plt.close()
        print(f'  Saved: {out}')


# ── Plot 3: MAST L3 nominal vs our L3 ext ─────────────────────────────────────
def plot_mast_vs_l3(grating, wav_grid, ks, als, bts):
    k_itp  = interp1d(wav_grid, ks,  bounds_error=False, fill_value=np.nan)
    al_itp = interp1d(wav_grid, als, bounds_error=False, fill_value=0.0)
    bt_itp = interp1d(wav_grid, bts, bounds_error=False, fill_value=0.0)

    fig, axes = plt.subplots(len(ALL_SOURCES), 1, figsize=(12, 4 * len(ALL_SOURCES)),
                             sharex=True)
    fig.suptitle(
        f'IFU v3 — {grating}: MAST Level-3 nominal vs our Level-3 ext (stage3_ext)',
        fontsize=13)

    for ax, src in zip(axes, ALL_SOURCES):
        name   = src['name']
        f_cs   = calspec_jy(src['calspec'])
        nrs2_lo = NRS2_LO[grating]

        mast_path = os.path.join(src['ifu_dir'], 'stage3', NOM_TAG[grating])
        ext_path  = os.path.join(src['ifu_dir'], 'stage3_ext', EXT_TAG[grating])

        w_mast, f_mast = load_spec(mast_path,  f'{name} MAST stage3 {grating}')
        w_ext,  f_ext  = load_spec(ext_path,   f'{name} stage3_ext {grating}')

        wav_fine = np.linspace(0.6, 5.6, 2000)
        ax.plot(wav_fine, f_cs(wav_fine), 'k:', lw=0.8, alpha=0.5, label='CALSPEC')

        if len(w_mast) > 5:
            ax.plot(w_mast, f_mast, color='steelblue', lw=1.5,
                    label=f'MAST L3 {grating} (stage3, NRS1 only)')

        if len(w_ext) > 5:
            # Full NRS1+NRS2 from our stage3_ext
            m_nrs1 = w_ext < nrs2_lo
            if m_nrs1.any():
                ax.plot(w_ext[m_nrs1], f_ext[m_nrs1], color='cornflowerblue', lw=1.2,
                        alpha=0.8, label=f'Our L3 ext {grating} NRS1')
            m_nrs2 = w_ext >= nrs2_lo
            w_e, f_e = w_ext[m_nrs2], f_ext[m_nrs2]
            if len(w_e) > 10:
                ax.plot(w_e, f_e, color='darkorange', lw=1.0, alpha=0.7,
                        label=f'Our L3 ext {grating} NRS2 (raw)')
                k_v  = k_itp(w_e)
                al_v = al_itp(w_e)
                bt_v = bt_itp(w_e)
                f_sc = f_cs(w_e / 2)
                f_tc = f_cs(w_e / 3)
                f_c  = (f_e - al_v * np.where(np.isfinite(f_sc), f_sc, 0)
                            - bt_v * np.where(np.isfinite(f_tc), f_tc, 0)) / k_v
                good = np.isfinite(f_c) & (k_v > 0.01)
                ax.plot(w_e[good], f_c[good], 'r-', lw=2.0,
                        label=f'Our L3 ext {grating} NRS2 (corrected v3)')

        ax.set_ylabel('Flux [Jy]'); ax.set_title(name, fontsize=11)
        ax.legend(fontsize=8); ax.grid(True, alpha=0.2); ax.set_ylim(bottom=0)

    axes[-1].set_xlabel('Wavelength [µm]', fontsize=11)
    plt.tight_layout()
    out = os.path.join(OUTPUT_DIR, f'ifu_v3_mast_vs_l3_{grating}.png')
    plt.savefig(out, dpi=150, bbox_inches='tight'); plt.close()
    print(f'  Saved: {out}')


# ── Plot 4: full spectrum ──────────────────────────────────────────────────────
def plot_full_spectrum(src, results):
    name  = src['name']
    ifu_d = src['ifu_dir']
    f_cs  = calspec_jy(src['calspec'])

    fig, ax = plt.subplots(figsize=(14, 5))
    ax.set_title(
        f'IFU v3 — {name}: full spectrum  (MAST L3 nominal + IFU Level-3 ext corrected)',
        fontsize=12)

    wav_all = np.linspace(0.6, 5.6, 3000)
    ax.plot(wav_all, f_cs(wav_all), 'k:', lw=0.8, alpha=0.5, label='CALSPEC')

    # MAST L3 nominal (stage3)
    for grating, tag, color, label in [
        ('G140M', 'f100lp_g140m-f100lp_x1d.fits', 'steelblue',  'MAST L3 G140M IFU (stage3, NRS1)'),
        ('G235M', 'f170lp_g235m-f170lp_x1d.fits', 'mediumblue', 'MAST L3 G235M IFU (stage3, NRS1)'),
    ]:
        p = os.path.join(ifu_d, 'stage3', tag)
        w, f = load_spec(p, f'{name} stage3 {grating}')
        if len(w) > 5:
            ax.plot(w, f, color=color, lw=1.3, alpha=0.8, label=label)

    # MAST L3 G395M from FS (truth for G235M NRS2)
    w_g395, f_g395 = load_spec(src['mast_g395m_fs'], f'{name} G395M FS')
    if len(w_g395) > 5:
        ax.plot(w_g395, f_g395, color='navy', lw=1.3, alpha=0.8, label='MAST L3 G395M FS')

    # Our L3 ext corrected (NRS2 portion only)
    for grating in ('G140M', 'G235M'):
        wav_grid, ks, als, bts, _, _ = results[grating]
        k_itp  = interp1d(wav_grid, ks,  bounds_error=False, fill_value=np.nan)
        al_itp = interp1d(wav_grid, als, bounds_error=False, fill_value=0.0)
        bt_itp = interp1d(wav_grid, bts, bounds_error=False, fill_value=0.0)

        nrs2_lo = NRS2_LO[grating]
        ext_path = os.path.join(ifu_d, 'stage3_ext', EXT_TAG[grating])
        w_all, f_all = load_spec(ext_path, f'{name} stage3_ext {grating}')
        if len(w_all) > 5:
            m = w_all >= nrs2_lo
            w_e, f_e = w_all[m], f_all[m]
            if len(w_e) > 10:
                k_v  = k_itp(w_e)
                al_v = al_itp(w_e)
                bt_v = bt_itp(w_e)
                f_sc = f_cs(w_e / 2)
                f_tc = f_cs(w_e / 3)
                f_c  = (f_e - al_v * np.where(np.isfinite(f_sc), f_sc, 0)
                            - bt_v * np.where(np.isfinite(f_tc), f_tc, 0)) / k_v
                good = np.isfinite(f_c) & (k_v > 0.01)
                color = 'crimson' if grating == 'G140M' else 'orangered'
                ax.plot(w_e[good], f_c[good], color=color, lw=2.0,
                        label=f'IFU L3 ext {grating} NRS2 corrected (v3)')

    ax.set_xlabel('Wavelength [µm]', fontsize=11)
    ax.set_ylabel('Flux [Jy]', fontsize=11)
    ax.legend(fontsize=8, loc='upper right')
    ax.grid(True, alpha=0.2); ax.set_ylim(bottom=0)
    plt.tight_layout()
    out = os.path.join(OUTPUT_DIR, f'ifu_v3_full_spectrum_{name}.png')
    plt.savefig(out, dpi=150, bbox_inches='tight'); plt.close()
    print(f'  Saved: {out}')


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    print('=== IFU Calibration v3 Solver (Level-3 stage3_ext; Parlanti model) ===')
    print(f'Output: {OUTPUT_DIR}')

    results = {}
    for grating in ('G140M', 'G235M'):
        wav_grid, ks, als, bts, k_raw, active = solve_grating(grating)
        results[grating] = (wav_grid, ks, als, bts, k_raw, active)

        plot_coeffs(grating, wav_grid, k_raw, ks, als, bts, active)
        plot_source_spectra(grating, wav_grid, ks, als, bts)
        plot_mast_vs_l3(grating, wav_grid, ks, als, bts)

    for src in ALL_SOURCES:
        plot_full_spectrum(src, results)

    print('\n=== IFU v3 Done ===')
    return results


if __name__ == '__main__':
    main()
