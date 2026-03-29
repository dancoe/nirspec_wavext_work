"""
Spectral decomposition diagnostic for the Parlanti IFU extended-wavelength exercise.

Compares IFU stage3_ext spectra (NRS1+NRS2, Jy) against IFU stage3 nominal spectra
of the NEXT grating (the ground truth for the NRS2-extended overlap region).

Goal: verify that the extended G140M NRS2 data (1.87–3.60 µm) matches the
      independently-observed G235M nominal data (1.66–3.17 µm) after k/α/β correction.

Data used (all in Jy from the standard JWST pipeline):
  stage3      = nominal IFU cube-extracted x1d  (NRS1 only)
  stage3_ext  = extended IFU cube-extracted x1d  (NRS1 + NRS2)
  G395M FS    = Fixed-Slit stage3 x1d used as truth for G235M NRS2 overlap

Key comparisons:
  G140M NRS2 extension (1.87–3.60 µm) ↔ G235M stage3 nominal (1.66–3.17 µm)
  G235M NRS2 extension (3.15–5.50 µm) ↔ G395M FS nominal    (2.87–5.14 µm)
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy.io import fits
from scipy.interpolate import interp1d
from scipy.ndimage import median_filter

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE         = '/Users/dcoe/NIRSpec/wavext'
IFU_DIR      = f'{BASE}/data/IFU'
FS_DIR       = f'{BASE}/data'          # FS per-PID directories live here
CALSPEC_DIR  = f'{BASE}/data/CALSPEC'
COEFF_DIR    = f'{BASE}/nirspec_wavext_work/plots/Parlanti/cal/153678_v3'
OUTPUT_DIR   = COEFF_DIR
os.makedirs(OUTPUT_DIR, exist_ok=True)

C_ANG_S = 2.99792458e18  # Å/s

# ── Source table ──────────────────────────────────────────────────────────────
# (pid, label, calspec_file, ifu_subdir, fs_subdir, obs_id_g140m, obs_id_g235m, obs_id_g395m)
SOURCES = {
    '1538_P330E': {
        'label':    'P330E',
        'calspec':  'p330e_mod_008.fits',
        'ifu_dir':  f'{IFU_DIR}/PID1538_P330E',
        'g395m_fs': f'{FS_DIR}/PID1538_P330E/jw01538-o160_t002-s000000001_nirspec_f290lp-g395m-s1600a1-sub2048_x1d.fits',
        'color':    '#f39c12',   # orange
    },
    '1537_G191B2B': {
        'label':    'G191-B2B',
        'calspec':  'g191b2b_mod_012.fits',
        'ifu_dir':  f'{IFU_DIR}/PID1537_G191-B2B',
        'g395m_fs': f'{FS_DIR}/PID1537_G191-B2B/jw01537-o007_t001-s000000001_nirspec_f290lp-g395m-s1600a1-sub2048_x1d.fits',
        'color':    '#2980b9',   # blue
    },
    '1536_J1743': {
        'label':    'J1743045',
        'calspec':  '1743045_mod_007.fits',
        'ifu_dir':  f'{IFU_DIR}/PID1536_J1743045',
        'g395m_fs': f'{FS_DIR}/PID1536_J1743045/jw01536-o002_t004-s000000001_nirspec_f290lp-g395m-s1600a1-sub2048_x1d.fits',
        'color':    '#27ae60',   # green
    },
}

# ── Helpers ────────────────────────────────────────────────────────────────────
def load_spec(path, label=''):
    if not os.path.exists(path):
        print(f"  WARNING: missing {os.path.basename(path)}" + (f' [{label}]' if label else ''))
        return np.array([]), np.array([])
    with fits.open(path) as h:
        d   = h[1].data
        wl  = d['WAVELENGTH'].astype(float)
        fl  = d['FLUX'].astype(float)
        dq  = d['DQ'].astype(int) if 'DQ' in d.dtype.names else np.zeros(len(wl), int)
    good = (wl > 0) & np.isfinite(fl) & (dq == 0)
    return wl[good], fl[good]

def load_calspec_jy(fname):
    with fits.open(f'{CALSPEC_DIR}/{fname}') as h:
        d     = h[1].data
        wl_a  = d['WAVELENGTH'].astype(float)
        flam  = d['FLUX'].astype(float)
    fnu_jy = flam * (wl_a**2) / C_ANG_S * 1e23
    order  = np.argsort(wl_a)
    return wl_a[order] / 1e4, fnu_jy[order]

def load_coeffs(path):
    if not os.path.exists(path):
        return None
    data = np.loadtxt(path, delimiter=',', skiprows=1)
    return {'wav': data[:, 0], 'k': data[:, 1], 'alpha': data[:, 2], 'beta': data[:, 3]}

def smooth(y, size=21):
    ok = np.isfinite(y)
    if ok.sum() < size:
        return y
    out = np.full_like(y, np.nan)
    fill = np.where(ok, y, np.nanmedian(y[ok]))
    out[ok] = median_filter(fill, size=size)[ok]
    return out

def apply_parlanti(w_obs, f_obs, w_coeff, k, alpha, beta, f_ref_interp):
    """Apply Parlanti correction: (f_obs - α·f(λ/2) - β·f(λ/3)) / k"""
    ki     = interp1d(w_coeff, k,     bounds_error=False, fill_value=np.nan)(w_obs)
    ai     = interp1d(w_coeff, alpha, bounds_error=False, fill_value=np.nan)(w_obs)
    bi     = interp1d(w_coeff, beta,  bounds_error=False, fill_value=np.nan)(w_obs)
    ki     = np.where(np.isfinite(ki), ki, np.nanmedian(k))
    ai     = np.nan_to_num(ai)
    bi     = np.nan_to_num(bi)
    f_hl   = f_ref_interp(w_obs / 2)
    f_trd  = f_ref_interp(w_obs / 3)
    return (f_obs - ai * f_hl - bi * f_trd) / ki


# ══════════════════════════════════════════════════════════════════════════════
# MAIN PLOT: one figure per source
# ══════════════════════════════════════════════════════════════════════════════
def plot_components(src_key):
    src       = SOURCES[src_key]
    label     = src['label']
    ifu_dir   = src['ifu_dir']
    color     = src['color']

    print(f"\n=== {label} ===")

    # ── Load IFU data ──────────────────────────────────────────────────────────
    # stage3 nominal (NRS1 only, Jy) — used as ground truth
    w_nom140, f_nom140 = load_spec(f'{ifu_dir}/stage3/f100lp_g140m-f100lp_x1d.fits',   'G140M nom')
    w_nom235, f_nom235 = load_spec(f'{ifu_dir}/stage3/f170lp_g235m-f170lp_x1d.fits',   'G235M nom')

    # stage3_ext (NRS1+NRS2, Jy) — NRS2 portion has incomplete flat/throughput model
    w_ext140, f_ext140 = load_spec(f'{ifu_dir}/stage3_ext/f100lp_g140m-f100lp_x1d.fits', 'G140M ext')
    w_ext235, f_ext235 = load_spec(f'{ifu_dir}/stage3_ext/f170lp_g235m-f170lp_x1d.fits', 'G235M ext')

    # G395M Fixed-Slit nominal (Jy) — ground truth for G235M NRS2 overlap
    w_395, f_395 = load_spec(src['g395m_fs'], 'G395M FS')

    # ── CALSPEC reference ──────────────────────────────────────────────────────
    wl_cs, fj_cs  = load_calspec_jy(src['calspec'])
    f_calspec     = interp1d(wl_cs, fj_cs, bounds_error=False, fill_value=np.nan)
    wgrid_full    = np.linspace(0.9, 5.6, 4000)

    # ── Parlanti coefficients ──────────────────────────────────────────────────
    c140 = load_coeffs(f'{COEFF_DIR}/coeffs_G140M.csv')
    c235 = load_coeffs(f'{COEFF_DIR}/coeffs_G235M.csv')

    # ── Figure layout ──────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(14, 10), layout='constrained')
    gs  = gridspec.GridSpec(2, 2, figure=fig, hspace=0.35, wspace=0.30)
    ax_top   = fig.add_subplot(gs[0, :])   # full spectrum comparison
    ax_k140  = fig.add_subplot(gs[1, 0])   # k(λ) for G140M NRS2
    ax_k235  = fig.add_subplot(gs[1, 1])   # k(λ) for G235M NRS2

    fig.suptitle(f'NIRSpec IFU Extended-Wavelength Calibration Check — {label}\n'
                 f'stage3_ext (Jy) vs stage3 nominal (ground truth)',
                 fontsize=13)

    # ── Top panel: full spectra ────────────────────────────────────────────────
    ax = ax_top

    # CALSPEC reference
    ax.plot(wgrid_full, f_calspec(wgrid_full),
            color='0.5', lw=1.5, ls='--', label='CALSPEC reference', zorder=0)

    # G140M: stage3 nominal (NRS1 truth)
    if len(w_nom140):
        ax.plot(w_nom140, smooth(f_nom140, 15),
                color='#1a5276', lw=1.5, alpha=0.9, label='G140M stage3 nominal (NRS1 truth)')

    # G140M: stage3_ext (NRS1+NRS2 — NRS2 under-calibrated)
    if len(w_ext140):
        ax.plot(w_ext140, smooth(f_ext140, 15),
                color='#5dade2', lw=1.0, alpha=0.7, label='G140M stage3_ext (NRS1+NRS2, uncorrected)')

    # G235M: stage3 nominal (ground truth for G140M NRS2 overlap)
    if len(w_nom235):
        ax.plot(w_nom235, smooth(f_nom235, 15),
                color='#784212', lw=1.8, alpha=0.9, label='G235M stage3 nominal (truth @ 1.7–3.2 µm)')

    # G235M: stage3_ext (NRS2 under-calibrated)
    if len(w_ext235):
        ax.plot(w_ext235, smooth(f_ext235, 15),
                color='#f0b27a', lw=1.0, alpha=0.7, label='G235M stage3_ext (NRS1+NRS2, uncorrected)')

    # G395M FS nominal (truth for G235M NRS2 overlap)
    if len(w_395):
        ax.plot(w_395, smooth(f_395, 15),
                color='#922b21', lw=1.8, alpha=0.9, label='G395M FS nominal (truth @ 2.9–5.1 µm)')

    # G140M corrected (after Parlanti k/α/β applied to NRS2 portion)
    if len(w_ext140) and c140 is not None:
        # Only correct NRS2 portion; leave NRS1 as-is
        nrs1_mask = w_ext140 < 1.87
        nrs2_mask = w_ext140 >= 1.87
        f_corr_140 = f_ext140.copy()
        if nrs2_mask.sum() > 5:
            f_corr_140[nrs2_mask] = apply_parlanti(
                w_ext140[nrs2_mask], f_ext140[nrs2_mask],
                c140['wav'], c140['k'], c140['alpha'], c140['beta'],
                f_calspec
            )
        ax.plot(w_ext140, smooth(f_corr_140, 15),
                color='#1a5276', lw=1.5, ls=':', alpha=0.8,
                label='G140M stage3_ext corrected (V3 k/α/β)')

    # G235M corrected  
    if len(w_ext235) and c235 is not None:
        nrs2_mask = w_ext235 >= 3.15
        f_corr_235 = f_ext235.copy()
        if nrs2_mask.sum() > 5:
            f_corr_235[nrs2_mask] = apply_parlanti(
                w_ext235[nrs2_mask], f_ext235[nrs2_mask],
                c235['wav'], c235['k'], c235['alpha'], c235['beta'],
                f_calspec
            )
        ax.plot(w_ext235, smooth(f_corr_235, 15),
                color='#784212', lw=1.5, ls=':', alpha=0.8,
                label='G235M stage3_ext corrected (V3 k/α/β)')

    ax.set_yscale('log')
    ax.set_xlim(0.9, 5.6)
    cs_vals = f_calspec(wgrid_full[(wgrid_full > 0.9) & (wgrid_full < 5.6)])
    ok_cs = np.isfinite(cs_vals) & (cs_vals > 0)
    if ok_cs.any():
        ylo = cs_vals[ok_cs].min() / 8
        yhi = cs_vals[ok_cs].max() * 8
        ax.set_ylim(ylo, yhi)
    else:
        ylo, yhi = 1e-4, 1e-1
        ax.set_ylim(ylo, yhi)
    ax.set_xlabel('Wavelength (µm)', fontsize=11)
    ax.set_ylabel('Flux (Jy)', fontsize=11)
    ax.legend(fontsize=8, loc='lower left', ncol=2, frameon=True, framealpha=0.9)
    ax.grid(True, alpha=0.2, which='both')
    # Mark the two overlap comparison windows
    ax.axvspan(1.87, 3.17, color='#5dade2', alpha=0.07, label='_')
    ax.axvspan(3.15, 5.14, color='#f0b27a', alpha=0.07, label='_')
    # Annotations at top of shaded regions
    ax.text(2.3,  yhi * 0.55, 'G140M NRS2\nvs G235M NOM', fontsize=7,
            color='#1a5276', ha='center', va='top')
    ax.text(4.1,  yhi * 0.55, 'G235M NRS2\nvs G395M NOM', fontsize=7,
            color='#784212', ha='center', va='top')

    # ── Bottom-left: k(λ) for G140M NRS2 ─────────────────────────────────────
    ax = ax_k140
    ax.axhline(1.0, color='k', lw=0.8, ls='--', alpha=0.5)

    if len(w_ext140) and len(w_nom235):
        # Both are in Jy; build k = G235M_NOM / G140M_EXT in overlap
        f_235_interp = interp1d(w_nom235, f_nom235, bounds_error=False, fill_value=np.nan)
        ovl = (w_ext140 >= 1.95) & (w_ext140 <= 3.17) & np.isfinite(f_ext140) & (f_ext140 > 0)
        if ovl.sum() > 10:
            k_meas = f_235_interp(w_ext140[ovl]) / f_ext140[ovl]
            ax.plot(w_ext140[ovl], smooth(k_meas, 31),
                    color='#5dade2', lw=1.5, label='k = G235M NOM / G140M EXT (measured)')
    if c140 is not None:
        ax.plot(c140['wav'], c140['k'], color='#1a5276', lw=1.5, ls='--', label='k (V3 coefficients)')

    ax.set_xlim(1.87, 3.60)
    ax.set_ylim(0, 3)
    ax.set_xlabel('Wavelength (µm)', fontsize=10)
    ax.set_ylabel('k(λ)  [G235M NOM / G140M EXT]', fontsize=10)
    ax.set_title('G140M NRS2 throughput factor', fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.2)

    # ── Bottom-right: k(λ) for G235M NRS2 ────────────────────────────────────
    ax = ax_k235
    ax.axhline(1.0, color='k', lw=0.8, ls='--', alpha=0.5)

    if len(w_ext235) and len(w_395):
        f_395_interp = interp1d(w_395, f_395, bounds_error=False, fill_value=np.nan)
        ovl = (w_ext235 >= 3.20) & (w_ext235 <= 5.14) & np.isfinite(f_ext235) & (f_ext235 > 0)
        if ovl.sum() > 10:
            k_meas2 = f_395_interp(w_ext235[ovl]) / f_ext235[ovl]
            ax.plot(w_ext235[ovl], smooth(k_meas2, 31),
                    color='#f0b27a', lw=1.5, label='k = G395M NOM / G235M EXT (measured)')
    if c235 is not None:
        ax.plot(c235['wav'], c235['k'], color='#784212', lw=1.5, ls='--', label='k (V3 coefficients)')

    ax.set_xlim(3.15, 5.50)
    ax.set_ylim(0, 5)
    ax.set_xlabel('Wavelength (µm)', fontsize=10)
    ax.set_ylabel('k(λ)  [G395M NOM / G235M EXT]', fontsize=10)
    ax.set_title('G235M NRS2 throughput factor', fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.2)

    outpath = f'{OUTPUT_DIR}/CAL_COMPONENTS_{label}.png'
    plt.savefig(outpath, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  SAVED: {outpath}")
    return outpath


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description='IFU extended-wavelength calibration diagnostic')
    ap.add_argument('--source', default='all',
                    help='Source key (1538_P330E, 1537_G191B2B, 1536_J1743, or "all")')
    args = ap.parse_args()

    keys = list(SOURCES.keys()) if args.source == 'all' else [args.source]
    for k in keys:
        if k not in SOURCES:
            print(f"Unknown source '{k}'. Choose from: {list(SOURCES.keys())}")
            continue
        plot_components(k)
    print("\nDone.")
