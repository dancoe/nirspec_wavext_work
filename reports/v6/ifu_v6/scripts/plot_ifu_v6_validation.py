"""
IFU v6 — Validation Plots

Generates three plot types:
1. Coefficient plots: k, α, β vs wavelength (v6 vs v5 vs Parlanti)
2. Per-source NRS2 spectra: obs vs corrected vs CALSPEC truth (3 standards)
3. UGC-5101 G235M cross-validation: corrected G235M ext vs G395M nominal

Output directory: reports/v6/ifu_v6/plots/
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import interp1d

BASE         = '/Users/dcoe/NIRSpec/wavext'
IFU_DIR      = f'{BASE}/data/IFU'
DATA_DIR     = f'{BASE}/data'
CALSPEC_DIR  = f'{BASE}/data/CALSPEC'
PARLANTI_CAL = f'{BASE}/data/parlanti_repo/calibration_files'
V5_DIR       = f'{BASE}/results/v5'
V6_DIR       = f'{BASE}/results/v6'
PLOTDIR      = f'{BASE}/nirspec_wavext_work/reports/v6/ifu_v6/plots'
os.makedirs(PLOTDIR, exist_ok=True)

C_ANG_S = 2.99792458e18

CALSPEC_SOURCES = [
    {'name': 'P330E',    'sptype': 'G2V',    'calspec': 'p330e_mod_008.fits',
     'color': '#1f77b4', 'ifu_dir': f'{IFU_DIR}/PID1538_P330E'},
    {'name': 'G191-B2B', 'sptype': 'WD',     'calspec': 'g191b2b_mod_012.fits',
     'color': '#ff7f0e', 'ifu_dir': f'{IFU_DIR}/PID1537_G191-B2B'},
    {'name': 'J1743045', 'sptype': 'A8III',  'calspec': '1743045_mod_007.fits',
     'color': '#2ca02c', 'ifu_dir': f'{IFU_DIR}/PID1536_J1743045'},
]

EXT_TAG = {'G140M': 'f100lp_g140m-f100lp_x1d.fits',
           'G235M': 'f170lp_g235m-f170lp_x1d.fits'}
NRS2_LO = {'G140M': 1.87, 'G235M': 3.15}


# ── Helpers ───────────────────────────────────────────────────────────────────
def load_spec(path):
    if not os.path.exists(path):
        return np.array([]), np.array([])
    with fits.open(path) as h:
        d  = h[1].data
        wl = d['WAVELENGTH'].astype(float)
        fl = d['FLUX'].astype(float)
        dq = (d['DQ'].astype(int) if 'DQ' in d.dtype.names
              else np.zeros(len(wl), int))
    ok = (wl > 0.3) & np.isfinite(fl) & (fl > 0) & (dq == 0)
    return wl[ok], fl[ok]


def calspec_jy(fname):
    with fits.open(os.path.join(CALSPEC_DIR, fname)) as h:
        d    = h[1].data
        wl_a = d['WAVELENGTH'].astype(float)
        flam = d['FLUX'].astype(float)
    fnu = flam * wl_a**2 / C_ANG_S * 1e23
    idx = np.argsort(wl_a)
    return interp1d(wl_a[idx] / 1e4, fnu[idx],
                    bounds_error=False, fill_value=np.nan)


def load_coeffs(path):
    with fits.open(path) as h:
        d = h[1].data
        wl, k, a, b = d['WAVELENGTH'], d['K'], d['ALPHA'], d['BETA']
    return (wl,
            interp1d(wl, k, bounds_error=False, fill_value=np.nan),
            interp1d(wl, a, bounds_error=False, fill_value=0.0),
            interp1d(wl, b, bounds_error=False, fill_value=0.0))


def load_parlanti(grating):
    fname = {'G140M': 'calibration_functions_g140m_f100lp.fits',
             'G235M': 'calibration_functions_g235m_f170lp.fits'}[grating]
    path = os.path.join(PARLANTI_CAL, fname)
    if not os.path.exists(path):
        return None
    with fits.open(path) as h:
        d = h[1].data
        wlp, kp, ap, bp = d['wavelength'], d['k'], d['alpha'], d['beta']
    return (wlp,
            interp1d(wlp, kp, bounds_error=False, fill_value=np.nan),
            interp1d(wlp, np.clip(ap, 0, None), bounds_error=False, fill_value=0.0),
            interp1d(wlp, np.clip(bp, 0, None), bounds_error=False, fill_value=0.0))


def correct_spectrum(wl_nrs, fl_nrs, f_truth, k_itp, a_itp, b_itp):
    """Apply: S_corr(λ) = (S_obs(λ) - α(λ)·T(λ/2) - β(λ)·T(λ/3)) / k(λ)"""
    k_v = k_itp(wl_nrs)
    a_v = a_itp(wl_nrs)
    b_v = b_itp(wl_nrs)
    t1  = f_truth(wl_nrs / 2.0)
    t3  = f_truth(wl_nrs / 3.0)
    corr = (fl_nrs - a_v * t1 - b_v * t3) / np.where(k_v > 0.01, k_v, np.nan)
    return corr


# ── Plot 1: Coefficient comparison (v6 vs v5 vs Parlanti) ─────────────────────
def plot_coefficients():
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('IFU v6 — k, α, β Coefficients\n'
                 'v6 (CALSPEC + UGC5101 cross-grating) vs v5 (FS-derived) vs Parlanti (2025)',
                 fontsize=12)

    for col, grating in enumerate(['G140M', 'G235M']):
        key = 'g140m_f100lp' if grating == 'G140M' else 'g235m_f170lp'

        v6_path = os.path.join(V6_DIR, f'calib_v6_ifu_{key}.fits')
        v5_path = os.path.join(V5_DIR, f'calib_v5_{key}.fits')

        if not os.path.exists(v6_path):
            print(f'  MISSING: {v6_path}')
            continue

        wl6, k6, a6, b6 = load_coeffs(v6_path)
        wv6 = np.linspace(wl6.min(), wl6.max(), 500)

        parl = load_parlanti(grating)
        nrs2_lo = NRS2_LO[grating]

        for row, (coeff_fn, label, ylim, ylab) in enumerate([
            (k6, 'k(λ)', (0, 1.3), 'k(λ)'),
            (a6, 'α(λ)', (0, 0.25), 'α(λ)'),
        ]):
            ax = axes[row, col]
            cv6 = coeff_fn(wv6)
            ax.plot(wv6, cv6, 'k-', lw=2, label='v6 (IFU, smoothed)')

            if os.path.exists(v5_path):
                _, kv5, av5, bv5 = load_coeffs(v5_path)
                cv5 = (kv5 if row == 0 else av5)(wv6)
                ax.plot(wv6, cv5, 'b--', lw=1.5, alpha=0.7, label='v5 (FS-derived)')

            if parl is not None:
                cp = (parl[1] if row == 0 else parl[2])(wv6)
                ax.plot(wv6, cp, 'r:', lw=1.5, alpha=0.7, label='Parlanti (2025)')

            ax.axvline(nrs2_lo, color='orange', lw=1, ls='--', alpha=0.5,
                       label=f'NRS2 lower boundary {nrs2_lo} µm')
            ax.axhline(1.0 if row == 0 else 0.0, color='gray', lw=0.6, ls=':')
            ax.set_xlim(wv6.min(), wv6.max())
            ax.set_ylim(*ylim)
            ax.set_xlabel('Wavelength (µm)')
            ax.set_ylabel(ylab)
            ax.set_title(f'{grating} — {label}')
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.2)

    plt.tight_layout()
    out = os.path.join(PLOTDIR, 'ifu_v6_coefficients.png')
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  SAVED: {out}')


# ── Plot 2: Per-source NRS2 validation (CALSPEC standards) ───────────────────
def plot_calspec_validation():
    for grating in ('G140M', 'G235M'):
        key     = 'g140m_f100lp' if grating == 'G140M' else 'g235m_f170lp'
        nrs2_lo = NRS2_LO[grating]

        v6_path = os.path.join(V6_DIR, f'calib_v6_ifu_{key}.fits')
        if not os.path.exists(v6_path):
            continue
        wl6, k6, a6, b6 = load_coeffs(v6_path)

        fig, axes = plt.subplots(len(CALSPEC_SOURCES), 2,
                                 figsize=(14, 4 * len(CALSPEC_SOURCES)))
        fig.suptitle(f'IFU v6 — {grating} NRS2 Extended: Obs / Corrected / CALSPEC Truth',
                     fontsize=12)

        for i, src in enumerate(CALSPEC_SOURCES):
            path = os.path.join(src['ifu_dir'], 'stage3_ext', EXT_TAG[grating])
            wl_all, fl_all = load_spec(path)
            if len(wl_all) < 10:
                continue
            mask = wl_all >= nrs2_lo
            wl_nrs, fl_nrs = wl_all[mask], fl_all[mask]

            f_cs    = calspec_jy(src['calspec'])
            fl_corr = correct_spectrum(wl_nrs, fl_nrs, f_cs, k6, a6, b6)
            fl_true = f_cs(wl_nrs)

            # Spectrum panel
            ax = axes[i, 0]
            ax.plot(wl_nrs, fl_nrs,  color='silver',         lw=1,   label='Observed (NRS2 ext)')
            ax.plot(wl_nrs, fl_corr, color=src['color'],     lw=1.5, label='Corrected (v6)')
            ax.plot(wl_nrs, fl_true, color='k', ls='--',     lw=1.5, label='CALSPEC truth')
            ax.set_ylabel('Flux (Jy)')
            ax.set_title(f'{src["name"]} ({src["sptype"]})')
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.2)
            ax.set_xlabel('Wavelength (µm)')

            # Ratio panel
            ax2 = axes[i, 1]
            ratio = fl_corr / fl_true
            ax2.plot(wl_nrs, ratio, color=src['color'], lw=1.5)
            ax2.axhline(1.0, color='k', ls='--', lw=1)
            ax2.fill_between(wl_nrs, 0.9, 1.1, alpha=0.12, color='green',
                             label='±10%')
            ax2.set_ylim(0.3, 1.7)
            ax2.set_ylabel('Corrected / Truth')
            ax2.set_title(f'{src["name"]} — residual ratio')
            ax2.legend(fontsize=8)
            ax2.grid(True, alpha=0.2)
            ax2.set_xlabel('Wavelength (µm)')

        plt.tight_layout()
        out = os.path.join(PLOTDIR, f'ifu_v6_{grating.lower()}_calspec_validation.png')
        plt.savefig(out, dpi=150, bbox_inches='tight')
        plt.close()
        print(f'  SAVED: {out}')


# ── Plot 3: UGC-5101 G235M cross-validation ────────────────────────────────────
def plot_ugc5101_xval():
    key     = 'g235m_f170lp'
    nrs2_lo = NRS2_LO['G235M']
    v6_path = os.path.join(V6_DIR, f'calib_v6_ifu_{key}.fits')
    if not os.path.exists(v6_path):
        print('  MISSING v6 G235M coefficients')
        return

    _, k6, a6, b6 = load_coeffs(v6_path)

    path_g235m = f'{DATA_DIR}/PID2186_UGC-5101/stage3_ext/g235m_f170lp_g235m-f170lp_x1d.fits'
    path_g395m = f'{DATA_DIR}/PID2186_UGC5101/stage3_nom/g395m_f290lp_g395m-f290lp_x1d.fits'

    wl235, fl235 = load_spec(path_g235m)
    wl395, fl395 = load_spec(path_g395m)

    if len(wl235) < 10 or len(wl395) < 10:
        print('  UGC5101 data missing — skip xval plot')
        return

    mask = wl235 >= nrs2_lo
    wl_nrs, fl_nrs = wl235[mask], fl235[mask]

    # Build stitched truth (same as solver)
    STITCH = 3.0
    wl_lo = wl235[wl235 < STITCH]; fl_lo = fl235[wl235 < STITCH]
    # Scale G395M to G235M in overlap
    mask_ov = (wl395 >= 2.87) & (wl395 < STITCH)
    mask_235_ov = (wl235 >= 2.87) & (wl235 < STITCH)
    if mask_ov.sum() > 5 and mask_235_ov.sum() > 5:
        f235_med = np.nanmedian(np.interp(wl395[mask_ov], wl235, fl235))
        f395_med = np.nanmedian(fl395[mask_ov])
        scale = f235_med / f395_med if f395_med > 0 else 1.0
    else:
        scale = 1.0
    wl_hi = wl395[wl395 >= STITCH]; fl_hi = fl395[wl395 >= STITCH] * scale
    wl_comb = np.concatenate([wl_lo, wl_hi])
    fl_comb = np.concatenate([fl_lo, fl_hi])
    idx = np.argsort(wl_comb)
    f_truth = interp1d(wl_comb[idx], fl_comb[idx],
                       bounds_error=False, fill_value='extrapolate')

    fl_corr = correct_spectrum(wl_nrs, fl_nrs, f_truth, k6, a6, b6)
    fl_ref  = np.interp(wl_nrs, wl395, fl395,
                        left=np.nan, right=np.nan)

    fig, axes = plt.subplots(2, 1, figsize=(12, 9))
    fig.suptitle('IFU v6 — UGC-5101 G235M Cross-Validation\n'
                 'Corrected NRS2 ext vs G395M Nominal (independent reference)',
                 fontsize=12)

    ax = axes[0]
    ax.plot(wl235[wl235 < nrs2_lo], fl235[wl235 < nrs2_lo],
            color='steelblue', lw=1.2, label='G235M NRS1 nominal')
    ax.plot(wl_nrs, fl_nrs,  color='silver',      lw=1,   label='G235M NRS2 observed')
    ax.plot(wl_nrs, fl_corr, color='darkorange',  lw=1.5, label='G235M NRS2 corrected (v6)')
    ax.plot(wl395,  fl395,   color='k', ls='--',  lw=1.5, label='G395M nominal (truth)')
    ax.axvline(nrs2_lo, color='gray', ls=':', lw=1)
    ax.set_ylabel('Flux (Jy)')
    ax.set_xlabel('Wavelength (µm)')
    ax.set_title('UGC-5101 — G235M Extended vs G395M Nominal')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.2)
    ax.set_xlim(3.0, 5.4)

    ax2 = axes[1]
    ratio = fl_corr / fl_ref
    ax2.plot(wl_nrs, ratio, color='darkorange', lw=1.5)
    ax2.axhline(1.0, color='k', ls='--', lw=1)
    ax2.fill_between(wl_nrs, 0.8, 1.2, alpha=0.12, color='green', label='±20%')
    ax2.fill_between(wl_nrs, 0.9, 1.1, alpha=0.12, color='green', label='±10%')
    ax2.set_ylim(0.2, 2.0)
    ax2.axvline(nrs2_lo, color='gray', ls=':', lw=1)
    ax2.set_ylabel('Corrected NRS2 / G395M Nominal')
    ax2.set_xlabel('Wavelength (µm)')
    ax2.set_title('Flux Ratio: G235M corrected (v6) / G395M nominal')
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.2)
    ax2.set_xlim(3.0, 5.4)

    plt.tight_layout()
    out = os.path.join(PLOTDIR, 'ifu_v6_ugc5101_g235m_xval.png')
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  SAVED: {out}')


# ── Plot 4: Full merged spectrum for each standard ────────────────────────────
def plot_full_merged():
    for grating in ('G140M', 'G235M'):
        key     = 'g140m_f100lp' if grating == 'G140M' else 'g235m_f170lp'
        nrs2_lo = NRS2_LO[grating]
        v6_path = os.path.join(V6_DIR, f'calib_v6_ifu_{key}.fits')
        if not os.path.exists(v6_path):
            continue
        _, k6, a6, b6 = load_coeffs(v6_path)

        for src in CALSPEC_SOURCES:
            path = os.path.join(src['ifu_dir'], 'stage3_ext', EXT_TAG[grating])
            wl_all, fl_all = load_spec(path)
            if len(wl_all) < 10:
                continue

            f_cs = calspec_jy(src['calspec'])

            # NRS1 nominal (λ < NRS2_LO) and NRS2 extended (λ >= NRS2_LO)
            mask1 = wl_all < nrs2_lo
            mask2 = wl_all >= nrs2_lo
            wl1, fl1 = wl_all[mask1], fl_all[mask1]
            wl2, fl2 = wl_all[mask2], fl_all[mask2]
            fl2_corr = correct_spectrum(wl2, fl2, f_cs, k6, a6, b6)

            wl_cs = np.linspace(wl_all.min(), wl_all.max(), 3000)
            fl_cs = f_cs(wl_cs)

            fig, ax = plt.subplots(figsize=(12, 5))
            ax.plot(wl_cs, fl_cs,   'k--',          lw=1.5, alpha=0.8, label='CALSPEC truth')
            ax.plot(wl1,   fl1,     color='steelblue', lw=1.2, label='NRS1 nominal')
            ax.plot(wl2,   fl2,     color='silver',  lw=1,   alpha=0.7, label='NRS2 observed')
            ax.plot(wl2,   fl2_corr, color=src['color'], lw=1.5, label='NRS2 corrected (v6)')
            ax.axvline(nrs2_lo, color='orange', ls='--', lw=1, alpha=0.7,
                       label=f'NRS2 boundary {nrs2_lo} µm')
            ax.set_xlabel('Wavelength (µm)')
            ax.set_ylabel('Flux (Jy)')
            ax.set_title(f'IFU v6 — {src["name"]} ({src["sptype"]}) {grating} Full Spectrum')
            ax.legend(fontsize=9)
            ax.grid(True, alpha=0.2)
            y_max = float(np.nanpercentile(fl_cs, 99)) * 1.3
            ax.set_ylim(-0.001, max(y_max, 0.001))
            plt.tight_layout()
            name_safe = src['name'].replace('-', '')
            out = os.path.join(PLOTDIR,
                               f'ifu_v6_full_spectrum_{grating}_{name_safe}.png')
            plt.savefig(out, dpi=150, bbox_inches='tight')
            plt.close()
            print(f'  SAVED: {out}')


# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    print('Generating IFU v6 validation plots...')
    plot_coefficients()
    plot_calspec_validation()
    plot_ugc5101_xval()
    plot_full_merged()
    print('Done.')
