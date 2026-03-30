"""
IFU v5 Validation Plots — Science Targets

Generates extended wavelength validation plots for PID 2186 (UGC-5101 ULIRG)
and PID 2654 (SDSSJ0749/0841 AGN, when available).

For UGC-5101 (z=0.039, G235M/F170LP):
  - Raw extended spectrum (NRS2 extended: 3.15 – 5.5 µm beyond nominal)
  - Ghost-corrected spectrum using v5 k/α/β coefficients
  - Comparison: nominal vs extended regions

For SDSSJ0841 (z~2 AGN, G140M/F100LP):
  - Extended spectrum if available (shows H-alpha at expected λ)
  - Ghost-corrected spectrum
"""
import os
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from astropy.io import fits
from scipy.interpolate import interp1d
from scipy.ndimage import uniform_filter1d

BASE      = '/Users/dcoe/NIRSpec/wavext'
DATA_DIR  = f'{BASE}/data'
COEFF_DIR = f'{BASE}/results/v5'
OUT_DIR   = f'{BASE}/nirspec_wavext_work/reports/v5/ifu_v5/plots'
os.makedirs(OUT_DIR, exist_ok=True)

# Nominal (non-extended) wavelength limits for each grating
NOMINAL_LIMITS = {
    'G140M': (0.97, 1.87),  # NRS2 extended starts above ~1.87 µm
    'G235M': (1.66, 3.17),  # NRS2 extended starts above ~3.17 µm
}


def load_v5_coeffs(grating):
    """Load v5 k/α/β arrays from results/v5/."""
    filt = 'f100lp' if grating == 'G140M' else 'f170lp'
    path = f'{COEFF_DIR}/calib_v5_{grating.lower()}_{filt}.fits'
    if not os.path.exists(path):
        return None, None, None, None
    with fits.open(path) as h:
        d = h[1].data
        return d['WAVELENGTH'], d['K'], d['ALPHA'], d['BETA']


def apply_v5_correction(wl_obs, fl_obs, grating):
    """
    Apply v5 ghost correction:
        S_corr(λ) = [S_obs(λ) - α(λ)·S_obs(λ/2) - β(λ)·S_obs(λ/3)] / k(λ)

    Uses the observed spectrum itself for the ghost estimates (self-correction).
    """
    wv, kv, av, bv = load_v5_coeffs(grating)
    if wv is None:
        return fl_obs

    # Interpolator for obs spectrum → needed to evaluate S_obs(λ/2), S_obs(λ/3)
    good = np.isfinite(fl_obs)
    if good.sum() < 10:
        return fl_obs
    interp_obs = interp1d(wl_obs[good], fl_obs[good],
                          bounds_error=False, fill_value=0.0)

    # Interpolate v5 coefficients to the observed wavelength grid
    interp_k = interp1d(wv, kv, bounds_error=False, fill_value=np.nan)
    interp_a = interp1d(wv, av, bounds_error=False, fill_value=0.0)
    interp_b = interp1d(wv, bv, bounds_error=False, fill_value=0.0)

    k_obs = interp_k(wl_obs)
    a_obs = interp_a(wl_obs)
    b_obs = interp_b(wl_obs)

    # Ghost terms
    ghost_2nd = a_obs * interp_obs(wl_obs / 2.0)
    ghost_3rd = b_obs * interp_obs(wl_obs / 3.0)

    # Corrected flux (avoid division by near-zero k)
    k_safe = np.where(k_obs > 0.05, k_obs, np.nan)
    fl_corr = (fl_obs - ghost_2nd - ghost_3rd) / k_safe

    return fl_corr


def smooth(arr, width=20):
    """Apply a simple box-car smooth."""
    good = np.isfinite(arr)
    if good.sum() < width:
        return arr
    out = arr.copy()
    out[good] = uniform_filter1d(arr[good], size=width)
    return out


# ── UGC-5101 G235M Extended ───────────────────────────────────────────────────

def plot_ugc5101_g235m():
    """Extended G235M spectrum of UGC-5101 with v5 ghost correction."""
    x1d_path = f'{DATA_DIR}/PID2186_UGC-5101/stage3_ext/g235m_f170lp_g235m-f170lp_x1d.fits'
    if not os.path.exists(x1d_path):
        print(f'  MISSING: {x1d_path}')
        return

    with fits.open(x1d_path) as h:
        d = h[1].data
        wl = d['WAVELENGTH']
        fl = d['FLUX']

    # v5 corrected
    fl_corr = apply_v5_correction(wl, fl, 'G235M')
    fl_smooth = smooth(fl, width=15)
    fl_corr_smooth = smooth(fl_corr, width=15)

    # Nominal G235M upper limit
    wl_nom_max = NOMINAL_LIMITS['G235M'][1]

    fig, axes = plt.subplots(2, 1, figsize=(14, 10), sharex=True)
    fig.suptitle(
        'UGC-5101 (ULIRG, z=0.039) — G235M/F170LP Extended Wavelength\n'
        'PID 2186 | v5 Ghost Correction Applied',
        fontsize=14, fontweight='bold')

    ax, ax_r = axes

    # Top panel: raw + corrected
    mask_nom = wl <= wl_nom_max
    mask_ext = wl > wl_nom_max

    ax.plot(wl[mask_nom], fl_smooth[mask_nom], color='steelblue', lw=1.8,
            label='Extended pipeline (nominal region, λ < 3.17 µm)', zorder=5)
    ax.plot(wl[mask_ext], fl_smooth[mask_ext], color='darkorange', lw=1.8,
            label='Extended pipeline (extended region, λ > 3.17 µm)', zorder=5)
    ax.plot(wl, fl_corr_smooth, color='green', lw=1.5, ls='--', alpha=0.85,
            label='v5 Ghost-corrected', zorder=6)

    ax.axvline(wl_nom_max, color='gray', ls=':', lw=1.5, label='Nominal NRS2 limit')
    ax.set_ylabel('Flux (Jy)', fontsize=12)
    ax.legend(fontsize=10, loc='upper right')
    ax.grid(True, alpha=0.25)
    ax.set_ylim(bottom=0)

    # Annotate extended region
    ax.fill_betweenx([0, ax.get_ylim()[1] if ax.get_ylim()[1] > 0 else 1],
                     wl_nom_max, wl.max(), alpha=0.05, color='orange',
                     label='_Extended range')
    ymax = np.nanpercentile(fl_smooth[np.isfinite(fl_smooth)], 98)
    ax.set_ylim(0, 1.3 * ymax)
    ax.fill_betweenx([0, 1.3 * ymax], wl_nom_max, wl.max(),
                     alpha=0.06, color='orange')
    ax.text(0.73, 0.92, 'Extended NRS2\n(beyond nominal)',
            transform=ax.transAxes, ha='center', fontsize=10, color='darkorange',
            style='italic')

    # Bottom panel: ratio corrected/raw
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = fl_corr / fl
        ratio = np.where(np.abs(fl) > 1e-5, ratio, np.nan)
    ax_r.plot(wl, smooth(ratio, width=20), color='green', lw=1.5)
    ax_r.axhline(1.0, color='black', lw=0.8, ls='--')
    ax_r.axvline(wl_nom_max, color='gray', ls=':', lw=1.5)
    ax_r.set_ylabel('Corrected / Raw', fontsize=11)
    ax_r.set_xlabel('Wavelength (µm)', fontsize=12)
    ax_r.set_ylim(0.5, 1.5)
    ax_r.grid(True, alpha=0.25)

    plt.tight_layout()
    outpath = os.path.join(OUT_DIR, 'ifu_v5_ugc5101_g235m_extended.png')
    plt.savefig(outpath, dpi=200, bbox_inches='tight')
    plt.close()
    print(f'Saved: {outpath}')


# ── SDSSJ0841 G140M Extended ──────────────────────────────────────────────────

def plot_sdssj0841_g140m(redshift_guess=2.0):
    """Extended G140M spectrum of SDSSJ0841 with H-alpha annotation."""
    x1d_path = f'{DATA_DIR}/PID2654_SDSSJ0841/stage3_ext/g140m_f100lp_g140m-f100lp_x1d.fits'
    if not os.path.exists(x1d_path):
        # Try alternate naming
        patterns = glob.glob(f'{DATA_DIR}/PID2654_SDSSJ0841/stage3_ext/*g140m*x1d*.fits')
        if not patterns:
            print(f'  SDSSJ0841 G140M x1d not yet available (Stage 3 still running)')
            return False
        x1d_path = patterns[0]

    with fits.open(x1d_path) as h:
        d = h[1].data
        wl = d['WAVELENGTH']
        fl = d['FLUX']

    # v5 corrected
    fl_corr = apply_v5_correction(wl, fl, 'G140M')
    fl_smooth   = smooth(fl, width=20)
    fl_corr_smooth = smooth(fl_corr, width=20)

    wl_nom_max = NOMINAL_LIMITS['G140M'][1]
    ha_rest   = 0.6563  # µm
    ha_obs    = ha_rest * (1 + redshift_guess)

    fig, axes = plt.subplots(2, 1, figsize=(14, 10), sharex=True)
    fig.suptitle(
        f'SDSSJ0841 (AGN, z~{redshift_guess}) — G140M/F100LP Extended Wavelength\n'
        'PID 2654 | v5 Ghost Correction | H-α Expected Location Marked',
        fontsize=13, fontweight='bold')

    ax, ax_r = axes

    mask_nom = (wl > 0) & (wl <= wl_nom_max)
    mask_ext = wl > wl_nom_max

    if mask_nom.any():
        ax.plot(wl[mask_nom], fl_smooth[mask_nom], color='steelblue', lw=1.8,
                label='Extended pipeline (nominal, λ < 1.87 µm)', zorder=5)
    if mask_ext.any():
        ax.plot(wl[mask_ext], fl_smooth[mask_ext], color='darkorange', lw=1.8,
                label='Extended pipeline (extended, λ > 1.87 µm)', zorder=5)
    ax.plot(wl, fl_corr_smooth, color='green', lw=1.5, ls='--', alpha=0.85,
            label='v5 Ghost-corrected', zorder=6)

    ax.axvline(wl_nom_max, color='gray', ls=':', lw=1.5, label='Nominal NRS2 limit')

    # Mark expected H-alpha
    if 1.0 < ha_obs < 4.0:
        ymax = np.nanpercentile(fl_smooth[np.isfinite(fl_smooth)], 98) if np.any(np.isfinite(fl_smooth)) else 1.0
        ax.axvline(ha_obs, color='red', ls='--', lw=1.8, alpha=0.85,
                   label=f'H-α expected (z={redshift_guess}: {ha_obs:.3f} µm)')
        ax.annotate(f'H-α\nz={redshift_guess}',
                    xy=(ha_obs, 0.6 * ymax), fontsize=10, color='red',
                    ha='center', style='italic',
                    arrowprops=None)

    ax.set_ylabel('Flux (Jy)', fontsize=12)
    ax.legend(fontsize=10, loc='upper right')
    ax.grid(True, alpha=0.25)

    ymax_val = np.nanpercentile(fl[np.isfinite(fl)], 99) if np.any(np.isfinite(fl)) else 1
    ax.set_ylim(max(-0.2 * ymax_val, -1e-2), 1.5 * ymax_val)

    # Ratio panel
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = fl_corr / fl
        ratio = np.where(np.abs(fl) > 1e-6, ratio, np.nan)
    ax_r.plot(wl, smooth(ratio, width=25), color='green', lw=1.5)
    ax_r.axhline(1.0, color='black', lw=0.8, ls='--')
    ax_r.axvline(wl_nom_max, color='gray', ls=':', lw=1.5)
    if 1.0 < ha_obs < 4.0:
        ax_r.axvline(ha_obs, color='red', ls='--', lw=1.5, alpha=0.7)
    ax_r.set_ylabel('Corrected / Raw', fontsize=11)
    ax_r.set_xlabel('Wavelength (µm)', fontsize=12)
    ax_r.set_ylim(0.5, 1.5)
    ax_r.grid(True, alpha=0.25)

    plt.tight_layout()
    outpath = os.path.join(OUT_DIR, 'ifu_v5_sdssj0841_g140m_extended.png')
    plt.savefig(outpath, dpi=200, bbox_inches='tight')
    plt.close()
    print(f'Saved: {outpath}')
    return True


# ── Coefficient overview for IFU v5 context ───────────────────────────────────

def plot_v5_coeffs_ifu_context():
    """Simple 2-panel k(λ) plot for the IFU v5 report."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle(
        'v5 Calibration Coefficients Applied to IFU Science Targets\n'
        r'$S_\mathrm{corr}(\lambda) = [S_\mathrm{obs}(\lambda) - '
        r'\tilde{\alpha}\,S_\mathrm{obs}(\lambda/2) - '
        r'\tilde{\beta}\,S_\mathrm{obs}(\lambda/3)] / k(\lambda)$',
        fontsize=12)

    for ax, grating in zip(axes, ['G140M', 'G235M']):
        wv, kv, av, bv = load_v5_coeffs(grating)
        if wv is None:
            continue
        av_plot = np.where(av > 1e-7, av, 1e-7)
        bv_plot = np.where(bv > 1e-7, bv, 1e-7)
        ax.plot(wv, kv, color='black',       lw=2.2, label=r'$k(\lambda)$')
        ax.plot(wv, av_plot, color='darkorange', lw=1.5,
                label=r'$\tilde{\alpha}(\lambda)$')
        ax.plot(wv, bv_plot, color='royalblue',  lw=1.5,
                label=r'$\tilde{\beta}(\lambda)$')

        # Shade extended region
        nom_max = NOMINAL_LIMITS[grating][1]
        ax.axvspan(nom_max, wv.max(), alpha=0.07, color='orange',
                   label='Extended NRS2 region')
        ax.axvline(nom_max, color='gray', ls=':', lw=1.2)
        ax.set_title(grating, fontsize=13)
        ax.set_xlabel('Wavelength (µm)', fontsize=11)
        ax.set_ylabel('Coefficient', fontsize=11)
        ax.set_ylim(0, 2.0)
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.25)

    plt.tight_layout()
    out = os.path.join(OUT_DIR, 'ifu_v5_coefficients.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out}')


# ── UGC-5101 G235M extended vs G395M nominal cross-validation ─────────────────

def plot_ugc5101_g395m_crossval():
    """
    Cross-validation: compare v5-corrected G235M extended (NRS2 extended region)
    against nominal G395M pipeline in the overlap 3.17–5.25 µm.

    Both are from PID 2186 (UGC-5101 ULIRG, z=0.039).
    G395M is run with standard JWST pipeline (no Parlanti overrides) to serve
    as an independent reference for the G235M extended region.
    """
    g235m_path = f'{DATA_DIR}/PID2186_UGC-5101/stage3_ext/g235m_f170lp_g235m-f170lp_x1d.fits'
    g395m_path = f'{DATA_DIR}/PID2186_UGC5101/stage3_nom/g395m_f290lp_g395m-f290lp_x1d.fits'

    if not os.path.exists(g235m_path):
        print(f'  MISSING G235M ext: {g235m_path}')
        return
    if not os.path.exists(g395m_path):
        print(f'  MISSING G395M nom: {g395m_path}')
        return

    with fits.open(g235m_path) as h:
        d = h[1].data
        wl_g235 = d['WAVELENGTH']
        fl_g235 = d['FLUX']

    with fits.open(g395m_path) as h:
        d = h[1].data
        wl_g395 = d['WAVELENGTH']
        fl_g395 = d['FLUX']

    # Apply v5 ghost correction to G235M
    fl_g235_corr = apply_v5_correction(wl_g235, fl_g235, 'G235M')

    # Smoothed versions
    fl_g235_s      = smooth(fl_g235, width=15)
    fl_g235_corr_s = smooth(fl_g235_corr, width=15)
    fl_g395_s      = smooth(fl_g395, width=10)

    # Overlap region: above G235M nominal limit (3.17 µm), below G395M max
    wl_nom_max = NOMINAL_LIMITS['G235M'][1]   # 3.17 µm
    wl_overlap_max = min(wl_g235.max(), wl_g395.max()) - 0.05  # ~5.25 µm
    wl_plot_min = wl_nom_max - 0.3  # Show a bit of the nominal G235M for context

    # Compute flux ratio in overlap (interpolate G395M onto G235M grid)
    interp_g395 = interp1d(wl_g395, fl_g395_s, bounds_error=False, fill_value=np.nan)
    g395_on_g235 = interp_g395(wl_g235)
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio_raw  = np.where(np.abs(g395_on_g235) > 1e-6,
                               fl_g235_s / g395_on_g235, np.nan)
        ratio_corr = np.where(np.abs(g395_on_g235) > 1e-6,
                               fl_g235_corr_s / g395_on_g235, np.nan)

    # --- Figure ---
    fig, axes = plt.subplots(2, 1, figsize=(14, 10), sharex=True)
    fig.suptitle(
        'UGC-5101 (z=0.039) — G235M Extended vs G395M Nominal Cross-Validation\n'
        'PID 2186 | Overlap Region: 3.17–5.25 µm | v5 Ghost Correction Applied',
        fontsize=13, fontweight='bold')

    ax, ax_r = axes

    # Top: both spectra
    ext_mask = wl_g235 >= wl_nom_max
    nom_mask = (wl_g235 >= wl_plot_min) & (wl_g235 < wl_nom_max)

    ax.plot(wl_g235[nom_mask], fl_g235_s[nom_mask],
            color='steelblue', lw=1.5, ls='-', alpha=0.6,
            label='G235M ext (nominal region, raw)')
    ax.plot(wl_g235[ext_mask], fl_g235_s[ext_mask],
            color='darkorange', lw=2.0,
            label='G235M ext (extended NRS2, raw)')
    ax.plot(wl_g235[ext_mask], fl_g235_corr_s[ext_mask],
            color='orangered', lw=1.8, ls='--',
            label='G235M ext (v5 ghost-corrected)')
    ax.plot(wl_g395[wl_g395 >= wl_nom_max - 0.3],
            fl_g395_s[wl_g395 >= wl_nom_max - 0.3],
            color='green', lw=2.0,
            label='G395M nom (standard pipeline, independent reference)')

    ax.axvline(wl_nom_max, color='gray', ls=':', lw=1.5,
               label='G235M nominal NRS2 limit (3.17 µm)')
    ax.fill_betweenx([0, 1], wl_nom_max, wl_overlap_max + 0.1,
                     transform=ax.get_xaxis_transform(),
                     alpha=0.05, color='orange', label='_overlap region')

    ax.set_ylabel('Flux (Jy)', fontsize=12)
    ax.legend(fontsize=10, loc='upper right')
    ax.grid(True, alpha=0.25)

    # Set y limits from data
    yvals = np.concatenate([fl_g235_s[np.isfinite(fl_g235_s) & (wl_g235 > wl_plot_min)],
                             fl_g395_s[np.isfinite(fl_g395_s) & (wl_g395 > wl_plot_min)]])
    ymax = np.nanpercentile(yvals, 98)
    ax.set_ylim(0, 1.35 * ymax)

    ax.text(0.75, 0.92, 'Cross-validation\nregion',
            transform=ax.transAxes, ha='center', fontsize=10, color='darkorange',
            style='italic')

    # Bottom panel: G235M / G395M ratio
    mask_ratio = (wl_g235 >= wl_nom_max) & (wl_g235 <= wl_overlap_max)
    ax_r.plot(wl_g235[mask_ratio], smooth(ratio_raw[mask_ratio], width=20),
              color='darkorange', lw=1.8, label='G235M raw / G395M nom')
    ax_r.plot(wl_g235[mask_ratio], smooth(ratio_corr[mask_ratio], width=20),
              color='orangered', lw=1.8, ls='--', label='G235M corrected / G395M nom')
    ax_r.axhline(1.0, color='black', lw=1.0, ls='--', label='Perfect agreement')
    ax_r.axvline(wl_nom_max, color='gray', ls=':', lw=1.5)

    # Shading for ±20% agreement band
    ax_r.fill_between([wl_nom_max, wl_overlap_max], [0.8, 0.8], [1.2, 1.2],
                      alpha=0.1, color='green', label='±20% band')

    ax_r.set_xlabel('Wavelength (µm)', fontsize=12)
    ax_r.set_ylabel('G235M / G395M ratio', fontsize=11)
    ax_r.set_ylim(0.3, 2.0)
    ax_r.legend(fontsize=10)
    ax_r.grid(True, alpha=0.25)
    ax_r.set_xlim(wl_plot_min, wl_overlap_max + 0.1)

    plt.tight_layout()
    outpath = os.path.join(OUT_DIR, 'ifu_v5_ugc5101_g235m_vs_g395m_xval.png')
    plt.savefig(outpath, dpi=200, bbox_inches='tight')
    plt.close()
    print(f'Saved: {outpath}')
    return True


# ── Main ──────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    print('Generating IFU v5 validation plots ...')

    print('\n[1] UGC-5101 G235M extended ...')
    plot_ugc5101_g235m()

    print('\n[2] Coefficient overview ...')
    plot_v5_coeffs_ifu_context()

    print('\n[3] SDSSJ0841 G140M extended ...')
    available = plot_sdssj0841_g140m(redshift_guess=2.0)
    if not available:
        print('    Will try again once Stage 3 completes.')

    print('\n[4] G395M cross-validation ...')
    plot_ugc5101_g395m_crossval()

    print('\nDone.')
