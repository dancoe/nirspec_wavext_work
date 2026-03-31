"""
NIRSpec Wavelength Extension — IFU v7 Validation Plots

Generates:
  1. k/α/β coefficient plots for all IFU sources (all + hold-outs)
  2. Corrected IFU spectra vs CALSPEC truth (NRS1 + corrected NRS2)
  3. Hold-out residuals (fractional)
  4. IFU vs FS coefficient comparison

Outputs: reports/v7/ifu_v7/plots/
"""
import os
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy.io import fits
from scipy.interpolate import interp1d

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE        = '/Users/dcoe/NIRSpec/wavext'
IFU_DIR     = f'{BASE}/data/IFU'
DATA_DIR    = f'{BASE}/data'
CALSPEC_DIR = f'{BASE}/data/CALSPEC'
RESULTS_DIR = f'{BASE}/results/v7'
PLOT_DIR    = f'{BASE}/nirspec_wavext_work/reports/v7/ifu_v7/plots'
os.makedirs(PLOT_DIR, exist_ok=True)

C_ANG_S = 2.99792458e18  # Å/s

SOURCES = {
    '1537': {'name': 'G191-B2B',  'sptype': 'WD-DA.8', 'color': '#1f77b4', 'calspec': 'g191b2b_mod_012.fits',
             'ifu_dir': 'PID1537_G191-B2B'},
    '1538': {'name': 'P330E',     'sptype': 'G2V',     'color': '#ff7f0e', 'calspec': 'p330e_mod_008.fits',
             'ifu_dir': 'PID1538_P330E'},
    '1536': {'name': 'J1743045',  'sptype': 'A8III',   'color': '#2ca02c', 'calspec': '1743045_mod_007.fits',
             'ifu_dir': 'PID1536_J1743045'},
}

NRS2_LO  = {'g140m': 1.87,  'g235m': 3.15}
NRS2_HI  = {'g140m': 3.60,  'g235m': 5.50}
NRS1_LO  = {'g140m': 0.97,  'g235m': 1.66}
NRS1_HI  = {'g140m': 1.88,  'g235m': 3.16}

X1D_PATTERNS = {
    'g140m': ['*g140m*x1d.fits', '*G140M*x1d.fits', '*f100lp*x1d.fits'],
    'g235m': ['*g235m*x1d.fits', '*G235M*x1d.fits', '*f170lp*x1d.fits'],
}
X1D_NOM_PATTERNS = {
    'g140m': ['*g140m*x1d.fits', '*G140M*x1d.fits', '*f100lp*x1d.fits'],
    'g235m': ['*g235m*x1d.fits', '*G235M*x1d.fits', '*f170lp*x1d.fits'],
}


# ── Helpers ────────────────────────────────────────────────────────────────────
def load_spec(path, lo=None, hi=None):
    if not path or not os.path.exists(path):
        return np.array([]), np.array([])
    with fits.open(path) as h:
        ext = 'EXTRACT1D' if 'EXTRACT1D' in [e.name for e in h] else 1
        d  = h[ext].data
        wl = d['WAVELENGTH'].astype(float)
        fl = d['FLUX'].astype(float)
        dq = (d['DQ'].astype(int) if 'DQ' in d.dtype.names else np.zeros(len(wl), int))
    good = (wl > 0.3) & np.isfinite(fl) & (fl > 0) & (dq == 0)
    if lo is not None:
        good &= wl >= lo
    if hi is not None:
        good &= wl <= hi
    return wl[good], fl[good]


def calspec_jy(fname, wl_out):
    path = os.path.join(CALSPEC_DIR, fname)
    if not os.path.exists(path):
        return np.full_like(wl_out, np.nan)
    with fits.open(path) as h:
        d    = h[1].data
        wl_a = d['WAVELENGTH'].astype(float)
        flam = d['FLUX'].astype(float)
    fnu = flam * wl_a**2 / C_ANG_S * 1e23
    idx = np.argsort(wl_a)
    interp = interp1d(wl_a[idx] / 1e4, fnu[idx], bounds_error=False, fill_value=np.nan)
    return interp(wl_out)


def load_calib(grating, tag='all'):
    fname = os.path.join(RESULTS_DIR, f'calib_v7_ifu_{grating}_{tag}.fits')
    if not os.path.exists(fname):
        return None, None, None, None
    with fits.open(fname) as h:
        d = h[1].data
        return d['WAVELENGTH'], d['K'], d['ALPHA'], d['BETA']


def find_ifu_x1d(pid, grating, stage='stage3_ext'):
    info = SOURCES.get(pid)
    if info is None:
        return None
    base = os.path.join(IFU_DIR, info['ifu_dir'], stage)
    for pat in X1D_PATTERNS[grating]:
        hits = sorted(glob.glob(os.path.join(base, pat)))
        if hits:
            return hits[0]
    return None


def apply_correction(wl_obs, fl_obs, wl_cal, k, alpha, beta, calspec_fname):
    """Apply the Parlanti k/α/β correction to an IFU x1d spectrum."""
    from scipy.interpolate import interp1d as itp1d
    k_interp     = itp1d(wl_cal, k,     bounds_error=False, fill_value=(k[0],     k[-1]))
    alpha_interp = itp1d(wl_cal, alpha, bounds_error=False, fill_value=(alpha[0], alpha[-1]))
    beta_interp  = itp1d(wl_cal, beta,  bounds_error=False, fill_value=(beta[0],  beta[-1]))

    cs_interp    = itp1d(
        *_calspec_arrays(calspec_fname), bounds_error=False, fill_value=np.nan)

    fl_corr = np.full_like(fl_obs, np.nan)
    for i, lam in enumerate(wl_obs):
        ki = k_interp(lam)
        ai = alpha_interp(lam)
        bi = beta_interp(lam)
        f2 = float(cs_interp(lam / 2.0))
        f3 = float(cs_interp(lam / 3.0))
        if not (np.isfinite(ki) and ki > 0 and np.isfinite(f2) and np.isfinite(f3)):
            continue
        fl_corr[i] = (fl_obs[i] - ai * f2 - bi * f3) / ki
    return fl_corr


def _calspec_arrays(fname):
    """Return (wl_um, flux_jy) arrays for a CALSPEC model."""
    path = os.path.join(CALSPEC_DIR, fname)
    with fits.open(path) as h:
        d    = h[1].data
        wl_a = d['WAVELENGTH'].astype(float)
        flam = d['FLUX'].astype(float)
    fnu = flam * wl_a**2 / C_ANG_S * 1e23
    idx = np.argsort(wl_a)
    return wl_a[idx] / 1e4, fnu[idx]


# ── Plot 1: k, α, β coefficients ──────────────────────────────────────────────
def plot_coefficients(grating):
    """Plot k/α/β from all-source and held-out IFU solutions."""
    fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
    fig.suptitle(f'IFU v7 Parlanti Coefficients — {grating.upper()}', fontsize=13)

    wl_all, k_all, a_all, b_all = load_calib(grating, 'all')
    if wl_all is None:
        print(f'  No all-source result for {grating}, skipping coefficient plot.')
        return

    for ax, arr, ylab in zip(axes, [k_all, a_all, b_all], ['k(λ)', 'α(λ)', 'β(λ)']):
        ax.plot(wl_all, arr, 'k-', lw=2, label='All sources', zorder=5)
        ax.set_ylabel(ylab, fontsize=11)
        ax.axhline(0, color='gray', lw=0.5, ls='--')

    colors = [s['color'] for s in SOURCES.values()]
    for (pid, s), col in zip(SOURCES.items(), colors):
        tag = f'holdout_{pid}'
        wl, k, a, b = load_calib(grating, tag)
        if wl is None:
            continue
        lbl = f'Hold-out {s["name"]}'
        for ax, arr in zip(axes, [k, a, b]):
            ax.plot(wl, arr, '--', color=col, lw=1.2, alpha=0.7, label=lbl)

    axes[0].legend(fontsize=8, ncol=2, loc='upper right')
    axes[-1].set_xlabel('Wavelength (µm)', fontsize=11)
    axes[-1].set_xlim(NRS2_LO[grating] - 0.05, NRS2_HI[grating] + 0.05)

    for ax in axes:
        ax.grid(True, alpha=0.3)

    fname = os.path.join(PLOT_DIR, f'ifu_v7_coefficients_{grating}.png')
    fig.tight_layout()
    fig.savefig(fname, dpi=150)
    plt.close(fig)
    print(f'  Saved: {os.path.relpath(fname, BASE)}')


# ── Plot 2: Hold-out residuals ─────────────────────────────────────────────────
def plot_holdout_residuals(grating):
    """For each source: apply hold-out coefficients to the held-out spectrum and
    compare to CALSPEC truth, showing fractional residual."""
    wl_all, k_all, a_all, b_all = load_calib(grating, 'all')
    if wl_all is None:
        return

    fig, axes = plt.subplots(len(SOURCES), 1, figsize=(10, 4 * len(SOURCES)), sharex=True)
    if len(SOURCES) == 1:
        axes = [axes]
    fig.suptitle(f'IFU v7 Hold-Out Residuals — {grating.upper()}', fontsize=13)

    for ax, (pid, s) in zip(axes, SOURCES.items()):
        wl_cal, k, alpha, beta = load_calib(grating, f'holdout_{pid}')
        if wl_cal is None:
            ax.set_title(f'{s["name"]}: no hold-out result')
            continue

        x1d_path = find_ifu_x1d(pid, grating)
        wl_obs, fl_obs = load_spec(x1d_path, lo=NRS2_LO[grating], hi=NRS2_HI[grating])
        if len(wl_obs) < 10:
            ax.set_title(f'{s["name"]}: no data')
            continue

        # Apply correction
        fl_corr = apply_correction(wl_obs, fl_obs, wl_cal, k, alpha, beta, s['calspec'])
        fl_truth = calspec_jy(s['calspec'], wl_obs)

        good = np.isfinite(fl_corr) & np.isfinite(fl_truth) & (fl_truth > 0)
        if good.sum() < 5:
            ax.set_title(f'{s["name"]}: correction failed')
            continue

        resid = (fl_corr[good] - fl_truth[good]) / fl_truth[good] * 100

        ax.plot(wl_obs[good], resid, color=s['color'], lw=1, alpha=0.7)
        ax.axhline(0, color='k', lw=1)
        ax.axhspan(-10, 10, alpha=0.1, color='green', label='±10%')
        med_resid = np.nanmedian(resid)
        ax.set_title(f'{s["name"]} ({s["sptype"]}) — median residual: {med_resid:.1f}%',
                     fontsize=10)
        ax.set_ylabel('(corrected − truth) / truth (%)', fontsize=9)
        ax.set_ylim(-50, 50)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)

    axes[-1].set_xlabel('Wavelength (µm)', fontsize=11)
    axes[-1].set_xlim(NRS2_LO[grating] - 0.05, NRS2_HI[grating] + 0.05)

    fname = os.path.join(PLOT_DIR, f'ifu_v7_holdout_residuals_{grating}.png')
    fig.tight_layout()
    fig.savefig(fname, dpi=150)
    plt.close(fig)
    print(f'  Saved: {os.path.relpath(fname, BASE)}')


# ── Plot 3: Full spectrum (NRS1 nom + NRS2 corrected vs CALSPEC) ───────────────
def plot_full_spectrum(grating):
    """Show full IFU spectrum (NRS1 nominal + corrected NRS2) vs CALSPEC model."""
    wl_all, k_all, a_all, b_all = load_calib(grating, 'all')
    if wl_all is None:
        return

    fig, axes = plt.subplots(len(SOURCES), 1, figsize=(12, 4 * len(SOURCES)), sharex=True)
    if len(SOURCES) == 1:
        axes = [axes]
    fig.suptitle(f'IFU v7 Full Spectrum — {grating.upper()}', fontsize=13)

    for ax, (pid, s) in zip(axes, SOURCES.items()):
        # NRS1 nominal spectrum (from stage3, not stage3_ext)
        wl_n1, fl_n1 = load_spec(find_ifu_x1d(pid, grating, 'stage3'),
                                  lo=NRS1_LO[grating], hi=NRS1_HI[grating])
        # NRS2 extended spectrum (from stage3_ext)
        wl_ext, fl_ext = load_spec(find_ifu_x1d(pid, grating),
                                   lo=NRS2_LO[grating], hi=NRS2_HI[grating])

        cs_wl = np.linspace(NRS1_LO[grating], NRS2_HI[grating], 1000)
        cs_fl = calspec_jy(s['calspec'], cs_wl)

        # Apply correction to NRS2
        fl_corr = apply_correction(wl_ext, fl_ext, wl_all, k_all, a_all, b_all, s['calspec'])

        ax.plot(cs_wl, cs_fl * 1e3, 'gray', lw=1.5, ls='--', label='CALSPEC truth', alpha=0.8)
        if len(wl_n1) > 0:
            ax.plot(wl_n1, fl_n1 * 1e3, color=s['color'], lw=1.5, label='NRS1 nominal')
        ax.plot(wl_ext, fl_corr * 1e3, color=s['color'], lw=1.5, ls=':', label='NRS2 corr (v7)')
        # Also show raw (uncorrected)
        ax.plot(wl_ext, fl_ext * 1e3, color='red', lw=0.8, alpha=0.5, label='NRS2 raw')

        ax.set_ylabel('Flux (mJy)', fontsize=9)
        ax.set_title(f'{s["name"]} ({s["sptype"]}, PID {pid})', fontsize=10)
        ax.legend(fontsize=8, ncol=2)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(NRS1_LO[grating] - 0.05, NRS2_HI[grating] + 0.05)

    axes[-1].set_xlabel('Wavelength (µm)', fontsize=11)

    fname = os.path.join(PLOT_DIR, f'ifu_v7_full_spectrum_{grating}.png')
    fig.tight_layout()
    fig.savefig(fname, dpi=150)
    plt.close(fig)
    print(f'  Saved: {os.path.relpath(fname, BASE)}')


# ── Plot 4: IFU vs FS coefficient comparison ──────────────────────────────────
def plot_ifu_vs_fs(grating):
    """Overplot IFU and FS v7 all-source k/α/β on the same axes."""
    wl_ifu, k_ifu, a_ifu, b_ifu = load_calib(grating, 'all')  # IFU
    # Load FS result from results/v7/
    fs_fname = os.path.join(RESULTS_DIR, f'calib_v7_fs_{grating}_all.fits')
    if not os.path.exists(fs_fname):
        print(f'  No FS v7 result for {grating}, skipping IFU vs FS plot.')
        return
    with fits.open(fs_fname) as h:
        d = h[1].data
        wl_fs, k_fs, a_fs, b_fs = d['WAVELENGTH'], d['K'], d['ALPHA'], d['BETA']

    if wl_ifu is None:
        print(f'  No IFU v7 result for {grating}, skipping IFU vs FS plot.')
        return

    fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
    fig.suptitle(f'IFU v7 vs FS v7 Coefficients — {grating.upper()}', fontsize=13)

    for ax, ifu_arr, fs_arr, ylab in zip(
            axes,
            [k_ifu, a_ifu, b_ifu],
            [k_fs,  a_fs,  b_fs],
            ['k(λ)', 'α(λ)', 'β(λ)']):
        ax.plot(wl_ifu, ifu_arr, 'b-',  lw=2, label='IFU v7 (all sources)')
        ax.plot(wl_fs,  fs_arr,  'r--', lw=2, label='FS v7 (all sources)')
        ax.set_ylabel(ylab, fontsize=11)
        ax.axhline(0, color='gray', lw=0.5, ls='--')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

    axes[-1].set_xlabel('Wavelength (µm)', fontsize=11)
    axes[-1].set_xlim(NRS2_LO[grating] - 0.05, NRS2_HI[grating] + 0.05)

    fname = os.path.join(PLOT_DIR, f'ifu_v7_vs_fs_{grating}.png')
    fig.tight_layout()
    fig.savefig(fname, dpi=150)
    plt.close(fig)
    print(f'  Saved: {os.path.relpath(fname, BASE)}')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='IFU v7 validation plots')
    parser.add_argument('--grating', default='all',
                        help='g140m, g235m, or all')
    args = parser.parse_args()

    gratings = ['g140m', 'g235m'] if args.grating == 'all' else [args.grating.lower()]

    for g in gratings:
        print(f'\n--- Generating plots for {g.upper()} ---')
        plot_coefficients(g)
        plot_holdout_residuals(g)
        plot_full_spectrum(g)
        plot_ifu_vs_fs(g)

    print(f'\nAll IFU v7 plots saved to: {os.path.relpath(PLOT_DIR, BASE)}')
