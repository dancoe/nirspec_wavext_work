"""
Parlanti comparison v5 — coefficient comparison plots.

Compares:
  - Parlanti et al. (2025) published calibration (black dashed)
  - FS v4 (3-source, degenerate — for context)
  - FS v5 (4-source joint-solve, degeneracy broken — primary result)

The v5 addition of NGC2506-G31 (G1V) as a 4th calibration source resolves
the k/α degeneracy seen in v4 (hot-star only). This plot shows the recovery
of k(λ) back to ~0.96 and α settling to Parlanti-consistent values.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from astropy.io import fits
from scipy.interpolate import interp1d

BASE       = '/Users/dcoe/NIRSpec/wavext'
PAR_DIR    = f'{BASE}/data/parlanti_repo/calibration_files'
FS_V4_DIR  = f'{BASE}/nirspec_wavext_work/plots/Parlanti/cal/fs_v4'
FS_V5_DIR  = f'{BASE}/results/v5'
OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))
os.makedirs(OUTPUT_DIR, exist_ok=True)


def load_parlanti(grating):
    filt = 'f100lp' if grating == 'G140M' else 'f170lp'
    fname = f'calibration_functions_{grating.lower()}_{filt}.fits'
    with fits.open(f'{PAR_DIR}/{fname}') as h:
        d = h[1].data
        return d['wavelength'], d['k'], d['alpha'], d['beta']


def load_v4_csv(grating, mode='fs'):
    """Load v4 FS CSV coefficients (if present)."""
    path = f'{FS_V4_DIR}/coeffs_fs_v4_{grating}.csv'
    if not os.path.exists(path):
        return None, None, None, None
    d = np.loadtxt(path, delimiter=',', skiprows=1)
    return d[:, 0], d[:, 1], d[:, 2], d[:, 3]


def load_v5_fits(grating):
    """Load v5 FITS calibration coefficients."""
    filt = 'f100lp' if grating == 'G140M' else 'f170lp'
    path = f'{FS_V5_DIR}/calib_v5_{grating.lower()}_{filt}.fits'
    if not os.path.exists(path):
        print(f'  WARNING: missing v5 FITS {path}')
        return None, None, None, None
    with fits.open(path) as h:
        d = h[1].data
        return d['WAVELENGTH'], d['K'], d['ALPHA'], d['BETA']


def plot_panel(ax, rax, grating, coeff_idx, coeff_name, ylims, ylabel):
    """
    Plot one panel: top=coefficient curves, bottom=ratio to Parlanti.
    coeff_idx: 0=k, 1=alpha, 2=beta
    """
    wp, kp, ap, bp = load_parlanti(grating)
    par_vals = [kp, ap, bp]
    pv = par_vals[coeff_idx]

    ax.plot(wp, pv, color='black', lw=1.0, ls='--', alpha=0.9,
            label='Parlanti et al. (2025)', zorder=5)

    # v4 FS (gray thin, for context)
    wv4, kv4, av4, bv4 = load_v4_csv(grating)
    if wv4 is not None:
        v4_vals = [kv4, av4, bv4][coeff_idx]
        ax.plot(wv4, v4_vals, color='tomato', lw=1.2, ls='-', alpha=0.7,
                label='FS v4 (3-star, degenerate)', zorder=4)
        fi4 = interp1d(wv4, v4_vals, bounds_error=False, fill_value=np.nan)
        ratio4 = fi4(wp) / np.where(np.abs(pv) > 1e-7, pv, np.nan)
        rax.plot(wp, ratio4, color='tomato', lw=1.2, ls='-', alpha=0.7)

    # v5 FS (primary result, thick blue)
    wv5, kv5, av5, bv5 = load_v5_fits(grating)
    if wv5 is not None:
        v5_vals = [kv5, av5, bv5][coeff_idx]
        # Small positive floor for log-scale safety
        v5_plot = np.where(v5_vals > 1e-7, v5_vals, 1e-7)
        ax.plot(wv5, v5_plot, color='royalblue', lw=2.0, ls='-',
                label='FS v5 (4-star, degeneracy broken)', zorder=6)
        fi5 = interp1d(wv5, v5_plot, bounds_error=False, fill_value=np.nan)
        ratio5 = fi5(wp) / np.where(np.abs(pv) > 1e-7, pv, np.nan)
        rax.plot(wp, ratio5, color='royalblue', lw=1.8, ls='-')

    ax.set_ylim(*ylims)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(grating, fontsize=13)
    ax.legend(fontsize=10, loc='upper right')
    ax.grid(True, alpha=0.25)

    rax.axhline(1.0, color='black', lw=0.9, ls=':')
    rax.set_yscale('log')
    rax.set_ylim(0.05, 20)
    for y_tick in [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]:
        rax.axhline(y_tick, color='gray', lw=0.3, ls=':', alpha=0.5)
    rax.yaxis.set_major_formatter(ticker.ScalarFormatter())
    rax.yaxis.set_major_locator(ticker.FixedLocator([0.1, 0.3, 1.0, 3.0, 10.0]))
    rax.yaxis.set_minor_locator(ticker.NullLocator())
    rax.set_ylabel('Ratio\n(v/Parlanti)', fontsize=9)
    rax.set_xlabel('Wavelength (µm)', fontsize=11)
    rax.grid(True, alpha=0.2)


def make_comparison_plot(coeff_name, coeff_idx, ylims_g140m, ylims_g235m, ylabel, outname):
    """Generate a 2×2 panel comparison (2 gratings × top/bottom)."""
    fig = plt.figure(figsize=(16, 10))
    hr  = [3, 1]
    gs  = fig.add_gridspec(2, 2, height_ratios=hr, hspace=0.07, wspace=0.25)
    ax_g140 = fig.add_subplot(gs[0, 0])
    rx_g140 = fig.add_subplot(gs[1, 0], sharex=ax_g140)
    ax_g235 = fig.add_subplot(gs[0, 1])
    rx_g235 = fig.add_subplot(gs[1, 1], sharex=ax_g235)

    coeff_labels = {
        'k':     r'Throughput $k(\lambda)$',
        'alpha': r'2nd-Order Ghost $\tilde{\alpha}(\lambda)$',
        'beta':  r'3rd-Order Ghost $\tilde{\beta}(\lambda)$',
    }
    fig.suptitle(
        f'Parlanti Comparison v5 — {coeff_labels[coeff_name]}\n'
        r'FS v5 (4-star NNLS: G191-B2B + P330E + J1743045 + NGC2506-G31) '
        'vs Parlanti et al. (2025) vs FS v4 (3-star degenerate)',
        fontsize=13, fontweight='bold')

    plot_panel(ax_g140, rx_g140, 'G140M', coeff_idx, coeff_name, ylims_g140m, ylabel)
    plot_panel(ax_g235, rx_g235, 'G235M', coeff_idx, coeff_name, ylims_g235m, ylabel)

    plt.tight_layout(rect=[0, 0, 1, 0.93])
    outpath = os.path.join(OUTPUT_DIR, outname)
    plt.savefig(outpath, dpi=200, bbox_inches='tight')
    plt.close()
    print(f'Saved: {outpath}')


def make_overview_plot():
    """3-row summary: k, α, β all together in one 3x2 figure."""
    fig = plt.figure(figsize=(16, 18))
    gs  = fig.add_gridspec(3, 2, hspace=0.40, wspace=0.25)

    gratings  = ['G140M', 'G235M']
    coeff_info = [
        (0, 'k',     r'$k(\lambda)$',     (0.0, 2.5), (0.0, 2.5)),
        (1, 'alpha', r'$\tilde{\alpha}$',  (-0.005, 0.30), (-0.005, 0.30)),
        (2, 'beta',  r'$\tilde{\beta}$',   (-0.002, 0.02), (-0.002, 0.02)),
    ]

    for row, (cidx, cname, ylabel, yl_140, yl_235) in enumerate(coeff_info):
        for col, grating in enumerate(gratings):
            ax = fig.add_subplot(gs[row, col])
            wp, kp, ap, bp = load_parlanti(grating)
            pv = [kp, ap, bp][cidx]
            ax.plot(wp, pv, color='black', lw=1.0, ls='--', alpha=0.9,
                    label='Parlanti (2025)', zorder=5)

            wv4, kv4, av4, bv4 = load_v4_csv(grating)
            if wv4 is not None:
                v4v = [kv4, av4, bv4][cidx]
                ax.plot(wv4, v4v, color='tomato', lw=1.2, ls='-', alpha=0.7,
                        label='FS v4 (3-star)', zorder=4)

            wv5, kv5, av5, bv5 = load_v5_fits(grating)
            if wv5 is not None:
                v5v = [kv5, av5, bv5][cidx]
                v5v_plot = np.where(v5v > 1e-7, v5v, 1e-7)
                ax.plot(wv5, v5v_plot, color='royalblue', lw=2.0, ls='-',
                        label='FS v5 (4-star)', zorder=6)

            ylims = yl_140 if grating == 'G140M' else yl_235
            ax.set_ylim(*ylims)
            ax.set_ylabel(ylabel, fontsize=12)
            ax.set_title(f'{grating} — {coeff_labels_short[cname]}', fontsize=12)
            if row == 0 and col == 0:
                ax.legend(fontsize=9, loc='upper right')
            ax.grid(True, alpha=0.25)
            if row == 2:
                ax.set_xlabel('Wavelength (µm)', fontsize=11)

    fig.suptitle(
        'Parlanti Comparison v5 — All Coefficients\n'
        r'(FS v4: degenerate 3-star  $\rightarrow$  FS v5: degeneracy broken with G1V)',
        fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    out = os.path.join(OUTPUT_DIR, 'parlanti_v5_comparison_overview.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out}')


coeff_labels_short = {'k': 'Throughput k', 'alpha': 'Ghost α', 'beta': 'Ghost β'}


if __name__ == '__main__':
    # Individual coefficient comparisons
    make_comparison_plot('k',     0, (0.0, 2.5), (0.0, 2.5),
                         r'$k(\lambda)$',
                         'comp_kappa_v5.png')
    make_comparison_plot('alpha', 1, (-0.005, 0.30), (-0.005, 0.30),
                         r'$\tilde{\alpha}(\lambda)$',
                         'comp_alpha_v5.png')
    make_comparison_plot('beta',  2, (-0.002, 0.02), (-0.002, 0.02),
                         r'$\tilde{\beta}(\lambda)$',
                         'comp_beta_v5.png')
    # Overview 3x2 summary
    make_overview_plot()
    print('\nAll comparison plots done.')
