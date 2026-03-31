"""
Parlanti comparison v4 — coefficient comparison plots.

Compares:
  - Parlanti et al. (2025) published calibration (black thin)
  - IFU v4 (cyan)
  - FS v4 (red)

Unlike v3, both IFU v4 and FS v4 now solve for k, alpha, and beta independently.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from astropy.io import fits
from scipy.interpolate import interp1d

BASE       = '/Users/dcoe/NIRSpec/wavext'
PAR_DIR    = f'{BASE}/data/parlanti_repo/calibration_files'
FS_V4_DIR  = f'{BASE}/nirspec_wavext_work/plots/Parlanti/cal/fs_v4'
IFU_V4_DIR = f'{BASE}/nirspec_wavext_work/plots/Parlanti/cal/ifu_v4'
OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))
os.makedirs(OUTPUT_DIR, exist_ok=True)


def load_par_fits(grating):
    fname = (f'calibration_functions_{grating.lower()}_'
             f'{"f100lp" if grating == "G140M" else "f170lp"}.fits')
    with fits.open(f'{PAR_DIR}/{fname}') as h:
        d = h[1].data
        return d['wavelength'], d['k'], d['alpha'], d['beta']


def load_csv(path):
    if not os.path.exists(path):
        print(f'  WARNING: missing {path}')
        return None, None, None, None
    d = np.loadtxt(path, delimiter=',', skiprows=1)
    return d[:, 0], d[:, 1], d[:, 2], d[:, 3]


def plot_comparison(coeff_type):
    titles = {
        'kappa': r'Throughput $k(\lambda)$',
        'alpha': r'2nd Order Ghost $\tilde{\alpha}(\lambda)$',
        'beta':  r'3rd Order Ghost $\tilde{\beta}(\lambda)$',
    }
    ylabels = {
        'kappa': 'k',
        'alpha': r'$\tilde{\alpha}$',
        'beta':  r'$\tilde{\beta}$',
    }

    fig = plt.figure(figsize=(15, 12))
    gs  = fig.add_gridspec(2, 2, height_ratios=[3, 1], hspace=0.12)
    axes_top = [fig.add_subplot(gs[0, 0]),
                fig.add_subplot(gs[0, 1])]
    axes_bot = [fig.add_subplot(gs[1, 0], sharex=axes_top[0]),
                fig.add_subplot(gs[1, 1], sharex=axes_top[1])]

    fig.suptitle(
        f'329 v4 Comparison — {titles[coeff_type]}\n'
        f'IFU v4 (NNLS 3-source) vs FS v4 (NNLS 3-source) vs Parlanti (2025)',
        fontsize=15)

    for i, g in enumerate(('G140M', 'G235M')):
        ax  = axes_top[i]
        rax = axes_bot[i]
        idx = {'kappa': 1, 'alpha': 2, 'beta': 3}[coeff_type]

        # Parlanti
        wp, kp, ap, bp = load_par_fits(g)
        vals_p = [kp, ap, bp][idx - 1]
        ax.plot(wp, vals_p, color='black', lw=0.9, ls='--', alpha=0.9,
                label='Parlanti et al. (2025)')

        # FS v4 (red)
        wf, kf, af, bf = load_csv(f'{FS_V4_DIR}/coeffs_fs_v4_{g}.csv')
        if wf is not None:
            vals_f = [kf, af, bf][idx - 1]
            ax.plot(wf, vals_f, color='red', lw=1.5, label='FS v4 (NNLS)')
            fi = interp1d(wf, vals_f, bounds_error=False, fill_value=np.nan)
            ratio_f = fi(wp) / vals_p
            rax.plot(wp, ratio_f, color='red', lw=1.2, alpha=0.8)

        # IFU v4 (cyan)
        wi, ki, ai, bi = load_csv(f'{IFU_V4_DIR}/coeffs_ifu_v4_{g}.csv')
        if wi is not None:
            vals_i = [ki, ai, bi][idx - 1]
            ax.plot(wi, vals_i, color='cyan', lw=1.5, label='IFU v4 (NNLS)')
            ii = interp1d(wi, vals_i, bounds_error=False, fill_value=np.nan)
            ratio_i = ii(wp) / vals_p
            rax.plot(wp, ratio_i, color='cyan', lw=1.2, alpha=0.8)

        if coeff_type in ('alpha', 'beta'):
            ax.set_yscale('linear') # Linear for alpha/beta comparison in v4 as offsets are larger
            if coeff_type == 'alpha':
                ax.set_ylim(-0.01, 0.4)
            else:
                ax.set_ylim(-0.005, 0.1)
        else:
            ax.set_ylim(0, 2.5)

        ax.set_title(g, fontsize=14)
        if i == 0:
            ax.set_ylabel(ylabels[coeff_type], fontsize=12)
        ax.grid(True, alpha=0.2, which='both')
        ax.legend(fontsize=11)

        rax.axhline(1.0, color='black', lw=0.8, ls=':')
        rax.set_yscale('log')
        rax.set_ylim(0.1, 10.0) # Larger range for v4 residuals as alpha can be 6x higher
        ticks = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]
        rax.yaxis.set_major_formatter(ticker.ScalarFormatter())
        rax.yaxis.set_major_locator(ticker.FixedLocator(ticks))
        rax.yaxis.set_minor_locator(ticker.NullLocator())
        if i == 0:
            rax.set_ylabel('Ratio to Parlanti', fontsize=11)
        rax.set_xlabel('Wavelength [µm]', fontsize=12)
        rax.grid(True, alpha=0.2, which='both')

    plt.tight_layout()
    out = os.path.join(OUTPUT_DIR, f'comp_{coeff_type}_v4.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out}')


if __name__ == '__main__':
    for ct in ('kappa', 'alpha', 'beta'):
        plot_comparison(ct)
    print('Done.')
