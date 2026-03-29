import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import interp1d
import matplotlib.ticker as ticker

# Paths
BASE        = '/Users/dcoe/NIRSpec/wavext'
PAR_DIR     = f'{BASE}/data/parlanti_repo/calibration_files'
FS_V2_DIR   = f'{BASE}/nirspec_wavext_work/plots/Parlanti/cal/fs_v2'
IFU_V2_DIR  = f'{BASE}/nirspec_wavext_work/plots/Parlanti/cal/ifu_v2'
OUTPUT_DIR  = f'{BASE}/nirspec_wavext_work/reports/329_parlanti-comparison_v2'
os.makedirs(OUTPUT_DIR, exist_ok=True)

def load_par_fits(grating):
    fname = f'calibration_functions_{grating.lower()}_{"f100lp" if grating=="G140M" else "f170lp"}.fits'
    with fits.open(f'{PAR_DIR}/{fname}') as h:
        d = h[1].data
        return d['wavelength'], d['k'], d['alpha'], d['beta']

def load_csv(path):
    if not os.path.exists(path): return None, None, None, None
    d = np.loadtxt(path, delimiter=',', skiprows=1)
    return d[:, 0], d[:, 1], d[:, 2], d[:, 3]

def plot_comparison(coeff_type):
    # coeff_type: 'kappa', 'alpha', 'beta'
    # Use GridSpec to have different height ratios (3:1) for values vs residuals
    fig = plt.figure(figsize=(15, 12))
    gs = fig.add_gridspec(2, 2, height_ratios=[3, 1], hspace=0.08)
    
    axes_top = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])]
    axes_bot = [fig.add_subplot(gs[1, 0], sharex=axes_top[0]), fig.add_subplot(gs[1, 1], sharex=axes_top[1])]
    
    titles = {'kappa': r'Throughput $k(\lambda)$', 
              'alpha': r'2nd Order Ghost $\tilde{\alpha}(\lambda)$', 
              'beta':  r'3rd Order Ghost $\tilde{\beta}(\lambda)$'}
    ylabel = {'kappa': 'k', 'alpha': r'$\tilde{\alpha}$', 'beta': r'$\tilde{\beta}$'}
    
    fig.suptitle(f'329 v2 Comparison — {titles[coeff_type]}', fontsize=18)

    gratings = ['G140M', 'G235M']
    
    for i, g in enumerate(gratings):
        ax = axes_top[i]
        rax = axes_bot[i]
        
        # Parlanti (Black Thin)
        wp, kp, ap, bp = load_par_fits(g)
        vals_p = {'kappa': kp, 'alpha': ap, 'beta': bp}[coeff_type]
        ax.plot(wp, vals_p, color='black', lw=0.8, alpha=0.9, label='Parlanti et al. (2025)')
        
        # FS v2 (Red)
        wf, kf, af, bf = load_csv(f'{FS_V2_DIR}/coeffs_fs_v2_{g}.csv')
        if wf is not None:
            vals_f = {'kappa': kf, 'alpha': af, 'beta': bf}[coeff_type]
            ax.plot(wf, vals_f, color='red', lw=1.5, label='FS v2 (329)')
            
            # Interpolate onto Parlanti grid to get residuals
            f_interp = interp1d(wf, vals_f, bounds_error=False, fill_value=np.nan)
            vals_f_p = f_interp(wp)
            ratio_f  = vals_f_p / vals_p
            rax.plot(wp, ratio_f, color='red', lw=1.2, alpha=0.8)
            
        # IFU v2 (Cyan)
        wi, ki, ai, bi = load_csv(f'{IFU_V2_DIR}/coeffs_ifu_v2_{g}.csv')
        if wi is not None:
            vals_i = {'kappa': ki, 'alpha': ai, 'beta': bi}[coeff_type]
            ax.plot(wi, vals_i, color='cyan', lw=1.5, label='IFU v2 (329)')
            
            # Interpolate onto Parlanti grid to get residuals
            i_interp = interp1d(wi, vals_i, bounds_error=False, fill_value=np.nan)
            vals_i_p = i_interp(wp)
            ratio_i  = vals_i_p / vals_p
            rax.plot(wp, ratio_i, color='cyan', lw=1.2, alpha=0.8)

        if coeff_type in ['alpha', 'beta']:
            ax.set_yscale('log')
            ax.set_ylim(1e-4, 0.5)
        else:
            ax.set_ylim(0, 2.5)

        ax.set_title(f'{g}', fontsize=14)
        if i == 0:
            ax.set_ylabel(ylabel[coeff_type], fontsize=12)
        ax.grid(True, alpha=0.2, which='both')
        ax.legend(fontsize=11)
        
        # Residual Styling
        rax.axhline(1.0, color='black', lw=0.8, ls=':')
        rax.set_yscale('log')
        
        # Fixed Y-range [0.5, 2.0]
        rax.set_ylim(0.5, 2.0)
            
        # Custom log ticks in linear style
        ticks = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
            
        rax.yaxis.set_major_formatter(ticker.ScalarFormatter())
        rax.yaxis.set_major_locator(ticker.FixedLocator(ticks))
        rax.yaxis.set_minor_locator(ticker.NullLocator()) # No minor ticks
        
        if i == 0:
            rax.set_ylabel('Ratio to Parlanti', fontsize=11)
        rax.set_xlabel('Wavelength [µm]', fontsize=12)
        rax.grid(True, alpha=0.2, which='both')
        
    plt.tight_layout()
    outpath = f'{OUTPUT_DIR}/comp_{coeff_type}_v2.png'
    plt.savefig(outpath, dpi=200, bbox_inches='tight')
    plt.close()
    print(f'Saved: {outpath}')

if __name__ == '__main__':
    for ct in ['kappa', 'alpha', 'beta']:
        plot_comparison(ct)
