"""
Calibration comparison plots for Parlanti et al. (2025) NIRSpec IFU
higher-order contamination correction.

Generates:
  1. k(λ), α(λ), β(λ) vs wavelength — our derivation vs Parlanti reference
  2. Residuals S_obs/S_cal before and after correction, per source and grating
  3. Summary residual panel across all three standards
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import astropy.io.fits as fits
from scipy.interpolate import interp1d

# ── Paths ────────────────────────────────────────────────────────────────────
BASE         = '/Users/dcoe/NIRSpec/wavext'
IFU_DIR      = f'{BASE}/data/IFU'
CALSPEC_DIR  = f'{BASE}/data/CALSPEC'
PARLANTI_DIR = f'{BASE}/data/parlanti_repo/calibration_files'
CALIB_DIR    = f'{BASE}/results/calibration'
PLOT_DIR     = f'{BASE}/results/plots'
os.makedirs(PLOT_DIR, exist_ok=True)

C_ANG_S = 2.99792458e18  # Å/s

# ── Source configuration ─────────────────────────────────────────────────────
SOURCES = {
    '1537': ('G191-B2B', 'g191b2b_mod_012.fits', 'PID1537_G191-B2B', 'blue'),
    '1538': ('P330E',    'p330e_mod_008.fits',    'PID1538_P330E',    'goldenrod'),
    '1536': ('J1743045', '1743045_mod_007.fits',  'PID1536_J1743045', 'forestgreen'),
}

GRATING_CONFIGS = {
    'g140m_f100lp': {
        'label':          'G140M/F100LP',
        'short':          'G140M',
        'color':          'blue',
        'color_ext':      'cornflowerblue',
        'x1d_stem':       'f100lp_g140m-f100lp_x1d.fits',
        'stage':          'stage3_ext',
        'wl_range':       (0.97, 3.60),
        'wl_range_core':  (1.00, 1.80),  # NRS1 core, for residual zoom
        'calib_ours':     f'{CALIB_DIR}/calib_g140m_f100lp.fits',
        'calib_parlanti': f'{PARLANTI_DIR}/calibration_functions_g140m_f100lp.fits',
        'pids': ['1537', '1538', '1536'],
    },
    'g235m_f170lp': {
        'label':          'G235M/F170LP',
        'short':          'G235M',
        'color':          'goldenrod',
        'color_ext':      'palegoldenrod',
        'x1d_stem':       'f170lp_g235m-f170lp_x1d.fits',
        'stage':          'stage3_ext',
        'wl_range':       (1.66, 5.50),
        'wl_range_core':  (1.70, 3.20),
        'calib_ours':     f'{CALIB_DIR}/calib_g235m_f170lp.fits',
        'calib_parlanti': f'{PARLANTI_DIR}/calibration_functions_g235m_f170lp.fits',
        'pids': ['1537', '1538', '1536'],
    },
    'g395m_f290lp': {
        'label':          'G395M/F290LP',
        'short':          'G395M',
        'color':          'red',
        'color_ext':      'lightcoral',
        'x1d_stem':       'f290lp_g395m-f290lp_x1d.fits',
        'stage':          'stage3_ext',
        'wl_range':       (2.87, 8.50),
        'wl_range_core':  (2.90, 5.10),
        'calib_ours':     f'{CALIB_DIR}/calib_g395m_f290lp.fits',
        'calib_parlanti': f'{PARLANTI_DIR}/calibration_functions_g395m_f290lp.fits',
        'pids': ['1537', '1538', '1536'],
    },
}

# ── Utility functions ─────────────────────────────────────────────────────────
def flam_to_jy(wl_ang, f_flam):
    fnu_cgs = f_flam * wl_ang**2 / C_ANG_S
    return fnu_cgs * 1e23

def load_calspec(fname):
    h = fits.open(f'{CALSPEC_DIR}/{fname}')
    d = h[1].data
    wl_ang = d['WAVELENGTH'].astype(float)
    fj     = flam_to_jy(wl_ang, d['FLUX'].astype(float))
    h.close()
    order   = np.argsort(wl_ang / 1e4)
    wl_um   = wl_ang[order] / 1e4
    return wl_um, fj[order]

def load_x1d(path):
    h   = fits.open(path)
    d   = h[1].data
    wl  = d['WAVELENGTH'].astype(float)
    fl  = d['FLUX'].astype(float)
    dq  = d['DQ'].astype(int)
    h.close()
    return wl, fl, dq

def load_calib(path):
    """Load k/alpha/beta FITS. Returns (wl_um, k, alpha, beta)."""
    h = fits.open(path)
    d = h[1].data
    # Handle both our format (WAVELENGTH, K, ALPHA, BETA) and Parlanti's
    wl = d['WAVELENGTH'].astype(float)
    k   = d['K'].astype(float)
    # Parlanti uses lowercase column names, try both
    try:
        alpha = d['ALPHA'].astype(float)
        beta  = d['BETA'].astype(float)
    except KeyError:
        alpha = d['alpha'].astype(float)
        beta  = d['beta'].astype(float)
    h.close()
    # Ensure wavelength is in µm (Parlanti uses µm, we use µm too)
    if wl.max() > 100:   # probably in Angstroms
        wl = wl / 1e4
    return wl, k, alpha, beta

def interp_safe(wl_src, vals_src, wl_tgt):
    ok = np.isfinite(vals_src)
    if ok.sum() < 2:
        return np.full(len(wl_tgt), np.nan)
    f = interp1d(wl_src[ok], vals_src[ok],
                 kind='linear', bounds_error=False, fill_value=np.nan)
    return f(wl_tgt)

def apply_calib(wl_obs, flux_obs, wl_cal, k, alpha, beta):
    """Apply Parlanti correction: (flux - alpha*flux_half - beta*flux_third) / k."""
    k_terp     = interp_safe(wl_cal, k,     wl_obs)
    alpha_terp = interp_safe(wl_cal, alpha, wl_obs)
    beta_terp  = interp_safe(wl_cal, beta,  wl_obs)

    # Fill NaNs in coefficients (default to 1.0 for k, 0.0 for contamination)
    k_terp     = np.nan_to_num(k_terp, nan=np.nanmedian(k))
    alpha_terp = np.nan_to_num(alpha_terp, nan=0.0)
    beta_terp  = np.nan_to_num(beta_terp, nan=0.0)

    # Interpolate observed flux at λ/2 and λ/3
    f_obs_interp = interp1d(wl_obs, flux_obs, kind='linear',
                            bounds_error=False, fill_value=np.nan)
    flux_half  = f_obs_interp(wl_obs / 2)
    flux_third = f_obs_interp(wl_obs / 3)

    # Crucial: if we don't have coverage at half/third wavelength, assume 0 contamination
    flux_half  = np.nan_to_num(flux_half, nan=0.0)
    flux_third = np.nan_to_num(flux_third, nan=0.0)

    flux_corr = (flux_obs - alpha_terp * flux_half - beta_terp * flux_third) / k_terp
    return flux_corr


# ═══════════════════════════════════════════════════════════════════════════════
# PLOT 1: k, α, β vs wavelength — comparison with Parlanti
# ═══════════════════════════════════════════════════════════════════════════════
def plot_coefficients():
    num_gratings = len(GRATING_CONFIGS)
    fig, axes = plt.subplots(3, num_gratings, figsize=(7 * num_gratings, 10))
    fig.suptitle('Parlanti k(λ), α(λ), β(λ) calibration coefficients\n'
                 'Ours (blue) vs Parlanti reference (orange)', fontsize=13)

    coeff_names = ['k', 'alpha', 'beta']
    coeff_labels = [r'$k(\lambda)$', r'$\alpha(\lambda)$', r'$\beta(\lambda)$']
    coeff_ylims = [(0, 1.3), (-0.01, 0.25), (-0.001, 0.025)]

    for col, (gkey, gcfg) in enumerate(GRATING_CONFIGS.items()):
        try:
            wl_ours, k_ours, alpha_ours, beta_ours = load_calib(gcfg['calib_ours'])
        except Exception as e:
            print(f"  SKIPPING coefficients for {gkey}: {e}")
            for row in range(3):
                ax = axes[row, col] if num_gratings > 1 else axes[row]
                ax.text(0.5, 0.5, f"Missing Data\n{gkey}", ha='center', transform=ax.transAxes)
                ax.set_title(gcfg['label'], fontsize=12)
            continue

        try:
            wl_par, k_par, alpha_par, beta_par = load_calib(gcfg['calib_parlanti'])
        except Exception:
            wl_par = k_par = alpha_par = beta_par = None

        vals_ours = [k_ours, alpha_ours, beta_ours]
        vals_par  = [k_par,  alpha_par,  beta_par] if wl_par is not None else [None]*3

        for row, (coeff, label, ylim) in enumerate(zip(coeff_names, coeff_labels, coeff_ylims)):
            ax = axes[row, col] if num_gratings > 1 else axes[row]
            vo = vals_ours[row]
            vp = vals_par[row]

            # Mask NaN and clip to plot range
            ok = np.isfinite(vo)
            ax.plot(wl_ours[ok], vo[ok], lw=1.2, color=gcfg['color'],
                    alpha=0.8, label='Ours (JWST 1.20.2)')

            if vp is not None:
                ok_p = np.isfinite(vp)
                ax.plot(wl_par[ok_p], vp[ok_p], lw=1.2, color=gcfg['color'],
                        alpha=0.4, ls='--', label='Parlanti 2025')

            ax.set_ylabel(label, fontsize=11)
            ax.set_xlim(*gcfg['wl_range'])
            ax.set_ylim(*ylim)
            ax.axhline(0, color='k', lw=0.5, ls=':')
            ax.grid(True, alpha=0.3)
            if row == 0:
                ax.set_title(gcfg['label'], fontsize=12)
                ax.legend(fontsize=8, loc='upper right')
            if row == 2:
                ax.set_xlabel('Wavelength (µm)', fontsize=10)

    plt.tight_layout()
    outpath = f'{PLOT_DIR}/fig1_calibration_coefficients.png'
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved: {outpath}')


# ═══════════════════════════════════════════════════════════════════════════════
# PLOT 2: S_obs/S_cal ratio before and after calibration per PID (all gratings)
# ═══════════════════════════════════════════════════════════════════════════════
def plot_residuals():
    """Plot residuals for each PID, showing all available gratings in one figure."""
    print("\n--- Plot 2: Residuals before/after correction (one figure per PID) ---")
    
    for pid in SOURCES.keys():
        name, cs_fname, pid_dir, source_color = SOURCES[pid]
        
        fig, axes = plt.subplots(2, 1, figsize=(12, 10), sharex=True)
        fig.suptitle(f'PID {pid} ({name}) — Flux calibration quality\n'
                     f'S_obs / S_CALSPEC before and after Parlanti correction',
                     fontsize=14)

        ax_before = axes[0]
        ax_after  = axes[1]
        
        # Determine global wavelength range for this PID's plot
        all_wl_min, all_wl_max = 100, 0

        for gkey, gcfg in GRATING_CONFIGS.items():
            if pid not in gcfg['pids']:
                continue
                
            try:
                wl_cal, k_cal, alpha_cal, beta_cal = load_calib(gcfg['calib_ours'])
            except Exception as e:
                print(f'  WARNING: cannot load calib for {gkey}: {e}')
                continue

            x1d_path = f'{IFU_DIR}/{pid_dir}/{gcfg["stage"]}/{gcfg["x1d_stem"]}'
            if not os.path.exists(x1d_path):
                print(f'  WARNING: missing x1d for PID {pid} grating {gkey}')
                continue

            wl_obs, fl_obs, dq_obs = load_x1d(x1d_path)
            cs_wl, cs_fj = load_calspec(cs_fname)
            cs_interp = interp1d(cs_wl, cs_fj, bounds_error=False, fill_value=np.nan)

            # Mask bad pixels
            bad = (dq_obs > 0) | ~np.isfinite(fl_obs) | (fl_obs <= 0)
            fl_masked = fl_obs.copy()
            fl_masked[bad] = np.nan

            # Ratio before correction
            cs_at_obs = cs_interp(wl_obs)
            ratio_before = fl_masked / cs_at_obs

            # Apply calibration correction
            fl_corr = apply_calib(wl_obs, fl_masked, wl_cal, k_cal, alpha_cal, beta_cal)
            ratio_after = fl_corr / cs_at_obs

            # Median filter for display
            from scipy.ndimage import median_filter
            def smooth_ratio(r):
                ok = np.isfinite(r)
                if ok.sum() < 10: return r
                tmp = r.copy()
                tmp[~ok] = np.nan
                # Fill NaNs with median of valid for filtering, then mask again
                fill_val = np.nanmedian(tmp[ok])
                sm = median_filter(np.where(ok, tmp, fill_val), size=11)
                sm[~ok] = np.nan
                return sm

            ratio_before_sm = smooth_ratio(ratio_before)
            ratio_after_sm = smooth_ratio(ratio_after)

            label = f'{gcfg["short"]}'
            ax_before.plot(wl_obs, ratio_before_sm, lw=1.5, color=gcfg['color'],
                           alpha=0.6, label=label)
            ax_after.plot(wl_obs, ratio_after_sm, lw=1.8, color=gcfg['color'],
                          alpha=0.9, label=label)
            
            all_wl_min = min(all_wl_min, gcfg['wl_range'][0])
            all_wl_max = max(all_wl_max, gcfg['wl_range'][1])

        for ax, title in [(ax_before, 'Before calibration correction'),
                          (ax_after,  'After calibration correction')]:
            ax.axhline(1.0, color='k', lw=1.0, ls='-', alpha=0.8)
            ax.axhspan(0.95, 1.05, color='green', alpha=0.08)
            ax.set_ylabel('S_obs / S_CALSPEC', fontsize=11)
            ax.set_title(title, fontsize=12)
            ax.legend(fontsize=8, loc='upper right', ncol=3)
            ax.set_ylim(0.4, 1.8)
            ax.set_xlim(min(all_wl_min, 0.95), max(all_wl_max, 5.6))
            ax.grid(True, alpha=0.3)

        ax_after.set_xlabel('Wavelength (µm)', fontsize=12)
        plt.tight_layout()
        outpath = f'{PLOT_DIR}/fig2_PID{pid}_residuals.png'
        plt.savefig(outpath, dpi=150, bbox_inches='tight')
        plt.close()
        print(f'Saved: {outpath}')


# ═══════════════════════════════════════════════════════════════════════════════
# PLOT 3: Calibrated spectra comparison (standard range vs extended range)
# ═══════════════════════════════════════════════════════════════════════════════
def plot_spectra_comparison():
    """Compare standard-range x1d to extended-range x1d to CALSPEC."""
    print("\n--- Plot 3: Spectra comparison (one figure per PID, gratings overlaid) ---")

    for pid in SOURCES.keys():
        name, cs_fname, pid_dir, _ = SOURCES[pid]

        # 2 main panels: Top for spectra (log), Bottom for residuals (linear)
        fig = plt.figure(figsize=(12, 10))
        gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1], hspace=0.1)
        
        ax = fig.add_subplot(gs[0])
        ax_res = fig.add_subplot(gs[1], sharex=ax)

        fig.suptitle(f'PID {pid} ({name}) — Observed vs CALSPEC spectra\n'
                     f'Comparing Nominal and Extended/Corrected Pipelines', fontsize=15, y=0.96)

        all_wl_min, all_wl_max = 100.0, 0.0

        # Load CALSPEC once for the background
        cs_wl, cs_fj = load_calspec(cs_fname)
        cs_interp = interp1d(cs_wl, cs_fj, bounds_error=False, fill_value=np.nan)
        
        for gkey, gcfg in GRATING_CONFIGS.items():
            if pid not in gcfg['pids']:
                continue
            
            # Calibration for this specific grating
            try:
                wl_cal, k_cal, alpha_cal, beta_cal = load_calib(gcfg['calib_ours'])
            except Exception as e:
                print(f'  WARNING: cannot load calib for {gkey}: {e}')
                continue

            # Load extended x1d and standard x1d
            x1d_ext = f'{IFU_DIR}/{pid_dir}/{gcfg["stage"]}/{gcfg["x1d_stem"]}'
            x1d_std = f'{IFU_DIR}/{pid_dir}/stage3/{gcfg["x1d_stem"]}'

            if not os.path.exists(x1d_ext):
                print(f'  WARNING: missing extended x1d for {gkey}')
                continue

            wl_e, fl_e, dq_e = load_x1d(x1d_ext)
            bad_e = (dq_e > 0) | ~np.isfinite(fl_e) | (fl_e <= 0)
            fl_e[bad_e] = np.nan

            # Apply calibration
            fl_corr = apply_calib(wl_e, fl_e, wl_cal, k_cal, alpha_cal, beta_cal)
            cs_at_obs = cs_interp(wl_e)

            # --- Top Panel: Spectra ---
            # Uncorrected (lighter)
            ax.plot(wl_e, fl_e, lw=0.8, color=gcfg['color_ext'], alpha=0.3, label='_nolegend_')
            # Corrected (bolder)
            ax.plot(wl_e, fl_corr, lw=1.2, color=gcfg['color'], alpha=0.9, 
                    label=f'{gcfg["short"]} Corrected')

            if os.path.exists(x1d_std):
                wl_s, fl_s, dq_s = load_x1d(x1d_std)
                bad_s = (dq_s > 0) | ~np.isfinite(fl_s) | (fl_s <= 0)
                fl_s[bad_s] = np.nan
                ax.plot(wl_s, fl_s, lw=1.0, color='gray', alpha=0.4, ls='--', 
                        label=f'{gcfg["short"]} Nominal' if gkey == 'g140m_f100lp' else '_nolegend_')

            # --- Bottom Panel: Residuals ---
            ratio_corr = fl_corr / cs_at_obs
            ratio_uncorr = fl_e / cs_at_obs
            
            # Median filter for visualization
            from scipy.ndimage import median_filter
            def smooth(y):
                ok = np.isfinite(y)
                if ok.sum() < 10: return y
                y_sm = np.full_like(y, np.nan)
                y_sm[ok] = median_filter(y[ok], size=11)
                return y_sm

            ax_res.plot(wl_e, smooth(ratio_uncorr), lw=1.0, color=gcfg['color_ext'], alpha=0.3, ls=':')
            ax_res.plot(wl_e, smooth(ratio_corr), lw=1.2, color=gcfg['color'], alpha=1.0)
            
            # Update global limits
            valid_wl = wl_e[np.isfinite(fl_e)]
            if len(valid_wl) > 0:
                all_wl_min = min(all_wl_min, valid_wl.min())
                all_wl_max = max(all_wl_max, valid_wl.max())

        # Finalize panels
        wl_grid = np.linspace(all_wl_min, all_wl_max, 2000)
        cs_grid = cs_interp(wl_grid)
        ax.plot(wl_grid, cs_grid, lw=1.5, color='black', ls=(0, (5, 5)), 
                alpha=0.8, label='CALSPEC', zorder=0)

        ax.set_ylabel('Flux (Jy)', fontsize=12)
        ax.set_yscale('log')
        if all_wl_max > all_wl_min:
            ax.set_xlim(all_wl_min*0.98, all_wl_max*1.02)
        ax.grid(True, alpha=0.15, which='both')
        ax.legend(fontsize=9, loc='upper right', ncol=3, frameon=True)

        # Dynamic limits based on CALSPEC
        ok_cs = np.isfinite(cs_grid) & (cs_grid > 0)
        if ok_cs.any():
            cmin, cmax = cs_grid[ok_cs].min(), cs_grid[ok_cs].max()
            ax.set_ylim(cmin / 5, cmax * 5)

        ax_res.axhline(1.0, color='black', lw=1.0, ls='-', alpha=0.8)
        ax_res.axhspan(0.98, 1.02, color='green', alpha=0.05)
        ax_res.set_ylabel('Ratio', fontsize=11)
        ax_res.set_xlabel('Wavelength (µm)', fontsize=12)
        ax_res.set_ylim(0.75, 1.25)
        ax_res.grid(True, alpha=0.15)

        plt.tight_layout()
        plt.subplots_adjust(top=0.92)
        outpath = f'{PLOT_DIR}/fig3_spectra_PID{pid}.png'
        plt.savefig(outpath, dpi=200, bbox_inches='tight')
        plt.close()
        print(f'Saved: {outpath}')




# ═══════════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════════
def main():
    print('Generating calibration plots ...')
    print(f'Output directory: {PLOT_DIR}')

    print('\n--- Plot 1: Calibration coefficients vs wavelength ---')
    plot_coefficients()

    print('\n--- Plot 2: Residuals before/after correction ---')
    plot_residuals()

    print('\n--- Plot 3: Spectra comparison (extended vs standard vs CALSPEC) ---')
    plot_spectra_comparison()

    print('\nAll plots saved.')


if __name__ == '__main__':
    main()
