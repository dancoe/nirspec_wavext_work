"""
FS v6 — Validation Plots

Generates:
1. Coefficient plots: k, α, β vs wavelength (v6 vs v5 vs Parlanti)
2. Per-source NRS2 spectra: obs vs corrected vs CALSPEC truth (4 standards)
3. PID 1492 cross-validation: corrected G140M NRS2 ext vs G235M nominal,
                               corrected G235M NRS2 ext vs G395M nominal
4. Full merged spectra for CALSPEC standards

Output directory: reports/v6/fs_v6/plots/
"""
import os
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import interp1d

BASE         = '/Users/dcoe/NIRSpec/wavext'
DATA_DIR     = f'{BASE}/data'
CALSPEC_DIR  = f'{BASE}/data/CALSPEC'
PARLANTI_CAL = f'{BASE}/data/parlanti_repo/calibration_files'
V5_DIR       = f'{BASE}/results/v5'
V6_DIR       = f'{BASE}/results/v6'
PLOTDIR      = f'{BASE}/nirspec_wavext_work/reports/v6/fs_v6/plots'
os.makedirs(PLOTDIR, exist_ok=True)

C_ANG_S = 2.99792458e18

CALSPEC_SOURCES = {
    '1537': ('G191-B2B',    'g191b2b_mod_012.fits',  '#ff7f0e', ['PID1537_G191-B2B']),
    '1538': ('P330E',       'p330e_mod_008.fits',    '#1f77b4', ['PID1538_P330E']),
    '1536': ('J1743045',    '1743045_mod_007.fits',  '#2ca02c', ['PID1536_J1743045']),
    '6644': ('NGC2506-G31', 'ngc2506g31_mod_003.fits','#9467bd', ['PID6644_NGC2506G31','PID6644_NGC2506-G31']),
}

NRS2_LO = {'g140m': 1.87, 'g235m': 3.15}

PID1492_DIR = f'{DATA_DIR}/PID1492'
PID1492 = {
    'g140m': {
        'nrs2_ext':   f'{PID1492_DIR}/jw01492003001_03102_00005_nrs2_g140m_extract_1d.fits',
        'truth':      f'{PID1492_DIR}/jw01492-o001_t001-s000000001_nirspec_f170lp-g235m-s200a1-subs200a1_x1d.fits',
        'nrs1_truth': f'{PID1492_DIR}/jw01492-o001_t001-s000000001_nirspec_f100lp-g140m-s200a1-subs200a1_x1d.fits',
        'truth_label': 'G235M nominal (MAST)',
    },
    'g235m': {
        'nrs2_ext':   f'{PID1492_DIR}/jw01492003001_03104_00004_nrs2_g235m_extract_1d.fits',
        'truth':      f'{PID1492_DIR}/jw01492-o001_t001-s000000001_nirspec_f290lp-g395m-s200a1-subs200a1_x1d.fits',
        'nrs1_truth': f'{PID1492_DIR}/jw01492-o001_t001-s000000001_nirspec_f170lp-g235m-s200a1-subs200a1_x1d.fits',
        'truth_label': 'G395M nominal (MAST)',
    },
}


# ── Helpers ────────────────────────────────────────────────────────────────────
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
    fname = {'g140m': 'calibration_functions_g140m_f100lp.fits',
             'g235m': 'calibration_functions_g235m_f170lp.fits'}[grating]
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


def correct_spectrum(wl, fl, f_truth, k_itp, a_itp, b_itp):
    kv = k_itp(wl)
    av = a_itp(wl)
    bv = b_itp(wl)
    return (fl - av * f_truth(wl / 2.0) - bv * f_truth(wl / 3.0)) / np.where(kv > 0.01, kv, np.nan)


def find_nrs2_ext(pid, grating, root_dirs):
    for root in root_dirs:
        rp = os.path.join(DATA_DIR, root)
        if not os.path.isdir(rp):
            continue
        for f in glob.glob(os.path.join(rp, '**', f'*{grating}*x1d.fits'), recursive=True):
            if 'nrs2' in os.path.basename(f).lower():
                return f
    return None


def find_nrs1_nom(pid, grating, root_dirs):
    for root in root_dirs:
        rp = os.path.join(DATA_DIR, root)
        if not os.path.isdir(rp):
            continue
        hits = glob.glob(os.path.join(rp, f'*{grating}*x1d.fits'))
        hits = [f for f in hits if 'nrs2' not in os.path.basename(f).lower()]
        if hits:
            return hits[0]
        for f in glob.glob(os.path.join(rp, '**', f'*{grating}*x1d.fits'), recursive=True):
            if 'nrs2' not in os.path.basename(f).lower():
                return f
    return None


# ── Plot 1: Coefficient comparison ───────────────────────────────────────────
def plot_coefficients():
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('FS v6 — k, α, β Coefficients\n'
                 'v6 (4 stds + PID1492 cross-grating) vs v5 (4 stds) vs Parlanti (2025)',
                 fontsize=12)

    for col, grating in enumerate(['g140m', 'g235m']):
        key = f'{grating}_{"f100lp" if grating == "g140m" else "f170lp"}'
        v6_path = os.path.join(V6_DIR, f'calib_v6_fs_{key}.fits')
        v5_path = os.path.join(V5_DIR, f'calib_v5_{key}.fits')
        nrs2_lo = NRS2_LO[grating]

        if not os.path.exists(v6_path):
            print(f'  MISSING: {v6_path}')
            continue

        wl6, k6, a6, b6 = load_coeffs(v6_path)
        wv = np.linspace(wl6.min(), wl6.max(), 500)
        parl = load_parlanti(grating)

        for row, (fn6, label, ylim) in enumerate([
            (k6, 'k(λ)', (0, 1.3)),
            (a6, 'α(λ)', (0, 0.5)),
        ]):
            ax = axes[row, col]
            ax.plot(wv, fn6(wv), 'k-', lw=2, label='v6 (FS, smoothed)')

            if os.path.exists(v5_path):
                _, kv5, av5, _ = load_coeffs(v5_path)
                fn5 = kv5 if row == 0 else av5
                ax.plot(wv, fn5(wv), 'b--', lw=1.5, alpha=0.7, label='v5 (FS, 4 stds)')

            if parl is not None:
                fp = parl[1] if row == 0 else parl[2]
                ax.plot(wv, fp(wv), 'r:', lw=1.5, alpha=0.7, label='Parlanti (2025)')

            ax.axvline(nrs2_lo, color='orange', ls='--', lw=1, alpha=0.5)
            ax.axhline(1.0 if row == 0 else 0.0, color='gray', lw=0.6, ls=':')
            ax.set_ylim(*ylim)
            ax.set_xlabel('Wavelength (µm)')
            ax.set_ylabel(label)
            ax.set_title(f'{grating.upper()} — {label}')
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.2)

    plt.tight_layout()
    out = os.path.join(PLOTDIR, 'fs_v6_coefficients.png')
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  SAVED: {out}')


# ── Plot 2: CALSPEC standard validation ──────────────────────────────────────
def plot_calspec_validation():
    for grating in ('g140m', 'g235m'):
        key     = f'{grating}_{"f100lp" if grating == "g140m" else "f170lp"}'
        nrs2_lo = NRS2_LO[grating]
        v6_path = os.path.join(V6_DIR, f'calib_v6_fs_{key}.fits')
        if not os.path.exists(v6_path):
            continue
        _, k6, a6, b6 = load_coeffs(v6_path)

        srcs = list(CALSPEC_SOURCES.items())
        fig, axes = plt.subplots(len(srcs), 2, figsize=(14, 4 * len(srcs)))
        fig.suptitle(f'FS v6 — {grating.upper()} NRS2 Extended Validation (CALSPEC standards)',
                     fontsize=12)

        for i, (pid, (name, calspec_fname, color, root_dirs)) in enumerate(srcs):
            nrs2_path = find_nrs2_ext(pid, grating, root_dirs)
            if nrs2_path is None:
                continue
            wl2, fl2 = load_spec(nrs2_path)
            mask = wl2 >= nrs2_lo
            wl_nrs, fl_nrs = wl2[mask], fl2[mask]
            if len(wl_nrs) < 10:
                continue

            f_cs    = calspec_jy(calspec_fname)
            fl_corr = correct_spectrum(wl_nrs, fl_nrs, f_cs, k6, a6, b6)
            fl_true = f_cs(wl_nrs)

            ax = axes[i, 0]
            ax.plot(wl_nrs, fl_nrs,  color='silver', lw=1,   label='Observed (NRS2 ext)')
            ax.plot(wl_nrs, fl_corr, color=color,    lw=1.5, label='Corrected (v6)')
            ax.plot(wl_nrs, fl_true, 'k--',          lw=1.5, label='CALSPEC truth')
            ax.set_ylabel('Flux (Jy)')
            ax.set_title(f'{name}')
            ax.legend(fontsize=8); ax.grid(True, alpha=0.2)
            ax.set_xlabel('Wavelength (µm)')

            ax2 = axes[i, 1]
            ratio = fl_corr / fl_true
            ax2.plot(wl_nrs, ratio, color=color, lw=1.5)
            ax2.axhline(1.0, color='k', ls='--', lw=1)
            ax2.fill_between(wl_nrs, 0.9, 1.1, alpha=0.12, color='green', label='±10%')
            ax2.set_ylim(0.3, 1.7)
            ax2.set_ylabel('Corrected / Truth')
            ax2.set_title(f'{name} — ratio')
            ax2.legend(fontsize=8); ax2.grid(True, alpha=0.2)
            ax2.set_xlabel('Wavelength (µm)')

        plt.tight_layout()
        out = os.path.join(PLOTDIR, f'fs_v6_{grating}_calspec_validation.png')
        plt.savefig(out, dpi=150, bbox_inches='tight')
        plt.close()
        print(f'  SAVED: {out}')


# ── Plot 3: PID 1492 cross-validation ─────────────────────────────────────────
def plot_pid1492_xval():
    for grating in ('g140m', 'g235m'):
        key     = f'{grating}_{"f100lp" if grating == "g140m" else "f170lp"}'
        nrs2_lo = NRS2_LO[grating]
        v6_path = os.path.join(V6_DIR, f'calib_v6_fs_{key}.fits')
        if not os.path.exists(v6_path):
            continue
        _, k6, a6, b6 = load_coeffs(v6_path)

        fmap = PID1492[grating]
        wl_ext, fl_ext = load_spec(fmap['nrs2_ext'])
        wl_truth, fl_truth = load_spec(fmap['truth'])
        wl_nrs1, fl_nrs1 = load_spec(fmap['nrs1_truth'])

        if len(wl_ext) < 10 or len(wl_truth) < 10:
            print(f'  PID1492 {grating.upper()} data missing')
            continue

        mask = wl_ext >= nrs2_lo
        wl_nrs, fl_nrs = wl_ext[mask], fl_ext[mask]

        # Compute normalisation scale (DN/s → pseudo-Jy)
        lo = max(wl_nrs.min(), wl_truth.min())
        hi = min(wl_nrs.max(), wl_truth.max())
        wl_ov = np.linspace(lo, hi, 50)
        f_truth_itp = interp1d(wl_truth, fl_truth, bounds_error=False, fill_value=np.nan)
        f_nrs2_itp  = interp1d(wl_nrs,   fl_nrs,   bounds_error=False, fill_value=np.nan)
        ratio_ov = f_truth_itp(wl_ov) / f_nrs2_itp(wl_ov)
        good = np.isfinite(ratio_ov) & (ratio_ov > 0)
        scale = float(np.nanmedian(ratio_ov[good])) if good.sum() >= 5 else 1.0
        fl_nrs_jy = fl_nrs * scale

        # Build stitched truth for ghost terms
        stitch = wl_truth.min()
        if len(wl_nrs1) > 10:
            lo_part = wl_nrs1 < stitch
            hi_part = wl_truth >= stitch
            wl_c = np.concatenate([wl_nrs1[lo_part], wl_truth[hi_part]])
            fl_c = np.concatenate([fl_nrs1[lo_part], fl_truth[hi_part]])
        else:
            wl_c, fl_c = wl_truth, fl_truth
        idx = np.argsort(wl_c)
        f_stitched = interp1d(wl_c[idx], fl_c[idx],
                              bounds_error=False, fill_value='extrapolate')

        fl_corr = correct_spectrum(wl_nrs, fl_nrs_jy, f_stitched, k6, a6, b6)
        fl_ref  = f_truth_itp(wl_nrs)

        fig, axes = plt.subplots(2, 1, figsize=(12, 9))
        fig.suptitle(f'FS v6 — PID 1492 {grating.upper()} Cross-Validation\n'
                     f'Corrected NRS2 ext vs {fmap["truth_label"]}',
                     fontsize=12)

        ax = axes[0]
        ax.plot(wl_nrs,   fl_nrs_jy,  color='silver',     lw=1,   label=f'{grating.upper()} NRS2 ext (norm to Jy)')
        ax.plot(wl_nrs,   fl_corr,    color='darkorange',  lw=1.5, label=f'{grating.upper()} NRS2 corrected (v6)')
        ax.plot(wl_truth, fl_truth,   'k--',               lw=1.5, label=fmap['truth_label'])
        ax.axvline(nrs2_lo, color='gray', ls=':', lw=1)
        ax.set_ylabel('Flux (Jy)')
        ax.set_xlabel('Wavelength (µm)')
        ax.set_title(f'PID 1492 — {grating.upper()} Extended vs {fmap["truth_label"]}')
        ax.legend(fontsize=9); ax.grid(True, alpha=0.2)

        ax2 = axes[1]
        ratio = fl_corr / fl_ref
        ax2.plot(wl_nrs, ratio, color='darkorange', lw=1.5)
        ax2.axhline(1.0, color='k', ls='--', lw=1)
        ax2.fill_between(wl_nrs, 0.8, 1.2, alpha=0.12, color='green', label='±20%')
        ax2.fill_between(wl_nrs, 0.9, 1.1, alpha=0.12, color='green')
        ax2.set_ylim(0.2, 2.0)
        ax2.set_ylabel('Corrected / Truth')
        ax2.set_xlabel('Wavelength (µm)')
        ax2.set_title(f'Ratio: PID1492 corrected / {fmap["truth_label"]}')
        ax2.legend(fontsize=9); ax2.grid(True, alpha=0.2)

        plt.tight_layout()
        out = os.path.join(PLOTDIR, f'fs_v6_pid1492_{grating}_xval.png')
        plt.savefig(out, dpi=150, bbox_inches='tight')
        plt.close()
        print(f'  SAVED: {out}')


# ── Plot 4: Full merged spectra for CALSPEC standards ────────────────────────
def plot_full_merged():
    for grating in ('g140m', 'g235m'):
        key     = f'{grating}_{"f100lp" if grating == "g140m" else "f170lp"}'
        nrs2_lo = NRS2_LO[grating]
        v6_path = os.path.join(V6_DIR, f'calib_v6_fs_{key}.fits')
        if not os.path.exists(v6_path):
            continue
        _, k6, a6, b6 = load_coeffs(v6_path)

        for pid, (name, calspec_fname, color, root_dirs) in CALSPEC_SOURCES.items():
            nrs2_path = find_nrs2_ext(pid, grating, root_dirs)
            nrs1_path = find_nrs1_nom(pid, grating, root_dirs)
            if nrs2_path is None:
                continue

            wl2, fl2 = load_spec(nrs2_path)
            mask = wl2 >= nrs2_lo
            wl_nrs, fl_nrs = wl2[mask], fl2[mask]
            if len(wl_nrs) < 10:
                continue

            f_cs    = calspec_jy(calspec_fname)
            fl_corr = correct_spectrum(wl_nrs, fl_nrs, f_cs, k6, a6, b6)

            fig, ax = plt.subplots(figsize=(12, 5))
            if nrs1_path is not None:
                wl1, fl1 = load_spec(nrs1_path)
                if len(wl1) > 10:
                    ax.plot(wl1, fl1, color='steelblue', lw=1.2, label='NRS1 nominal')

            wl_cs = np.linspace(
                min(nrs2_lo * 0.5, 0.8),
                max(wl_nrs.max() + 0.2, 5.5), 3000)
            ax.plot(wl_cs, f_cs(wl_cs), 'k--', lw=1.5, alpha=0.8, label='CALSPEC truth')
            ax.plot(wl_nrs, fl_nrs,  color='silver', lw=1, label='NRS2 observed')
            ax.plot(wl_nrs, fl_corr, color=color,    lw=1.5, label='NRS2 corrected (v6)')
            ax.axvline(nrs2_lo, color='orange', ls='--', lw=1, alpha=0.7)
            ax.set_xlabel('Wavelength (µm)')
            ax.set_ylabel('Flux (Jy)')
            ax.set_title(f'FS v6 — {name} {grating.upper()} Full Spectrum')
            ax.legend(fontsize=9); ax.grid(True, alpha=0.2)
            y_max = float(np.nanpercentile(f_cs(wl_cs[np.isfinite(f_cs(wl_cs))]), 99)) * 1.3
            ax.set_ylim(-0.001, max(y_max, 0.001))
            plt.tight_layout()
            name_safe = name.replace('-', '')
            out = os.path.join(PLOTDIR, f'fs_v6_full_{grating}_{name_safe}.png')
            plt.savefig(out, dpi=150, bbox_inches='tight')
            plt.close()
            print(f'  SAVED: {out}')


if __name__ == '__main__':
    print('Generating FS v6 validation plots...')
    plot_coefficients()
    plot_calspec_validation()
    plot_pid1492_xval()
    plot_full_merged()
    print('Done.')
