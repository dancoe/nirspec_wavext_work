"""
Comparison plot: IFU v2 vs FS v2 vs Parlanti published k(λ) coefficients.

Shows k(λ) for both G140M and G235M NRS2 gratings across three calibration
derivations:
  - IFU v2   (from IFU stage3_ext data, P330E + J1743045)
  - FS v2    (from FS per-exposure NRS2 x1d data, P330E + J1743045)
  - Parlanti (published calibration from Parlanti et al. 2025 FITS files)

Also overlays α̃(λ) from Parlanti for scale reference.

Output saved to both:
  - reports/329_fs_v2/ifu_fs_parlanti_comparison.png
  - reports/329_fs_v2/ifu_fs_parlanti_comparison_alpha.png  (k + α, side-by-side)
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import interp1d

BASE         = '/Users/dcoe/NIRSpec/wavext'
IFU_V2_DIR   = f'{BASE}/nirspec_wavext_work/plots/Parlanti/cal/ifu_v2'
FS_V2_DIR    = f'{BASE}/nirspec_wavext_work/plots/Parlanti/cal/fs_v2'
PARLANTI_CAL = f'{BASE}/data/parlanti_repo/calibration_files'
OUTPUT_DIR   = f'{BASE}/nirspec_wavext_work/reports/329_fs_v2'
os.makedirs(OUTPUT_DIR, exist_ok=True)


def load_csv(path):
    if not os.path.exists(path):
        print(f'  WARNING: not found: {path}')
        return None
    data = np.loadtxt(path, delimiter=',', skiprows=1)
    return data[:, 0], data[:, 1], data[:, 2], data[:, 3]  # wl, k, alpha, beta


def load_parlanti(grating):
    fname = {
        'G140M': 'calibration_functions_g140m_f100lp.fits',
        'G235M': 'calibration_functions_g235m_f170lp.fits',
    }[grating]
    path = os.path.join(PARLANTI_CAL, fname)
    if not os.path.exists(path):
        return None, None, None
    with fits.open(path) as h:
        d = h[1].data
        wl = d['wavelength'].astype(float)
        k  = d['k'].astype(float)
        a  = d['alpha'].astype(float)
        b  = d['beta'].astype(float)
    return wl, k, a


# ── Figure 1: k(λ) comparison only ───────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle('NIRSpec NRS2 Throughput Correction k(λ)\n'
             'IFU v2 vs FS v2 vs Parlanti et al. (2025)',
             fontsize=13)

for ax, grating in zip(axes, ['G140M', 'G235M']):
    ifu_d = load_csv(f'{IFU_V2_DIR}/coeffs_ifu_v2_{grating}.csv')
    fs_d  = load_csv(f'{FS_V2_DIR}/coeffs_fs_v2_{grating}.csv')
    p_wl, p_k, p_a = load_parlanti(grating)

    if ifu_d is not None:
        ax.plot(ifu_d[0], ifu_d[1], color='#e74c3c', lw=2.2, label='IFU v2 k(λ)')
    if fs_d is not None:
        ax.plot(fs_d[0],  fs_d[1],  color='#2980b9', lw=2.2, label='FS v2 k(λ)')
    if p_wl is not None:
        # Clip Parlanti k to NRS2 range
        nrs2_lo = {'G140M': 1.87, 'G235M': 3.15}[grating]
        nrs2_hi = {'G140M': 3.55, 'G235M': 5.40}[grating]
        mask = (p_wl >= nrs2_lo) & (p_wl <= nrs2_hi) & np.isfinite(p_k)
        ax.plot(p_wl[mask], p_k[mask],
                color='#27ae60', lw=1.8, ls='--', label='Parlanti k(λ)')

    ax.axhline(1.0, color='0.5', lw=0.8, ls=':', alpha=0.7)
    ax.set_xlabel('Wavelength [µm]', fontsize=11)
    ax.set_ylabel('k(λ)',           fontsize=11)
    ax.set_title(f'{grating} NRS2', fontsize=12)
    ax.set_ylim(0, 2.2)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.2)

plt.tight_layout()
outpath = f'{OUTPUT_DIR}/ifu_fs_parlanti_comparison.png'
plt.savefig(outpath, dpi=200, bbox_inches='tight')
plt.close()
print(f'SAVING: {outpath}')

# ── Figure 2: k(λ) + α̃(λ) log scale, 2×2 grid ────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('NIRSpec NRS2 Calibration Coefficients — IFU v2 vs FS v2 vs Parlanti',
             fontsize=13)

grating_info = [
    ('G140M', 'G140M/F100LP (1.95–3.5 µm)', 1.87, 3.55),
    ('G235M', 'G235M/F170LP (3.3–5.3 µm)',  3.15, 5.40),
]

for col, (grating, title, lo, hi) in enumerate(grating_info):
    ifu_d = load_csv(f'{IFU_V2_DIR}/coeffs_ifu_v2_{grating}.csv')
    fs_d  = load_csv(f'{FS_V2_DIR}/coeffs_fs_v2_{grating}.csv')
    p_wl, p_k, p_a = load_parlanti(grating)

    # Top row: k(λ)
    ax = axes[0, col]
    if ifu_d is not None:
        ax.plot(ifu_d[0], ifu_d[1], color='#e74c3c', lw=2.0, label='IFU v2 k')
    if fs_d is not None:
        ax.plot(fs_d[0],  fs_d[1],  color='#2980b9', lw=2.0, label='FS v2 k')
    if p_wl is not None:
        mask = (p_wl >= lo) & (p_wl <= hi) & np.isfinite(p_k)
        ax.plot(p_wl[mask], p_k[mask], color='#27ae60', lw=1.8, ls='--', label='Parlanti k')
    ax.axhline(1.0, color='0.5', lw=0.8, ls=':', alpha=0.7)
    ax.set_ylabel('k(λ)', fontsize=11)
    ax.set_title(title, fontsize=10)
    ax.set_ylim(0, 2.2)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.2)
    ax.set_xlabel('Wavelength [µm]', fontsize=10)

    # Bottom row: α̃(λ) log scale
    ax = axes[1, col]
    if ifu_d is not None:
        a_ifu = np.where(ifu_d[2] > 0, ifu_d[2], 1e-5)
        ax.plot(ifu_d[0], a_ifu, color='#e74c3c', lw=2.0,
                label=r'IFU v2 $\tilde{\alpha}$')
    if fs_d is not None:
        a_fs = np.where(fs_d[2] > 0, fs_d[2], 1e-5)
        ax.plot(fs_d[0], a_fs, color='#2980b9', lw=2.0,
                label=r'FS v2 $\tilde{\alpha}$')
    if p_wl is not None:
        mask = (p_wl >= lo) & (p_wl <= hi) & np.isfinite(p_a) & (p_a > 0)
        ax.plot(p_wl[mask], p_a[mask], color='#27ae60', lw=1.8, ls='--',
                label=r'Parlanti $\tilde{\alpha}$')
    ax.set_yscale('log')
    ax.set_ylim(1e-4, 0.3)
    ax.set_ylabel(r'$\tilde{\alpha}(\lambda)$', fontsize=11)
    ax.set_xlabel('Wavelength [µm]', fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(True, which='both', alpha=0.2)

plt.tight_layout()
outpath = f'{OUTPUT_DIR}/ifu_fs_parlanti_comparison_full.png'
plt.savefig(outpath, dpi=200, bbox_inches='tight')
plt.close()
print(f'SAVING: {outpath}')
