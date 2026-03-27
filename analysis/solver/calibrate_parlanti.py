#!/usr/bin/env python
"""
Parlanti et al. 2025 Flux Correction: Coefficient Determination + Calibrated Plot
==================================================================================

Determines k(λ), a(λ), b(λ) coefficients using the Parlanti model:

    S(λ) = k(λ)·f(λ) + a(λ)·f(λ/2) + b(λ)·f(λ/3)

Where:
    S(λ) = extended NRS2 extraction (obs003, DN/s, scaled to Jy)
    f(λ) = PRISM baseline (obs001 Level 3, Jy)
    k(λ) = 1st-order throughput (drops from ~1 at nominal boundary)
    a(λ) = 2nd-order contamination fraction (< 5%)
    b(λ) = 3rd-order contamination fraction (< 1%)

Generates:
    plots/Parlanti/cal/FS_1492_cal.png  — 2-panel calibration plot (cf. Parlanti Fig. 5)
    plots/Parlanti/cal/parlanti_coefficients_v2.txt — improved coefficient table

Notes on units
--------------
The extended NRS2 extractions were run through assign_wcs + extract_2d + extract_1d
but NOT the photom step, so flux is in DN/s rather than Jy.  We recover the Jy scale
by matching each extended spectrum to the PRISM continuum at the first ~50 wavelength
points of the extended region, where k ≈ 1 and a, b ≈ 0 (still close to the nominal
boundary).  Emission-line peaks are excluded by 3-σ clipping of the DN/s→Jy ratios.
"""

import os, sys
import warnings
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
from numpy.polynomial.legendre import legvander

sys.path.insert(0, '/Users/dcoe/NIRSpec/wavext/jwst_nirspec_wavext')
os.environ.setdefault('CRDS_PATH', os.path.expanduser('~/crds_cache'))
os.environ.setdefault('CRDS_SERVER_URL', 'https://jwst-crds.stsci.edu')
warnings.filterwarnings('ignore')
from stdatamodels.jwst import datamodels

# ── Paths ─────────────────────────────────────────────────────────────────────
DATA      = '/Users/dcoe/NIRSpec/wavext/data/PID1492'
PLOTS_CAL = '/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/plots/Parlanti/cal'
os.makedirs(PLOTS_CAL, exist_ok=True)

# ── Helper: read 1-D extracted spectrum ───────────────────────────────────────
def read_x1d(fpath, col='flux'):
    """Return (wavelength_µm, flux) arrays from a JWST MultiSpecModel x1d.fits."""
    with datamodels.open(fpath) as m:
        all_w, all_f = [], []
        for s in getattr(m, 'spec', []):
            w = np.array(s.spec_table['wavelength'], dtype=float)
            f = np.array(s.spec_table[col], dtype=float)
            msk = (w > 0) & np.isfinite(f) & (f > 0)
            all_w.append(w[msk])
            all_f.append(f[msk])
    if not all_w:
        return np.array([]), np.array([])
    w = np.concatenate(all_w)
    f = np.concatenate(all_f)
    idx = np.argsort(w)
    return w[idx], f[idx]


# ── 1. Load all spectra ────────────────────────────────────────────────────────
print("Loading spectra …")
prism_w,    prism_f    = read_x1d(f'{DATA}/jw01492-o001_t001-s000000001_nirspec_clear-prism-s200a1-subs200a1_x1d.fits')
g140m_nom_w, g140m_nom_f = read_x1d(f'{DATA}/jw01492-o001_t001-s000000001_nirspec_f100lp-g140m-s200a1-subs200a1_x1d.fits')
g235m_nom_w, g235m_nom_f = read_x1d(f'{DATA}/jw01492-o001_t001-s000000001_nirspec_f170lp-g235m-s200a1-subs200a1_x1d.fits')
g395m_nom_w, g395m_nom_f = read_x1d(f'{DATA}/jw01492-o001_t001-s000000001_nirspec_f290lp-g395m-s200a1-subs200a1_x1d.fits')

g140m_ext_w, g140m_ext_raw = read_x1d(f'{DATA}/jw01492003001_03102_00005_nrs2_g140m_extract_1d.fits')
g235m_ext_w, g235m_ext_raw = read_x1d(f'{DATA}/jw01492003001_03104_00004_nrs2_g235m_extract_1d.fits')
g395m_ext_w, g395m_ext_raw = read_x1d(f'{DATA}/jw01492003001_03106_00002_nrs2_g395m_extract_1d.fits')

for name, w, f in [('PRISM', prism_w, prism_f),
                   ('G140M_nom', g140m_nom_w, g140m_nom_f),
                   ('G235M_nom', g235m_nom_w, g235m_nom_f),
                   ('G395M_nom', g395m_nom_w, g395m_nom_f),
                   ('G140M_ext', g140m_ext_w, g140m_ext_raw),
                   ('G235M_ext', g235m_ext_w, g235m_ext_raw),
                   ('G395M_ext', g395m_ext_w, g395m_ext_raw)]:
    if len(w):
        print(f"  {name:12s}: {len(w):5d} pts  "
              f"λ=[{w[0]:.3f},{w[-1]:.3f}] µm  "
              f"flux_med={np.median(f):.3e}")

# ── 2. PRISM interpolator ──────────────────────────────────────────────────────
prism_interp = interp1d(prism_w, prism_f, kind='linear',
                        bounds_error=False, fill_value=np.nan)

# ── 3. Scale extended DN/s spectra to Jy ──────────────────────────────────────
def boundary_scale(ext_w, ext_raw, prism_interp, n_bdy=80, sigma_clip=2.5):
    """
    Compute multiplicative scale factor to convert ext_raw (DN/s) → Jy by
    matching the PRISM continuum in the first n_bdy wavelength points
    (i.e., near the nominal–extended boundary where k ≈ 1 and a, b ≈ 0).

    Sigma-clipping on the ratio removes emission-line peaks.
    """
    w_b = ext_w[:n_bdy]
    f_b = ext_raw[:n_bdy]
    p_b = prism_interp(w_b)

    valid = np.isfinite(p_b) & (p_b > 0) & (f_b > 0)
    ratios = p_b[valid] / f_b[valid]          # [Jy / (DN/s)]

    # Iterative sigma clipping
    for _ in range(5):
        med = np.nanmedian(ratios)
        mad = np.nanmedian(np.abs(ratios - med))
        keep = np.abs(ratios - med) < sigma_clip * 1.4826 * mad
        if keep.sum() < 5:
            break
        ratios = ratios[keep]

    scale = np.nanmedian(ratios)
    return scale


scale_140m = boundary_scale(g140m_ext_w, g140m_ext_raw, prism_interp)
scale_235m = boundary_scale(g235m_ext_w, g235m_ext_raw, prism_interp)
scale_395m = boundary_scale(g395m_ext_w, g395m_ext_raw, prism_interp)

g140m_ext_jy = g140m_ext_raw * scale_140m
g235m_ext_jy = g235m_ext_raw * scale_235m
g395m_ext_jy = g395m_ext_raw * scale_395m

print(f"\nDN/s → Jy scale factors:")
print(f"  G140M: {scale_140m:.4e} Jy/(DN/s)  [boundary λ ≈ {g140m_ext_w[:80].mean():.2f} µm]")
print(f"  G235M: {scale_235m:.4e} Jy/(DN/s)  [boundary λ ≈ {g235m_ext_w[:80].mean():.2f} µm]")
print(f"  G395M: {scale_395m:.4e} Jy/(DN/s)  [boundary λ ≈ {g395m_ext_w[:20].mean():.2f} µm]")


# ── 4. Parlanti coefficient fitting ───────────────────────────────────────────
def fit_parlanti(s_wav, s_jy, prism_interp,
                 deg=4, lambda_reg=0.05,
                 clip_sigma=3.0, smooth_win=0):
    """
    Fit k(λ), a(λ), b(λ) as Legendre polynomials (degree=deg) in
    normalised wavelength using Tikhonov-regularised least squares.

    Physical priors are imposed after fitting:
        k clipped to [0.01, 2.0]   (allows for flat-field over/under-corrections)
        a clipped to [0, 0.20]     (2nd-order contamination ≤ 20%)
        b clipped to [0, 0.05]     (3rd-order contamination ≤ 5%)

    Parameters
    ----------
    s_wav      : wavelength array of the extended spectrum (µm)
    s_jy       : flux array in Jy (already boundary-scaled)
    deg        : polynomial degree for each of k, a, b
    lambda_reg : Tikhonov regularisation strength (fraction of max eigenvalue)
    clip_sigma : σ-clip on residuals (0 to skip)
    smooth_win : Savitzky-Golay window for pre-smoothing (0 to skip)

    Returns
    -------
    k, a, b : arrays evaluated at s_wav
    residuals: (s_jy – model) / s_jy  for diagnostic reporting
    """
    if smooth_win > 0:
        s_jy = savgol_filter(s_jy, window_length=smooth_win, polyorder=3)

    # Normalise λ to [-1, 1] for Legendre basis
    lam_min, lam_max = s_wav.min(), s_wav.max()
    lam_n = 2 * (s_wav - lam_min) / (lam_max - lam_min) - 1

    # Evaluate PRISM at 1st, 2nd, 3rd order wavelengths
    f1 = prism_interp(s_wav)          # f(λ)
    f2 = prism_interp(s_wav / 2.0)   # f(λ/2)  — 2nd-order contamination
    f3 = prism_interp(s_wav / 3.0)   # f(λ/3)  — 3rd-order contamination

    valid = (np.isfinite(f1) & np.isfinite(f2) & np.isfinite(f3) &
             np.isfinite(s_jy) & (f1 > 0) & (s_jy > 0))

    w_v = lam_n[valid];  S_v = s_jy[valid]
    F1  = f1[valid];     F2  = f2[valid];  F3 = f3[valid]

    def build_A(w_arr, F1, F2, F3, deg):
        n = deg + 1
        L = legvander(w_arr, deg)          # (N, n)
        A = np.zeros((len(w_arr), 3 * n))
        A[:, :n]      = L * F1[:, None]   # k-block
        A[:, n:2*n]   = L * F2[:, None]   # a-block
        A[:, 2*n:3*n] = L * F3[:, None]   # b-block
        return A

    A = build_A(w_v, F1, F2, F3, deg)

    # Tikhonov: solve (A^T A + λ I) x = A^T S
    ATA = A.T @ A
    lam = lambda_reg * np.linalg.norm(ATA, 'fro') / ATA.shape[0]
    ATA_reg = ATA + lam * np.eye(ATA.shape[0])
    ATs = A.T @ S_v
    x, _, _, _ = np.linalg.lstsq(ATA_reg, ATs, rcond=None)

    # Optional: iterative sigma-clip on residuals and re-fit
    if clip_sigma > 0:
        for _iter in range(3):
            model_v = A @ x
            resid = S_v - model_v
            std = np.std(resid)
            keep = np.abs(resid) < clip_sigma * std
            if keep.sum() < 20:
                break
            A_k  = build_A(w_v[keep], F1[keep], F2[keep], F3[keep], deg)
            S_k  = S_v[keep]
            ATA_k = A_k.T @ A_k + lam * np.eye(ATA.shape[0])
            x, _, _, _ = np.linalg.lstsq(ATA_k, A_k.T @ S_k, rcond=None)

    # Evaluate k, a, b at ALL wavelength points
    n = deg + 1
    L_all = legvander(lam_n, deg)
    k_raw = L_all @ x[:n]
    a_raw = L_all @ x[n:2*n]
    b_raw = L_all @ x[2*n:3*n]

    # Apply physical clipping
    k_phys = np.clip(k_raw, 0.01, 2.0)
    a_phys = np.clip(a_raw, 0.00, 0.20)
    b_phys = np.clip(b_raw, 0.00, 0.05)

    # Residuals for diagnostics
    model_all = k_phys * f1 + a_phys * f2 + b_phys * f3
    rel_resid  = (s_jy - model_all) / np.where(s_jy > 0, s_jy, np.nan)

    return k_phys, a_phys, b_phys, rel_resid


print("\nFitting G140M (λ = 1.98–3.27 µm) …")
k140, a140, b140, resid140 = fit_parlanti(
    g140m_ext_w, g140m_ext_jy, prism_interp, deg=4, lambda_reg=0.1,
    clip_sigma=3.0, smooth_win=15)
rms140 = np.nanstd(resid140) * 100
print(f"  G140M fit RMS residual: {rms140:.1f}%  "
      f"k=[{k140.min():.3f},{k140.max():.3f}]  "
      f"a=[{a140.min():.4f},{a140.max():.4f}]  "
      f"b=[{b140.min():.4f},{b140.max():.4f}]")

print("Fitting G235M (λ = 3.31–5.31 µm) …")
k235, a235, b235, resid235 = fit_parlanti(
    g235m_ext_w, g235m_ext_jy, prism_interp, deg=4, lambda_reg=0.1,
    clip_sigma=3.0, smooth_win=15)
rms235 = np.nanstd(resid235) * 100
print(f"  G235M fit RMS residual: {rms235:.1f}%  "
      f"k=[{k235.min():.3f},{k235.max():.3f}]  "
      f"a=[{a235.min():.4f},{a235.max():.4f}]  "
      f"b=[{b235.min():.4f},{b235.max():.4f}]")


# ── 5.  Build a combined high-resolution reference interpolator ───────────────
# For the CALIBRATION APPLICATION step, we subtract a·f(λ/2) + b·f(λ/3) using
# the highest-available spectral resolution estimate of the source.  This better
# removes sharp emission-line contamination than the low-res PRISM alone.
#
# Priority order per wavelength:
#   G140M_nom → G235M_nom → G395M_nom → PRISM (as fallback)
#
# To prevent double-counting we join with a smooth taper; here we simply prefer
# the grating nominal wherever available and fall back to PRISM elsewhere.

def build_hires_ref(nom_list, prism_w, prism_f):
    """
    Merge [(w_array, f_array), ...] nominal spectra and PRISM into one
    continuous interpolator, preferring the grating spectra where available.
    """
    # Collect all (w, f) pairs, nominals first
    segments = list(nom_list) + [(prism_w, prism_f)]
    # Pool all points; in overlapping regions the first-listed (grating) wins
    # because we build a dict keyed on rounded wavelength
    pool = {}
    for w_arr, f_arr in reversed(segments):   # reverse so earlier segs overwrite
        for wv, fv in zip(w_arr, f_arr):
            pool[round(float(wv), 5)] = float(fv)
    keys = sorted(pool.keys())
    w_hr = np.array(keys)
    f_hr = np.array([pool[k] for k in keys])
    return interp1d(w_hr, f_hr, kind='linear', bounds_error=False, fill_value=np.nan)

hires_interp = build_hires_ref(
    [(g140m_nom_w, g140m_nom_f),
     (g235m_nom_w, g235m_nom_f),
     (g395m_nom_w, g395m_nom_f)],
    prism_w, prism_f)


# ── 6. Apply calibration ───────────────────────────────────────────────────────
def apply_parlanti(s_wav, s_jy, k, a, b, ref_interp):
    """
    Recover the calibrated 1st-order spectrum from the contaminated one:
        f_cal = (S − a·f(λ/2) − b·f(λ/3)) / k
    Using a high-res reference for the contamination terms.
    """
    f2 = ref_interp(s_wav / 2.0)
    f3 = ref_interp(s_wav / 3.0)
    # Guard against division by near-zero k
    k_safe = np.where(k > 0.05, k, np.nan)
    f_cal  = (s_jy - a * np.where(np.isfinite(f2), f2, 0)
                   - b * np.where(np.isfinite(f3), f3, 0)) / k_safe
    return f_cal


g140m_cal = apply_parlanti(g140m_ext_w, g140m_ext_jy, k140, a140, b140, hires_interp)
g235m_cal = apply_parlanti(g235m_ext_w, g235m_ext_jy, k235, a235, b235, hires_interp)


# ── 6. Diagnostic coefficient report ──────────────────────────────────────────
coeff_path = os.path.join(PLOTS_CAL, 'parlanti_coefficients_v2.txt')
with open(coeff_path, 'w') as fh:
    fh.write("# Parlanti et al. 2025 Flux Correction Coefficients — Version 2\n")
    fh.write("# Method: Legendre polynomial fit (deg=4) with Tikhonov regularisation\n")
    fh.write("# Units: k, a, b are dimensionless (throughput fractions)\n")
    fh.write("# Physical bounds: k ∈ [0.01, 2.0], a ∈ [0, 0.20], b ∈ [0, 0.05]\n\n")
    for gname, w, k_ar, a_ar, b_ar, rms in [
            ('G140M', g140m_ext_w, k140, a140, b140, rms140),
            ('G235M', g235m_ext_w, k235, a235, b235, rms235)]:
        fh.write("=" * 60 + "\n")
        fh.write(f"Grating: {gname}\n")
        fh.write(f"λ range: {w.min():.3f} – {w.max():.3f} µm\n")
        fh.write(f"RMS residual: {rms:.1f}%\n")
        fh.write(f"k(λ): min={k_ar.min():.4f}  max={k_ar.max():.4f}  "
                 f"mean={k_ar.mean():.4f}\n")
        fh.write(f"a(λ): min={a_ar.min():.4f}  max={a_ar.max():.4f}  "
                 f"mean={a_ar.mean():.4f}\n")
        fh.write(f"b(λ): min={b_ar.min():.4f}  max={b_ar.max():.4f}  "
                 f"mean={b_ar.mean():.4f}\n\n")
print(f"\nCoefficients saved to {coeff_path}")


# ── 7. Plot FS_1492_cal.png ────────────────────────────────────────────────────
print("Generating FS_1492_cal.png …")

fig, axes = plt.subplots(2, 1, figsize=(10, 10))
fig.suptitle(
    'NIRSpec Extended Wavelength Calibration — PID 1492 FS\n'
    'IRAS 05248-7007  (Parlanti et al. 2025 model)',
    fontsize=12, y=0.99)

GRATING_COLORS = {'G140M': '#d98ad6', 'G235M': '#7878d6', 'G395M': '#78c3d6'}

def draw_grating_arrows(ax, grat_dict, xlim):
    """Draw <-> arrows at 94% of y-max, clipped to xlim."""
    y1 = ax.get_ylim()[1]
    y_arr = y1 * 0.940
    y_lbl = y1 * 0.972
    for g, (lo, hi, col) in grat_dict.items():
        lo_c = max(lo, xlim[0])
        hi_c = min(hi, xlim[1])
        if hi_c <= lo_c:
            continue
        ax.annotate('', xy=(hi_c, y_arr), xytext=(lo_c, y_arr),
                    arrowprops=dict(arrowstyle='<->', color=col, lw=2.0,
                                   mutation_scale=12))
        ax.text((lo_c + hi_c) / 2, y_lbl, g, ha='center', va='bottom',
                color=col, fontsize=9, fontweight='bold')

# ── Panel 1: G140M extension (1.98–3.27 µm) ──────────────────────────────────
ax = axes[0]
XLIM1 = (0.95, 3.30)

g235_mask_p1 = g235m_nom_w <= XLIM1[1]
ax.plot(g140m_nom_w, g140m_nom_f * 1e3, 'k-', lw=1.5, zorder=6,
        label='Nominal G140M+G235M')
ax.plot(g235m_nom_w[g235_mask_p1], g235m_nom_f[g235_mask_p1] * 1e3,
        'k-', lw=1.5, zorder=6)

# Extended: uncalibrated (green) and calibrated (red)
ax.plot(g140m_ext_w, g140m_ext_jy * 1e3, color='green',
        lw=0.7, alpha=0.55, zorder=3, label='Extended G140M — not calibrated')
cal_smooth = savgol_filter(np.where(np.isfinite(g140m_cal), g140m_cal, 0),
                           window_length=25, polyorder=3)
ax.plot(g140m_ext_w, cal_smooth * 1e3, color='red',
        lw=1.2, alpha=0.85, zorder=5, label='Extended G140M — calibrated')
ax.plot(prism_w, prism_f * 1e3, color='gray', lw=0.8, alpha=0.35,
        zorder=2, label='PRISM f(λ)')

ax.set_xlim(*XLIM1)
ref_p1 = np.concatenate([g140m_nom_f * 1e3,
                          g235m_nom_f[g235_mask_p1] * 1e3,
                          cal_smooth[np.isfinite(cal_smooth)] * 1e3,
                          prism_f * 1e3])
p1_ymax = np.percentile(ref_p1[ref_p1 > 0], 99.8) * 1.40
ax.set_ylim(-0.04 * p1_ymax, p1_ymax)

ax.set_xlabel('Wavelength [µm]', fontsize=11)
ax.set_ylabel('Flux [mJy]', fontsize=11)
ax.legend(fontsize=9, loc='upper left', framealpha=0.85)
ax.axvline(1.84, ls='--', color='gray', lw=0.8, alpha=0.5)
ax.axvline(1.98, ls=':',  color='gray', lw=0.8, alpha=0.5)
ax.set_title('G140M extension — NRS2 (1.98–3.27 µm)', fontsize=10)
draw_grating_arrows(axes[0],
    {'G140M': (0.97, 1.84, GRATING_COLORS['G140M']),
     'G235M': (1.66, 3.17, GRATING_COLORS['G235M'])}, XLIM1)


# ── Panel 2: G235M extension (3.31–5.31 µm) ──────────────────────────────────
ax = axes[1]
XLIM2 = (1.65, 5.35)

ax.plot(g235m_nom_w, g235m_nom_f * 1e3, 'k-', lw=1.5, zorder=6,
        label='Nominal G235M+G395M')
ax.plot(g395m_nom_w, g395m_nom_f * 1e3, 'k-', lw=1.5, zorder=6)

ax.plot(g235m_ext_w, g235m_ext_jy * 1e3, color='green',
        lw=0.7, alpha=0.55, zorder=3, label='Extended G235M — not calibrated')
cal235_smooth = savgol_filter(np.where(np.isfinite(g235m_cal), g235m_cal, 0),
                              window_length=25, polyorder=3)
ax.plot(g235m_ext_w, cal235_smooth * 1e3, color='red',
        lw=1.2, alpha=0.85, zorder=5, label='Extended G235M — calibrated')
ax.plot(prism_w, prism_f * 1e3, color='gray', lw=0.8, alpha=0.35,
        zorder=2, label='PRISM f(λ)')

ax.set_xlim(*XLIM2)
ref_p2 = np.concatenate([g235m_nom_f * 1e3, g395m_nom_f * 1e3,
                          cal235_smooth[np.isfinite(cal235_smooth)] * 1e3,
                          prism_f * 1e3])
p2_ymax = np.percentile(ref_p2[ref_p2 > 0], 99.8) * 1.40
ax.set_ylim(-0.04 * p2_ymax, p2_ymax)

ax.set_xlabel('Wavelength [µm]', fontsize=11)
ax.set_ylabel('Flux [mJy]', fontsize=11)
ax.legend(fontsize=9, loc='upper left', framealpha=0.85)
ax.axvline(3.07, ls='--', color='gray', lw=0.8, alpha=0.5)
ax.axvline(3.31, ls=':',  color='gray', lw=0.8, alpha=0.5)
ax.set_title('G235M extension — NRS2 (3.31–5.31 µm)', fontsize=10)
draw_grating_arrows(axes[1],
    {'G235M': (1.66, 3.17, GRATING_COLORS['G235M']),
     'G395M': (2.87, 5.27, GRATING_COLORS['G395M'])}, XLIM2)


plt.tight_layout()
out_path = os.path.join(PLOTS_CAL, 'FS_1492_cal.png')
plt.savefig(out_path, dpi=150, bbox_inches='tight')
print(f"Saved: {out_path}")

# ── 8. Diagnostic: coefficient profiles ───────────────────────────────────────
fig2, axes2 = plt.subplots(2, 3, figsize=(13, 7))
fig2.suptitle('Parlanti Coefficients k(λ), a(λ), b(λ) — Diagnostic', fontsize=11)

for row, (gname, w, k_ar, a_ar, b_ar, s_jy, cal) in enumerate([
        ('G140M', g140m_ext_w, k140, a140, b140, g140m_ext_jy, g140m_cal),
        ('G235M', g235m_ext_w, k235, a235, b235, g235m_ext_jy, g235m_cal)]):

    ax = axes2[row, 0]
    ax.plot(w, k_ar, 'navy', lw=1.5, label='k(λ)')
    ax.set_ylabel('k(λ)', color='navy')
    ax.set_title(f'{gname} — 1st-order throughput k(λ)')
    ax.set_xlabel('λ (µm)');  ax.set_ylim(0, 2.1)
    ax.axhline(1.0, ls='--', color='k', lw=0.8, alpha=0.4)
    ax.set_xlim(w[0], w[-1])

    ax = axes2[row, 1]
    ax.plot(w, a_ar * 100, 'darkorange', lw=1.5, label='a(λ) %')
    ax.set_ylabel('a(λ) [%]', color='darkorange')
    ax.set_title(f'{gname} — 2nd-order fraction a(λ)')
    ax.set_xlabel('λ (µm)');  ax.set_ylim(0, 22)
    ax.set_xlim(w[0], w[-1])

    ax = axes2[row, 2]
    ax.plot(w, b_ar * 100, 'darkgreen', lw=1.5, label='b(λ) %')
    ax.set_ylabel('b(λ) [%]', color='darkgreen')
    ax.set_title(f'{gname} — 3rd-order fraction b(λ)')
    ax.set_xlabel('λ (µm)');  ax.set_ylim(0, 6)
    ax.set_xlim(w[0], w[-1])

plt.tight_layout()
diag_path = os.path.join(PLOTS_CAL, 'parlanti_coefficients_v2_diagnostic.png')
plt.savefig(diag_path, dpi=150, bbox_inches='tight')
print(f"Saved: {diag_path}")

plt.close('all')
print("\nDone.")
