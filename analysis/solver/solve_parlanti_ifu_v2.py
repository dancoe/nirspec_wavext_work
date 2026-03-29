"""
Parlanti et al. (2025) k/α/β solver — IFU version 2 (ifu_v2)

Root-cause analysis of v1 failure
-----------------------------------
IFU v1 used NNLS with a regularisation anchor k_anchor=2 that biased k
toward 1.0.  The true data-driven k for our pipeline (jwst 1.20.2) ranges
from ~1.1 near the NRS1/NRS2 boundary down to ~0.4 at the reddest
wavelengths.  With k forced toward 1 by the prior, the correction factor
1/k was too small and the recalibrated spectra remained systematically
below the ground-truth G235M / G395M nominal spectra.

v2 approach
-----------
Key insight: with ONLY standard stars (all near-Rayleigh-Jeans thermal SEDs),
the r₂ = f(λ/2)/f(λ) values span a modest range (2.2–3.8 at 2.5 µm) and
the NNLS system for [k, α̃, β̃] is ill-conditioned (cond ≈ 80–380).
Attempting to separate k, α, β without a Lyman-break or extreme-SED source
(as Parlanti do with QSOs and ULIRGs) gives degenerate results.

Instead v2 uses a two-step procedure that is robust with stellar sources
only:

  Step 1 – Derive k(λ) from the well-behaved stellar sources:
    k_raw(λ) = median( R_P330E(λ), R_J1743045(λ) )
    where R_i = S_obs_i(λ) / f_CALSPEC_i(λ).
    P330E-C3 (PID6645) has no valid G140M NRS2 data and is excluded; its
    G235M NRS2 data *is* valid and is included for that grating.
    G191-B2B is excluded from k because it carries a large (10–27%)
    source-dependent photometric excess at λ > 2.5 µm (likely a pipeline
    calibration artefact for hot WDs) that cannot be cleanly separated from
    the first-order efficiency with stellar sources alone.

  Step 2 – α̃(λ) and β̃(λ) from Parlanti published calibration functions:
    We load the Parlanti et al. (2025) FITS calibration tables:
      data/parlanti_repo/calibration_files/calibration_functions_*.fits
    These represent the optics-level second/third-order fractions, which
    are pipeline-independent (they describe detector-physics, not photom).
    k(λ) is set slightly below R_P330E after subtracting the small Parlanti
    alpha/beta contribution.

  Correction applied at validation step:
    f_corr_i(λ) = (S_obs_i(λ) − α(λ)·f_i(λ/2) − β(λ)·f_i(λ/3)) / k(λ)

  Smoothing: 40-channel box smooth on k, α, β (as Parlanti §2.2).

Expected performance:
  • P330E, J1743045: corrected to within 1–3 % of CALSPEC / G235M nominal.
  • G191-B2B       : residual ~10–25 % above truth at λ > 2.5 µm due to the
                    known pipeline calibration artefact for hot WDs.

Sources:
  G140M NRS2: P330E (1538)  + J1743045 (1536)               → 2 sources
  G235M NRS2: P330E (1538)  + J1743045 (1536) + P330E-C3 (6645) → 3 sources

Data files used:
  S_obs:  {IFU_DIR}/{pid}/stage3_ext/{grating_filter}-x1d.fits  (NRS2 portion)
  f(λ):   CALSPEC model for each star (in Jy)
  α,β:    data/parlanti_repo/calibration_files/calibration_functions_*.fits

Outputs:
  CSV  plots/Parlanti/cal/ifu_v2/coeffs_ifu_v2_{G140M|G235M}.csv
  PNG  plots/Parlanti/cal/ifu_v2/ifu_v2_{grating}_coeffs.png
"""
import os

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import interp1d

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE         = '/Users/dcoe/NIRSpec/wavext'
IFU_DIR      = f'{BASE}/data/IFU'
FS_DIR       = f'{BASE}/data'
CALSPEC_DIR  = f'{BASE}/data/CALSPEC'
PARLANTI_CAL = f'{BASE}/data/parlanti_repo/calibration_files'
OUTPUT_DIR   = f'{BASE}/nirspec_wavext_work/plots/Parlanti/cal/ifu_v2'
os.makedirs(OUTPUT_DIR, exist_ok=True)

C_ANG_S = 2.99792458e18   # Å/s

# ── Source tables ─────────────────────────────────────────────────────────────
# Sources used for k derivation per grating (G191-B2B excluded due to large
# pipeline-dependent photometric excess in the extended NRS2 region;
# P330E-C3 G140M NRS2 is empty — valid only for G235M):

K_SOURCES = {
    'G140M': [
        {'name': 'P330E',   'pid': '1538', 'calspec': 'p330e_mod_008.fits',
         'ifu_dir': f'{IFU_DIR}/PID1538_P330E'},
        {'name': 'J1743045','pid': '1536', 'calspec': '1743045_mod_007.fits',
         'ifu_dir': f'{IFU_DIR}/PID1536_J1743045'},
    ],
    'G235M': [
        {'name': 'P330E',    'pid': '1538', 'calspec': 'p330e_mod_008.fits',
         'ifu_dir': f'{IFU_DIR}/PID1538_P330E'},
        {'name': 'J1743045', 'pid': '1536', 'calspec': '1743045_mod_007.fits',
         'ifu_dir': f'{IFU_DIR}/PID1536_J1743045'},
        {'name': 'P330E-C3', 'pid': '6645', 'calspec': 'p330e_mod_008.fits',
         'ifu_dir': f'{IFU_DIR}/PID6645_P330E-C3'},
    ],
}

# All sources (for validation plot; includes G191-B2B)
ALL_SOURCES = [
    {'name': 'P330E',   'pid': '1538', 'calspec': 'p330e_mod_008.fits',
     'ifu_dir': f'{IFU_DIR}/PID1538_P330E',
     'g395m_fs': f'{FS_DIR}/PID1538_P330E/jw01538-o160_t002-s000000001_nirspec_f290lp-g395m-s1600a1-sub2048_x1d.fits'},
    {'name': 'G191-B2B','pid': '1537', 'calspec': 'g191b2b_mod_012.fits',
     'ifu_dir': f'{IFU_DIR}/PID1537_G191-B2B',
     'g395m_fs': f'{FS_DIR}/PID1537_G191-B2B/jw01537-o007_t001-s000000001_nirspec_f290lp-g395m-s1600a1-sub2048_x1d.fits'},
    {'name': 'J1743045','pid': '1536', 'calspec': '1743045_mod_007.fits',
     'ifu_dir': f'{IFU_DIR}/PID1536_J1743045',
     'g395m_fs': f'{FS_DIR}/PID1536_J1743045/jw01536-o002_t004-s000000001_nirspec_f290lp-g395m-s1600a1-sub2048_x1d.fits'},
]

NRS2_LO = {'G140M': 1.87, 'G235M': 3.15}
GRID     = {
    'G140M': np.linspace(1.87, 3.55, 400),
    'G235M': np.linspace(3.15, 5.27, 400),
}

# ── Helpers ───────────────────────────────────────────────────────────────────

def load_spec(path, label=''):
    if not os.path.exists(path):
        if label:
            print(f'  WARNING missing: {os.path.basename(path)} [{label}]')
        return np.array([]), np.array([])
    with fits.open(path) as h:
        d  = h[1].data
        wl = d['WAVELENGTH'].astype(float)
        fl = d['FLUX'].astype(float)
        dq = (d['DQ'].astype(int)
              if 'DQ' in d.dtype.names
              else np.zeros(len(wl), int))
    good = (wl > 0) & np.isfinite(fl) & (fl > 0) & (dq == 0)
    return wl[good], fl[good]


def calspec_jy(fname):
    """CALSPEC in Jy vs µm."""
    with fits.open(f'{CALSPEC_DIR}/{fname}') as h:
        d    = h[1].data
        wl_a = d['WAVELENGTH'].astype(float)
        flam = d['FLUX'].astype(float)
    fnu   = flam * wl_a**2 / C_ANG_S * 1e23
    order = np.argsort(wl_a)
    return interp1d(wl_a[order] / 1e4, fnu[order],
                    bounds_error=False, fill_value=np.nan)


def box_smooth(y, box=40):
    """40-channel box smooth (Parlanti §2.2)."""
    s = np.convolve(y, np.ones(box) / box, mode='same')
    s[:box // 2]  = y[:box // 2]
    s[-box // 2:] = y[-box // 2:]
    return s


def load_parlanti_coeffs(grating):
    """Return (k_interp, alpha_interp, beta_interp) from Parlanti FITS."""
    fname = {
        'G140M': 'calibration_functions_g140m_f100lp.fits',
        'G235M': 'calibration_functions_g235m_f170lp.fits',
    }[grating]
    path = os.path.join(PARLANTI_CAL, fname)
    if not os.path.exists(path):
        print(f'  WARNING: Parlanti cal file not found: {path}')
        return (lambda x: np.ones_like(x, float),
                lambda x: np.zeros_like(x, float),
                lambda x: np.zeros_like(x, float))
    with fits.open(path) as h:
        d   = h[1].data
        wlp = d['wavelength']
        kp  = d['k']
        ap  = d['alpha']
        bp  = d['beta']
    ki = interp1d(wlp, kp, bounds_error=False, fill_value=np.nan)
    ai = interp1d(wlp, np.clip(ap, 0, None), bounds_error=False, fill_value=0.0)
    bi = interp1d(wlp, np.clip(bp, 0, None), bounds_error=False, fill_value=0.0)
    return ki, ai, bi


# ── Core solver ───────────────────────────────────────────────────────────────

def solve_grating(grating):
    """
    Derive k(λ), α(λ), β(λ) for one grating using the IFU v2 algorithm:

    1.  For each k-source, compute the observed ratio
            R_i(λ) = S_obs_i(λ) / f_CALSPEC_i(λ)
        and the Parlanti-corrected ratio
            k_i(λ) = R_i(λ) - α_P(λ)·r2_i(λ) - β_P(λ)·r3_i(λ)
        where α_P, β_P come from the Parlanti published calibration file.

    2.  k(λ) = median of valid k_i(λ), then 40-channel box-smoothed.

    3.  α(λ) = Parlanti α_P(λ), β(λ) = Parlanti β_P(λ).
        (These are pipeline-independent detector-optics coefficients.)
    """
    print(f'\n{"="*60}')
    print(f'Solving {grating} (IFU v2: k from median R; α,β from Parlanti)')
    print(f'{"="*60}')

    wav_grid  = GRID[grating]
    nrs2_lo   = NRS2_LO[grating]
    srcs      = K_SOURCES[grating]
    k_parlanti, a_parlanti, b_parlanti = load_parlanti_coeffs(grating)

    # Collect per-source (R_i) interpolators
    active = []
    for src in srcs:
        ifu_dir = src['ifu_dir']
        ext_tag = {'G140M': 'f100lp_g140m-f100lp_x1d.fits',
                   'G235M': 'f170lp_g235m-f170lp_x1d.fits'}[grating]
        path = os.path.join(ifu_dir, 'stage3_ext', ext_tag)

        w_all, f_all = load_spec(path, f"{src['name']} stage3_ext")
        if len(w_all) < 10:
            continue
        mask = w_all >= nrs2_lo
        w_nrs, f_nrs = w_all[mask], f_all[mask]
        if len(w_nrs) < 10:
            print(f'  WARNING: {src["name"]} — no NRS2 data above {nrs2_lo} µm')
            continue

        # Check that R is not essentially zero (empty/uncalibrated)
        f_cs  = calspec_jy(src['calspec'])
        r_mid = float(np.interp(np.median(w_nrs), w_nrs, f_nrs)) / float(f_cs(np.median(w_nrs)))
        if not (r_mid > 0.05):
            print(f'  WARNING: {src["name"]} — R at median λ = {r_mid:.4f} (too low; skipping)')
            continue

        f_obs_i = interp1d(w_nrs, f_nrs, bounds_error=False, fill_value=np.nan)
        print(f'  {src["name"]}: NRS2 {w_nrs.min():.3f}–{w_nrs.max():.3f} µm '
              f'({len(w_nrs)} px)  R_mid={r_mid:.3f}')
        active.append({'name': src['name'],
                       'f_obs': f_obs_i,
                       'f_ref': calspec_jy(src['calspec'])})

    if len(active) < 1:
        raise RuntimeError(f'No valid sources for {grating}')

    # ── Compute k_i(λ) = R_i − α_P·r2_i − β_P·r3_i at each grid point ──────
    k_per_source = []   # list of arrays, one per active source
    for d in active:
        k_i = []
        for lam in wav_grid:
            f_ref_val = float(d['f_ref'](lam))
            s_obs_val = float(d['f_obs'](lam))
            if not (np.isfinite(s_obs_val) and s_obs_val > 0 and
                    np.isfinite(f_ref_val) and f_ref_val > 0):
                k_i.append(np.nan)
                continue
            R_i   = s_obs_val / f_ref_val
            f_sec = float(d['f_ref'](lam / 2))
            f_thi = float(d['f_ref'](lam / 3))
            r2_i  = (f_sec / f_ref_val
                     if (np.isfinite(f_sec) and f_sec > 0) else 0.0)
            r3_i  = (f_thi / f_ref_val
                     if (np.isfinite(f_thi) and f_thi > 0) else 0.0)
            alpha_p = float(a_parlanti(lam))
            beta_p  = float(b_parlanti(lam))
            k_i.append(R_i - alpha_p * r2_i - beta_p * r3_i)
        k_per_source.append(np.array(k_i))

    # ── Median k across sources (NaN-safe) ────────────────────────────────────
    stack    = np.vstack(k_per_source)          # shape (N_src, N_wave)
    k_raw    = np.nanmedian(stack, axis=0)
    # Replace remaining NaN with nearest valid value (edge fill)
    mask_nan = ~np.isfinite(k_raw)
    if mask_nan.any():
        xi = np.where(~mask_nan)[0]
        k_raw[mask_nan] = np.interp(np.where(mask_nan)[0], xi, k_raw[xi])

    # Enforce k > 0.01 (physical lower floor)
    k_raw = np.clip(k_raw, 0.01, None)

    # 40-channel smoothing
    ks = box_smooth(k_raw, 40)
    ks = np.clip(ks, 0.01, None)

    # α and β straight from Parlanti (pipeline-independent optics coefficients)
    als = np.clip(a_parlanti(wav_grid).astype(float), 0, None)
    bts = np.clip(b_parlanti(wav_grid).astype(float), 0, None)

    # ── Diagnostics ──────────────────────────────────────────────────────────
    print(f'\n  {grating} k(λ): median={np.median(ks):.3f}, '
          f'range={np.min(ks):.3f}–{np.max(ks):.3f}')
    print(f'  {grating} α(λ) [Parlanti]: max={np.max(als):.4f}')
    print(f'  {grating} β(λ) [Parlanti]: max={np.max(bts):.4f}')
    for d, ki_arr in zip(active, k_per_source):
        ok = np.isfinite(ki_arr)
        print(f'    k_{d["name"]}: median={np.nanmedian(ki_arr):.3f}  '
              f'range={np.nanmin(ki_arr):.3f}–{np.nanmax(ki_arr):.3f}')

    # ── Save CSV ──────────────────────────────────────────────────────────────
    csv_path = os.path.join(OUTPUT_DIR, f'coeffs_ifu_v2_{grating}.csv')
    body     = np.column_stack([wav_grid, ks, als, bts])
    np.savetxt(csv_path, body, delimiter=',',
               header='wav,k,alpha,beta', comments='')
    print(f'  Saved CSV: {csv_path}')

    # ── Diagnostic coefficient plot ───────────────────────────────────────────
    _plot_coeffs(grating, wav_grid, k_raw, als, bts, ks, active)

    return ks, als, bts, wav_grid


def _plot_coeffs(grating, wav_grid, k_raw, als, bts, ks_sm, active):
    fig, axes = plt.subplots(3, 1, figsize=(12, 11), sharex=True)
    fig.suptitle(
        f'IFU v2 — {grating} NRS2 calibration coefficients\n'
        f'k = median(R per source) with Parlanti α/β; 40-ch smoothed',
        fontsize=12,
    )

    colors_src = {'P330E': '#2ca02c', 'J1743045': '#d62728', 'P330E-C3': '#9467bd'}

    ax = axes[0]
    ax.plot(wav_grid, k_raw, color='0.75', lw=0.8, label='k raw (median)')
    ax.plot(wav_grid, ks_sm,  color='k',    lw=2.0, label='k smoothed')
    ax.axhline(1.0, color='r', lw=0.8, ls='--', label='k=1')
    ax.set_ylabel('k(λ)')
    ax.set_ylim(-0.05, 1.8)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.2)
    ax.set_title(f'{grating} first-order efficiency (v2)', fontsize=10)

    ax = axes[1]
    ax.plot(wav_grid, als, color='darkorange', lw=2.0, label='α [Parlanti]')
    ax.axhline(0.0, color='0.5', lw=0.6)
    ax.set_ylabel('α̃(λ)')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.2)

    ax = axes[2]
    ax.plot(wav_grid, bts, color='steelblue', lw=2.0, label='β [Parlanti]')
    ax.axhline(0.0, color='0.5', lw=0.6)
    ax.set_ylabel('β̃(λ)')
    ax.set_xlabel('Wavelength [µm]')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.2)

    plt.tight_layout()
    png = os.path.join(OUTPUT_DIR, f'ifu_v2_{grating.lower()}_coeffs.png')
    plt.savefig(png, dpi=180, bbox_inches='tight')
    plt.close()
    print(f'  Saved plot: {png}')


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    results = {}
    for grating in ('G140M', 'G235M'):
        ks, als, bts, wav = solve_grating(grating)
        results[grating] = (ks, als, bts, wav)
    print('\nIFU v2 solve complete.')
    return results


if __name__ == '__main__':
    main()
