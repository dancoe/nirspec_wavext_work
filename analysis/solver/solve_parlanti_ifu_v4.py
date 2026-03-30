"""
Parlanti et al. (2025) k/α/β solver — IFU version 4 (ifu_v4)

Upgrade over IFU v3
--------------------
v3 used Parlanti's published α̃ and β̃ and only solved for k.
v4 solves for all three coefficients simultaneously using Non-Negative Least
Squares (NNLS) at each wavelength grid point. This gives a fully independent
empirical calibration — no Parlanti published α/β assumed.

Key change
-----------
Model at each λ:
  S_obs,i(λ) = k(λ)·f_i(λ) + α(λ)·f_i(λ/2) + β(λ)·f_i(λ/3)

Matrix form at wavelength λ:
  A = [[f1(λ), f1(λ/2), f1(λ/3)],
       [f2(λ), f2(λ/2), f2(λ/3)],
       [f3(λ), f3(λ/2), f3(λ/3)]]
  b = [S_obs,1(λ), S_obs,2(λ), S_obs,3(λ)]
  x = [k, α, β]   (all ≥ 0)

Solved with scipy.optimize.nnls. Each coefficient vector subsequently
smoothed with a 40-channel boxcar (consistent with Parlanti 2025).

Sources (3; all IFU stage3_ext)
---------------------------------
  P330E    (G2V,  PID 1538) — solar SED, primary anchor
  G191-B2B (WD,   PID 1537) — hot featureless SED
  J1743045 (A8III,PID 1536) — intermediate warm SED

Spectral diversity: G2V vs WD vs A8III ensures the r2 = f(λ/2)/f(λ) and
r3 = f(λ/3)/f(λ) columns of A are sufficiently linearly independent to
separate k from α at most wavelengths.

Note: NGC2506-G31 (G1V, PID 6644 — a 4th Cool-star source) is not yet
downloaded. When available it will make the system over-determined and
further stabilise the solution (especially for FS; see FS v4).

Pipeline provenance (identical to v3)
---------------------------------------
  stage3_ext x1d = Spec3Pipeline output from NRS1+NRS2  combined cube
  (Parlanti cubepar, extended wavelengthrange ASDF, patched nirspec.py)

Outputs
--------
  CSV   plots/Parlanti/cal/ifu_v4/coeffs_ifu_v4_{G140M|G235M}.csv
  PNG   plots/Parlanti/cal/ifu_v4/  (coeffs, spectra, MAST vs L3)
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import interp1d
from scipy.optimize import nnls

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE         = '/Users/dcoe/NIRSpec/wavext'
IFU_DIR      = f'{BASE}/data/IFU'
FS_DIR       = f'{BASE}/data'
CALSPEC_DIR  = f'{BASE}/data/CALSPEC'
PARLANTI_CAL = f'{BASE}/data/parlanti_repo/calibration_files'
OUTPUT_DIR   = f'{BASE}/nirspec_wavext_work/plots/Parlanti/cal/ifu_v4'
os.makedirs(OUTPUT_DIR, exist_ok=True)

C_ANG_S = 2.99792458e18  # Å/s

# ── Source tables ──────────────────────────────────────────────────────────────
# All 3 sources used in NNLS solver
ALL_SOURCES = [
    {
        'name':          'P330E',
        'pid':           '1538',
        'calspec':       'p330e_mod_008.fits',
        'ifu_dir':       f'{IFU_DIR}/PID1538_P330E',
        'mast_g395m_fs': f'{FS_DIR}/PID1538_P330E/jw01538-o160_t002-s000000001_nirspec_f290lp-g395m-s1600a1-sub2048_x1d.fits',
        'sptype':        'G2V',
    },
    {
        'name':          'G191-B2B',
        'pid':           '1537',
        'calspec':       'g191b2b_mod_012.fits',
        'ifu_dir':       f'{IFU_DIR}/PID1537_G191-B2B',
        'mast_g395m_fs': f'{FS_DIR}/PID1537_G191-B2B/jw01537-o007_t001-s000000001_nirspec_f290lp-g395m-s1600a1-sub2048_x1d.fits',
        'sptype':        'WD-DA.8',
    },
    {
        'name':          'J1743045',
        'pid':           '1536',
        'calspec':       '1743045_mod_007.fits',
        'ifu_dir':       f'{IFU_DIR}/PID1536_J1743045',
        'mast_g395m_fs': f'{FS_DIR}/PID1536_J1743045/jw01536-o002_t004-s000000001_nirspec_f290lp-g395m-s1600a1-sub2048_x1d.fits',
        'sptype':        'A8III',
    },
]

NRS2_LO = {'G140M': 1.87, 'G235M': 3.15}
GRID = {
    'G140M': np.linspace(1.87, 3.55, 400),
    'G235M': np.linspace(3.15, 5.27, 400),
}

EXT_TAG = {
    'G140M': 'f100lp_g140m-f100lp_x1d.fits',
    'G235M': 'f170lp_g235m-f170lp_x1d.fits',
}
NOM_TAG = EXT_TAG  # same tag; stage3 vs stage3_ext differ by directory

# ── Helpers ────────────────────────────────────────────────────────────────────
def load_spec(path, label=''):
    if not os.path.exists(path):
        if label:
            print(f'  WARNING missing: {os.path.basename(path)} [{label}]')
        return np.array([]), np.array([])
    with fits.open(path) as h:
        d  = h[1].data
        wl = d['WAVELENGTH'].astype(float)
        fl = d['FLUX'].astype(float)
        dq = (d['DQ'].astype(int) if 'DQ' in d.dtype.names
              else np.zeros(len(wl), int))
    good = (wl > 0) & np.isfinite(fl) & (fl > 0) & (dq == 0)
    return wl[good], fl[good]


def calspec_jy(fname):
    with fits.open(f'{CALSPEC_DIR}/{fname}') as h:
        d    = h[1].data
        wl_a = d['WAVELENGTH'].astype(float)
        flam = d['FLUX'].astype(float)
    fnu   = flam * wl_a**2 / C_ANG_S * 1e23
    order = np.argsort(wl_a)
    return interp1d(wl_a[order] / 1e4, fnu[order],
                    bounds_error=False, fill_value=np.nan)


def box_smooth(y, box=40):
    s = np.convolve(y, np.ones(box) / box, mode='same')
    s[:box // 2]  = y[:box // 2]
    s[-box // 2:] = y[-box // 2:]
    return s


def load_parlanti_coeffs(grating):
    """Load Parlanti published coefficients for comparison (not used in solver)."""
    fname = {
        'G140M': 'calibration_functions_g140m_f100lp.fits',
        'G235M': 'calibration_functions_g235m_f170lp.fits',
    }[grating]
    path = os.path.join(PARLANTI_CAL, fname)
    if not os.path.exists(path):
        print(f'  WARNING: Parlanti file {path} not found — comparison plots will be empty')
        return None
    with fits.open(path) as h:
        d = h[1].data
        wlp, kp, ap, bp = d['wavelength'], d['k'], d['alpha'], d['beta']
    return (wlp,
            interp1d(wlp, kp,                   bounds_error=False, fill_value=np.nan),
            interp1d(wlp, np.clip(ap, 0, None), bounds_error=False, fill_value=0.0),
            interp1d(wlp, np.clip(bp, 0, None), bounds_error=False, fill_value=0.0))


# ── Core NNLS solver ───────────────────────────────────────────────────────────
def solve_grating(grating):
    print(f'\n{"="*60}')
    print(f'IFU v4 NNLS Solving {grating} — solving k, α, β simultaneously')
    print(f'Sources: P330E (G2V) + G191-B2B (WD) + J1743045 (A8III)')
    print(f'{"="*60}')

    wav_grid = GRID[grating]
    nrs2_lo  = NRS2_LO[grating]

    # Load and interpolate all sources onto the grid
    active = []
    for src in ALL_SOURCES:
        path = os.path.join(src['ifu_dir'], 'stage3_ext', EXT_TAG[grating])
        w_all, f_all = load_spec(path, f"{src['name']} stage3_ext")
        if len(w_all) < 10:
            print(f'  SKIP {src["name"]}: no data in {os.path.basename(path)}')
            continue
        mask = w_all >= nrs2_lo
        w_nrs, f_nrs = w_all[mask], f_all[mask]
        if len(w_nrs) < 10:
            print(f'  SKIP {src["name"]}: <10 NRS2 pixels above {nrs2_lo} µm')
            continue

        f_cs  = calspec_jy(src['calspec'])
        # Sanity: median ratio should be reasonable (>0.05)
        r_mid = float(np.interp(np.median(w_nrs), w_nrs, f_nrs)) / float(f_cs(np.median(w_nrs)))
        if not (r_mid > 0.01):
            print(f'  SKIP {src["name"]}: R_mid={r_mid:.4f} implausibly low')
            continue

        f_obs_i = interp1d(w_nrs, f_nrs, bounds_error=False, fill_value=np.nan)
        print(f'  {src["name"]} ({src["sptype"]}): NRS2 {w_nrs.min():.3f}–'
              f'{w_nrs.max():.3f} µm ({len(w_nrs)} px)  R_mid={r_mid:.3f}')
        active.append({
            'name':  src['name'],
            'sptype': src['sptype'],
            'f_obs': f_obs_i,
            'f_ref': f_cs,
        })

    if len(active) < 2:
        raise RuntimeError(f'Need ≥2 sources for {grating}; only {len(active)} valid')
    if len(active) < 3:
        print(f'  WARNING: only {len(active)} sources — system may be under-determined')

    # ── NNLS at each wavelength point ─────────────────────────────────────────
    n_grid = len(wav_grid)
    k_raw  = np.full(n_grid, np.nan)
    a_raw  = np.full(n_grid, np.nan)
    b_raw  = np.full(n_grid, np.nan)
    cond_numbers = np.full(n_grid, np.nan)

    for j, lam in enumerate(wav_grid):
        rows_A = []
        rows_b = []
        for d in active:
            s_obs = float(d['f_obs'](lam))
            f1    = float(d['f_ref'](lam))
            f2    = float(d['f_ref'](lam / 2.0))
            f3    = float(d['f_ref'](lam / 3.0))
            if not (np.isfinite(s_obs) and s_obs > 0 and
                    np.isfinite(f1)    and f1    > 0 and
                    np.isfinite(f2)    and f2    > 0 and
                    np.isfinite(f3)    and f3    > 0):
                continue
            rows_A.append([f1, f2, f3])
            rows_b.append(s_obs)

        if len(rows_b) < 2:
            continue

        A = np.array(rows_A, dtype=float)
        bv = np.array(rows_b, dtype=float)

        # Condition-number estimate (ratio of singular values)
        try:
            sv = np.linalg.svd(A, compute_uv=False)
            cond_numbers[j] = sv[0] / sv[-1] if sv[-1] > 0 else np.inf
        except np.linalg.LinAlgError:
            pass

        x, _ = nnls(A, bv)
        k_raw[j], a_raw[j], b_raw[j] = x[0], x[1], x[2]

    # Fill isolated NaN gaps by linear interpolation, then clip
    for arr in (k_raw, a_raw, b_raw):
        bad = ~np.isfinite(arr)
        if bad.any():
            xi = np.where(~bad)[0]
            if len(xi) > 1:
                arr[bad] = np.interp(np.where(bad)[0], xi, arr[xi])
    k_raw = np.clip(k_raw, 0.01, None)
    a_raw = np.clip(a_raw, 0.0,  None)
    b_raw = np.clip(b_raw, 0.0,  None)

    # 40-channel smoothing (identical to Parlanti 2025 and v2/v3)
    ks  = box_smooth(k_raw, 40)
    als = box_smooth(a_raw, 40)
    bts = box_smooth(b_raw, 40)
    ks  = np.clip(ks,  0.01, None)
    als = np.clip(als, 0.0,  None)
    bts = np.clip(bts, 0.0,  None)

    # Condition number summary
    med_cond = np.nanmedian(cond_numbers)
    print(f'\n  {grating} — NNLS condition number: median={med_cond:.1f}, '
          f'max={np.nanmax(cond_numbers):.1f}')
    print(f'  {grating} k(λ):  median={np.median(ks):.3f}, '
          f'range={np.min(ks):.3f}–{np.max(ks):.3f}')
    print(f'  {grating} α(λ):  median={np.median(als):.4f}, max={np.max(als):.4f}')
    print(f'  {grating} β(λ):  median={np.median(bts):.4f}, max={np.max(bts):.4f}')

    # Save CSV
    csv_path = os.path.join(OUTPUT_DIR, f'coeffs_ifu_v4_{grating}.csv')
    np.savetxt(csv_path,
               np.column_stack([wav_grid, ks, als, bts, k_raw, a_raw, b_raw]),
               delimiter=',',
               header='wav,k_smooth,alpha_smooth,beta_smooth,k_raw,alpha_raw,beta_raw',
               comments='')
    print(f'  CSV saved: {csv_path}')

    return wav_grid, ks, als, bts, k_raw, a_raw, b_raw, active


# ── Plot 1: coefficients vs Parlanti ─────────────────────────────────────────
def plot_coeffs(grating, wav_grid, k_raw, ks, a_raw, als, b_raw, bts):
    parl = load_parlanti_coeffs(grating)

    fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
    fig.suptitle(
        f'IFU v4 — {grating} NRS2 NNLS coefficients (k, α, β solved simultaneously)\n'
        f'3 sources: P330E (G2V) + G191-B2B (WD) + J1743045 (A8III); 40-ch smoothed',
        fontsize=11)

    # k
    ax = axes[0]
    ax.plot(wav_grid, k_raw, color='0.75', lw=0.8, label='k raw (NNLS)')
    ax.plot(wav_grid, ks,    color='k',    lw=2.0, label='k smoothed (v4)')
    if parl is not None:
        ax.plot(wav_grid, parl[1](wav_grid), 'r--', lw=1.5, alpha=0.7,
                label='k Parlanti (published)')
    ax.axhline(1.0, color='gray', lw=0.7, ls=':', alpha=0.5)
    ax.set_ylabel('k(λ)'); ax.set_ylim(-0.05, 2.5)
    ax.legend(fontsize=8); ax.grid(True, alpha=0.2)
    ax.set_title(f'{grating} throughput k(λ) — IFU v4 NNLS vs Parlanti', fontsize=10)

    # α
    ax = axes[1]
    ax.plot(wav_grid, a_raw, color='wheat',      lw=0.8, label='α raw (NNLS)')
    ax.plot(wav_grid, als,   color='darkorange', lw=2.0, label='α smoothed (v4)')
    if parl is not None:
        ax.plot(wav_grid, parl[2](wav_grid), 'r--', lw=1.5, alpha=0.7,
                label='α̃ Parlanti (published)')
    ax.axhline(0, color='0.5', lw=0.6)
    ax.set_ylabel('α(λ)'); ax.legend(fontsize=8); ax.grid(True, alpha=0.2)
    ax.set_title(f'{grating} 2nd-order α(λ)', fontsize=10)

    # β
    ax = axes[2]
    ax.plot(wav_grid, b_raw, color='lightblue', lw=0.8, label='β raw (NNLS)')
    ax.plot(wav_grid, bts,   color='steelblue', lw=2.0, label='β smoothed (v4)')
    if parl is not None:
        ax.plot(wav_grid, parl[3](wav_grid), 'r--', lw=1.5, alpha=0.7,
                label='β̃ Parlanti (published)')
    ax.axhline(0, color='0.5', lw=0.6)
    ax.set_ylabel('β(λ)'); ax.set_xlabel('Wavelength [µm]')
    ax.legend(fontsize=8); ax.grid(True, alpha=0.2)
    ax.set_title(f'{grating} 3rd-order β(λ)', fontsize=10)

    plt.tight_layout()
    out = os.path.join(OUTPUT_DIR, f'ifu_v4_{grating.lower()}_coeffs.png')
    plt.savefig(out, dpi=150, bbox_inches='tight'); plt.close()
    print(f'  Saved: {out}')


# ── Plot 2: per-source NRS2 validation ────────────────────────────────────────
def plot_source_spectra(grating, wav_grid, ks, als, bts):
    k_itp  = interp1d(wav_grid, ks,  bounds_error=False, fill_value=np.nan)
    al_itp = interp1d(wav_grid, als, bounds_error=False, fill_value=0.0)
    bt_itp = interp1d(wav_grid, bts, bounds_error=False, fill_value=0.0)

    for src in ALL_SOURCES:
        name    = src['name']
        f_cs    = calspec_jy(src['calspec'])
        nrs2_lo = NRS2_LO[grating]

        ext_path = os.path.join(src['ifu_dir'], 'stage3_ext', EXT_TAG[grating])
        # Truth: G235M for G140M NRS2 check; G395M FS for G235M NRS2 check
        truth_path = (src['mast_g395m_fs'] if grating == 'G235M'
                      else os.path.join(src['ifu_dir'], 'stage3_ext',
                                        'f170lp_g235m-f170lp_x1d.fits'))

        w_all, f_all = load_spec(ext_path,   f'{name} stage3_ext')
        w_tru, f_tru = load_spec(truth_path, f'{name} truth')

        fig, ax = plt.subplots(figsize=(10, 5))
        ax.set_title(f'IFU v4 — {name} ({src["sptype"]}) / {grating}  [v4 NNLS coefficients]',
                     fontsize=12)

        wav_fine = np.linspace(min(wav_grid), max(wav_grid), 800)
        ax.plot(wav_fine, f_cs(wav_fine), 'k:', lw=1.0, alpha=0.6, label='CALSPEC model')

        if len(w_tru) > 5:
            truth_label = ('MAST L3 G395M FS (truth)' if grating == 'G235M'
                           else 'IFU stage3_ext G235M (truth)')
            ax.plot(w_tru, f_tru, 'b--', lw=1.5, alpha=0.8, label=truth_label)

        if len(w_all) > 5:
            m = w_all >= nrs2_lo
            w_e, f_e = w_all[m], f_all[m]
            ax.plot(w_e, f_e, color='0.6', lw=1.0, alpha=0.7,
                    label='IFU stage3_ext NRS2 (raw)')
            if len(w_e) > 10:
                k_v  = k_itp(w_e)
                al_v = al_itp(w_e)
                bt_v = bt_itp(w_e)
                f_sc = f_cs(w_e / 2)
                f_tc = f_cs(w_e / 3)
                f_c  = (f_e - al_v * np.where(np.isfinite(f_sc), f_sc, 0.0)
                            - bt_v * np.where(np.isfinite(f_tc), f_tc, 0.0)) / k_v
                good = np.isfinite(f_c) & (k_v > 0.01)
                ax.plot(w_e[good], f_c[good], 'r-', lw=2.0,
                        label='IFU stage3_ext NRS2 corrected (v4 NNLS)')

        ax.set_xlabel('Wavelength [µm]', fontsize=11)
        ax.set_ylabel('Flux [Jy]', fontsize=11)
        ax.legend(fontsize=9); ax.grid(True, alpha=0.2); ax.set_ylim(bottom=0)
        plt.tight_layout()
        out = os.path.join(OUTPUT_DIR, f'ifu_v4_spectra_{name}_{grating}.png')
        plt.savefig(out, dpi=150, bbox_inches='tight'); plt.close()
        print(f'  Saved: {out}')


# ── Plot 3: MAST L3 nominal vs our L3 ext ─────────────────────────────────────
def plot_mast_vs_l3(grating, wav_grid, ks, als, bts):
    k_itp  = interp1d(wav_grid, ks,  bounds_error=False, fill_value=np.nan)
    al_itp = interp1d(wav_grid, als, bounds_error=False, fill_value=0.0)
    bt_itp = interp1d(wav_grid, bts, bounds_error=False, fill_value=0.0)

    fig, axes = plt.subplots(len(ALL_SOURCES), 1,
                              figsize=(12, 4 * len(ALL_SOURCES)), sharex=True)
    fig.suptitle(
        f'IFU v4 — {grating}: MAST Level-3 nominal vs our Level-3 ext (v4 NNLS)',
        fontsize=13)

    for ax, src in zip(axes, ALL_SOURCES):
        name    = src['name']
        f_cs    = calspec_jy(src['calspec'])
        nrs2_lo = NRS2_LO[grating]

        mast_path = os.path.join(src['ifu_dir'], 'stage3',     NOM_TAG[grating])
        ext_path  = os.path.join(src['ifu_dir'], 'stage3_ext', EXT_TAG[grating])

        w_mast, f_mast = load_spec(mast_path, f'{name} stage3 {grating}')
        w_ext,  f_ext  = load_spec(ext_path,  f'{name} stage3_ext {grating}')

        wav_fine = np.linspace(0.6, 5.6, 2000)
        ax.plot(wav_fine, f_cs(wav_fine), 'k:', lw=0.8, alpha=0.5, label='CALSPEC')

        if len(w_mast) > 5:
            ax.plot(w_mast, f_mast, color='steelblue', lw=1.5,
                    label=f'MAST L3 {grating} (stage3, NRS1 only)')

        if len(w_ext) > 5:
            m_nrs1 = w_ext < nrs2_lo
            if m_nrs1.any():
                ax.plot(w_ext[m_nrs1], f_ext[m_nrs1], color='maroon', lw=1.2,
                        alpha=0.9, label='Our L3 ext (NRS1 portion — comparison)')
            m_nrs2 = w_ext >= nrs2_lo
            w_e, f_e = w_ext[m_nrs2], f_ext[m_nrs2]
            if len(w_e) > 10:
                ax.plot(w_e, f_e, color='darkorange', lw=1.0, alpha=0.7,
                        label='Our L3 ext NRS2 (raw)')
                k_v  = k_itp(w_e)
                al_v = al_itp(w_e)
                bt_v = bt_itp(w_e)
                f_sc = f_cs(w_e / 2)
                f_tc = f_cs(w_e / 3)
                f_c  = (f_e - al_v * np.where(np.isfinite(f_sc), f_sc, 0.0)
                            - bt_v * np.where(np.isfinite(f_tc), f_tc, 0.0)) / k_v
                good = np.isfinite(f_c) & (k_v > 0.01)
                ax.plot(w_e[good], f_c[good], 'r-', lw=2.0,
                        label='Our L3 ext NRS2 corrected (v4 NNLS)')

        ax.set_ylabel('Flux [Jy]')
        ax.set_title(f'{name} ({src["sptype"]})', fontsize=11)
        ax.legend(fontsize=8); ax.grid(True, alpha=0.2); ax.set_ylim(bottom=0)

    axes[-1].set_xlabel('Wavelength [µm]', fontsize=11)
    plt.tight_layout()
    out = os.path.join(OUTPUT_DIR, f'ifu_v4_mast_vs_l3_{grating}.png')
    plt.savefig(out, dpi=150, bbox_inches='tight'); plt.close()
    print(f'  Saved: {out}')


# ── Plot 4: IFU v4 vs v3 coefficient comparison ───────────────────────────────
def plot_v4_vs_v3(grating, wav_grid, ks, als, bts):
    """Compare v4 NNLS coefficients to v3 (Parlanti-fixed α/β)."""
    v3_csv = os.path.join(
        BASE, 'nirspec_wavext_work', 'plots', 'Parlanti', 'cal', 'ifu_v3',
        f'coeffs_ifu_v3_{grating}.csv')
    if not os.path.exists(v3_csv):
        print(f'  NOTE: v3 CSV not found ({v3_csv}); skipping v4 vs v3 comparison')
        return

    dat3  = np.loadtxt(v3_csv, delimiter=',', skiprows=1)
    wv3   = dat3[:, 0]; kv3 = dat3[:, 1]; av3 = dat3[:, 2]; bv3 = dat3[:, 3]
    parl  = load_parlanti_coeffs(grating)

    fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
    fig.suptitle(f'IFU v4 vs v3 — {grating} coefficient comparison', fontsize=12)

    # k
    axes[0].plot(wv3,      kv3, color='steelblue', lw=1.5, ls='--', label='k v3 (Parlanti α/β fixed)')
    axes[0].plot(wav_grid, ks,  color='k',         lw=2.0,          label='k v4 (NNLS)')
    if parl:
        axes[0].plot(wav_grid, parl[1](wav_grid), 'r:', lw=1.2, alpha=0.7,
                     label='k Parlanti')
    axes[0].axhline(1.0, color='gray', lw=0.6, ls=':')
    axes[0].set_ylabel('k(λ)'); axes[0].set_ylim(-0.05, 2.5)
    axes[0].legend(fontsize=8); axes[0].grid(True, alpha=0.2)

    # α
    axes[1].plot(wv3, av3, color='steelblue', lw=1.5, ls='--', label='α v3 (Parlanti)')
    axes[1].plot(wav_grid, als, color='darkorange', lw=2.0,    label='α v4 (NNLS)')
    if parl:
        axes[1].plot(wav_grid, parl[2](wav_grid), 'r:', lw=1.2, alpha=0.7,
                     label='α̃ Parlanti')
    axes[1].axhline(0, color='0.5', lw=0.6)
    axes[1].set_ylabel('α(λ)'); axes[1].legend(fontsize=8); axes[1].grid(True, alpha=0.2)

    # β
    axes[2].plot(wv3, bv3, color='steelblue', lw=1.5, ls='--', label='β v3 (Parlanti)')
    axes[2].plot(wav_grid, bts, color='steelblue', lw=2.0, label='β v4 (NNLS)')
    if parl:
        axes[2].plot(wav_grid, parl[3](wav_grid), 'r:', lw=1.2, alpha=0.7,
                     label='β̃ Parlanti')
    axes[2].axhline(0, color='0.5', lw=0.6)
    axes[2].set_ylabel('β(λ)'); axes[2].set_xlabel('Wavelength [µm]')
    axes[2].legend(fontsize=8); axes[2].grid(True, alpha=0.2)

    plt.tight_layout()
    out = os.path.join(OUTPUT_DIR, f'ifu_v4_vs_v3_{grating}.png')
    plt.savefig(out, dpi=150, bbox_inches='tight'); plt.close()
    print(f'  Saved: {out}')


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    print('=== IFU Calibration v4 Solver (NNLS: simultaneous k, α, β) ===')
    print(f'Output: {OUTPUT_DIR}')
    print(f'Sources: P330E (G2V) + G191-B2B (WD) + J1743045 (A8III)')

    results = {}
    for grating in ('G140M', 'G235M'):
        wav_grid, ks, als, bts, k_raw, a_raw, b_raw, active = solve_grating(grating)
        results[grating] = (wav_grid, ks, als, bts, k_raw, a_raw, b_raw, active)

        plot_coeffs(grating, wav_grid, k_raw, ks, a_raw, als, b_raw, bts)
        plot_source_spectra(grating, wav_grid, ks, als, bts)
        plot_mast_vs_l3(grating, wav_grid, ks, als, bts)
        plot_v4_vs_v3(grating, wav_grid, ks, als, bts)

    print('\n=== IFU v4 Done ===\n')

    # Print coefficient summary table
    print(f'{"Grating":<8} {"k med":>7} {"k range":>16} {"α med":>7} {"α max":>7} '
          f'{"β med":>7} {"β max":>7}')
    print('-' * 62)
    for grating in ('G140M', 'G235M'):
        wav_grid, ks, als, bts = results[grating][:4]
        print(f'{grating:<8} {np.median(ks):>7.3f} '
              f'{np.min(ks):>6.3f}–{np.max(ks):<6.3f}  '
              f'{np.median(als):>7.4f} {np.max(als):>7.4f} '
              f'{np.median(bts):>7.4f} {np.max(bts):>7.4f}')

    return results


if __name__ == '__main__':
    main()
