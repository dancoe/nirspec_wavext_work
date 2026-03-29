import os, sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import interp1d
from scipy.optimize import least_squares

# Paths
DATA_DIR = '/Users/dcoe/NIRSpec/wavext/data'
OUTPUT_DIR = '/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/plots/Parlanti/cal/low_order'
os.makedirs(OUTPUT_DIR, exist_ok=True)

def load_spec(path):
    if not os.path.exists(path): return np.array([]), np.array([])
    try:
        with fits.open(path) as hdul:
            data = hdul[1].data
            w, f = data['wavelength'], data['flux']
            msk = (w > 0) & np.isfinite(f) & (f > 0)
            return w[msk], f[msk]
    except: return np.array([]), np.array([])

def get_star_data(grating, name):
    if name == 'G191-B2B':
        d = f'{DATA_DIR}/PID1537_G191-B2B'
        ext = f'{d}/jw01537007001_{"07101_00003" if grating == "G140M" else "09101_00003"}_nrs2_extract_1d.fits'
        t3 = f'{d}/jw01537-o007_t001-s000000001_nirspec_clear-prism-s1600a1-sub2048_x1d.fits'
        nom = f'{d}/jw01537-o007_t001-s000000001_nirspec_f100lp-g140m-s1600a1-sub2048_x1d.fits'
    elif name == 'P330E':
        d = f'{DATA_DIR}/PID1538_P330E'
        ext = f'{d}/jw01538160001_{"06101_00001" if grating == "G140M" else "08101_00005"}_nrs2_extract_1d.fits'
        t3 = f'{d}/jw01538-o107_t002-s000000001_nirspec_clear-prism-s1600a1-sub512_x1d.fits'
        nom = f'{d}/jw01538-o160_t002-s000000001_nirspec_f100lp-g140m-s1600a1-sub2048_x1d.fits'
    
    w_e, f_e = load_spec(ext)
    w_t, f_t = load_spec(t3)
    w_n, f_n = load_spec(nom)
    
    if len(w_e) == 0 or len(w_t) == 0: return None
    
    if len(w_n) > 0:
        msk = (w_n > 1.2) & (w_n < 1.7)
        scale = np.nanmedian(interp1d(w_t, f_t, bounds_error=False)(w_n[msk])) / np.nanmedian(f_n[msk])
    else:
        scale = 1.0
        
    return w_e, f_e * scale, w_t, f_t

def solve_polynomial(grating='G140M', deg_a=0, deg_b=0, deg_k=4):
    print(f"\n Solving {grating} with Degrees: a={deg_a}, b={deg_b}, k={deg_k}")
    data = get_star_data(grating, 'P330E')
    if data is None: return None
    w_obs, f_obs, w_ref, f_ref = data
    f_interp = interp1d(w_ref, f_ref, bounds_error=False, fill_value=0)
    
    if grating == 'G140M':
        fit_msk = (w_obs > 2.0) & (w_obs < 3.1)
    else:
        fit_msk = (w_obs > 3.0) & (w_obs < 5.3)

    w_fit = w_obs[fit_msk]
    y_fit = f_obs[fit_msk]
    f_L = f_interp(w_fit)
    f_L2 = f_interp(w_fit/2)
    f_L3 = f_interp(w_fit/3)
    
    finite_msk = np.isfinite(y_fit) & np.isfinite(f_L) & np.isfinite(f_L2) & np.isfinite(f_L3)
    w_fit, y_fit, f_L, f_L2, f_L3 = w_fit[finite_msk], y_fit[finite_msk], f_L[finite_msk], f_L2[finite_msk], f_L3[finite_msk]

    if len(w_fit) < 10: return None
    w_0 = np.mean(w_fit)
    w_norm = w_fit - w_0
    
    def model(p):
        if deg_k >= 0:
            nk, na = deg_k+1, deg_a+1
            pk, pa, pb = p[:nk], p[nk:nk+na], p[nk+na:]
            k_val, a_val, b_val = np.polyval(pk, w_norm), np.polyval(pa, w_norm), np.polyval(pb, w_norm)
        else:
            na = deg_a + 1
            pa, pb = p[:na], p[na:]
            k_val, a_val, b_val = 1.0, np.polyval(pa, w_norm), np.polyval(pb, w_norm)
        return k_val * f_L + a_val * f_L2 + b_val * f_L3

    def residuals(p):
        res = model(p) - y_fit
        res[~np.isfinite(res)] = 0
        return res
    
    if deg_k >= 0:
        p0 = np.zeros( (deg_k+1) + (deg_a+1) + (deg_b+1) )
        p0[deg_k], p0[deg_k+1+deg_a], p0[deg_k+1+deg_a+1+deg_b] = 1.0, 0.2, 0.05
    else:
        p0 = np.zeros( (deg_a+1) + (deg_b+1) )
        p0[deg_a], p0[deg_a+1+deg_b] = 0.2, 0.05
    
    res = least_squares(residuals, p0, bounds=(0, np.inf))
    p_opt = res.x
    
    if deg_k >= 0:
        nk, na = deg_k+1, deg_a+1
        pk, pa, pb = p_opt[:nk], p_opt[nk:nk+na], p_opt[nk+na:]
        k_val = np.polyval(pk, w_norm)
    else:
        na = deg_a + 1
        pa, pb = p_opt[:na], p_opt[na:]
        k_val = 1.0

    a_val, b_val = np.polyval(pa, w_norm), np.polyval(pb, w_norm)
    print(f"  Median Alpha: {np.nanmedian(a_val):.4f}, Median Beta: {np.nanmedian(b_val):.4f}")
    
    return {'w_fit': w_fit, 'y_fit': y_fit, 'model': model(p_opt), 'k_val': k_val, 'a_val': a_val, 'b_val': b_val}

resolutions = [
    (0, 0, -1, 'Fixed k=1, Constant Alpha/Beta'),
    (1, 1, -1, 'Fixed k=1, Linear Alpha/Beta'),
    (0, 0, 4, 'Varying k, Constant Alpha/Beta'),
]

results = {}
for da, db, dk, label in resolutions:
    results[label] = solve_polynomial(grating='G140M', deg_a=da, deg_b=db, deg_k=dk)

plt.style.use('dark_background')
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), sharex=True)
for label, res in results.items():
    if res is None: continue
    ax1.plot(res['w_fit'], res['a_val'], label=f'α(λ) - {label}', lw=2)
    ax2.plot(res['w_fit'], res['b_val'], label=f'β(λ) - {label}', lw=2)
ax1.set_ylabel('Alpha'), ax1.legend(), ax1.grid(alpha=0.2)
ax2.set_ylabel('Beta'), ax2.legend(), ax2.grid(alpha=0.2)
plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/POLY_COMPARISON_G140M.png', dpi=300)
print(f"Analysis complete. Results saved to {OUTPUT_DIR}")
