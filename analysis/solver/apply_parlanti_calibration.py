#!/usr/bin/env python
"""
Apply Parlanti coefficients to calibrate extended M-grating spectra and
produce a combined calibration plot for PID 1492 FS (saved to plots/Parlanti/cal/FS_1492_cal.png).

Usage:
    python analysis/apply_parlanti_calibration.py
"""
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from stdatamodels.jwst import datamodels
from scipy.interpolate import interp1d

WORKDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
DATA_DIR = os.path.join(WORKDIR, 'data', 'PID1492')
PLOTS_DIR = os.path.join(WORKDIR, 'plots', 'Parlanti', 'cal')
COEFF_FILE = os.path.join(WORKDIR, 'plots', 'parlanti_coefficients.txt')
os.makedirs(PLOTS_DIR, exist_ok=True)


def read_level3(file_path):
    try:
        model = datamodels.open(file_path)
        specs = getattr(model, 'spec', [])
        if not specs and hasattr(model, 'spec_table'):
            specs = [model]
        all_w, all_f = [], []
        for s in specs:
            wav = s.spec_table['wavelength']
            flx = s.spec_table['flux']
            m = (wav > 0) & np.isfinite(flx) & (flx > 0)
            all_w.append(wav[m])
            all_f.append(flx[m])
        if len(all_w):
            w = np.concatenate(all_w)
            f = np.concatenate(all_f)
            idx = np.argsort(w)
            return w[idx], f[idx]
    except Exception as e:
        print('Error reading', file_path, e)
    return np.array([]), np.array([])


def parse_coefficients(coeff_file):
    # Parse the plain-text file produced earlier
    coeffs = {}
    with open(coeff_file) as f:
        lines = [l.rstrip('\n') for l in f]
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith('Grating:'):
            gr = line.split(':',1)[1].strip()
            # Find k,a,b blocks
            k, a, b = [], [], []
            # advance to "k(\u03bb)" header
            while i < len(lines) and 'k(\u03bb)' not in lines[i] and 'k(λ)' not in lines[i]:
                i += 1
            i += 1
            # read k lines
            while i < len(lines) and lines[i].strip():
                parts = lines[i].split(':')
                if len(parts) == 2 and 'λ' in parts[0]:
                    val = float(parts[1].strip())
                    k.append(val)
                i += 1
            # find a header
            while i < len(lines) and 'a(\u03bb)' not in lines[i] and 'a(λ)' not in lines[i]:
                i += 1
            i += 1
            while i < len(lines) and lines[i].strip():
                parts = lines[i].split(':')
                if len(parts) == 2 and 'λ' in parts[0]:
                    val = float(parts[1].strip())
                    a.append(val)
                i += 1
            # find b
            while i < len(lines) and 'b(\u03bb)' not in lines[i] and 'b(λ)' not in lines[i]:
                i += 1
            i += 1
            while i < len(lines) and lines[i].strip():
                parts = lines[i].split(':')
                if len(parts) == 2 and 'λ' in parts[0]:
                    val = float(parts[1].strip())
                    b.append(val)
                i += 1
            # store
            if k or a or b:
                # ensure numpy arrays in highest-to-lowest order
                coeffs[gr] = {'k': np.array(k), 'a': np.array(a), 'b': np.array(b)}
        else:
            i += 1
    return coeffs


def make_interp(wav, flux):
    return interp1d(wav, flux, kind='linear', bounds_error=False, fill_value=np.nan)


def apply_calibration(s_wav, s_flux, prism_interp, coeffs):
    # Evaluate polynomial coeffs (already highest->lowest)
    k = np.polyval(coeffs['k'], s_wav)
    a = np.polyval(coeffs['a'], s_wav)
    b = np.polyval(coeffs['b'], s_wav)
    f1 = prism_interp(s_wav)
    f2 = prism_interp(s_wav / 2.0)
    f3 = prism_interp(s_wav / 3.0)
    # Compute calibrated baseline estimate: f_cal = (S - a f2 - b f3) / k
    denom_mask = np.abs(k) > 1e-8
    f_cal = np.full_like(s_flux, np.nan)
    f_cal[denom_mask] = (s_flux[denom_mask] - a[denom_mask] * f2[denom_mask] - b[denom_mask] * f3[denom_mask]) / k[denom_mask]
    return f_cal, k, a, b


def main():
    # locate files
    prism_file = glob.glob(os.path.join(DATA_DIR, '*clear-prism*_x1d.fits'))
    prism_file = prism_file[0] if prism_file else None
    if prism_file is None:
        raise FileNotFoundError('PRISM level3 not found')

    prism_w, prism_f = read_level3(prism_file)
    prism_interp = make_interp(prism_w, prism_f)

    # nominal M-grating x1d files (top-level nominal x1d products)
    nominal_files = {
        'G140M': glob.glob(os.path.join(DATA_DIR, '*f100lp-g140m*_x1d.fits')),
        'G235M': glob.glob(os.path.join(DATA_DIR, '*f170lp-g235m*_x1d.fits')),
        'G395M': glob.glob(os.path.join(DATA_DIR, '*f290lp-g395m*_x1d.fits')),
    }
    nominal = {}
    for g, files in nominal_files.items():
        if files:
            w, f = read_level3(files[0])
            nominal[g] = (w, f)

    # extended (extracted) S(λ) files
    ext_files = {
        'G140M': glob.glob(os.path.join(DATA_DIR, '*g140m*extract_1d.fits')),
        'G235M': glob.glob(os.path.join(DATA_DIR, '*g235m*extract_1d.fits')),
        'G395M': glob.glob(os.path.join(DATA_DIR, '*g395m*extract_1d.fits')),
    }
    extended = {}
    for g, files in ext_files.items():
        if files:
            w, f = read_level3(files[0])
            extended[g] = (w, f)

    coeffs = parse_coefficients(COEFF_FILE)

    calibrated = {}
    for g in ['G140M', 'G235M', 'G395M']:
        if g in extended and g in coeffs:
            s_w, s_f = extended[g]
            fcal, k, a, b = apply_calibration(s_w, s_f, prism_interp, coeffs[g])
            calibrated[g] = (s_w, s_f, fcal, k, a, b)
            print(f'Calibrated {g}: {len(s_w)} points')

    # Create combined plot similar to pre-cal but overlay calibrated (red)
    plt.figure(figsize=(14,6))

    # Plot nominal M-gratings (if available) as darker lines
    for g, (w,f) in nominal.items():
        if g == 'G140M': plt.plot(w, f, color='black', lw=1.0, label='Nominal G140M (x1d)')
        if g == 'G235M': plt.plot(w, f, color='black', lw=1.0, label='Nominal G235M (x1d)')
        if g == 'G395M': plt.plot(w, f, color='black', lw=1.0, label='Nominal G395M (x1d)')

    # Plot PRISM baseline
    plt.plot(prism_w, prism_f, color='k', linewidth=1.2, label='PRISM (Reference f(λ))')

    color_map = {'G140M': 'tab:blue', 'G235M': 'tab:orange', 'G395M': 'tab:red'}

    # plot extended uncalibrated (light colors) and calibrated (red)
    for g, data in extended.items():
        s_w, s_f = data
        plt.plot(s_w, s_f, color=color_map.get(g,'gray'), alpha=0.15, linewidth=0.8, label=f'{g} Extended (uncal)')
    for g, vals in calibrated.items():
        s_w, s_f, fcal, k, a, b = vals
        # fcal is calibrated baseline; overlay as red "corrected" spectrum
        plt.plot(s_w, fcal, color='red', linewidth=1.0, label=f'{g} Extended (calibrated)')

    plt.yscale('log')
    plt.ylim(1e-4, 1e2)
    plt.xlim(0.6, 5.6)
    plt.xlabel('Wavelength (µm)')
    plt.ylabel('Flux')
    plt.title('NIRSpec Calibrated Extended Spectra (PID 1492 FS)')
    plt.legend(loc='upper right', fontsize='small')
    out = os.path.join(PLOTS_DIR, 'FS_1492_cal.png')
    plt.tight_layout()
    plt.savefig(out, dpi=150)
    print('Saved calibrated plot to', out)

if __name__ == '__main__':
    main()
