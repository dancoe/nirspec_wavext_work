"""Diagnose the pixel->wavelength conversion in show_MOS_rate_files.

Examines the WCS for the PID1537 G395M NRS2 extract_2d files to understand
why the 5.3 µm limit marker appears at the wrong position.
"""
import sys
import os
sys.path.insert(0, '/Users/dcoe/NIRSpec/wavext/jwst_nirspec_wavext')

from stdatamodels.jwst import datamodels
import numpy as np
from scipy.interpolate import interp1d

files = [
    '/Users/dcoe/NIRSpec/wavext/data/PID1537_G191-B2B/jw01537007001_07101_00003_nrs2_extract_2d.fits',
    '/Users/dcoe/NIRSpec/wavext/data/PID1537_G191-B2B/jw01537007001_09101_00003_nrs2_extract_2d.fits',
    '/Users/dcoe/NIRSpec/wavext/data/PID1537_G191-B2B/jw01537007001_11101_00002_nrs2_extract_2d.fits',
]

for fpath in files:
    if not os.path.exists(fpath):
        print(f"MISSING: {fpath}")
        continue
    print(f"\n=== {os.path.basename(fpath)} ===")
    model = datamodels.open(fpath)
    slit_sizes = [s.xsize for s in model.slits]
    best_idx = int(np.argmax(slit_sizes))
    slit = model.slits[best_idx]
    print(f"  Slit: {slit.name}  xstart={slit.xstart}, xsize={slit.xsize}, ysize={slit.ysize}")
    x_start = slit.xstart
    y_center_local = slit.ysize / 2.0

    wcs = slit.meta.wcs
    det_to_sky = wcs.get_transform('detector', 'world')

    # Reproduce pix2wave exactly as in mos.py
    def pix2wave(x):
        x = np.atleast_1d(np.asarray(x, dtype=float))
        x_local = x - x_start
        res = det_to_sky(x_local, np.full_like(x_local, y_center_local))
        lam = np.array(res[2], dtype=float)
        if np.nanmean(lam) < 1e-4:
            lam *= 1e6
        if np.all(np.isnan(lam)):
            return np.nan * x
        mask = ~np.isnan(lam)
        if np.any(mask) and np.any(~mask):
            lam[~mask] = np.interp(x[~mask], x[mask], lam[mask])
        return lam

    # Print a few key pixels
    print("  Sample wavelengths via pix2wave(x_raw):")
    for x_raw in [0, 200, 500, 1000, 1500, 1700, 1800, 1900, 1950, 2000, 2047]:
        w = pix2wave(np.array([float(x_raw)]))[0]
        print(f"    x_raw={x_raw:5d}  -> {w:.4f} um")

    # Build inv_func as in mos.py
    x_range = np.linspace(-50, 2100, 200)
    w_range = pix2wave(x_range)
    mask = ~np.isnan(w_range)
    print(f"  w_range: min={np.nanmin(w_range):.4f}, max={np.nanmax(w_range):.4f}, NaN count={np.sum(~mask)}")
    if np.sum(mask) > 10:
        idx = np.argsort(w_range[mask])
        uw, uidx = np.unique(w_range[mask][idx], return_index=True)
        ux = x_range[mask][idx][uidx]
        print(f"  Unique wavelength range for inv_func: {uw[0]:.4f} to {uw[-1]:.4f} um ({len(uw)} points)")
        if len(uw) > 1:
            inv_func = interp1d(uw, ux, fill_value="extrapolate", bounds_error=False)
            limit_wave = 5.3
            x_limit = float(inv_func(limit_wave))
            print(f"  inv_func(5.3 um) = x_raw = {x_limit:.1f}")
            # Check: what does pix2wave give at x_limit?
            w_at_limit = pix2wave(np.array([x_limit]))[0]
            print(f"  pix2wave({x_limit:.1f}) = {w_at_limit:.4f} um  (should be 5.3)")
            print(f"  Condition 5.3 in range: min({uw[0]:.3f}) < 5.3 < max({uw[-1]:.3f}) = {uw[0] < 5.3 < uw[-1]}")

print("\nDone.")
