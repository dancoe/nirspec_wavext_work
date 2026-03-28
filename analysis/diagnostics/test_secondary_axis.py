"""Test secondary_xaxis with the new pix2wave/inv_func."""
import sys
import os
sys.path.insert(0, '/Users/dcoe/NIRSpec/wavext/jwst_nirspec_wavext')
os.environ.setdefault('MPLBACKEND', 'Agg')

from stdatamodels.jwst import datamodels
from astropy.io import fits
import numpy as np
from scipy.interpolate import interp1d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import traceback

files = [
    ('/Users/dcoe/NIRSpec/wavext/data/PID1537_G191-B2B/jw01537007001_07101_00003_nrs2_rate.fits',
     '/Users/dcoe/NIRSpec/wavext/data/PID1537_G191-B2B/jw01537007001_07101_00003_nrs2_extract_2d.fits', 'G140M'),
    ('/Users/dcoe/NIRSpec/wavext/data/PID1537_G191-B2B/jw01537007001_09101_00003_nrs2_rate.fits',
     '/Users/dcoe/NIRSpec/wavext/data/PID1537_G191-B2B/jw01537007001_09101_00003_nrs2_extract_2d.fits', 'G235M'),
    ('/Users/dcoe/NIRSpec/wavext/data/PID1537_G191-B2B/jw01537007001_11101_00002_nrs2_rate.fits',
     '/Users/dcoe/NIRSpec/wavext/data/PID1537_G191-B2B/jw01537007001_11101_00002_nrs2_extract_2d.fits', 'G395M'),
]

fig, axs = plt.subplots(3, 1, figsize=(12, 9), layout='constrained')

for ax, (rate_path, ext_path, grating) in zip(axs, files):
    with fits.open(rate_path) as h:
        data = h['SCI'].data
    ax.imshow(data, origin='lower', aspect='auto', interpolation='nearest')
    ax.set_title(os.path.basename(rate_path))

    # Build wavelength axis
    model = datamodels.open(ext_path)
    slit = model.slits[0]
    x_start = slit.xstart
    y_center_local = slit.ysize / 2
    det_to_sky = slit.meta.wcs.get_transform('detector', 'world')

    x_lo = x_start
    x_hi = x_start + slit.xsize - 1
    _x_samp = np.linspace(x_lo, x_hi, max(50, slit.xsize))
    _x_local = _x_samp - x_start
    _res = det_to_sky(_x_local, np.full_like(_x_local, y_center_local))
    _lam = np.array(_res[2], dtype=float)
    if np.nanmean(_lam) < 1e-4:
        _lam *= 1e6
    _ok = ~np.isnan(_lam)

    _sidx = np.argsort(_x_samp[_ok])
    _ux = _x_samp[_ok][_sidx]
    _uw_sorted = _lam[_ok][_sidx]
    _, _uniq = np.unique(_ux, return_index=True)
    _ux = _ux[_uniq]
    _uw_sorted = _uw_sorted[_uniq]

    _p2w = interp1d(_ux, _uw_sorted, fill_value="extrapolate", bounds_error=False)
    _widx = np.argsort(_uw_sorted)
    _uw_uniq_sorted = _uw_sorted[_widx]
    _ux_for_inv, _u2 = np.unique(_uw_uniq_sorted, return_index=True)
    _ux_inv = _ux[_widx][_u2]
    _w2p = interp1d(_ux_for_inv, _ux_inv, fill_value="extrapolate", bounds_error=False)

    def pix2wave(x, _f=_p2w):
        return _f(np.asarray(x, dtype=float))
    def inv_func(w, _f=_w2p):
        return _f(np.asarray(w, dtype=float))

    try:
        ax_top = ax.secondary_xaxis('top', functions=(pix2wave, inv_func))
        ax_top.set_xlabel(r'Wavelength ($\mu$m)')
        print(f'{grating}: secondary axis OK')

        # Mark 5.3 um
        w_min, w_max = float(_uw_sorted.min()), float(_uw_sorted.max())
        if w_min < 5.3 < w_max:
            x_limit = float(inv_func(5.3))
            ax.axvline(x_limit, color='k', ls='--', lw=1, alpha=0.7)
            ax.text(x_limit, 0.95, ' 5.3 µm limit', transform=ax.get_xaxis_transform(),
                    ha='left', va='top', fontsize=7, color='k', alpha=0.7)
            print(f'{grating}: 5.3um marker at x={x_limit:.1f}')
    except Exception as e:
        print(f'{grating}: axis error: {e}')
        traceback.print_exc()

try:
    plt.savefig('/tmp/test_secondary_axis.png', dpi=100)
    print('Figure saved to /tmp/test_secondary_axis.png')
except Exception as e:
    print(f'Save error: {e}')
    traceback.print_exc()
