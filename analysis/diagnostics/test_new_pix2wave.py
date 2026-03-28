"""Test the new pix2wave/inv_func implementation."""
import sys
sys.path.insert(0, '/Users/dcoe/NIRSpec/wavext/jwst_nirspec_wavext')
from stdatamodels.jwst import datamodels
import numpy as np
from scipy.interpolate import interp1d
import traceback

for fpath, grating in [
    ('/Users/dcoe/NIRSpec/wavext/data/PID1537_G191-B2B/jw01537007001_07101_00003_nrs2_extract_2d.fits', 'G140M'),
    ('/Users/dcoe/NIRSpec/wavext/data/PID1537_G191-B2B/jw01537007001_09101_00003_nrs2_extract_2d.fits', 'G235M'),
    ('/Users/dcoe/NIRSpec/wavext/data/PID1537_G191-B2B/jw01537007001_11101_00002_nrs2_extract_2d.fits', 'G395M'),
]:
    print(f'\n=== {grating} ===')
    try:
        model = datamodels.open(fpath)
        slit = model.slits[0]
        x_start = slit.xstart
        y_center_local = slit.ysize / 2
        wcs_obj = slit.meta.wcs
        det_to_sky = wcs_obj.get_transform('detector', 'world')

        x_lo = x_start
        x_hi = x_start + slit.xsize - 1
        _x_samp = np.linspace(x_lo, x_hi, max(50, slit.xsize))
        _x_local = _x_samp - x_start
        _res = det_to_sky(_x_local, np.full_like(_x_local, y_center_local))
        _lam = np.array(_res[2], dtype=float)
        if np.nanmean(_lam) < 1e-4:
            _lam *= 1e6
        _ok = ~np.isnan(_lam)
        print(f'  Valid samples: {np.sum(_ok)}/{len(_lam)}')

        _sidx = np.argsort(_x_samp[_ok])
        _ux = _x_samp[_ok][_sidx]
        _uw_sorted = _lam[_ok][_sidx]
        _, _uniq = np.unique(_ux, return_index=True)
        _ux = _ux[_uniq]
        _uw_sorted = _uw_sorted[_uniq]

        print(f'  Pixel range: {_ux[0]:.1f}–{_ux[-1]:.1f}')
        print(f'  Wavelength range: {_uw_sorted[0]:.4f}–{_uw_sorted[-1]:.4f} um')

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

        # Test with scalar (as matplotlib calls it)
        print(f'  pix2wave(0.0) = {pix2wave(0.0)}')
        print(f'  pix2wave(1000.0) = {pix2wave(1000.0)}')
        print(f'  inv_func(5.3) = {inv_func(5.3)}')

        # Test with array
        test_x = np.array([0.0, 500.0, 1000.0, 1500.0, 2000.0])
        print(f'  pix2wave(array) = {pix2wave(test_x)}')

        # Test 5.3 limit condition
        w_min = float(_uw_sorted.min())
        w_max = float(_uw_sorted.max())
        limit_wave = 5.3
        print(f'  5.3 in range ({w_min:.4f}, {w_max:.4f})? {w_min < limit_wave < w_max}')
        if w_min < limit_wave < w_max:
            x_limit = float(inv_func(5.3))
            print(f'  x_limit for 5.3um = {x_limit:.1f}')

    except Exception as e:
        print(f'  ERROR: {e}')
        traceback.print_exc()
