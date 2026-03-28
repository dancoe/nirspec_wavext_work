"""Check the G395M WCS extrapolation for the full 0-2048 range."""
import sys
sys.path.insert(0, '/Users/dcoe/NIRSpec/wavext/jwst_nirspec_wavext')
from stdatamodels.jwst import datamodels
import numpy as np
from scipy.interpolate import interp1d

path = '/Users/dcoe/NIRSpec/wavext/data/PID1537_G191-B2B/jw01537007001_11101_00002_nrs2_extract_2d.fits'
model = datamodels.open(path)
slit = model.slits[0]
print(f'G395M: xstart={slit.xstart}, xsize={slit.xsize}, ysize={slit.ysize}')
x_start = slit.xstart
y_center_local = slit.ysize / 2
wcs_obj = slit.meta.wcs
det_to_sky = wcs_obj.get_transform('detector', 'world')

x_range = np.linspace(-50, 2100, 200)
x_local_range = x_range - x_start
res = det_to_sky(x_local_range, np.full_like(x_local_range, y_center_local))
lam = np.array(res[2], dtype=float)
print(f'lam mean={np.nanmean(lam):.6e}')
if np.nanmean(lam) < 1e-4:
    lam *= 1e6
nan_count = np.sum(np.isnan(lam))
print(f'NaN count: {nan_count}/{len(lam)}')
if not np.all(np.isnan(lam)):
    print(f'lam range: [{np.nanmin(lam):.4f}, {np.nanmax(lam):.4f}]')

# Also check what pix2wave gives for typical bottom-axis tick points
print('\nWhat pix2wave(x_raw) gives for G395M over full detector range:')
for x_raw in [0, 250, 500, 750, 1000, 1250, 1500, 1750, 2000]:
    xl = np.array([float(x_raw - x_start)])
    r = det_to_sky(xl, np.array([y_center_local]))
    lv = float(r[2][0]) if r[2][0] is not None else np.nan
    if not np.isnan(lv) and abs(lv) < 1e-4:
        lv *= 1e6
    print(f'  x_raw={x_raw}  x_local={x_raw - x_start}  lam={lv:.4f}')

# Build uw and ux as in mos.py
mask = ~np.isnan(lam)
print(f'\nValid samples: {np.sum(mask)}/{len(mask)}')
if np.sum(mask) > 10:
    idx = np.argsort(lam[mask])
    uw_raw, uidx = np.unique(lam[mask][idx], return_index=True)
    ux_raw = x_range[mask][idx][uidx]
    print(f'uw range: {uw_raw[0]:.4f} to {uw_raw[-1]:.4f} um ({len(uw_raw)} pts)')
    inv_func = interp1d(uw_raw, ux_raw, fill_value='extrapolate', bounds_error=False)
    x_lim = float(inv_func(5.3))
    print(f'inv_func(5.3) = {x_lim:.2f}')
    condition = (np.min(uw_raw) < 5.3 < np.max(uw_raw))
    print(f'Draw condition (5.3 in range): {condition}')
