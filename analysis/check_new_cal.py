"""Check G235M NRS2 flux quality vs CALSPEC and sflat headers."""
import numpy as np
from astropy.io import fits
from scipy.interpolate import interp1d

base = '/Users/dcoe/NIRSpec/wavext/data/PID1537_G191-B2B/nrs2_spec2_cal/'
calspec = '/Users/dcoe/NIRSpec/wavext/data/CALSPEC/g191b2b_mod_012.fits'
parlanti = '/Users/dcoe/NIRSpec/wavext/data/parlanti_repo/CRDS_1364/'

# Load CALSPEC
c = 3e18  # Ang/s
with fits.open(calspec) as h:
    cs_wl = h[1].data['WAVELENGTH'].astype(float)
    cs_fl = h[1].data['FLUX'].astype(float) * cs_wl**2 / c * 1e23
    cs_interp = interp1d(cs_wl/1e4, cs_fl, kind='linear', fill_value='extrapolate')

# Check G235M NRS2
with fits.open(base + 'jw01537007001_09101_00003_nrs2_x1d.fits') as h:
    wl = h[1].data['WAVELENGTH'].astype(float)
    fl = h[1].data['FLUX'].astype(float)
    dq = h[1].data['DQ'].astype(int)
    npix = h[1].data['NPIXELS'].astype(float)

good = (dq == 0) & np.isfinite(fl) & (fl != 0) & (npix > 0)
print(f'G235M NRS2: nok={good.sum()}/2048, WL=[{wl[good].min():.3f},{wl[good].max():.3f}]')
for w1, w2 in [(3.3, 3.8), (3.8, 4.5), (4.5, 5.0)]:
    m = good & (wl >= w1) & (wl <= w2)
    if m.sum() < 3:
        continue
    obs = np.nanmedian(fl[m])
    cs = np.nanmedian(cs_interp(wl[m]))
    print(f'  {w1}-{w2} µm: obs={obs:.4g} Jy, CS={cs:.4g} Jy, ratio={obs/cs:.3f}')

# Check sflat headers
print()
for fn in ['jwst_nirspec_sflat_0192.fits', 'jwst_nirspec_sflat_0211.fits']:
    with fits.open(parlanti + fn) as h:
        hdr = h[0].header
        det = hdr.get('DETECTOR', '?')
        grating = hdr.get('GRATING', '?')
        filt = hdr.get('FILTER', '?')
        exp = hdr.get('EXP_TYPE', '?')
        descrip = hdr.get('DESCRIP', '?')
        print(f'{fn}: GRATING={grating}, FILTER={filt}, DETECTOR={det}, EXP={exp}')
        print(f'  DESCRIP={descrip}')
