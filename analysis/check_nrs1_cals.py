"""Check photomjsr in NRS1 x1d header and calibration chain."""
from astropy.io import fits
import numpy as np

# Check NRS1 x1d keywords
f = '/Users/dcoe/NIRSpec/wavext/data/PID1537_G191-B2B/jw01537007001_07101_00003_nrs1_x1d.fits'
with fits.open(f) as h:
    print('Primary header (PHOTMJSR etc):')
    for k, v in h[0].header.items():
        if 'PHOTM' in k or 'PIXAR' in k or 'BUNIT' in k:
            print(f'  {k}: {v}')
    print('EXTRACT1D header (PHOTMJSR etc):')
    for k, v in h[1].header.items():
        if 'PHOTM' in k or 'PIXAR' in k or 'BUNIT' in k:
            print(f'  {k}: {v}')
    # Check flux values
    fl = h[1].data['FLUX'].astype(float)
    wl = h[1].data['WAVELENGTH'].astype(float)
    mask = (wl > 1.5) & (wl < 1.7) & np.isfinite(fl) & (fl > 0)
    print('NRS1 median flux at 1.5-1.7 um:', np.nanmedian(fl[mask]))

# Compare with CALSPEC G191-B2B
from astropy.io import fits as fitslib
calspec = '/Users/dcoe/NIRSpec/wavext/data/CALSPEC/g191b2b_mod_012.fits'
with fitslib.open(calspec) as h:
    d = h[1].data
    wl_a = d['WAVELENGTH'].astype(float)
    flam = d['FLUX'].astype(float)
    C_ANG_S = 2.99792458e18
    fnu_jy = flam * (wl_a**2) / C_ANG_S * 1e23
    mask2 = (wl_a > 1.5e4) & (wl_a < 1.7e4)
    print('CALSPEC G191-B2B at 1.5-1.7 um (Jy):', np.nanmedian(fnu_jy[mask2]))
