"""Check IFU stage2_ext NRS2 cal file for calibration quality."""
import numpy as np
from astropy.io import fits
import glob

# Find a G140M NRS2 IFU cal file
files = sorted(glob.glob('/Users/dcoe/NIRSpec/wavext/data/IFU/*/stage2_ext/*nrs2_cal.fits'))
print(f'Found {len(files)} IFU NRS2 cal files')

for f in files[:2]:
    with fits.open(f) as h:
        dq = h['DQ'].data.astype(int)
        sci = h['SCI'].data.astype(float)
        good = (dq & 1) == 0
        n_good = good.sum()
        # Read relevant header
        ph = h[0].header
        grating = ph.get('GRATING', '?')
        filt = ph.get('FILTER', '?')
        detector = ph.get('DETECTOR', '?')
        bunit = h['SCI'].header.get('BUNIT', '?')
        exts = [e.name for e in h]
        print(f'{f.split("/")[-1]}: {grating}/{filt}/{detector}')
        print(f'  SCI BUNIT={bunit}, good={n_good}/{dq.size}')
        print(f'  SCI[good] non-NaN={np.sum(np.isfinite(sci[good]))}, med={np.nanmedian(sci[good]):.4g}')
        print(f'  Extensions: {exts}')
        if 'PHOTOM_PS' in exts:
            phot = h['PHOTOM_PS'].data.astype(float)
            print(f'  PHOTOM_PS: min={np.nanmin(phot):.4g}, max={np.nanmax(phot):.4g}, med={np.nanmedian(phot):.4g}')
        print()
