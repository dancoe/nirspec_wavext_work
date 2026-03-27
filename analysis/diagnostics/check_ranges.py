import glob
from astropy.io import fits
import numpy as np

for pid in ['PID1537_G191-B2B', 'PID1538_P330E']:
    print(f"\n--- {pid} ---")
    files = glob.glob(f'/Users/dcoe/NIRSpec/wavext/data/{pid}/*extract_1d.fits')
    for f in files:
        with fits.open(f) as hdul:
            if 'EXTRACT1D' in hdul:
                w = hdul['EXTRACT1D'].data['wavelength']
                # filter out non-positive wavelengths
                w = w[w > 0]
                if len(w) > 0:
                    print(f"{f:70} | {w.min():.3f} - {w.max():.3f}")
