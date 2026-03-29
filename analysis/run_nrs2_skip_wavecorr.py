"""Definitive test: run full Spec2 WITHOUT wavecorr on NRS2 FS G140M rate file.
Checks if valid Jy output is produced.
"""
import os, warnings
import numpy as np
from astropy.io import fits
warnings.simplefilter('ignore')

os.environ.setdefault('CRDS_PATH', os.path.expanduser('~/crds_cache'))
os.environ.setdefault('CRDS_SERVER_URL', 'https://jwst-crds.stsci.edu')

from jwst.pipeline import Spec2Pipeline
import jwst
print('JWST:', jwst.__version__)

PARLANTI_DIR = '/Users/dcoe/NIRSpec/wavext/data/parlanti_repo/CRDS_1364'
WAVELENGTHRANGE = os.path.join(PARLANTI_DIR, 'jwst_nirspec_wavelengthrange_0008.asdf')

rate_g140m = '/Users/dcoe/NIRSpec/wavext/data/PID1537_G191-B2B/jw01537007001_07101_00003_nrs2_rate.fits'
rate_g235m = '/Users/dcoe/NIRSpec/wavext/data/PID1537_G191-B2B/jw01537007001_09101_00003_nrs2_rate.fits'

out = '/Users/dcoe/NIRSpec/wavext/data/PID1537_G191-B2B/nrs2_spec2_cal'
os.makedirs(out, exist_ok=True)

# Remove previous (all-NaN) outputs first
import glob
for old in glob.glob(os.path.join(out, '*')):
    os.remove(old)
    print('Removed:', os.path.basename(old))

cfg = {
    'assign_wcs':        {'override_wavelengthrange': WAVELENGTHRANGE},
    'bkg_subtract':      {'skip': True},
    'imprint_subtract':  {'skip': True},
    'msa_flagging':      {'skip': True},
    # Skip wavecorr — it causes all-NaN DQ flagging for extended NRS2
    'wavecorr':          {'skip': True},
    'cube_build':        {'skip': True},
    'extract_1d':        {'skip': False},
}

for rate, label in [(rate_g140m, 'G140M'), (rate_g235m, 'G235M')]:
    print(f'\nProcessing {label}...')
    Spec2Pipeline.call(rate, save_results=True, output_dir=out, steps=cfg)

# Check results
for f in sorted(glob.glob(os.path.join(out, '*_x1d.fits'))):
    with fits.open(f) as h:
        fl = h[1].data['FLUX'].astype(float)
        wl = h[1].data['WAVELENGTH'].astype(float)
        unit = str(h[1].columns['FLUX'].unit)
        n_fin = np.sum(np.isfinite(fl) & (fl > 0))
        med = float(np.nanmedian(fl[np.isfinite(fl) & (fl > 0)])) if n_fin > 0 else float('nan')
    print(f'{os.path.basename(f)}: unit={unit}, wl=[{wl.min():.3f},{wl.max():.3f}], '
          f'n_finite={n_fin}/{len(fl)}, med={med:.4g}')
