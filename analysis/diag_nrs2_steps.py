"""Diagnose which Spec2Pipeline step flags all NRS2 pixels as DO_NOT_USE."""
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

rate = '/Users/dcoe/NIRSpec/wavext/data/PID1537_G191-B2B/jw01537007001_07101_00003_nrs2_rate.fits'
out = '/tmp/diag_nrs2'
os.makedirs(out, exist_ok=True)

# Run steps one by one, checking DQ after each
steps_to_check = [
    ('wcs_only', {
        'assign_wcs':        {'override_wavelengthrange': WAVELENGTHRANGE},
        'flat_field':        {'skip': True},
        'bkg_subtract':      {'skip': True},
        'imprint_subtract':  {'skip': True},
        'msa_flagging':      {'skip': True},
        'wavecorr':          {'skip': True},
        'pathloss':          {'skip': True},
        'photom':            {'skip': True},
        'extract_1d':        {'skip': True},
        'cube_build':        {'skip': True},
    }),
    ('wcs_flat', {
        'assign_wcs':        {'override_wavelengthrange': WAVELENGTHRANGE},
        'bkg_subtract':      {'skip': True},
        'imprint_subtract':  {'skip': True},
        'msa_flagging':      {'skip': True},
        'wavecorr':          {'skip': True},
        'pathloss':          {'skip': True},
        'photom':            {'skip': True},
        'extract_1d':        {'skip': True},
        'cube_build':        {'skip': True},
    }),
    ('wcs_flat_photom', {
        'assign_wcs':        {'override_wavelengthrange': WAVELENGTHRANGE},
        'bkg_subtract':      {'skip': True},
        'imprint_subtract':  {'skip': True},
        'msa_flagging':      {'skip': True},
        'wavecorr':          {'skip': True},
        'pathloss':          {'skip': True},
        'extract_1d':        {'skip': True},
        'cube_build':        {'skip': True},
    }),
]

for label, cfg in steps_to_check:
    outf = os.path.join(out, f'{label}_cal.fits')
    if not os.path.exists(outf):
        Spec2Pipeline.call(rate, save_results=True, output_dir=out, steps=cfg)
        import glob
        produced = sorted(glob.glob(os.path.join(out, '*cal*.fits')))
        if produced:
            # rename last one
            import shutil
            shutil.move(produced[-1], outf)
    
    if os.path.exists(outf):
        with fits.open(outf) as h:
            dq = h['DQ'].data.astype(int)
            sci = h['SCI'].data.astype(float)
            n_bad_dq = np.sum((dq & 1) > 0)
            n_nan_sci = np.sum(~np.isfinite(sci))
            n_zero_sci = np.sum((sci == 0))
            print(f'{label}: DQ_bad={n_bad_dq}/{dq.size}, SCI_NaN={n_nan_sci}, SCI_zero={n_zero_sci}')
    else:
        print(f'{label}: output file not found')
