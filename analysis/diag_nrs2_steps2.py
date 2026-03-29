"""Diagnose which step (wavecorr vs pathloss) flags all NRS2 pixels."""
import os, warnings, glob, shutil
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

steps_to_check = [
    ('with_wavecorr', {
        'assign_wcs':        {'override_wavelengthrange': WAVELENGTHRANGE},
        'bkg_subtract':      {'skip': True},
        'imprint_subtract':  {'skip': True},
        'msa_flagging':      {'skip': True},
        'pathloss':          {'skip': True},
        'photom':            {'skip': True},
        'extract_1d':        {'skip': True},
        'cube_build':        {'skip': True},
    }),
    ('with_pathloss', {
        'assign_wcs':        {'override_wavelengthrange': WAVELENGTHRANGE},
        'bkg_subtract':      {'skip': True},
        'imprint_subtract':  {'skip': True},
        'msa_flagging':      {'skip': True},
        'wavecorr':          {'skip': True},
        'photom':            {'skip': True},
        'extract_1d':        {'skip': True},
        'cube_build':        {'skip': True},
    }),
    ('skip_wavecorr_run_photom', {
        'assign_wcs':        {'override_wavelengthrange': WAVELENGTHRANGE},
        'bkg_subtract':      {'skip': True},
        'imprint_subtract':  {'skip': True},
        'msa_flagging':      {'skip': True},
        'wavecorr':          {'skip': True},
        'cube_build':        {'skip': True},
        'extract_1d':        {'skip': False},
    }),
]

for label, cfg in steps_to_check:
    outf = os.path.join(out, f'{label}_cal.fits')
    out_x1d = os.path.join(out, f'{label}_x1d.fits')
    if not os.path.exists(outf):
        Spec2Pipeline.call(rate, save_results=True, output_dir=out, steps=cfg)
        produced_cal = sorted(glob.glob(os.path.join(out, '*_cal.fits')))
        produced_x1d = sorted(glob.glob(os.path.join(out, '*_x1d.fits')))
        # Only move files not already moved
        for pf in produced_cal:
            if not pf.endswith(f'{label}_cal.fits'):
                try:
                    shutil.move(pf, outf)
                except:
                    pass
        for pf in produced_x1d:
            if not pf.endswith(f'{label}_x1d.fits'):
                try:
                    shutil.move(pf, out_x1d)
                except:
                    pass

    if os.path.exists(outf):
        with fits.open(outf) as h:
            dq = h['DQ'].data.astype(int)
            sci = h['SCI'].data.astype(float)
            n_bad_dq = np.sum((dq & 1) > 0)
            n_nan_sci = np.sum(~np.isfinite(sci))
            n_zero_sci = np.sum((sci == 0))
            print(f'{label}: DQ_bad={n_bad_dq}/{dq.size}, SCI_NaN={n_nan_sci}, SCI_zero={n_zero_sci}')
    elif os.path.exists(out_x1d):
        with fits.open(out_x1d) as h:
            fl = h[1].data['FLUX'].astype(float)
            n_finite = np.sum(np.isfinite(fl) & (fl != 0))
            unit = str(h[1].columns['FLUX'].unit) if h[1].columns['FLUX'].unit else 'NONE'
            print(f'{label} (x1d): {n_finite}/{len(fl)} valid, unit={unit}, med={np.nanmedian(fl[np.isfinite(fl) & (fl>0)]):.4g}')
    else:
        print(f'{label}: NO output found')
