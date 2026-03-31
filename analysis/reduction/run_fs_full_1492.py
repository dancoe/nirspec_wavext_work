"""
Run Spec2Pipeline on ALL PID 1492 Fixed Slit rate files to produce flux-calibrated cal/x1d in Jy.

Following the WAVEXT.md v7 strategy:
- NRS1 (Nominal): Standard CRDS reference files.
- NRS2 (Extended): Custom Parlanti/v7 reference files for extended wavelength and Jy calibration.

Produces results in data/PID1492/spec2_full_cal/
"""
import os, glob, argparse, warnings
from astropy.io import fits
from jwst.pipeline import Spec2Pipeline
import jwst

warnings.simplefilter('ignore', RuntimeWarning)

os.environ['CRDS_PATH'] = os.path.expanduser('~/crds_cache')
os.environ['CRDS_SERVER_URL'] = 'https://jwst-crds.stsci.edu'

# Paths
BASE = '/Users/dcoe/NIRSpec/wavext/data'
PARLANTI_DIR = os.path.join(BASE, 'parlanti_repo/CRDS_1364')
WAVELENGTHRANGE = os.path.join(PARLANTI_DIR, 'jwst_nirspec_wavelengthrange_0008.asdf')
PHOTOM_NRS2_EXT = '/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/reference_files/jwst_nirspec_photom_fs_nrs2_ext.fits'

# Parlanti extended flat-field references
FFLAT_G140M = os.path.join(PARLANTI_DIR, 'jwst_nirspec_fflat_0105.fits')
FFLAT_G235M = os.path.join(PARLANTI_DIR, 'jwst_nirspec_fflat_0091.fits')
SFLAT_G140M = os.path.join(PARLANTI_DIR, 'jwst_nirspec_sflat_0191.fits')
SFLAT_G235M = os.path.join(PARLANTI_DIR, 'jwst_nirspec_sflat_0192.fits')

FLAT_OVERRIDES = {
    ('G140M', 'F100LP'): {'override_fflat': FFLAT_G140M, 'override_sflat': SFLAT_G140M},
    ('G235M', 'F170LP'): {'override_fflat': FFLAT_G235M, 'override_sflat': SFLAT_G235M},
}

def run_spec2(rate_file, out_dir):
    base = os.path.basename(rate_file).replace('_rate.fits', '')
    h = fits.getheader(rate_file)
    detector = h.get('DETECTOR')
    grating = h.get('GRATING')
    filt = h.get('FILTER')
    
    print(f'Processing {base} ({detector}, {grating}/{filt})...')
    
    # Base configuration
    steps = {
        'bkg_subtract':     {'skip': True},
        'imprint_subtract': {'skip': True},
        'msa_flagging':     {'skip': True},
        'wavecorr':         {'skip': True},
        'resample_spec':    {'skip': True},
        'cube_build':       {'skip': True},
        'extract_1d':       {'skip': False},
    }
    
    # Overrides for NRS2 (Extended Wavelength)
    if detector == 'NRS2':
        steps['assign_wcs'] = {'override_wavelengthrange': WAVELENGTHRANGE}
        steps['photom'] = {'override_photom': PHOTOM_NRS2_EXT}
        if (grating, filt) in FLAT_OVERRIDES:
            steps['flat_field'] = FLAT_OVERRIDES[(grating, filt)]
    
    try:
        Spec2Pipeline.call(
            rate_file,
            save_results=True,
            output_dir=out_dir,
            steps=steps,
        )
        print(f'  DONE: {base}')
    except BaseException as e:
        print(f'  ERROR was FATAL for {base}: {e}')

if __name__ == '__main__':
    data_dir = os.path.join(BASE, 'PID1492/MAST/rate')
    out_dir = os.path.join(BASE, 'PID1492/spec2_full_cal')
    os.makedirs(out_dir, exist_ok=True)
    
    rate_files = sorted(glob.glob(os.path.join(data_dir, '*.fits')))
    print(f'Found {len(rate_files)} rate files.')
    
    for f in rate_files:
        run_spec2(f, out_dir)
