"""
Run Spec2Pipeline on NRS2 Fixed Slit rate files to produce flux-calibrated x1d in Jy.

Uses Parlanti (2025) extended reference files to enable calibration at NRS2
wavelengths beyond the nominal detector cutoff:
  - wavelengthrange_0008.asdf : extended assign_wcs coverage
  - fflat_0105 / fflat_0091   : extended f-flat for G140M / G235M
  - sflat_0191 / sflat_0192   : extended s-flat for G140M / G235M NRS2
  - custom photom ext          : extends relresponse to NRS2 wavelengths

The resample_spec step is skipped because extract_1d on the s2d fails to
properly account for the extended-wavelength photom correction, whereas
extracting directly from the cal gives correct Jy flux.

Produces *_cal.fits and *_x1d.fits in <data_dir>/nrs2_spec2_cal/.

Usage:
    micromamba run -n jwst_1.20.2 python run_fs_nrs2_spec2.py
    micromamba run -n jwst_1.20.2 python run_fs_nrs2_spec2.py --pid 1537
"""
import os, sys, glob, argparse, warnings
import numpy as np
from astropy.io import fits
warnings.simplefilter('ignore', RuntimeWarning)

os.environ.setdefault('CRDS_PATH', os.path.expanduser('~/crds_cache'))
os.environ.setdefault('CRDS_SERVER_URL', 'https://jwst-crds.stsci.edu')

from jwst.pipeline import Spec2Pipeline
import jwst
print(f'JWST pipeline: {jwst.__version__}')

BASE = '/Users/dcoe/NIRSpec/wavext/data'
PARLANTI_DIR = '/Users/dcoe/NIRSpec/wavext/data/parlanti_repo/CRDS_1364'
WAVELENGTHRANGE = os.path.join(PARLANTI_DIR, 'jwst_nirspec_wavelengthrange_0008.asdf')

# Custom photom reference: extends relresponse wavelength coverage to NRS2
PHOTOM_NRS2_EXT = os.path.join(
    '/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/reference_files',
    'jwst_nirspec_photom_fs_nrs2_ext.fits'
)

# Parlanti extended flat-field references covering NRS2 wavelengths:
# fflat: fibre-optics flat (wavelength-dependent throughput correction)
# sflat: slit flat (slit-profile wavelength-dependent correction)
# Both have slit_name=ANY so they apply to all FS slits
FFLAT_G140M = os.path.join(PARLANTI_DIR, 'jwst_nirspec_fflat_0105.fits')  # 0.97-5.5 µm
FFLAT_G235M = os.path.join(PARLANTI_DIR, 'jwst_nirspec_fflat_0091.fits')  # 1.66-5.5 µm
SFLAT_G140M = os.path.join(PARLANTI_DIR, 'jwst_nirspec_sflat_0191.fits')  # NRS2, 0.97-5.0 µm
SFLAT_G235M = os.path.join(PARLANTI_DIR, 'jwst_nirspec_sflat_0192.fits')  # NRS2, 1.66-6.0 µm

# Map grating/filter to Parlanti flat overrides
FLAT_OVERRIDES = {
    ('G140M', 'F100LP'): {'override_fflat': FFLAT_G140M, 'override_sflat': SFLAT_G140M},
    ('G235M', 'F170LP'): {'override_fflat': FFLAT_G235M, 'override_sflat': SFLAT_G235M},
}

# PIDs and data directories for Fixed Slit targets
PIDS = {
    '1536': os.path.join(BASE, 'PID1536_J1743045'),
    '1537': os.path.join(BASE, 'PID1537_G191-B2B'),
    '1538': os.path.join(BASE, 'PID1538_P330E'),
    '1492': os.path.join(BASE, 'PID1492'),
    '6644': os.path.join(BASE, 'PID6644_NGC2506G31'),
}

# Only process these grating/filter combinations (extended wavelength range configs)
TARGET_CONFIGS = {('G140M', 'F100LP'), ('G235M', 'F170LP')}


def get_header_kw(path, *keys):
    h = fits.getheader(path)
    return tuple(h.get(k) for k in keys)


def run_spec2_on_rate(rate_file, out_dir, grating, filt):
    """Run Spec2Pipeline on one NRS2 rate file. Returns output x1d path."""
    base = os.path.basename(rate_file).replace('_rate.fits', '')
    out_x1d = os.path.join(out_dir, f'{base}_x1d.fits')

    if os.path.exists(out_x1d):
        print(f'    EXISTS: {os.path.basename(out_x1d)}')
        return out_x1d

    print(f'    Spec2: {base} ...')

    flat_cfg = FLAT_OVERRIDES.get((grating, filt), {})
    if not flat_cfg:
        print(f'      WARNING: no Parlanti flat overrides for {grating}/{filt}, using CRDS defaults')

    spec2_cfg = {
        'assign_wcs': {'override_wavelengthrange': WAVELENGTHRANGE},
        'bkg_subtract':     {'skip': True},
        'imprint_subtract': {'skip': True},
        'msa_flagging':     {'skip': True},
        # Use Parlanti extended flat references (cover NRS2 wavelengths)
        'flat_field': flat_cfg,
        # Skip wavecorr: wavelength correction reference lacks NRS2 coverage
        'wavecorr':    {'skip': True},
        # Custom photom: extends relresponse to NRS2 wavelengths
        'photom':      {'override_photom': PHOTOM_NRS2_EXT},
        # Skip resample_spec: extract_1d from cal produces correct Jy;
        # resample_spec→s2d→extract_1d fails for extended-wavelength photom
        'resample_spec': {'skip': True},
        'cube_build':    {'skip': True},
        'extract_1d':    {'skip': False},
    }

    try:
        Spec2Pipeline.call(
            rate_file,
            save_results=True,
            output_dir=out_dir,
            steps=spec2_cfg,
        )
    except Exception as e:
        print(f'      ERROR: {e}')
        return None

    # Find the output x1d (pipeline may add step suffixes)
    produced = sorted(glob.glob(os.path.join(out_dir, f'{base}*_x1d.fits')))
    if produced:
        print(f'      → {os.path.basename(produced[-1])}')
        return produced[-1]
    elif os.path.exists(out_x1d):
        print(f'      → {os.path.basename(out_x1d)}')
        return out_x1d
    else:
        print(f'      WARNING: no x1d produced for {base}')
        return None


def process_pid(pid):
    data_dir = PIDS[pid]
    out_dir = os.path.join(data_dir, 'nrs2_spec2_cal')
    os.makedirs(out_dir, exist_ok=True)

    rate_files = sorted(glob.glob(os.path.join(data_dir, '*_nrs2_rate.fits')))
    if not rate_files:
        print(f'  [{pid}] No NRS2 rate files found in {data_dir}')
        return []

    results = []
    for rate_file in rate_files:
        grating, filt, exp_type = get_header_kw(rate_file, 'GRATING', 'FILTER', 'EXP_TYPE')
        if exp_type != 'NRS_FIXEDSLIT':
            print(f'  [{pid}] Skipping {os.path.basename(rate_file)} (EXP_TYPE={exp_type})')
            continue
        if (grating, filt) not in TARGET_CONFIGS:
            print(f'  [{pid}] Skipping {os.path.basename(rate_file)} ({grating}/{filt} not in target configs)')
            continue
        print(f'  [{pid}] {grating}/{filt}: {os.path.basename(rate_file)}')
        x1d = run_spec2_on_rate(rate_file, out_dir, grating, filt)
        if x1d:
            results.append(x1d)

    print(f'  [{pid}] Done: {len(results)} x1d files in {out_dir}')
    return results


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pid', default='all', help='PID to process (1536, 1537, 1538, or "all")')
    args = parser.parse_args()

    if not os.path.exists(WAVELENGTHRANGE):
        raise FileNotFoundError(f'Parlanti wavelengthrange not found: {WAVELENGTHRANGE}')

    pids = list(PIDS.keys()) if args.pid == 'all' else [p.strip() for p in args.pid.split(',')]
    for pid in pids:
        print(f'\n=== PID {pid} ===')
        process_pid(pid)

    print('\nAll done.')
