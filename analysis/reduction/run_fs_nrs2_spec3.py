"""
Run FS NRS2 Extended Wavelength Pipeline: Spec2 → Spec3 (Level 3).

This script produces true Level-3 Fixed Slit x1d spectra for the NRS2
extended-wavelength range using Parlanti (2025) reference file overrides
and our custom pipeline patch (NRS2 hard-coded error removed).

Workflow per PID / grating:
  1. Spec2Pipeline on each NRS2 rate file: assign_wcs (ext range) + flat
     (Parlanti sflat/fflat) + photom (ext NRS2 ref) + extract_1d.
     Output: *_nrs2_cal.fits in nrs2_spec2_cal/
  2. Spec3Pipeline on the resulting cal files: resample_spec + extract_1d.
     Output: Level-3 *_nrs2_x1d.fits in nrs2_spec3_ext/

Differences from run_fs_nrs2_spec2.py (v2):
  - resample_spec is now ENABLED in Spec3 (applied to post-photom cal files,
    so wavelength-range issues in Spec2 do not apply here).
  - Spec3 uses all available NRS2 cal files per PID/grating, giving the
    pipeline's standard drizzle-combine and optimal extraction.
  - Even with a single exposure, the Spec3 product is properly rectified onto
    a uniform wavelength grid and carries Level-3 metadata.

Custom pipeline environment:
  export PYTHONPATH=/Users/dcoe/NIRSpec/wavext/jwst_nirspec_wavext
  micromamba run -n jwst_1.20.2 python run_fs_nrs2_spec3.py

Reference:
  Parlanti et al. 2025, arXiv:2512.14844
"""
import os
import sys
import glob
import json
import warnings
import argparse

import numpy as np
from astropy.io import fits

warnings.simplefilter('ignore', RuntimeWarning)

os.environ.setdefault('CRDS_PATH', os.path.expanduser('~/crds_cache'))
os.environ.setdefault('CRDS_SERVER_URL', 'https://jwst-crds.stsci.edu')

from jwst.pipeline import Spec2Pipeline, Spec3Pipeline
from jwst.associations import asn_from_list as afl
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base
import jwst

print(f'JWST pipeline: {jwst.__version__}')

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE         = '/Users/dcoe/NIRSpec/wavext/data'
PARLANTI_DIR = os.path.join(BASE, 'parlanti_repo', 'CRDS_1364')
WAVRANGE     = os.path.join(PARLANTI_DIR, 'jwst_nirspec_wavelengthrange_0008.asdf')
PHOTOM_EXT   = os.path.join(
    '/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/reference_files',
    'jwst_nirspec_photom_fs_nrs2_ext.fits',
)

# Parlanti extended flat references keyed by (GRATING, FILTER)
FFLAT = {
    'G140M': os.path.join(PARLANTI_DIR, 'jwst_nirspec_fflat_0105.fits'),
    'G235M': os.path.join(PARLANTI_DIR, 'jwst_nirspec_fflat_0091.fits'),
}
SFLAT = {
    'G140M': os.path.join(PARLANTI_DIR, 'jwst_nirspec_sflat_0191.fits'),
    'G235M': os.path.join(PARLANTI_DIR, 'jwst_nirspec_sflat_0192.fits'),
}

# PIDs and target configurations
PIDS = {
    '1536': {'dir': os.path.join(BASE, 'PID1536_J1743045'), 'name': 'J1743045'},
    '1537': {'dir': os.path.join(BASE, 'PID1537_G191-B2B'),  'name': 'G191-B2B'},
    '1538': {'dir': os.path.join(BASE, 'PID1538_P330E'),     'name': 'P330E'},
    '6644': {'dir': os.path.join(BASE, 'PID6644_NGC2506G31'), 'name': 'NGC2506-G31'},
}

TARGET_CONFIGS = {('G140M', 'F100LP'), ('G235M', 'F170LP')}

# ── Helpers ────────────────────────────────────────────────────────────────────

def get_hdr(path, *keys):
    h = fits.getheader(path)
    return tuple(h.get(k) for k in keys)


def build_spec3_asn(cal_files, product_name, asn_path):
    """Write a Level-3 ASN for a list of NRS2 cal files."""
    asn = afl.asn_from_list(cal_files, rule=DMS_Level3_Base,
                             product_name=product_name)
    _, serialized = asn.dump()
    with open(asn_path, 'w') as fh:
        fh.write(serialized)
    return asn_path


# ── Stage-2 (single exposure) ─────────────────────────────────────────────────

def run_spec2(rate_file, out_dir, grating, filt):
    """Run Spec2Pipeline on one NRS2 rate file. Returns cal path or None."""
    base    = os.path.basename(rate_file).replace('_rate.fits', '')
    out_cal = os.path.join(out_dir, f'{base}_cal.fits')

    if os.path.exists(out_cal):
        print(f'    Spec2 EXISTS: {os.path.basename(out_cal)}')
        return out_cal

    print(f'    Spec2: {base} …')

    spec2_cfg = {
        'assign_wcs':       {'override_wavelengthrange': WAVRANGE},
        'bkg_subtract':     {'skip': True},
        'imprint_subtract': {'skip': True},
        'msa_flagging':     {'skip': True},
        'flat_field': {
            'override_fflat': FFLAT[grating],
            'override_sflat': SFLAT[grating],
        },
        'wavecorr':     {'skip': True},   # reference lacks NRS2 coverage
        'photom':       {'override_photom': PHOTOM_EXT},
        'resample_spec': {'skip': True},  # skip in Spec2; Spec3 handles combine+resample
        'cube_build':   {'skip': True},
        'extract_1d':   {'skip': True},   # skip per-exposure extraction; use Spec3
    }

    try:
        Spec2Pipeline.call(
            rate_file,
            save_results=True,
            output_dir=out_dir,
            steps=spec2_cfg,
        )
    except Exception as exc:
        print(f'      Spec2 ERROR: {exc}')
        return None

    # Handle pipeline step-suffix naming variants
    produced = sorted(glob.glob(os.path.join(out_dir, f'{base}*_cal.fits')))
    if produced:
        print(f'      → {os.path.basename(produced[-1])}')
        return produced[-1]
    if os.path.exists(out_cal):
        return out_cal
    print(f'      WARNING: no cal.fits produced for {base}')
    return None


# ── Stage-3 (combine exposures) ───────────────────────────────────────────────

def run_spec3(cal_files, out_dir, product_name):
    """Run Spec3Pipeline on a list of cal files. Returns x1d path or None."""
    asn_path  = os.path.join(out_dir, f'{product_name}_l3asn.json')
    out_x1d   = os.path.join(out_dir, f'{product_name}_x1d.fits')

    if os.path.exists(out_x1d):
        print(f'    Spec3 EXISTS: {os.path.basename(out_x1d)}')
        return out_x1d

    build_spec3_asn(cal_files, product_name, asn_path)

    n = len(cal_files)
    print(f'    Spec3: {product_name}  ({n} cal file{"s" if n!=1 else ""})')
    for f in cal_files:
        print(f'      {os.path.basename(f)}')

    spec3_cfg = {
        'master_background':  {'skip': True},
        'outlier_detection':  {'skip': n < 2},   # skip if single exposure
        'resample_spec':      {'skip': False},    # ← enabled at Level 3
        'extract_1d':         {'skip': False},
    }

    try:
        Spec3Pipeline.call(
            asn_path,
            save_results=True,
            output_dir=out_dir,
            steps=spec3_cfg,
        )
    except Exception as exc:
        print(f'      Spec3 ERROR: {exc}')
        return None

    produced = sorted(glob.glob(os.path.join(out_dir, f'{product_name}*_x1d.fits')))
    if produced:
        print(f'      → {os.path.basename(produced[-1])}')
        return produced[-1]
    if os.path.exists(out_x1d):
        return out_x1d
    print(f'      WARNING: no x1d produced for {product_name}')
    return None


# ── Per-PID driver ────────────────────────────────────────────────────────────

def process_pid(pid):
    info     = PIDS[pid]
    data_dir = info['dir']
    s2_dir   = os.path.join(data_dir, 'nrs2_spec2_cal')
    s3_dir   = os.path.join(data_dir, 'nrs2_spec3_ext')
    os.makedirs(s2_dir, exist_ok=True)
    os.makedirs(s3_dir, exist_ok=True)

    rate_files = sorted(glob.glob(os.path.join(data_dir, '*_nrs2_rate.fits')))
    if not rate_files:
        print(f'  [{pid}] No NRS2 rate files found in {data_dir}')
        return

    # Group rate files by (GRATING, FILTER)
    from collections import defaultdict
    groups = defaultdict(list)
    for rf in rate_files:
        grat, filt, exp_type = get_hdr(rf, 'GRATING', 'FILTER', 'EXP_TYPE')
        if exp_type != 'NRS_FIXEDSLIT':
            print(f'  [{pid}] Skipping non-FS file: {os.path.basename(rf)} (EXP_TYPE={exp_type})')
            continue
        if (grat, filt) not in TARGET_CONFIGS:
            print(f'  [{pid}] Skipping {os.path.basename(rf)} ({grat}/{filt} not in target configs)')
            continue
        groups[(grat, filt)].append(rf)

    for (grat, filt), rate_list in sorted(groups.items()):
        print(f'\n  [{pid}] {grat}/{filt}: {len(rate_list)} NRS2 rate file(s)')

        # ── Spec2: produce cal.fits for each rate file ────────────────────────
        cal_files = []
        for rf in rate_list:
            cal = run_spec2(rf, s2_dir, grat, filt)
            if cal and os.path.exists(cal):
                cal_files.append(cal)

        if not cal_files:
            print(f'    No cal files produced; skipping Spec3 for {grat}/{filt}')
            continue

        # ── Spec3: combine cal files → Level-3 x1d ───────────────────────────
        product_name = f'nrs2_l3_{pid}_{grat.lower()}'
        run_spec3(cal_files, s3_dir, product_name)


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description='FS NRS2 Level-3 extended-wavelength pipeline (Spec2 + Spec3)'
    )
    p.add_argument('--pid', default='all',
                   help='PID or comma-separated list (default: all)')
    return p.parse_args()


if __name__ == '__main__':
    args = parse_args()
    pids = list(PIDS.keys()) if args.pid == 'all' else args.pid.split(',')

    print(f'Pipeline: {jwst.__version__}')
    print(f'PYTHONPATH: {os.environ.get("PYTHONPATH", "(not set)")}')
    print(f'Processing PIDs: {pids}')
    print()

    for pid in pids:
        if pid not in PIDS:
            print(f'Unknown PID: {pid}'); continue
        print(f'\n{"="*60}')
        print(f'PID {pid} — {PIDS[pid]["name"]}')
        print(f'{"="*60}')
        process_pid(pid)

    print('\nDone.')
