"""
NIRSpec IFU Reduction Script
============================

Reduces NIRSpec IFU data (standard stars G191-B2B from PID 1537 and P330E from PID 1538)
through JWST pipeline Stages 1-3, following the STScI IFU pipeline notebook:

  https://github.com/spacetelescope/jwst-pipeline-notebooks/blob/main/notebooks/NIRSPEC/IFU/JWPipeNB-NIRSpec-IFU.ipynb

These are the same calibration stars used by Parlanti et al. (arXiv:2512.14844) to
derive the extended-wavelength correction coefficients k(λ), α̃(λ), β̃(λ) for
G140M/F100LP and G235M/F170LP in IFU mode.

Targets G140M/F100LP and G235M/F170LP grating/filter combinations for extension analysis.

Usage (activate jwst_1.20.2 environment first):
    micromamba activate jwst_1.20.2
    export CRDS_PATH=~/crds_cache
    export CRDS_SERVER_URL=https://jwst-crds.stsci.edu
    export PYTHONPATH=/Users/dcoe/NIRSpec/wavext/jwst_nirspec_wavext
    python reduce_ifu.py [--pid 1537|1538|all] [--stage 1|2|3|all]
                         [--download] [--skip-existing] [--gratings G140M,G235M]

Author: Dan Coe (following Parlanti et al. 2025 methodology)
"""

import os
import sys
import glob
import copy
import json
import time
import argparse
import warnings
import numpy as np
from astropy.io import fits

warnings.simplefilter("ignore", RuntimeWarning)

# ---------------------------------------------------------------------------
# Environment setup – must be done before importing jwst/crds
# ---------------------------------------------------------------------------

if os.getenv('CRDS_PATH') is None:
    os.environ['CRDS_PATH'] = os.path.join(os.path.expanduser('~'), 'crds_cache')
if os.getenv('CRDS_SERVER_URL') is None:
    os.environ['CRDS_SERVER_URL'] = 'https://jwst-crds.stsci.edu'

print(f"CRDS local filepath: {os.environ['CRDS_PATH']}")
print(f"CRDS file server:    {os.environ['CRDS_SERVER_URL']}")

# ---------------------------------------------------------------------------
# JWST pipeline imports
# ---------------------------------------------------------------------------

import crds
from crds.client import api
from stpipe import crds_client
from jwst.pipeline import Detector1Pipeline   # calwebb_detector1
from jwst.pipeline import Spec2Pipeline       # calwebb_spec2
from jwst.pipeline import Spec3Pipeline       # calwebb_spec3
from jwst.extract_1d import Extract1dStep
from jwst import datamodels
from jwst.associations import asn_from_list as afl
from jwst.associations.lib.rules_level2_base import DMSLevel2bBase
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base
import jwst

print(f"JWST pipeline version: {jwst.__version__}")

# ---------------------------------------------------------------------------
# Global configuration
# ---------------------------------------------------------------------------

BASE_DATA_DIR = '/Users/dcoe/NIRSpec/wavext/data/IFU'

# Program → target name mapping
PROGRAMS = {
    '1537': 'G191-B2B',
    '1538': 'P330E',
}

# Observation IDs where IFU mode with G140M+G235M was observed.
# These are confirmed from MAST – adjust if needed.
# Format: PID -> list of (obs_id, filter, grating) tuples for IFU
IFU_OBS = {
    '1537': [
        # G191-B2B: multi-grating IFU observations
        # All grating configs in one observation or separate – check MAST
        ('007', 'F100LP', 'G140M'),
        ('007', 'F170LP', 'G235M'),
        ('007', 'F290LP', 'G395M'),
    ],
    '1538': [
        # P330E: multi-grating IFU observations
        ('160', 'F100LP', 'G140M'),
        ('160', 'F170LP', 'G235M'),
        ('160', 'F290LP', 'G395M'),
    ],
}

# Target gratings to process (for Parlanti-style extended-wavelength analysis)
TARGET_GRATINGS = ['G140M', 'G235M']

EXP_TYPE = 'NRS_IFU'

# ---------------------------------------------------------------------------
# Helper: get_matching (from STScI notebook)
# ---------------------------------------------------------------------------

def get_matching(files, detector, filt, grating, exp_type):
    """
    Filter a list of FITS files to find those matching detector/filter/grating
    and split into regular vs. imprint exposures.
    """
    files_regular, files_imprint = [], []
    for file in files:
        try:
            hdr = fits.getheader(file)
        except Exception:
            continue
        if hdr.get('EXP_TYPE') != exp_type:
            continue
        if (hdr.get('DETECTOR') == detector and
                hdr.get('FILTER') == filt and
                hdr.get('GRATING') == grating):
            is_imprt = hdr.get('IS_IMPRT', False)
            (files_imprint if is_imprt else files_regular).append(file)
    return files_regular, files_imprint


def match_gwa(file1, file2):
    """Check if GWA tilt values match closely enough to be associated."""
    hdr1 = fits.getheader(file1)
    hdr2 = fits.getheader(file2)
    return np.allclose(
        (hdr1.get('GWA_XTIL', 0), hdr1.get('GWA_YTIL', 0)),
        (hdr2.get('GWA_XTIL', 0), hdr2.get('GWA_YTIL', 0)),
        atol=1e-8, rtol=0,
    )


# ---------------------------------------------------------------------------
# Helper: write Level 2 ASN file (following STScI notebook)
# ---------------------------------------------------------------------------

def writel2asn(onescifile, allscifiles, bgfiles, asnfile, product_name):
    """
    Create a Level 2 association file for a single IFU science exposure.
    Handles both nodded and dithered patterns.
    Returns True if the ASN was written successfully.
    """
    try:
        hdr = fits.getheader(onescifile)
    except Exception as e:
        print(f"  [writel2asn] Cannot read header of {onescifile}: {e}")
        return False

    if hdr.get('EXP_TYPE') != EXP_TYPE:
        return False
    if hdr.get('IS_IMPRT', False):
        return False   # Do not process imprint exposures alone

    detector  = hdr['DETECTOR']
    filt      = hdr['FILTER']
    grating   = hdr['GRATING']
    patttype  = hdr.get('PATTTYPE', 'NONE')
    pattnum   = hdr.get('PATT_NUM', 1)

    asn = afl.asn_from_list([onescifile], rule=DMSLevel2bBase,
                             product_name=product_name)
    asn.data['program'] = hdr.get('PROGRAM', '')

    use_sci, use_sci_imprint = get_matching(allscifiles, detector, filt,
                                            grating, EXP_TYPE)
    use_bg,  use_bg_imprint  = (get_matching(bgfiles, detector, filt,
                                             grating, EXP_TYPE)
                                 if bgfiles else ([], []))

    members = asn['products'][0]['members']

    is_nod = 'NOD' in patttype.split('-')
    if is_nod:
        # Nodded: other dither positions used as pixel-by-pixel backgrounds
        for file in use_sci:
            if fits.getval(file, 'PATT_NUM') != pattnum:
                members.append({'expname': file, 'exptype': 'background'})
        for file in use_sci_imprint:
            if match_gwa(onescifile, file):
                members.append({'expname': file, 'exptype': 'imprint'})
        for file in use_sci + use_sci_imprint:
            members.append({'expname': file, 'exptype': 'selfcal'})
    else:
        # Dithered: dedicated backgrounds (if any)
        for file in use_bg:
            members.append({'expname': file, 'exptype': 'background'})
        for file in use_sci_imprint:
            if match_gwa(onescifile, file):
                members.append({'expname': file, 'exptype': 'imprint'})
        for file in use_bg_imprint:
            if match_gwa(use_bg[0] if use_bg else onescifile, file):
                members.append({'expname': file, 'exptype': 'imprint'})
        for file in use_sci + use_sci_imprint + use_bg + use_bg_imprint:
            members.append({'expname': file, 'exptype': 'selfcal'})

    _, serialized = asn.dump()
    with open(asnfile, 'w') as f:
        f.write(serialized)
    return True


# ---------------------------------------------------------------------------
# Helper: write Level 3 ASN file (following STScI notebook)
# ---------------------------------------------------------------------------

def writel3asn(scifiles, bgfiles, asn_dir):
    """
    Create Level 3 association files, one per GRATING/FILTER combination.
    scifiles: list of _cal.fits paths from Stage 2.
    bgfiles:  list of _x1d.fits background paths (empty list if none).
    Returns list of ASN file paths.
    """
    from collections import defaultdict
    grouped = defaultdict(lambda: {'sci': [], 'bg': []})

    for f in scifiles:
        try:
            k = (fits.getval(f, 'FILTER'), fits.getval(f, 'GRATING'))
            grouped[k]['sci'].append(f)
        except Exception:
            continue
    for f in bgfiles:
        try:
            k = (fits.getval(f, 'FILTER'), fits.getval(f, 'GRATING'))
            grouped[k]['bg'].append(f)
        except Exception:
            continue

    asn_files = []
    for (filt, grat), files in grouped.items():
        name = f"{filt}_{grat}".lower()
        asnfile = os.path.join(asn_dir, f"{name}_l3asn.json")
        asn = afl.asn_from_list(files['sci'], rule=DMS_Level3_Base,
                                 product_name=name)
        for bg in files['bg']:
            asn['products'][0]['members'].append({'expname': bg,
                                                   'exptype': 'background'})
        _, serialized = asn.dump()
        with open(asnfile, 'w') as f:
            f.write(serialized)
        asn_files.append(asnfile)
        print(f"  Written L3 ASN: {asnfile}")

    return asn_files


# ---------------------------------------------------------------------------
# Stage 1: Detector1Pipeline  (uncal → rate)
# ---------------------------------------------------------------------------

def run_stage1(uncal_dir, det1_dir, skip_existing=True):
    """Run calwebb_detector1 on all *_uncal.fits in uncal_dir."""
    uncal_files = sorted(glob.glob(os.path.join(uncal_dir, '*_uncal.fits')))
    if not uncal_files:
        print(f"  [Stage1] No uncal.fits in {uncal_dir}")
        return []

    # Stage 1 configuration (following STScI notebook recommendations for IFU)
    det1dict = {}
    det1dict['group_scale'] = {}
    det1dict['dq_init']     = {}
    det1dict['saturation']  = {}
    det1dict['superbias']   = {}
    det1dict['refpix']      = {}
    det1dict['linearity']   = {}
    det1dict['dark_current'] = {}
    det1dict['jump'] = {
        'maximum_cores': 'half',
        'expand_large_events': True,
        'expand_factor': 3,
    }
    det1dict['clean_flicker_noise'] = {}
    det1dict['ramp_fit']   = {}
    det1dict['gain_scale'] = {}

    rate_files = []
    for uncal_file in uncal_files:
        base = os.path.basename(uncal_file).replace('_uncal.fits', '_rate.fits')
        out  = os.path.join(det1_dir, base)
        if skip_existing and os.path.exists(out):
            print(f"  [Stage1] Already exists: {base}")
            rate_files.append(out)
            continue
        print(f"  [Stage1] Processing {os.path.basename(uncal_file)}")
        try:
            result = Detector1Pipeline.call(
                uncal_file,
                save_results=True,
                steps=det1dict,
                output_dir=det1_dir,
            )
        except Exception as e:
            print(f"    ERROR: {e}")
            continue
        # Collect output rate files
        produced = sorted(glob.glob(os.path.join(det1_dir, '*_rate.fits')))
        rate_files.extend(produced)

    return sorted(set(rate_files))


# ---------------------------------------------------------------------------
# Stage 2: Spec2Pipeline  (rate → cal)
# ---------------------------------------------------------------------------

def run_stage2(rate_files, asn_dir, spec2_dir, bg_rate_files=None,
               skip_existing=True):
    """
    Build per-exposure Level 2 ASN files and run calwebb_spec2.
    Returns list of *_cal.fits output files.
    """
    if bg_rate_files is None:
        bg_rate_files = []

    # Stage 2 configuration (following STScI notebook)
    spec2dict = {}
    spec2dict['assign_wcs']   = {}
    spec2dict['badpix_selfcal'] = {}
    spec2dict['msa_flagging'] = {}
    spec2dict['nsclean']      = {}
    spec2dict['imprint_subtract'] = {}
    spec2dict['bkg_subtract'] = {'skip': True}   # no dedicated bg for now
    spec2dict['srctype']      = {}
    spec2dict['wavecorr']     = {}
    spec2dict['flat_field']   = {}
    spec2dict['pathloss']     = {}
    spec2dict['photom']       = {}
    spec2dict['pixel_replace'] = {}
    spec2dict['cube_build']   = {'skip': True}   # skip, done in Stage 3
    spec2dict['extract_1d']   = {'skip': True}   # skip, done in Stage 3

    cal_files = []
    for rate_file in sorted(set(rate_files)):
        base = os.path.basename(rate_file)
        prod = base.replace('_rate.fits', '')
        asnfile = os.path.join(asn_dir, base.replace('_rate.fits', '_l2asn.json'))
        out_cal = os.path.join(spec2_dir, base.replace('_rate.fits', '_cal.fits'))

        if skip_existing and os.path.exists(out_cal):
            print(f"  [Stage2] Already exists: {os.path.basename(out_cal)}")
            cal_files.append(out_cal)
            continue

        print(f"  [Stage2] Processing {base}")
        try:
            written = writel2asn(rate_file, rate_files, bg_rate_files,
                                 asnfile, prod)
            if not written:
                print(f"    Skipped (not IFU or imprint): {base}")
                continue
            result = Spec2Pipeline.call(
                asnfile,
                save_results=True,
                steps=spec2dict,
                output_dir=spec2_dir,
            )
        except Exception as e:
            print(f"    ERROR: {e}")
            continue
        produced = sorted(glob.glob(os.path.join(spec2_dir, '*_cal.fits')))
        cal_files.extend(produced)

    return sorted(set(cal_files))


# ---------------------------------------------------------------------------
# Stage 3: Spec3Pipeline  (cal → s3d + x1d)
# ---------------------------------------------------------------------------

def run_stage3(cal_files, asn_dir, spec3_dir, bg_x1d_files=None,
               skip_existing=True):
    """
    Build Level 3 ASN files (per GRATING/FILTER) and run calwebb_spec3.
    Returns list of *_s3d.fits and *_x1d.fits output files.
    """
    if bg_x1d_files is None:
        bg_x1d_files = []
    if not cal_files:
        print("  [Stage3] No cal.fits files to process.")
        return [], []

    # Stage 3 configuration (following STScI notebook)
    spec3dict = {}
    spec3dict['assign_mtwcs']       = {}
    spec3dict['master_background']  = {'skip': True}
    spec3dict['outlier_detection']  = {'kernel_size': '3 3'}
    spec3dict['pixel_replace']      = {}
    spec3dict['cube_build']         = {}
    spec3dict['extract_1d']         = {}

    # Group cal files by filter/grating and check if already done
    from collections import defaultdict
    groups = defaultdict(list)
    for f in cal_files:
        try:
            k = (fits.getval(f, 'FILTER'), fits.getval(f, 'GRATING'))
            groups[k].append(f)
        except Exception:
            pass

    s3d_files, x1d_files = [], []
    asn_files = writel3asn(cal_files, bg_x1d_files, asn_dir)

    for asnf in sorted(asn_files):
        name = os.path.basename(asnf).replace('_l3asn.json', '')
        out_s3d = os.path.join(spec3_dir, f"{name}_s3d.fits")
        if skip_existing and os.path.exists(out_s3d):
            print(f"  [Stage3] Already exists: {os.path.basename(out_s3d)}")
            s3d_files.append(out_s3d)
            continue
        print(f"  [Stage3] Running on {os.path.basename(asnf)}")
        try:
            result = Spec3Pipeline.call(
                asnf,
                save_results=True,
                steps=spec3dict,
                output_dir=spec3_dir,
            )
        except Exception as e:
            print(f"    ERROR: {e}")
            continue

    s3d_files  = sorted(glob.glob(os.path.join(spec3_dir, '*_s3d.fits')))
    x1d_files  = sorted(glob.glob(os.path.join(spec3_dir, '*_x1d.fits')))
    print(f"  [Stage3] S3D cubes:  {[os.path.basename(f) for f in s3d_files]}")
    print(f"  [Stage3] X1D spectra:{[os.path.basename(f) for f in x1d_files]}")
    return s3d_files, x1d_files


# ---------------------------------------------------------------------------
# MAST download
# ---------------------------------------------------------------------------

def download_ifu_data(pid, target_name, base_dir, gratings_filter=None,
                      download_type='rate'):
    """
    Query MAST for IFU _rate.fits files from the given PID and download.
    download_type: 'rate' (Stage 1 products) or 'uncal' (raw).
    Returns list of downloaded file paths.
    """
    from astroquery.mast import Observations

    if gratings_filter is None:
        gratings_filter = ['F100LP;G140M', 'F170LP;G235M']

    uncal_dir = os.path.join(base_dir, 'uncal')
    det1_dir  = os.path.join(base_dir, 'stage1')
    os.makedirs(uncal_dir, exist_ok=True)
    os.makedirs(det1_dir,  exist_ok=True)

    target_dir = det1_dir if download_type == 'rate' else uncal_dir
    prod_subgroup = 'RATE' if download_type == 'rate' else 'UNCAL'
    calib_lvl = 2 if download_type == 'rate' else 1

    print(f"\nQuerying MAST for PID {pid} ({target_name}) IFU {download_type}...")

    obs_table = Observations.query_criteria(
        instrument_name=['NIRSPEC/IFU'],
        provenance_name=['CALJWST'],
        obs_id=[f'*{pid}*'],
    )

    if len(obs_table) == 0:
        print(f"  No IFU observations found for PID {pid}.")
        return []

    print(f"  Found {len(obs_table)} observations")

    downloads = []
    for exposure in obs_table:
        products = Observations.get_product_list(exposure)
        # Filter for target file type and grating
        filtered = Observations.filter_products(
            products,
            filters=gratings_filter,
            calib_level=[calib_lvl],
            productSubGroupDescription=prod_subgroup,
        )
        if len(filtered) == 0:
            continue
        # Exclude tiny confirmation images by keeping only large files
        avg_size = np.nanmean(filtered['size']) if len(filtered) > 1 else 0
        if avg_size > 0:
            filtered = filtered[filtered['size'] > avg_size * 0.5]
        for uri in filtered['dataURI']:
            # MAST URIs use zero-padded 5-digit PID (e.g. jw01537...)
            pid_padded = pid.zfill(5)
            if f'jw{pid_padded}' in str(uri) or f'jw{pid}' in str(uri):
                downloads.append(uri)

    downloads = list(set(downloads))
    print(f"  Files to download: {len(downloads)}")

    downloaded = []
    for uri in sorted(downloads):
        fname = os.path.join(target_dir, os.path.basename(uri))
        if os.path.exists(fname):
            print(f"  Already downloaded: {os.path.basename(fname)}")
            downloaded.append(fname)
            continue
        print(f"  Downloading: {os.path.basename(uri)}")
        try:
            manifest = Observations.download_file(uri, local_path=target_dir)
            if manifest[0] == 'COMPLETE':
                downloaded.append(fname)
        except Exception as e:
            print(f"    ERROR: {e}")

    return sorted(downloaded)


# ---------------------------------------------------------------------------
# Per-PID pipeline runner
# ---------------------------------------------------------------------------

def process_pid(pid, run_download=True, stages=(1, 2, 3), skip_existing=True,
                gratings_filter=None):
    """
    Run the full IFU pipeline for a given PID.
    stages: tuple of stages to run (subset of 1, 2, 3).
    """
    target_name = PROGRAMS.get(pid, pid)
    base_dir = os.path.join(BASE_DATA_DIR, f'PID{pid}_{target_name}')

    uncal_dir  = os.path.join(base_dir, 'uncal')
    det1_dir   = os.path.join(base_dir, 'stage1')
    asn_dir    = os.path.join(base_dir, 'asn')
    spec2_dir  = os.path.join(base_dir, 'stage2')
    spec3_dir  = os.path.join(base_dir, 'stage3')

    for d in [uncal_dir, det1_dir, asn_dir, spec2_dir, spec3_dir]:
        os.makedirs(d, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"Processing PID {pid} ({target_name})")
    print(f"{'='*60}")

    # ----- Download -----
    if run_download:
        rate_files = download_ifu_data(
            pid, target_name, base_dir,
            gratings_filter=gratings_filter,
            download_type='rate',
        )
        if not rate_files:
            print(f"  No rate files downloaded for PID {pid}. "
                  "Checking for uncal files...")
            uncal_files = download_ifu_data(
                pid, target_name, base_dir,
                gratings_filter=gratings_filter,
                download_type='uncal',
            )
    else:
        # Collect whatever already exists
        rate_files = sorted(glob.glob(os.path.join(det1_dir, '*_rate.fits')))
        uncal_files = sorted(glob.glob(os.path.join(uncal_dir, '*_uncal.fits')))
        if not rate_files:
            print(f"  No rate.fits yet in {det1_dir}; found {len(uncal_files)} uncal")

    # ----- Stage 1 -----
    if 1 in stages:
        uncal_files = sorted(glob.glob(os.path.join(uncal_dir, '*_uncal.fits')))
        if uncal_files:
            print(f"\n[Stage 1] Running Detector1Pipeline for {len(uncal_files)} files")
            t0 = time.perf_counter()
            new_rate = run_stage1(uncal_dir, det1_dir,
                                  skip_existing=skip_existing)
            print(f"[Stage 1] Done – {len(new_rate)} rate files "
                  f"({(time.perf_counter()-t0)/60:.1f} min)")
        else:
            print("[Stage 1] No uncal files; skipping")

    # ----- Gather rate files for Stage 2 -----
    rate_files = sorted(glob.glob(os.path.join(det1_dir, '*_rate.fits')))
    if not rate_files:
        print(f"  No rate files in {det1_dir}. Cannot run Stages 2/3.")
        return

    # Filter to target gratings if specified
    if gratings_filter:
        rate_files = [f for f in rate_files if any(
            g.split(';')[-1] in (fits.getheader(f).get('GRATING', ''))
            for g in gratings_filter
        )]
    print(f"  Rate files for Stage 2: {len(rate_files)}")

    # ----- Stage 2 -----
    if 2 in stages:
        print(f"\n[Stage 2] Running Spec2Pipeline for {len(rate_files)} files")
        t0 = time.perf_counter()
        cal_files = run_stage2(rate_files, asn_dir, spec2_dir,
                               skip_existing=skip_existing)
        print(f"[Stage 2] Done – {len(cal_files)} cal files "
              f"({(time.perf_counter()-t0)/60:.1f} min)")
    else:
        cal_files = sorted(glob.glob(os.path.join(spec2_dir, '*_cal.fits')))
        print(f"[Stage 2] Skipped; found {len(cal_files)} existing cal files")

    # ----- Stage 3 -----
    if 3 in stages:
        print(f"\n[Stage 3] Running Spec3Pipeline for {len(cal_files)} cal files")
        t0 = time.perf_counter()
        s3d_files, x1d_files = run_stage3(cal_files, asn_dir, spec3_dir,
                                           skip_existing=skip_existing)
        print(f"[Stage 3] Done – {len(s3d_files)} S3D + {len(x1d_files)} X1D "
              f"({(time.perf_counter()-t0)/60:.1f} min)")
    else:
        s3d_files = sorted(glob.glob(os.path.join(spec3_dir, '*_s3d.fits')))
        x1d_files = sorted(glob.glob(os.path.join(spec3_dir, '*_x1d.fits')))
        print(f"[Stage 3] Skipped; found {len(s3d_files)} S3D, {len(x1d_files)} X1D")

    print(f"\nPID {pid} complete.")
    return {
        'pid': pid,
        'target': target_name,
        'base_dir': base_dir,
        's3d': s3d_files,
        'x1d': x1d_files,
    }


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description='Reduce NIRSpec IFU data (Parlanti et al. calibration targets)')
    parser.add_argument('--pid', default='all',
                        help='Program ID to process: 1537|1538|all (default: all)')
    parser.add_argument('--stage', default='all',
                        help='Stages to run: 1|2|3|all|2,3 (default: all)')
    parser.add_argument('--download', action='store_true',
                        help='Download data from MAST before processing')
    parser.add_argument('--skip-existing', dest='skip_existing',
                        action='store_true', default=True,
                        help='Skip files that have already been processed')
    parser.add_argument('--no-skip', dest='skip_existing',
                        action='store_false')
    parser.add_argument('--gratings', default='G140M,G235M',
                        help='Comma-separated gratings to process (default: G140M,G235M)')
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    args = parse_args()

    # Resolve which PIDs to run
    if args.pid == 'all':
        pids = list(PROGRAMS.keys())
    else:
        pids = [p.strip() for p in args.pid.split(',')]

    # Resolve which stages to run
    if args.stage == 'all':
        stages = (1, 2, 3)
    else:
        stages = tuple(int(s) for s in args.stage.split(','))

    # Resolve target gratings → MAST filter strings
    target_grat = [g.strip() for g in args.gratings.split(',')]
    grat_to_filter = {
        'G140M': 'F100LP;G140M',
        'G235M': 'F170LP;G235M',
        'G395M': 'F290LP;G395M',
    }
    gratings_filter = [grat_to_filter[g] for g in target_grat
                       if g in grat_to_filter]

    print(f"\nIFU Reduction Plan:")
    print(f"  PIDs:     {pids}")
    print(f"  Stages:   {stages}")
    print(f"  Gratings: {gratings_filter}")
    print(f"  Download: {args.download}")
    print(f"  Skip existing: {args.skip_existing}\n")

    t_total = time.perf_counter()

    results = []
    for pid in pids:
        r = process_pid(
            pid,
            run_download=args.download,
            stages=stages,
            skip_existing=args.skip_existing,
            gratings_filter=gratings_filter,
        )
        if r:
            results.append(r)

    elapsed = (time.perf_counter() - t_total) / 60
    print(f"\n{'='*60}")
    print(f"All done in {elapsed:.1f} min")
    print(f"{'='*60}")
    for r in results:
        print(f"  PID {r['pid']} ({r['target']}): "
              f"{len(r['s3d'])} S3D cubes, {len(r['x1d'])} X1D spectra")
