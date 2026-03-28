"""
Run Stages 2 and 3 of the JWST IFU pipeline for downloaded rate files.
Uses only one set of exposures per grating (3 exposures for a 3-point dither)
to keep processing time manageable.

Usage:
    micromamba activate jwst_1.20.2
    export CRDS_PATH=~/crds_cache
    export CRDS_SERVER_URL=https://jwst-crds.stsci.edu
    python run_ifu_pipeline.py --pid 1538|1537|all [--stage 2|3|23]
"""
import os, sys, glob, copy, json, time, argparse, warnings
import numpy as np
from astropy.io import fits

warnings.simplefilter("ignore", RuntimeWarning)

os.environ.setdefault('CRDS_PATH', os.path.expanduser('~/crds_cache'))
os.environ.setdefault('CRDS_SERVER_URL', 'https://jwst-crds.stsci.edu')

from jwst.pipeline import Spec2Pipeline, Spec3Pipeline
from jwst import datamodels
from jwst.associations import asn_from_list as afl
from jwst.associations.lib.rules_level2_base import DMSLevel2bBase
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base
import jwst

print(f"JWST pipeline: {jwst.__version__}")

BASE_DIR = '/Users/dcoe/NIRSpec/wavext/data/IFU'
PROGRAMS = {
    '1536': 'J1743045',    # 3rd calibration source — G140M + G235M + G395M
    '1537': 'G191-B2B',    # white dwarf
    '1538': 'P330E',       # gray star Cycle 1
    '6645': 'P330E-C3',    # gray star Cycle 3 (4th source)
}
EXP_TYPE = 'NRS_IFU'

# Limit how many NRS1+NRS2 pairs to use per grating/PID (saves disk/time)
# For a nominal dither pattern PID 1537 has 4 dithers; use all to get more
# coverage, but limit NRS2 pairs since Spec2 errors when all slices miss NRS2.
MAX_EXPOSURES_PER_GRATING = 4   # per detector; use all available (≤4)


def get_header(f, key=None):
    h = fits.getheader(f)
    return h if key is None else h.get(key)


def writel2asn(onescifile, allscifiles, asnfile, product_name):
    """Create a Level 2 ASN for a single IFU science exposure."""
    hdr = get_header(onescifile)
    if hdr.get('EXP_TYPE') != EXP_TYPE or hdr.get('IS_IMPRT', False):
        return False

    detector = hdr['DETECTOR']
    filt     = hdr['FILTER']
    grating  = hdr['GRATING']
    patttype = hdr.get('PATTTYPE', 'NONE')
    pattnum  = hdr.get('PATT_NUM', 1)

    asn = afl.asn_from_list([onescifile], rule=DMSLevel2bBase,
                             product_name=product_name)
    asn.data['program'] = hdr.get('PROGRAM', '')

    members = asn['products'][0]['members']

    # For dithered IFU, add selfcal members and nodded backgrounds
    sci_match = [f for f in allscifiles if f != onescifile
                 and get_header(f, 'DETECTOR') == detector
                 and get_header(f, 'FILTER') == filt
                 and get_header(f, 'GRATING') == grating
                 and not get_header(f, 'IS_IMPRT')]

    is_nod = 'NOD' in patttype.split('-')
    if is_nod:
        for f in sci_match:
            if get_header(f, 'PATT_NUM') != pattnum:
                members.append({'expname': f, 'exptype': 'background'})
    # Selfcal: all exposures of same config
    for f in sci_match:
        members.append({'expname': f, 'exptype': 'selfcal'})

    _, serialized = asn.dump()
    with open(asnfile, 'w') as fh:
        fh.write(serialized)
    return True


def writel3asn(scifiles, asn_dir):
    """Create Level 3 ASNs grouped by FILTER/GRATING."""
    from collections import defaultdict
    grouped = defaultdict(list)
    for f in scifiles:
        try:
            k = (get_header(f, 'FILTER'), get_header(f, 'GRATING'))
            grouped[k].append(f)
        except Exception:
            pass

    asn_files = []
    for (filt, grat), files in grouped.items():
        name = f"{filt}_{grat}".lower()
        asnf = os.path.join(asn_dir, f"{name}_l3asn.json")
        asn = afl.asn_from_list(files, rule=DMS_Level3_Base, product_name=name)
        _, s = asn.dump()
        with open(asnf, 'w') as fh:
            fh.write(s)
        print(f"  L3 ASN: {os.path.basename(asnf)} ({len(files)} cal files)")
        asn_files.append(asnf)
    return asn_files


def run_stage2(pid):
    target = PROGRAMS[pid]
    base  = os.path.join(BASE_DIR, f'PID{pid}_{target}')
    det1  = os.path.join(base, 'stage1')
    asn   = os.path.join(base, 'asn')
    spec2 = os.path.join(base, 'stage2')
    for d in [asn, spec2]:
        os.makedirs(d, exist_ok=True)

    rate_files = sorted(glob.glob(os.path.join(det1, '*_rate.fits')))
    if not rate_files:
        print(f"[S2] No rate files for PID {pid}")
        return []

    # Group by (GRATING, DETECTOR) and limit to MAX_EXPOSURES_PER_GRATING
    from collections import defaultdict
    groups = defaultdict(list)
    for f in rate_files:
        try:
            h = get_header(f)
            if h.get('EXP_TYPE') != EXP_TYPE:
                continue
            k = (h['GRATING'], h['DETECTOR'])
            groups[k].append(f)
        except Exception as e:
            print(f"  Skipping {os.path.basename(f)}: {e}")

    # Spec2 pipeline config
    spec2dict = {}
    for step in ['assign_wcs', 'badpix_selfcal', 'msa_flagging', 'nsclean',
                 'imprint_subtract', 'srctype', 'wavecorr',
                 'flat_field', 'pathloss', 'photom', 'pixel_replace']:
        spec2dict[step] = {}
    spec2dict['bkg_subtract'] = {'skip': True}
    spec2dict['cube_build']   = {'skip': True}
    spec2dict['extract_1d']   = {'skip': True}

    all_cal = []
    for (grat, det), files in sorted(groups.items()):
        files = files[:MAX_EXPOSURES_PER_GRATING]
        print(f"\n  [{pid}] {grat}/{det}: {len(files)} exposures")
        for rate_file in files:
            base_name = os.path.basename(rate_file).replace('_rate.fits', '')
            out_cal = os.path.join(spec2, f'{base_name}_cal.fits')
            if os.path.exists(out_cal):
                print(f"    EXISTS: {os.path.basename(out_cal)}")
                all_cal.append(out_cal)
                continue
            asnf = os.path.join(asn, f'{base_name}_l2asn.json')
            print(f"    Spec2: {base_name}")
            ok = writel2asn(rate_file, rate_files, asnf, base_name)
            if not ok:
                print(f"      Skipped (not IFU)")
                continue
            try:
                Spec2Pipeline.call(asnf, save_results=True,
                                   steps=spec2dict, output_dir=spec2)
                if os.path.exists(out_cal):
                    all_cal.append(out_cal)
                    print(f"      → {os.path.basename(out_cal)}")
                else:
                    # cal file may have different name
                    produced = sorted(glob.glob(os.path.join(spec2, f'*{base_name}*_cal.fits')))
                    all_cal.extend(produced)
            except Exception as e:
                print(f"      ERROR: {e}")

    print(f"\n  [{pid}] Stage 2 complete: {len(all_cal)} cal files")
    return all_cal


def run_stage3(pid):
    target = PROGRAMS[pid]
    base  = os.path.join(BASE_DIR, f'PID{pid}_{target}')
    asn   = os.path.join(base, 'asn')
    spec2 = os.path.join(base, 'stage2')
    spec3 = os.path.join(base, 'stage3')
    os.makedirs(spec3, exist_ok=True)

    cal_files = sorted(glob.glob(os.path.join(spec2, '*_cal.fits')))
    # Filter to target gratings G140M and G235M
    cal_sel = []
    for f in cal_files:
        try:
            grat = get_header(f, 'GRATING')
            if grat in ('G140M', 'G235M'):
                cal_sel.append(f)
        except Exception:
            pass
    if not cal_sel:
        print(f"[S3] No G140M/G235M cal files for PID {pid}")
        # Fall back to all cal files
        cal_sel = cal_files

    spec3dict = {}
    for step in ['assign_mtwcs', 'outlier_detection', 'pixel_replace',
                 'cube_build', 'extract_1d']:
        spec3dict[step] = {}
    spec3dict['master_background'] = {'skip': True}
    spec3dict['outlier_detection']['kernel_size'] = '3 3'

    print(f"\n  [{pid}] Stage 3: {len(cal_sel)} cal files")
    asn_files = writel3asn(cal_sel, asn)

    for asnf in sorted(asn_files):
        name = os.path.basename(asnf).replace('_l3asn.json', '')
        out_s3d = os.path.join(spec3, f'{name}_s3d.fits')
        if os.path.exists(out_s3d):
            print(f"    EXISTS: {os.path.basename(out_s3d)}")
            continue
        print(f"    Spec3: {os.path.basename(asnf)}")
        try:
            Spec3Pipeline.call(asnf, save_results=True,
                               steps=spec3dict, output_dir=spec3)
        except Exception as e:
            print(f"      ERROR: {e}")

    s3d = sorted(glob.glob(os.path.join(spec3, '*_s3d.fits')))
    x1d = sorted(glob.glob(os.path.join(spec3, '*_x1d.fits')))
    print(f"  [{pid}] Stage 3 done: {len(s3d)} S3D cubes, {len(x1d)} X1D spectra")
    return s3d, x1d


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--pid', default='all')
    p.add_argument('--stage', default='23')
    return p.parse_args()


if __name__ == '__main__':
    args = parse_args()
    pids = list(PROGRAMS.keys()) if args.pid == 'all' else args.pid.split(',')
    stages = [int(s) for s in args.stage]

    t0 = time.perf_counter()
    for pid in pids:
        print(f"\n{'='*60}\nPID {pid} ({PROGRAMS[pid]})\n{'='*60}")
        if 2 in stages:
            cal = run_stage2(pid)
        if 3 in stages:
            s3d, x1d = run_stage3(pid)

    print(f"\nTotal time: {(time.perf_counter()-t0)/60:.1f} min")
