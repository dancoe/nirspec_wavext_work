"""
Extended-range IFU pipeline runner using Parlanti et al. (2025) reference files.

Applies modified S-flat, F-flat, wavelengthrange, and cubepar reference files
to recover the extended wavelength range of NIRSpec/IFU medium-grating data:
  G140M/F100LP:  0.97 – 3.57 µm  (nominal: 0.97 – 1.88 µm)
  G235M/F170LP:  1.66 – 5.25 µm  (nominal: 1.70 – 3.15 µm)

Prerequisites:
  1. Run modify_pipeline.py to comment out the NRS2 check in jwst/assign_wcs/nirspec.py
  2. Parlanti reference files must be present in data/parlanti_repo/CRDS_1364/
  3. Rate files must be downloaded (run download_ifu_rates.py first)

Usage:
    micromamba run -n jwst_1.20.2 python run_ifu_pipeline_ext.py --pid 1537
    micromamba run -n jwst_1.20.2 python run_ifu_pipeline_ext.py --pid all
    micromamba run -n jwst_1.20.2 python run_ifu_pipeline_ext.py --pid 1537,1538 --stage 2
    micromamba run -n jwst_1.20.2 python run_ifu_pipeline_ext.py --pid 1537 --stage 3

Reference:
    Parlanti et al. 2025, arXiv:2512.14844
    https://github.com/eleonoraparlanti/nirspecIFU-extended
"""
import os, sys, glob, json, time, argparse, warnings
import numpy as np
from astropy.io import fits

warnings.simplefilter("ignore", RuntimeWarning)

os.environ.setdefault('CRDS_PATH', os.path.expanduser('~/crds_cache'))
os.environ.setdefault('CRDS_SERVER_URL', 'https://jwst-crds.stsci.edu')

from jwst.pipeline import Spec2Pipeline, Spec3Pipeline
from jwst.associations import asn_from_list as afl
from jwst.associations.lib.rules_level2_base import DMSLevel2bBase
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base
import jwst

print(f"JWST pipeline: {jwst.__version__}")

BASE_DIR = '/Users/dcoe/NIRSpec/wavext/data/IFU'
PROGRAMS = {
    '1536': 'J1743045',
    '1537': 'G191-B2B',
    '1538': 'P330E',
    '6645': 'P330E-C3',
}
EXP_TYPE = 'NRS_IFU'
MAX_EXPOSURES_PER_GROUP = 4   # per (GRATING, FILTER, DETECTOR)

# Parlanti-modified reference file paths
PARLANTI_DIR = '/Users/dcoe/NIRSpec/wavext/data/parlanti_repo/CRDS_1364'
WAVELENGTHRANGE = os.path.join(PARLANTI_DIR, 'jwst_nirspec_wavelengthrange_0008.asdf')
CUBEPAR = os.path.join(PARLANTI_DIR, 'jwst_nirspec_cubepar_0009.fits')

# S-flat: keyed by (GRATING, FILTER, DETECTOR)
SFLAT = {
    ('G140M', 'F100LP', 'NRS1'): os.path.join(PARLANTI_DIR, 'jwst_nirspec_sflat_0208.fits'),
    ('G140M', 'F100LP', 'NRS2'): os.path.join(PARLANTI_DIR, 'jwst_nirspec_sflat_0191.fits'),
    ('G235M', 'F170LP', 'NRS1'): os.path.join(PARLANTI_DIR, 'jwst_nirspec_sflat_0211.fits'),
    ('G235M', 'F170LP', 'NRS2'): os.path.join(PARLANTI_DIR, 'jwst_nirspec_sflat_0192.fits'),
}

# F-flat: keyed by GRATING
FFLAT = {
    'G140M': os.path.join(PARLANTI_DIR, 'jwst_nirspec_fflat_0105.fits'),
    'G235M': os.path.join(PARLANTI_DIR, 'jwst_nirspec_fflat_0091.fits'),
}

# Gratings and filters to process (Parlanti extended configurations only)
TARGET_CONFIGS = {('G140M', 'F100LP'), ('G235M', 'F170LP')}


def get_header(f, key=None):
    h = fits.getheader(f)
    return h if key is None else h.get(key)


def check_reffiles():
    """Verify all required Parlanti reference files are present and non-empty."""
    missing = []
    for path in [WAVELENGTHRANGE, CUBEPAR] + list(SFLAT.values()) + list(FFLAT.values()):
        if not os.path.exists(path) or os.path.getsize(path) == 0:
            missing.append(os.path.basename(path))
    if missing:
        raise FileNotFoundError(
            f"Parlanti reference files missing or empty: {missing}\n"
            f"  Expected in: {PARLANTI_DIR}\n"
            f"  Run the download commands in PARLANTI_IFU.md first."
        )
    print(f"  All Parlanti reference files present in {PARLANTI_DIR}")


def writel2asn(onescifile, allscifiles, asnfile, product_name):
    """Create a Level 2 ASN for a single IFU science exposure."""
    hdr = get_header(onescifile)
    if hdr.get('EXP_TYPE') != EXP_TYPE:
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
    sci_match = [f for f in allscifiles if f != onescifile
                 and get_header(f, 'DETECTOR') == detector
                 and get_header(f, 'FILTER') == filt
                 and get_header(f, 'GRATING') == grating]

    is_nod = 'NOD' in patttype.split('-')
    if is_nod:
        for f in sci_match:
            if get_header(f, 'PATT_NUM') != pattnum:
                members.append({'expname': f, 'exptype': 'background'})
    for f in sci_match:
        members.append({'expname': f, 'exptype': 'selfcal'})

    _, serialized = asn.dump()
    with open(asnfile, 'w') as fh:
        fh.write(serialized)
    return True


def writel3asn(scifiles, asn_dir):
    """Create Level 3 ASNs grouped by FILTER/GRATING (combining NRS1+NRS2)."""
    from collections import defaultdict
    grouped = defaultdict(list)
    for f in scifiles:
        try:
            filt = get_header(f, 'FILTER')
            grat = get_header(f, 'GRATING')
            grouped[(filt, grat)].append(f)
        except Exception:
            pass

    asn_files = []
    for (filt, grat), files in sorted(grouped.items()):
        name = f"{filt}_{grat}".lower()
        asnf = os.path.join(asn_dir, f"{name}_ext_l3asn.json")
        asn = afl.asn_from_list(files, rule=DMS_Level3_Base, product_name=name)
        _, s = asn.dump()
        with open(asnf, 'w') as fh:
            fh.write(s)
        n_nrs1 = sum(1 for f in files if 'nrs1' in f)
        n_nrs2 = sum(1 for f in files if 'nrs2' in f)
        print(f"  L3 ASN: {os.path.basename(asnf)} "
              f"({len(files)} cal files: {n_nrs1} NRS1, {n_nrs2} NRS2)")
        asn_files.append(asnf)
    return asn_files


def run_stage2_ext(pid):
    """Run Stage 2 with Parlanti reference file overrides."""
    target = PROGRAMS[pid]
    base   = os.path.join(BASE_DIR, f'PID{pid}_{target}')
    det1   = os.path.join(base, 'stage1')
    asn_d  = os.path.join(base, 'asn_ext')
    spec2  = os.path.join(base, 'stage2_ext')
    for d in [asn_d, spec2]:
        os.makedirs(d, exist_ok=True)

    rate_files = sorted(glob.glob(os.path.join(det1, '*_rate.fits')))
    if not rate_files:
        print(f"[S2ext] No rate files for PID {pid}")
        return []

    # Group by (GRATING, FILTER, DETECTOR) — keep F100LP and F170LP separate
    from collections import defaultdict
    groups = defaultdict(list)
    for f in rate_files:
        try:
            h = get_header(f)
            if h.get('EXP_TYPE') != EXP_TYPE:
                continue
            key = (h['GRATING'], h['FILTER'], h['DETECTOR'])
            if (h['GRATING'], h['FILTER']) not in TARGET_CONFIGS:
                continue  # skip F070LP/G140M, G395M, etc.
            groups[key].append(f)
        except Exception as e:
            print(f"  Skipping {os.path.basename(f)}: {e}")

    all_cal = []
    for (grat, filt, det), files in sorted(groups.items()):
        files = files[:MAX_EXPOSURES_PER_GROUP]
        sflat_path = SFLAT.get((grat, filt, det))
        fflat_path = FFLAT.get(grat)
        if not sflat_path:
            print(f"\n  [{pid}] WARNING: No Parlanti sflat for {grat}/{filt}/{det}, skipping")
            continue
        print(f"\n  [{pid}] {grat}/{filt}/{det}: {len(files)} exposures (extended)")
        print(f"    sflat: {os.path.basename(sflat_path)}")
        print(f"    fflat: {os.path.basename(fflat_path)}")

        spec2_steps = {
            'assign_wcs': {
                'override_wavelengthrange': WAVELENGTHRANGE,
            },
            'flat_field': {
                'override_sflat': sflat_path,
                'override_fflat': fflat_path,
            },
            'bkg_subtract': {'skip': True},
            'cube_build':   {'skip': True},
            'extract_1d':   {'skip': True},
        }

        for rate_file in files:
            base_name = os.path.basename(rate_file).replace('_rate.fits', '')
            out_cal = os.path.join(spec2, f'{base_name}_cal.fits')
            if os.path.exists(out_cal):
                print(f"    EXISTS: {os.path.basename(out_cal)}")
                all_cal.append(out_cal)
                continue
            asnf = os.path.join(asn_d, f'{base_name}_l2asn.json')
            print(f"    Spec2: {base_name}")
            ok = writel2asn(rate_file, rate_files, asnf, base_name)
            if not ok:
                print(f"      Skipped (not IFU)")
                continue
            try:
                Spec2Pipeline.call(asnf, save_results=True,
                                   steps=spec2_steps, output_dir=spec2)
                if os.path.exists(out_cal):
                    all_cal.append(out_cal)
                    print(f"      → {os.path.basename(out_cal)}")
                else:
                    produced = sorted(glob.glob(os.path.join(spec2, f'*{base_name}*_cal.fits')))
                    all_cal.extend(produced)
                    if produced:
                        print(f"      → {[os.path.basename(p) for p in produced]}")
            except Exception as e:
                print(f"      ERROR: {e}")

    print(f"\n  [{pid}] Stage 2 ext complete: {len(all_cal)} cal files")
    return all_cal


def run_stage3_ext(pid):
    """Run Stage 3 with extended cubepar, combining NRS1+NRS2 cal files."""
    target = PROGRAMS[pid]
    base   = os.path.join(BASE_DIR, f'PID{pid}_{target}')
    asn_d  = os.path.join(base, 'asn_ext')
    spec2  = os.path.join(base, 'stage2_ext')
    spec3  = os.path.join(base, 'stage3_ext')
    os.makedirs(spec3, exist_ok=True)

    cal_files = sorted(glob.glob(os.path.join(spec2, '*_cal.fits')))
    # Filter to only the Parlanti target configs
    cal_sel = []
    for f in cal_files:
        try:
            grat = get_header(f, 'GRATING')
            filt = get_header(f, 'FILTER')
            if (grat, filt) in TARGET_CONFIGS:
                cal_sel.append(f)
        except Exception:
            pass
    if not cal_sel:
        print(f"[S3ext] No target cal files for PID {pid}")
        return [], []

    spec3_steps = {
        'master_background': {'skip': True},
        'outlier_detection': {'kernel_size': '3 3'},
        'cube_build': {
            'override_cubepar': CUBEPAR,
        },
        'extract_1d': {},
    }

    print(f"\n  [{pid}] Stage 3 ext: {len(cal_sel)} cal files")
    n_nrs1 = sum(1 for f in cal_sel if 'nrs1' in f)
    n_nrs2 = sum(1 for f in cal_sel if 'nrs2' in f)
    print(f"    ({n_nrs1} NRS1, {n_nrs2} NRS2)")
    print(f"    cubepar: {os.path.basename(CUBEPAR)}")
    asn_files = writel3asn(cal_sel, asn_d)

    for asnf in sorted(asn_files):
        name = os.path.basename(asnf).replace('_ext_l3asn.json', '')
        out_s3d = os.path.join(spec3, f'{name}_s3d.fits')
        if os.path.exists(out_s3d):
            print(f"    EXISTS: {os.path.basename(out_s3d)}")
            continue
        print(f"    Spec3: {os.path.basename(asnf)}")
        try:
            Spec3Pipeline.call(asnf, save_results=True,
                               steps=spec3_steps, output_dir=spec3)
        except Exception as e:
            print(f"      ERROR: {e}")

    s3d = sorted(glob.glob(os.path.join(spec3, '*_s3d.fits')))
    x1d = sorted(glob.glob(os.path.join(spec3, '*_x1d.fits')))
    print(f"  [{pid}] Stage 3 ext done: {len(s3d)} S3D cubes, {len(x1d)} X1D spectra")
    return s3d, x1d


def parse_args():
    p = argparse.ArgumentParser(description='Extended IFU pipeline with Parlanti reffiles')
    p.add_argument('--pid', default='all',
                   help='PID or comma-separated list, or "all"')
    p.add_argument('--stage', default='23',
                   help='Stage(s) to run: 2, 3, or 23 (default)')
    return p.parse_args()


if __name__ == '__main__':
    args = parse_args()
    pids = list(PROGRAMS.keys()) if args.pid == 'all' else args.pid.split(',')
    stages = [int(s) for s in args.stage]

    # Verify reference files before starting
    print("Checking Parlanti reference files ...")
    check_reffiles()
    print(f"Pipeline patch: NRS2 check commented out in jwst/assign_wcs/nirspec.py")
    print(f"Processing PIDs: {pids}, stages: {stages}")
    print()

    t0 = time.perf_counter()
    for pid in pids:
        if pid not in PROGRAMS:
            print(f"Unknown PID: {pid}, skipping"); continue
        target = PROGRAMS[pid]
        print(f"\n{'='*60}")
        print(f"PID {pid} — {target}")
        print(f"{'='*60}")

        if 2 in stages:
            cal = run_stage2_ext(pid)
        if 3 in stages:
            s3d, x1d = run_stage3_ext(pid)

    elapsed = time.perf_counter() - t0
    print(f"\nTotal time: {elapsed/60:.1f} min")
