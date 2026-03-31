"""
Run FS NRS1 Nominal Pipeline: Spec2 → Spec3 (Level 3).

This script produces "our" Level-3 Fixed Slit x1d spectra for the NRS1
detector using standard pipeline settings (no wavelength extension overrides).
This is for comparison with MAST Level-3 products.

Workflow per PID / grating:
  1. Spec2Pipeline on each NRS1 rate file.
     Output: *_nrs1_cal.fits in nrs1_spec2_cal/
  2. Spec3Pipeline on the resulting cal files.
     Output: Level-3 *_nrs1_x1d.fits in nrs1_spec3_ext/
"""
import os
import glob
import warnings
import argparse

from astropy.io import fits
from jwst.pipeline import Spec2Pipeline, Spec3Pipeline
from jwst.associations import asn_from_list as afl
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base
import jwst

warnings.simplefilter('ignore', RuntimeWarning)

os.environ.setdefault('CRDS_PATH', os.path.expanduser('~/crds_cache'))
os.environ.setdefault('CRDS_SERVER_URL', 'https://jwst-crds.stsci.edu')

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE = '/Users/dcoe/NIRSpec/wavext/data'

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
    asn = afl.asn_from_list(cal_files, rule=DMS_Level3_Base,
                             product_name=product_name)
    _, serialized = asn.dump()
    with open(asn_path, 'w') as fh:
        fh.write(serialized)
    return asn_path

# ── Stage-2 ───────────────────────────────────────────────────────────────────

def run_spec2(rate_file, out_dir):
    base    = os.path.basename(rate_file).replace('_rate.fits', '')
    out_cal = os.path.join(out_dir, f'{base}_cal.fits')

    if os.path.exists(out_cal):
        return out_cal

    print(f'    Spec2 (NRS1): {base} …')
    try:
        Spec2Pipeline.call(
            rate_file,
            save_results=True,
            output_dir=out_dir,
            steps={'wavecorr': {'skip': False}} # Standard
        )
    except Exception as exc:
        print(f'      Spec2 ERROR: {exc}')
        return None

    produced = sorted(glob.glob(os.path.join(out_dir, f'{base}*_cal.fits')))
    return produced[-1] if produced else None

# ── Stage-3 ───────────────────────────────────────────────────────────────────

def run_spec3(cal_files, out_dir, product_name):
    asn_path  = os.path.join(out_dir, f'{product_name}_l3asn.json')
    out_x1d   = os.path.join(out_dir, f'{product_name}_x1d.fits')

    if os.path.exists(out_x1d):
        return out_x1d

    build_spec3_asn(cal_files, product_name, asn_path)
    print(f'    Spec3 (NRS1): {product_name} ({len(cal_files)} cals)')
    
    try:
        Spec3Pipeline.call(
            asn_path,
            save_results=True,
            output_dir=out_dir,
        )
    except Exception as exc:
        print(f'      Spec3 ERROR: {exc}')
        return None

    produced = sorted(glob.glob(os.path.join(out_dir, f'{product_name}*_x1d.fits')))
    return produced[-1] if produced else None

# ── Driver ─────────────────────────────────────────────────────────────────────

def process_pid(pid):
    info     = PIDS[pid]
    data_dir = info['dir']
    s2_dir   = os.path.join(data_dir, 'nrs1_spec2_cal')
    s3_dir   = os.path.join(data_dir, 'nrs1_spec3_nom')
    os.makedirs(s2_dir, exist_ok=True)
    os.makedirs(s3_dir, exist_ok=True)

    rate_files = sorted(glob.glob(os.path.join(data_dir, '*_nrs1_rate.fits')))
    
    from collections import defaultdict
    groups = defaultdict(list)
    for rf in rate_files:
        grat, filt, exp_type = get_hdr(rf, 'GRATING', 'FILTER', 'EXP_TYPE')
        if exp_type == 'NRS_FIXEDSLIT' and (grat, filt) in TARGET_CONFIGS:
            groups[(grat, filt)].append(rf)

    for (grat, filt), rate_list in sorted(groups.items()):
        print(f'\n  [{pid}] NRS1 {grat}/{filt}: {len(rate_list)} files')
        cal_files = []
        for rf in rate_list:
            cal = run_spec2(rf, s2_dir)
            if cal: cal_files.append(cal)

        if cal_files:
            product_name = f'nrs1_l3_{pid}_{grat.lower()}'
            run_spec3(cal_files, s3_dir, product_name)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pid', default='all', help='PID to process (all or comma-separated list)')
    args = parser.parse_args()

    pids = list(PIDS.keys()) if args.pid == 'all' else [p.strip() for p in args.pid.split(',')]
    
    for pid in pids:
        if pid not in PIDS:
            print(f'Unknown PID: {pid}'); continue
        print(f'\n=== PID {pid} ===')
        process_pid(pid)

    print('\nAll done.')
