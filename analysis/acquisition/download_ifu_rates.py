"""
Download IFU rate files from MAST for PIDs 1537 (G191-B2B) and 1538 (P330E).
Gets G140M/F100LP and G235M/F170LP IFU observations.
Only downloads main science exposures (~84MB), skips tiny TA/confirmation images.
"""
import os
import sys

os.environ['CRDS_PATH'] = os.path.expanduser('~/crds_cache')
os.environ['CRDS_SERVER_URL'] = 'https://jwst-crds.stsci.edu'

from astroquery.mast import Observations
import numpy as np

BASE_DIR = '/Users/dcoe/NIRSpec/wavext/data/IFU'
TARGET_FILTERS = ['F100LP;G140M', 'F170LP;G235M', 'F290LP;G395M']
SIZE_THRESH_MB = 50  # Only download files larger than this (skips TA images)

def download_ifu_rates(pid, target_name):
    """Download IFU _rate.fits files for both NRS1 and NRS2."""
    out_dir = os.path.join(BASE_DIR, f'PID{pid}_{target_name}', 'stage1')
    os.makedirs(out_dir, exist_ok=True)

    obs_table = Observations.query_criteria(
        instrument_name=['NIRSPEC/IFU'],
        provenance_name=['CALJWST'],
        obs_id=[f'*{pid}*'],
        filters=TARGET_FILTERS,
    )
    print(f'PID {pid} ({target_name}): {len(obs_table)} IFU observations')

    all_uris = []
    for row in obs_table:
        products = Observations.get_product_list(row)
        rate_files = Observations.filter_products(
            products,
            calib_level=[2],
            productSubGroupDescription='RATE',
        )
        # Skip tiny TA/confirmation images
        rate_files = rate_files[rate_files['size'] > SIZE_THRESH_MB * 1e6]
        # All remaining URIs are already from this PID (obs_table filtered by PID)
        for uri in rate_files['dataURI']:
            all_uris.append(str(uri))

    all_uris = sorted(set(all_uris))
    print(f'  Large rate files found: {len(all_uris)}')

    downloaded = []
    for uri in all_uris:
        fname = os.path.join(out_dir, os.path.basename(uri))
        # Skip only if file is already fully downloaded (> threshold)
        if os.path.exists(fname) and os.path.getsize(fname) > SIZE_THRESH_MB * 1e6:
            sz = os.path.getsize(fname)/1e6
            print(f'  EXISTS: {os.path.basename(fname)} ({sz:.1f}MB)')
            downloaded.append(fname)
            continue
        # Remove any stub/partial file before re-downloading
        if os.path.exists(fname):
            os.remove(fname)
        print(f'  Downloading: {os.path.basename(uri)}')
        for attempt in range(3):
            try:
                manifest = Observations.download_file(uri, local_path=fname)
                status = manifest[0] if hasattr(manifest, '__getitem__') else str(manifest)
                if os.path.exists(fname) and os.path.getsize(fname) > SIZE_THRESH_MB * 1e6:
                    sz = os.path.getsize(fname)/1e6
                    print(f'    OK ({sz:.1f}MB)')
                    downloaded.append(fname)
                    break
                else:
                    sz = os.path.getsize(fname)/1e6 if os.path.exists(fname) else 0
                    print(f'    Attempt {attempt+1}: incomplete ({sz:.1f}MB), status={status}')
                    if os.path.exists(fname):
                        os.remove(fname)
            except Exception as e:
                print(f'    Attempt {attempt+1} ERROR: {e}')
                if os.path.exists(fname):
                    os.remove(fname)
        else:
            print(f'    FAILED after 3 attempts: {os.path.basename(uri)}')

    print(f'  Downloaded/cached: {len(downloaded)} files')
    return downloaded

# All calibration standard star programs
PROGRAMS = [
    ('1536', 'J1743045'),   # 3rd source — all gratings including G395M
    ('1537', 'G191-B2B'),   # white dwarf
    ('1538', 'P330E'),      # gray star Cycle 1
    ('6645', 'P330E-C3'),   # gray star Cycle 3 (4th source)
]

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--pid', default='all',
                        help='Comma-separated PIDs or "all"')
    args = parser.parse_args()

    pids_to_run = [p for p, _ in PROGRAMS]
    if args.pid != 'all':
        pids_to_run = args.pid.split(',')

    all_files = {}
    for pid, target in PROGRAMS:
        if pid not in pids_to_run:
            continue
        files = download_ifu_rates(pid, target)
        all_files[pid] = files

    print('\n=== Summary ===')
    for pid, files in all_files.items():
        nrs1 = [f for f in files if 'nrs1' in f]
        nrs2 = [f for f in files if 'nrs2' in f]
        print(f'PID {pid}: {len(nrs1)} NRS1, {len(nrs2)} NRS2 rate files')
