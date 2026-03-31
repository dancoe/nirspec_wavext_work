#!/usr/bin/env python3
"""
Watch for IFU rate file downloads and launch Stage 2+3 pipeline
automatically when G140M and G235M rate files are ready for each PID.

Usage:
    micromamba run -n jwst_1.20.2 python watch_and_reduce_ifu.py
"""
import os
import glob
import time
import subprocess
import sys

STAGE1_DIRS = {
    '1537': '/Users/dcoe/NIRSpec/wavext/data/IFU/PID1537_G191-B2B/stage1',
    '1538': '/Users/dcoe/NIRSpec/wavext/data/IFU/PID1538_P330E/stage1',
}
PIPELINE_SCRIPT = '/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/analysis/reduction/run_ifu_pipeline_ext.py'
LOG_DIR = '/tmp'
SIZE_THRESH = 50 * 1e6  # 50MB


def count_ready_configs(stage1_dir):
    """Count how many G140M and G235M NRS1+NRS2 rate file pairs are ready."""
    g140m_nrs1 = sorted(glob.glob(os.path.join(stage1_dir, '*04101*nrs1_rate.fits')))
    g140m_nrs2 = sorted(glob.glob(os.path.join(stage1_dir, '*04101*nrs2_rate.fits')))
    g235m_nrs1 = sorted(glob.glob(os.path.join(stage1_dir, '*06101*nrs1_rate.fits')))
    g235m_nrs2 = sorted(glob.glob(os.path.join(stage1_dir, '*06101*nrs2_rate.fits')))
    # Check all are fully downloaded (>50MB)
    def all_ready(files):
        return all(os.path.exists(f) and os.path.getsize(f) > SIZE_THRESH for f in files)
    g140m_ok = len(g140m_nrs1) >= 4 and len(g140m_nrs2) >= 4 and all_ready(g140m_nrs1 + g140m_nrs2)
    g235m_ok = len(g235m_nrs1) >= 4 and len(g235m_nrs2) >= 4 and all_ready(g235m_nrs1 + g235m_nrs2)
    return g140m_ok, g235m_ok


started = set()

print(f"Watching for rate file downloads for PIDs: {list(STAGE1_DIRS.keys())}")
print(f"Will launch pipeline when G140M + G235M files are all ready.")
print()

while len(started) < len(STAGE1_DIRS):
    for pid, stage1_dir in STAGE1_DIRS.items():
        if pid in started:
            continue
        os.makedirs(stage1_dir, exist_ok=True)
        g140m_ok, g235m_ok = count_ready_configs(stage1_dir)
        n = len(glob.glob(os.path.join(stage1_dir, '*_rate.fits')))
        print(f"  PID {pid}: {n} rate files | G140M {'✓' if g140m_ok else '...'} | G235M {'✓' if g235m_ok else '...'}", flush=True)

        if g140m_ok and g235m_ok:
            print(f"\n  *** PID {pid}: G140M + G235M ready — launching pipeline! ***", flush=True)
            log = os.path.join(LOG_DIR, f'ifu_{pid}_reduce.log')
            cmd = [
                'micromamba', 'run', '-n', 'jwst_1.20.2',
                'python', PIPELINE_SCRIPT,
                '--pid', pid, '--stage', '23',
            ]
            with open(log, 'w') as fh:
                proc = subprocess.Popen(cmd, stdout=fh, stderr=subprocess.STDOUT,
                                        cwd='/Users/dcoe/NIRSpec/wavext')
            print(f"  Pipeline PID={proc.pid}, log: {log}")
            started.add(pid)

    if len(started) < len(STAGE1_DIRS):
        time.sleep(30)

print("\nAll PIDs launched. Monitor logs in /tmp/ifu_*_reduce.log")
