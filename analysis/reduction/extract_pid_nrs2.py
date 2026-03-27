import os, glob, subprocess

ANALYSIS_ROOT = '/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/analysis'
DATA_DIR_TEMPLATE = '/Users/dcoe/NIRSpec/wavext/data/PID{pid}_*'

def run_step(script, input_file, output_file):
    cmd = ['python3', os.path.join(ANALYSIS_ROOT, script), input_file, output_file]
    print(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

def process_pid(pid):
    dirs = glob.glob(DATA_DIR_TEMPLATE.format(pid=pid))
    if not dirs:
        print(f"No directory found for PID {pid}")
        return
    data_dir = dirs[0]
    # We want NRS2 rate.fits files
    rates = glob.glob(os.path.join(data_dir, '*_nrs2_rate.fits'))
    for rate in rates:
        wcs = rate.replace('_rate.fits', '_nrs2_assign_wcs.fits')
        e2d = rate.replace('_rate.fits', '_nrs2_extract_2d.fits')
        e1d = rate.replace('_rate.fits', '_nrs2_extract_1d.fits')
        
        # Determine grating from header or file/dir context if needed
        # Actually, our extraction scripts (e.g. run_extract_1d.py) should handle the NRS2 correctly
        # if the WCS is correct.
        
        if not os.path.exists(wcs):
             run_step('reduction/run_assign_wcs.py', rate, wcs)
        if not os.path.exists(e2d):
             run_step('reduction/run_extract_2d.py', wcs, e2d)
        if not os.path.exists(e1d):
             run_step('reduction/run_extract_1d.py', e2d, e1d)

for pid in ['1536', '1537', '1538']:
    process_pid(pid)
