import os, glob
import subprocess

def run_plot(input_file):
    cmd = [
        "micromamba", "run", "-n", "jwst_1.20.2",
        "python", "/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/analysis/plot_extract_1d.py",
        input_file
    ]
    env = os.environ.copy()
    env["CRDS_SERVER_URL"] = "https://jwst-crds.stsci.edu"
    env["CRDS_PATH"] = os.path.expanduser("~/crds_cache")
    env["PYTHONPATH"] = "/Users/dcoe/NIRSpec/wavext/jwst_nirspec_wavext"
    
    print(f"Plotting {os.path.basename(input_file)}...")
    subprocess.run(cmd, env=env)

# Find all new extract_1d.fits
files = glob.glob("/Users/dcoe/NIRSpec/wavext/data/PID153[78]_*/*nrs2_extract_1d.fits")

for f in sorted(files):
    run_plot(f)
