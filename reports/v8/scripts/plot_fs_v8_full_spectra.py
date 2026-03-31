
import os
import glob
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

# Paths
FS_DIR = '/Users/dcoe/NIRSpec/wavext/data/FS'
OUTPUT_DIR = '/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/reports/v8/fs_v8'
os.makedirs(OUTPUT_DIR, exist_ok=True)

TARGETS = [
    {'pid': '1536', 'name': 'J1743045', 'pid_dir': f'{FS_DIR}/PID1536_J1743045'},
    {'pid': '1537', 'name': 'G191-B2B', 'pid_dir': f'{FS_DIR}/PID1537_G191-B2B'},
    {'pid': '1538', 'name': 'P330E', 'pid_dir': f'{FS_DIR}/PID1538_P330E'},
    {'pid': '6644', 'name': 'NGC2506-G31', 'pid_dir': f'{FS_DIR}/PID6644_NGC2506G31'},
    {'pid': '1492', 'name': 'IRAS-05248', 'pid_dir': f'{FS_DIR}/PID1492*'},
    {'pid': '1128', 'name': 'J1808347', 'pid_dir': f'{FS_DIR}/PID1128*'},
]

def plot_target(target):
    pid = target['pid']
    print(f"Plotting PID {pid}...")
    
    plt.figure(figsize=(12, 6))
    
    # Standard Gratings
    gratings = ['g140m', 'g235m', 'g395m', 'prism']
    colors = {'g140m': 'blue', 'g235m': 'green', 'g395m': 'red', 'prism': 'purple'}
    
    # Resolve the directory if it's a glob
    dirs = glob.glob(target['pid_dir'])
    if not dirs:
        print(f"Directory not found for {target['pid']}")
        plt.close()
        return
        
    pid_dir = dirs[0]
    
    found_any = False
    for grating in gratings:
        # Search recursively for x1d files within the PID folder
        files = glob.glob(f"{pid_dir}/**/*{grating}*x1d.fits", recursive=True)
        # Avoid duplicate gratings from different runs if possible?
        # Prefer 'ext' or 'v7' if present
        if not files:
            continue
            
        found_any = True
        
        # Simple preference: prefer 'ext' or 'v7' in the path
        best_file = files[0]
        for f in files:
            if 'ext' in f.lower() or 'v7' in f.lower():
                best_file = f
                break
        
        with fits.open(best_file) as h:
            wl = h[1].data['WAVELENGTH']
            spec = h[1].data['FLUX']
            
        # Gap aware plotting
        if grating == 'g140m':
            mask1 = wl < 2.17
            mask2 = wl > 2.28
            plt.plot(wl[mask1], spec[mask1], color=colors[grating], label=f'{grating.upper()} (v8)')
            plt.plot(wl[mask2], spec[mask2], color=colors[grating])
        elif grating == 'g235m':
            mask1 = wl < 3.66
            mask2 = wl > 3.82
            plt.plot(wl[mask1], spec[mask1], color=colors[grating], label=f'{grating.upper()} (v8)')
            plt.plot(wl[mask2], spec[mask2], color=colors[grating])
        else:
            plt.plot(wl, spec, color=colors[grating], label=f'{grating.upper()} (v8)')
            
    if not found_any:
        print(f"No files found for PID {pid} in {pid_dir}")
        plt.close()
        return

    plt.yscale('log')
    plt.xlabel("Wavelength (um)")
    plt.ylabel("Flux (Jy)")
    plt.title(f"{target['name']} – PID {pid} – FS wavext v8")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.xlim(0.6, 5.6)
    
    plt.savefig(f"{OUTPUT_DIR}/full_spectrum_v8_{pid}.png", dpi=150)
    plt.close()

if __name__ == "__main__":
    for target in TARGETS:
        try:
            plot_target(target)
        except Exception as e:
            print(f"Failed PID {target['pid']}: {e}")
    print("All FS plots generated.")
