import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import ZScaleInterval

def plot_rate_pair(nrs1_rate, nrs2_rate, output_path):
    fig, axes = plt.subplots(1, 2, figsize=(15, 8))
    zscale = ZScaleInterval()
    
    for ax, path, label in zip(axes, [nrs1_rate, nrs2_rate], ['NRS1', 'NRS2']):
        with fits.open(path) as hd:
            data = hd['SCI'].data
            if len(data.shape) == 3: data = data[0] # Take first int
            
            vmin, vmax = zscale.get_limits(data)
            im = ax.imshow(data, origin='lower', cmap='gray', vmin=vmin, vmax=vmax)
            ax.set_title(f'1492 Rate File: {label}\n{os.path.basename(path)}')
            ax.set_xticks([]); ax.set_yticks([])
            
    plt.suptitle("NIRSpec PID 1492 FS Rate Files (NRS1 vs NRS2)", fontsize=14)
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    print(f"Saved: {output_path}")

if __name__ == "__main__":
    base = "data/FS/PID1492/MAST/rate/"
    out = "reports/1492/plots/"
    os.makedirs(out, exist_ok=True)
    
    f1 = os.path.join(base, "jw01492001001_03102_00005_nrs1_rate.fits")
    f2 = os.path.join(base, "jw01492003001_03102_00005_nrs2_rate.fits")
    
    plot_rate_pair(f1, f2, os.path.join(out, "rate_1492_G140M_nrs1_nrs2.png"))
