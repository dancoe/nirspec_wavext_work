
import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import ImageNormalize, AsinhStretch

def get_config(f):
    try:
        with fits.open(f) as h:
            hdr = h[0].header
            return (hdr.get('GRATING'), hdr.get('FILTER'), hdr.get('DETECTOR'), hdr.get('EXP_TYPE'))
    except:
        return None

def main():
    rate_dir = '/Users/dcoe/NIRSpec/wavext/data/FS/PID1492/MAST/rate'
    plot_dir = '/Users/dcoe/NIRSpec/wavext/reports/1492/plots'
    os.makedirs(plot_dir, exist_ok=True)
    
    rate_files = sorted(glob.glob(os.path.join(rate_dir, '*_rate.fits')))
    
    # Representative file for each combination
    unique_configs = {}
    for f in rate_files:
        cfg = get_config(f)
        if cfg:
            # We want both NRS1 and NRS2 for each Grating/Filter if possible
            key = (cfg[0], cfg[1]) # Grating, Filter
            if key not in unique_configs:
                unique_configs[key] = {}
            if cfg[2] not in unique_configs[key]: # Detector
                unique_configs[key][cfg[2]] = f
    
    # Flatten to a list of files to plot
    files_to_plot = []
    for key in sorted(unique_configs.keys()):
        for det in sorted(unique_configs[key].keys()):
            files_to_plot.append(unique_configs[key][det])
            
    num_files = len(files_to_plot)
    if num_files == 0:
        print("No rate files found.")
        return
        
    ncols = 2
    nrows = (num_files + 1) // 2
    fig, axes = plt.subplots(nrows, ncols, figsize=(20, 4 * nrows), squeeze=False)
    axes = axes.flatten()
    
    for i, f in enumerate(files_to_plot):
        ax = axes[i]
        with fits.open(f) as h:
            data = h['SCI'].data
            hdr = h[0].header
            
        # Detect subarray positioning
        xstart = hdr.get('SUBSTRT1', 1)
        ystart = hdr.get('SUBSTRT2', 1)
        xsize = hdr.get('SUBSIZE1', data.shape[1])
        ysize = hdr.get('SUBSIZE2', data.shape[0])
        
        x_base, x_top = xstart - 1, xstart - 1 + xsize
        y_base, y_top = ystart - 1, ystart - 1 + ysize
        
        norm = ImageNormalize(vmin=-0.001, vmax=0.01, stretch=AsinhStretch())
        im = ax.imshow(data, origin='lower', cmap='viridis', norm=norm, 
                       interpolation='nearest', aspect='auto',
                       extent=[x_base, x_top, y_base, y_top])
        
        # Style as in mos.py
        ax.set_ylabel('Pixel Row')
        ax.set_xlabel('Pixel Column')
        ax.set_title(f"{hdr.get('GRATING')}/{hdr.get('FILTER')} - {hdr.get('DETECTOR')}\n{os.path.basename(f)}", fontsize=10)
        
        # Pixel ticks styling as requested
        ax.set_ylim(y_base - 10, y_top + 10)
        ax.set_xlim(0, 2048)
        # Enable grid as often seen in MOS plots for alignment check
        ax.grid(True, color='w', alpha=0.1, ls='-')
        
    # Hide unused axes
    for j in range(i + 1, len(axes)):
        axes[j].axis('off')
        
    plt.tight_layout()
    output_path = os.path.join(plot_dir, 'all_1492_rate_configs.png')
    plt.savefig(output_path, dpi=200, bbox_inches='tight')
    print(f"Saved plot to: {output_path}")

if __name__ == '__main__':
    main()
