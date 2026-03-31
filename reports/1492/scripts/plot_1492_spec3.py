
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table

def main():
    s3_dir = '/Users/dcoe/NIRSpec/wavext/data/FS/PID1492/spec3_full_cal'
    plot_dir = '/Users/dcoe/NIRSpec/wavext/reports/1492/plots'
    os.makedirs(plot_dir, exist_ok=True)
    
    gratings = [
        ('g140m', 'f100lp'),
        ('g235m', 'f170lp'),
        ('g395m', 'f290lp')
    ]
    
    fig, axes = plt.subplots(3, 1, figsize=(15, 12))
    
    for i, (grat, filt) in enumerate(gratings):
        ax = axes[i]
        
        nrs1_files = sorted(glob.glob(os.path.join(s3_dir, f'*_{grat}_{filt}_*nrs1*_x1d.fits')))
        nrs2_files = sorted(glob.glob(os.path.join(s3_dir, f'*_{grat}_{filt}_*nrs2*_x1d.fits')))
        
        if nrs1_files:
            with fits.open(nrs1_files[0]) as hdul:
                tab = Table(hdul['EXTRACT1D'].data)
                ax.step(tab['WAVELENGTH'], tab['FLUX'], label='NRS1', where='mid', color='C0', alpha=0.7)
                det1 = hdul[0].header.get('DETECTOR')
        
        if nrs2_files:
            with fits.open(nrs2_files[0]) as hdul:
                tab = Table(hdul['EXTRACT1D'].data)
                ax.step(tab['WAVELENGTH'], tab['FLUX'], label='NRS2', where='mid', color='C3', alpha=0.7)
                det2 = hdul[0].header.get('DETECTOR')
                    
        ax.set_title(f"PID 1492: {grat.upper()} / {filt.upper()}")
        ax.set_ylabel("Flux (Jy)")
        ax.legend()
        ax.grid(True, alpha=0.3)
        
    axes[2].set_xlabel("Wavelength ($\mu$m)")
    plt.tight_layout()
    
    output_path = os.path.join(plot_dir, 'spec3_1492_nrs1_nrs2_sidebyside.png')
    plt.savefig(output_path, dpi=200)
    print(f"Saved spectral plot to: {output_path}")

if __name__ == '__main__':
    main()
