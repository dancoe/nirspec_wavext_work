"""
NIRSpec Wavelength Extension — Create Custom Photom Reference File
Version: FS v7

This script creates a new photom reference file with LARGER column buffers
for wavelength etc. (500 elements instead of 150) to ensure we can extend relresponses
all the way to the NRS2 edge.
"""
import os
import numpy as np
from astropy.io import fits

BASE = '/Users/dcoe/NIRSpec/wavext'
SRC  = '/Users/dcoe/crds_cache/references/jwst/nirspec/jwst_nirspec_photom_0014.fits'
OUT  = os.path.join(BASE, 'nirspec_wavext_work/reference_files/jwst_nirspec_photom_fs_nrs2_ext.fits')

def make_extended_photom_robust():
    print(f"Reading {SRC}")
    with fits.open(SRC) as h:
        data = h[1].data
        cols = h[1].columns
        hdr  = h[1].header
        nrows = len(data)
        
        # Define new columns with larger shapes for the 1D arrays
        new_cols = []
        for c in cols:
            if c.name in ['wavelength', 'relresponse', 'reluncertainty']:
                # Change format from e.g. 150D to 500D
                new_fmt = '500D'
                # Initialize empty array for this column
                arr = np.zeros((nrows, 500), dtype='float64')
                new_cols.append(fits.Column(name=c.name, format=new_fmt, unit=c.unit, array=arr))
            else:
                # Keep other columns as they are
                new_cols.append(fits.Column(name=c.name, format=c.format, unit=c.unit, array=data[c.name].copy()))
        
        new_hdu = fits.BinTableHDU.from_columns(new_cols, header=hdr)
        new_data = new_hdu.data
        
        for i in range(nrows):
            row = data[i]
            filt = row['filter'].strip().upper()
            grat = row['grating'].strip().upper()
            ne   = int(row['nelem'])
            
            # Original arrays
            wl = row['wavelength'][:ne]
            rr = row['relresponse'][:ne]
            ru = row['reluncertainty'][:ne]
            
            if grat in ['G140M', 'G235M']:
                # Extend to 6.0 um
                ext_wl = np.linspace(wl[-1] + 0.01, 6.0, 50)
                wl = np.concatenate([wl, ext_wl])
                rr = np.concatenate([rr, np.ones_like(ext_wl)])
                ru = np.concatenate([ru, np.ones_like(ext_wl)])
                ne = len(wl)
            
            # Place into new_data
            new_data[i]['nelem'] = ne
            new_data[i]['wavelength'][:ne] = wl
            new_data[i]['relresponse'][:ne] = rr
            new_data[i]['reluncertainty'][:ne] = ru
            
        os.makedirs(os.path.dirname(OUT), exist_ok=True)
        fits.HDUList([h[0], new_hdu]).writeto(OUT, overwrite=True)
        print(f"Successfully saved robustly extended file: {OUT}")

if __name__ == "__main__":
    make_extended_photom_robust()
