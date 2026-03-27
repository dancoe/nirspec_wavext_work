#!/usr/bin/env python
"""Check units and flux levels from FITS files."""
import sys, os
sys.path.insert(0, '/Users/dcoe/NIRSpec/wavext/jwst_nirspec_wavext')
os.environ.setdefault('CRDS_PATH', os.path.expanduser('~/crds_cache'))
os.environ.setdefault('CRDS_SERVER_URL', 'https://jwst-crds.stsci.edu')
from astropy.io import fits
import numpy as np

files = {
    'PRISM':     '/Users/dcoe/NIRSpec/wavext/data/PID1492/jw01492-o001_t001-s000000001_nirspec_clear-prism-s200a1-subs200a1_x1d.fits',
    'G140M_nom': '/Users/dcoe/NIRSpec/wavext/data/PID1492/jw01492-o001_t001-s000000001_nirspec_f100lp-g140m-s200a1-subs200a1_x1d.fits',
    'G235M_nom': '/Users/dcoe/NIRSpec/wavext/data/PID1492/jw01492-o001_t001-s000000001_nirspec_f170lp-g235m-s200a1-subs200a1_x1d.fits',
    'G395M_nom': '/Users/dcoe/NIRSpec/wavext/data/PID1492/jw01492-o001_t001-s000000001_nirspec_f290lp-g395m-s200a1-subs200a1_x1d.fits',
    'G140M_ext': '/Users/dcoe/NIRSpec/wavext/data/PID1492/jw01492003001_03102_00005_nrs2_g140m_extract_1d.fits',
    'G235M_ext': '/Users/dcoe/NIRSpec/wavext/data/PID1492/jw01492003001_03104_00004_nrs2_g235m_extract_1d.fits',
    'G395M_ext': '/Users/dcoe/NIRSpec/wavext/data/PID1492/jw01492003001_03106_00002_nrs2_g395m_extract_1d.fits',
}

for name, fpath in files.items():
    with fits.open(fpath) as h:
        # Print extension list
        print(f'\n=== {name} ===')
        print(f'  HDU names: {[e.name for e in h]}')
        # Check primary header for target
        hdr = h[0].header
        print(f'  TARGET: {hdr.get("TARGNAME", "?")}  PROG: {hdr.get("PROGRAM", "?")}  OBS: {hdr.get("OBSERVTN", "?")}')
        # Check EXTRACT1D or SPEC binary table
        for ext in h:
            if ext.name in ('EXTRACT1D', 'BINTABLE', 'SPEC') or hasattr(ext, 'columns'):
                colnames = [c.name for c in ext.columns] if hasattr(ext, 'columns') else []
                if 'FLUX' in colnames or 'WAVELENGTH' in [c.upper() for c in colnames]:
                    print(f'  Ext {ext.name}: columns={colnames[:8]}')
                    # Check TUNIT for flux
                    for key in ext.header.keys():
                        if 'TUNIT' in key:
                            print(f'    {key} = {ext.header[key]}')
                    # Print sample values
                    if hasattr(ext, 'data') and ext.data is not None:
                        w_col = 'WAVELENGTH' if 'WAVELENGTH' in [c.name.upper() for c in ext.columns] else 'wavelength'
                        f_col = 'FLUX' if 'FLUX' in [c.name.upper() for c in ext.columns] else 'flux'
                        try:
                            w = ext.data[w_col.lower() if w_col.lower() in ext.data.names else w_col]
                            f = ext.data[f_col.lower() if f_col.lower() in ext.data.names else f_col]
                            msk = (w > 0) & np.isfinite(f)
                            if msk.sum() > 0:
                                print(f'    wav range: [{w[msk].min():.3f}, {w[msk].max():.3f}]')
                                print(f'    flux median/max: {np.median(f[msk]):.3e} / {f[msk].max():.3e}')
                        except Exception as e:
                            print(f'    Error reading data: {e}')
                    break
