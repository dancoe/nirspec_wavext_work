"""Quick check of flux units in NRS2 per-exposure vs MAST L3 extract_1d files."""
from astropy.io import fits
import numpy as np

base = '/Users/dcoe/NIRSpec/wavext'
files = [
    ('NRS2_1537_G140M', base+'/data/PID1537_G191-B2B/jw01537007001_07101_00003_nrs2_extract_1d.fits'),
    ('L3_1537_G140M',   base+'/data/PID1537_G191-B2B/jw01537-o007_t001-s000000001_nirspec_f100lp-g140m-s1600a1-sub2048_x1d.fits'),
    ('NRS2_1538_G140M', base+'/data/PID1538_P330E/jw01538160001_06101_00001_nrs2_extract_1d.fits'),
    ('L3_1538_G140M',   base+'/data/PID1538_P330E/jw01538-o160_t002-s000000001_nirspec_f100lp-g140m-s1600a1-sub2048_x1d.fits'),
]

for label, path in files:
    with fits.open(path) as h:
        ext = h[1]
        wl = ext.data['WAVELENGTH'].astype(float)
        fl = ext.data['FLUX'].astype(float)
        mask = (wl > 1.97) & (wl < 2.1) & np.isfinite(fl) & (fl > 0)
        med = float(np.nanmedian(fl[mask])) if mask.sum() > 0 else float('nan')
        colunit = str(ext.columns['FLUX'].unit) if ext.columns['FLUX'].unit else 'NONE'
        htunit  = ext.header.get('TUNIT2', 'NONE')
        bunit   = ext.header.get('BUNIT',  'NONE')
        print(f'{label}: med={med:.4g}, col_unit={colunit}, TUNIT2={htunit}, BUNIT={bunit}')
