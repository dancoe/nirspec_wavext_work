"""Check IFU photom reference photomj scalar and how it's applied."""
import numpy as np
from astropy.io import fits

ifu_photom = '/Users/dcoe/crds_cache/references/jwst/nirspec/jwst_nirspec_photom_0016.fits'
fs_photom = '/Users/dcoe/crds_cache/references/jwst/nirspec/jwst_nirspec_photom_0014.fits'

print('=== IFU photom (0016) ===')
with fits.open(ifu_photom) as h:
    data = h[1].data
    mask = (data['filter'] == 'F100LP') & (data['grating'] == 'G140M')
    if mask.sum() > 0:
        row = data[mask][0]
        print('photomj:', row['photmj'])
        print('nelem:', row['nelem'])
        print('relresponse (first 5):', row['relresponse'][:5])
        print('relresponse (last 5):', row['relresponse'][row['nelem']-5:row['nelem']])

print('\n=== FS photom (0014) ===')
with fits.open(fs_photom) as h:
    data = h[1].data
    mask = (data['filter'] == 'F100LP') & (data['grating'] == 'G140M') & (data['slit'] == 'S1600A1')
    if mask.sum() > 0:
        row = data[mask][0]
        print('photomj:', row['photmj'])
        print('nelem:', row['nelem'])
        wl = row['wavelength']
        rel = row['relresponse']
        nelem = row['nelem']
        print('wl (first 5):', wl[:5])
        print('wl (last 5):', wl[nelem-5:nelem])
        print('relresponse (first 5):', rel[:5])
        print('relresponse (last 5):', rel[nelem-5:nelem])

# The key difference: does FS photom code extrapolate relresponse to 0 for NRS2 wavelengths?
# (which would flag all NRS2 pixels as DO_NOT_USE)
# Check the photom step code
import inspect
from jwst.photom.photom import DataSet
src = inspect.getsource(DataSet.photom_io)
# Look for out-of-range wavelength handling
for i, line in enumerate(src.split('\n')):
    if 'fill' in line.lower() or 'bounds_error' in line.lower() or 'nan' in line.lower() or 'zero' in line.lower():
        print(f'Line ~{i}: {line}')
