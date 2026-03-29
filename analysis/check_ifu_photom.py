"""Check IFU photom reference wavelength coverage for NRS2."""
import numpy as np
from astropy.io import fits

ifu_photom = '/Users/dcoe/crds_cache/references/jwst/nirspec/jwst_nirspec_photom_0016.fits'
with fits.open(ifu_photom) as h:
    data = h[1].data
    print('Columns:', list(data.dtype.names))
    print('Filters (unique):', list(set(data['filter'])))
    print('Gratings (unique):', list(set(data['grating'])))
    mask = (data['filter'] == 'F100LP') & (data['grating'] == 'G140M')
    print('\nIFU F100LP/G140M entries:', mask.sum())
    for row in data[mask]:
        nelem = row['nelem']
        wl = row['wavelength']
        rel = row['relresponse']
        print(f'  nelem={nelem}, wl_range=[{wl[:nelem].min():.3f},{wl[:nelem].max():.3f}] um')
        print(f'  rel_range=[{rel[:nelem].min():.4g},{rel[:nelem].max():.4g}]')
        
    mask2 = (data['filter'] == 'F170LP') & (data['grating'] == 'G235M')
    print('\nIFU F170LP/G235M entries:', mask2.sum())
    for row in data[mask2]:
        nelem = row['nelem']
        wl = row['wavelength']
        print(f'  nelem={nelem}, wl_range=[{wl[:nelem].min():.3f},{wl[:nelem].max():.3f}] um')
