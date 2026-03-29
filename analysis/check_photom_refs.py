"""Check if IFU photom reference covers NRS2 wavelengths."""
import os
import numpy as np
from astropy.io import fits
os.environ.setdefault('CRDS_PATH', os.path.expanduser('~/crds_cache'))
os.environ.setdefault('CRDS_SERVER_URL', 'https://jwst-crds.stsci.edu')
from jwst.photom.photom_step import PhotomStep

step = PhotomStep()

# FS NRS2 rate
fs_nrs2 = '/Users/dcoe/NIRSpec/wavext/data/PID1537_G191-B2B/jw01537007001_07101_00003_nrs2_rate.fits'
fs_ref = step.get_reference_file(fs_nrs2, 'photom')
print('FS NRS2 photom ref:', fs_ref)

with fits.open(fs_ref) as h:
    data = h[1].data
    mask = (data['filter'] == 'F100LP') & (data['grating'] == 'G140M')
    print('FS entries for F100LP/G140M:', mask.sum())
    if mask.sum() > 0:
        entry = data[mask][0]
        nelem = entry['nelem']
        wl = entry['wavelength']
        print('  wl range: %.3f - %.3f um' % (wl[:nelem].min(), wl[:nelem].max()))
        if 'slit' in data.dtype.names:
            print('  slit:', entry['slit'])
        # print all 5 entries
        for row in data[mask]:
            nelem2 = row['nelem']
            wl2 = row['wavelength']
            print('  entry: wl=[%.3f,%.3f] slit=%s' % (wl2[:nelem2].min(), wl2[:nelem2].max(),
                  row['slit'] if 'slit' in data.dtype.names else '?'))

# IFU NRS2 rate example
ifu_nrs2 = '/Users/dcoe/NIRSpec/wavext/data/IFU/PID1537_G191-B2B/stage1/jw01537007001_07101_00004_nrs2_rate.fits'
if os.path.exists(ifu_nrs2):
    ifu_ref = step.get_reference_file(ifu_nrs2, 'photom')
    print('IFU NRS2 photom ref:', ifu_ref)
    if ifu_ref != fs_ref:
        with fits.open(ifu_ref) as h:
            data = h[1].data
            mask_ifu = (data['filter'] == 'F100LP') & (data['grating'] == 'G140M')
            print('IFU entries for F100LP/G140M:', mask_ifu.sum())
            if mask_ifu.sum() > 0:
                entry = data[mask_ifu][0]
                nelem = entry['nelem']
                wl = entry['wavelength']
                print('  wl range: %.3f - %.3f um' % (wl[:nelem].min(), wl[:nelem].max()))
    else:
        print('Same photom ref for IFU and FS!')
else:
    print('IFU NRS2 rate file not found')
    # Try to find any IFU rate file
    import glob
    ifu_rates = glob.glob('/Users/dcoe/NIRSpec/wavext/data/IFU/*/stage1/*nrs2*rate*.fits')[:1]
    if ifu_rates:
        ifu_ref = step.get_reference_file(ifu_rates[0], 'photom')
        print('IFU photom ref:', ifu_ref)
        if ifu_ref != fs_ref:
            with fits.open(ifu_ref) as h:
                data = h[1].data
                print('IFU photom file columns:', list(data.dtype.names))
