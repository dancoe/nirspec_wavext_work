"""
Create custom FS photom reference file that extends wavelength coverage to NRS2.

Copies jwst_nirspec_photom_0014.fits and extends the relresponse arrays
for the G140M/F100LP and G235M/F170LP entries to cover NRS2 wavelengths,
using relresponse = 1.0 (no change) for the extended range.

This prevents the photom step from issuing NaN fill_value for NRS2 pixels
and flagging them all as DO_NOT_USE.
"""
import os
import shutil
import numpy as np
from astropy.io import fits

SRC = '/Users/dcoe/crds_cache/references/jwst/nirspec/jwst_nirspec_photom_0014.fits'
DST = '/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/reference_files/jwst_nirspec_photom_fs_nrs2_ext.fits'
os.makedirs(os.path.dirname(DST), exist_ok=True)

# Extensions to add for NRS2 coverage
# key: (filter, grating), value: max NRS2 wavelength in microns
NRS2_EXT = {
    ('F100LP', 'G140M'): 3.20,  # G140M NRS2 goes to ~3.16 um
    ('F170LP', 'G235M'): 5.50,  # G235M NRS2 goes to ~5.31 um
    # Could also do G395M NRS2:
    # ('F290LP', 'G395M'): 6.00,
}

# New wavelength spacing to add for NRS2 (will be appended to the existing array)
NRS2_WSTEP = 0.02  # 20 nm steps, coarser than NRS1 is fine (relresponse = 1.0 everywhere)

shutil.copy(SRC, DST)
print(f'Copied: {os.path.basename(SRC)} → {os.path.basename(DST)}')

with fits.open(DST, mode='update') as h:
    data = h[1].data
    nrows = len(data)

    # Find the column shapes
    wl_col_len = data.dtype['wavelength'].shape[0]
    rel_col_len = data.dtype['relresponse'].shape[0]
    print(f'Ref file: {nrows} rows, wavelength array length={wl_col_len}')

    n_updated = 0
    for i in range(nrows):
        row = data[i]
        filt = row['filter']
        grat = row['grating']

        if (filt, grat) not in NRS2_EXT:
            continue

        nrs2_max = NRS2_EXT[(filt, grat)]
        nelem = row['nelem']
        wl = row['wavelength'].copy()
        rel = row['relresponse'].copy()

        nrs1_max = wl[nelem - 1]  # last NRS1 wavelength
        print(f'  Row {i}: {filt}/{grat}, nrs1_max={nrs1_max:.3f} um, nrs2_max={nrs2_max:.3f} um')

        if nrs1_max >= nrs2_max:
            print(f'    Already covers NRS2 range!')
            continue

        # Create new NRS2 wavelengths
        nrs2_wl = np.arange(nrs1_max + NRS2_WSTEP, nrs2_max + NRS2_WSTEP, NRS2_WSTEP)
        n_new = len(nrs2_wl)

        # New combined arrays
        new_wl = np.concatenate([wl[:nelem], nrs2_wl])
        new_rel = np.concatenate([rel[:nelem], np.ones(n_new, dtype=rel.dtype)])
        new_nelem = nelem + n_new

        if new_nelem > wl_col_len:
            print(f'    WARNING: need {new_nelem} elements but column only has {wl_col_len}!')
            # Pad to fit: use coarser spacing
            n_can_fit = wl_col_len - nelem
            nrs2_wl = np.linspace(nrs1_max + NRS2_WSTEP, nrs2_max, n_can_fit)
            n_new = len(nrs2_wl)
            new_wl = np.concatenate([wl[:nelem], nrs2_wl])
            new_rel = np.concatenate([rel[:nelem], np.ones(n_new, dtype=rel.dtype)])
            new_nelem = nelem + n_new

        # Update the arrays in-place (zero-padding at end)
        data[i]['nelem'] = new_nelem
        data[i]['wavelength'][:new_nelem] = new_wl[:new_nelem]
        data[i]['wavelength'][new_nelem:] = 0.0
        data[i]['relresponse'][:new_nelem] = new_rel[:new_nelem]
        data[i]['relresponse'][new_nelem:] = 0.0

        # Also update reluncertainty if present
        if 'reluncertainty' in data.dtype.names:
            data[i]['reluncertainty'][:new_nelem] = 1.0  # no uncertainty
            data[i]['reluncertainty'][new_nelem:] = 0.0

        n_updated += 1
        print(f'    Updated: nelem {nelem} → {new_nelem}, wl=[{new_wl.min():.3f},{new_wl.max():.3f}]')

    print(f'\nTotal rows updated: {n_updated}')
    h.flush()

print(f'\nSaved: {DST}')

# Verify
with fits.open(DST) as h:
    data = h[1].data
    for (filt, grat) in NRS2_EXT.keys():
        mask = (data['filter'] == filt) & (data['grating'] == grat) & (data['slit'] == 'S1600A1')
        if mask.sum() > 0:
            row = data[mask][0]
            nelem = row['nelem']
            wl = row['wavelength']
            print(f'Verify {filt}/{grat}/S1600A1: nelem={nelem}, wl=[{wl[:nelem].min():.3f},{wl[:nelem].max():.3f}]')
