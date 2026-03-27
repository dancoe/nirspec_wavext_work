import os
from stdatamodels.jwst import datamodels
import numpy as np

f = '/Users/dcoe/NIRSpec/wavext/data/PID1492/jw01492001001_17101_00006_nrs2_g395h_extract_1d.fits'
model = datamodels.open(f)
for spec in model.spec:
    wav = spec.spec_table['wavelength']
    valid_wav = wav[wav > 0]
    if len(valid_wav) > 0:
        print(f"Slit: {spec.name}, Wavelength range: {np.nanmin(valid_wav):.3f} to {np.nanmax(valid_wav):.3f} um")
