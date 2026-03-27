import sys
import os
from stdatamodels.jwst import datamodels
import numpy as np

if len(sys.argv) < 2:
    print("Usage: python inspect_wcs_range.py <fits_file>")
    sys.exit(1)

f = os.path.abspath(sys.argv[1])
print(f"Inspecting WCS and wavelength range for {f}...")

model = datamodels.open(f)
grating = model.meta.instrument.grating
filt = model.meta.instrument.filter
detector = model.meta.instrument.detector
print(f"Instrument: {grating} / {filt} on {detector}")

# Try to find the slitlets and their wavelengths
if hasattr(model, 'slits'):
    for slit in model.slits:
        wcs = slit.meta.wcs
        # Evaluate WCS on the central line
        ny, nx = slit.data.shape
        x = np.arange(nx)
        y = np.ones(nx) * (ny / 2.0)
        _, _, wav = wcs(x, y)
        valid_wav = wav[wav > 0]
        if len(valid_wav) > 0:
            print(f"Slit: {slit.name}, Shape: {slit.data.shape}, WCS range: {np.nanmin(valid_wav):.3f} to {np.nanmax(valid_wav):.3f} um")
        else:
             print(f"Slit: {slit.name}, no valid wavelengths found.")
elif hasattr(model, 'spec'):
    # Multi-slit extracted 1D
    for spec in model.spec:
        wav = spec.spec_table['wavelength']
        valid_wav = wav[wav > 0]
        if len(valid_wav) > 0:
            print(f"1D Slit: {spec.name}, range: {np.nanmin(valid_wav):.3f} to {np.nanmax(valid_wav):.3f} um")
else:
     print("Data model doesn't have standard slits or spec structures.")
