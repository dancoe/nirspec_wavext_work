import sys
import os
from stdatamodels.jwst import datamodels
import numpy as np

# Use the extracted file
path = '/Users/dcoe/NIRSpec/wavext/data/PID1492/jw01492001001_17101_00006_nrs2_g395h_extract_2d.fits'
model = datamodels.open(path)

print(f"Number of slits: {len(model.slits)}")
slit = model.slits[0]
print(f"Slit name: {slit.name}")
print(f"Wavelength array shape: {slit.wavelength.shape}")
print(f"Wavelength finite values: {np.count_nonzero(np.isfinite(slit.wavelength))}")
if np.count_nonzero(np.isfinite(slit.wavelength)) > 0:
    print(f"Min: {np.nanmin(slit.wavelength)}, Max: {np.nanmax(slit.wavelength)}")
else:
    print("All wavelength are NaN")

print(f"Extracted WCS available: {hasattr(slit.meta, 'wcs')}")
if hasattr(slit.meta, 'wcs'):
    wcs = slit.meta.wcs
    print(f"Extracted WCS bounding box: {wcs.bounding_box}")
    # Try evaluate at 0,0 (center of a pixel in the extraction)
    print(f"Evaluating at x=1024, y=28 (middle of extraction area)")
    try:
        res = wcs(1024, 28)
        print(f"Result: {res}")
    except Exception as e:
        print(f"Evaluation failed: {e}")
