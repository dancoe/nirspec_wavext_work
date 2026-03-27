import os
import sys
import numpy as np
from jwst import datamodels

input_path = os.path.abspath('/Users/dcoe/NIRSpec/wavext/data/PID1492/jw01492001001_17101_00006_nrs2_wcs.fits')

print(f"Inspecting WCS in {input_path}")
model = datamodels.open(input_path)
wcs = model.meta.wcs

# Slit name
slit_name = 'S200A1'

# Try a full grid of points to find where wavelength is NOT nan
x = np.linspace(0, 2047, 100)
y = np.linspace(0, 2047, 100)
xv, yv = np.meshgrid(x, y)
names = np.array([slit_name] * len(xv.flatten()))

try:
    ra, dec, lam = wcs(xv.flatten(), yv.flatten(), names)
    print(f"Result shape: {lam.shape}")
    valid_lam = lam[~np.isnan(lam)]
    print(f"Found {len(valid_lam)} valid wavelength points out of {len(xv.flatten())}")
    if len(valid_lam) > 0:
        print(f"Valid wavelength range: {valid_lam.min():.4f} - {valid_lam.max():.4f} um")
        # Find some coordinates
        idx = np.where(~np.isnan(lam))[0]
        for i in idx[:5]:
             print(f"  ({xv.flatten()[i]:.1f}, {yv.flatten()[i]:.1f}) -> {lam[i]:.4f}")
    else:
        print("No valid wavelengths found in the grid!")

    # Check the bounding box again
    print(f"Bounding box: {wcs.bounding_box}")
except Exception as e:
    print(f"WCS evaluation failed: {e}")
    
# Let's inspect the wavelengthrange attribute of the model if it's there
# JwstDataModel might have different meta for spectral range
print(f"Meta wavelength range start: {model.meta.wcsinfo.waverange_start}")
print(f"Meta wavelength range end: {model.meta.wcsinfo.waverange_end}")
print(f"Meta spectral order: {model.meta.wcsinfo.spectral_order}")
