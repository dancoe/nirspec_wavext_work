import os
import sys
import numpy as np
from jwst import datamodels

input_path = os.path.abspath('/Users/dcoe/NIRSpec/wavext/data/PID1492/jw01492001001_17101_00006_nrs2_wcs.fits')

print(f"Centroid search in {input_path}")
model = datamodels.open(input_path)
wcs = model.meta.wcs

# Try smaller grid to avoid truncation
x = np.linspace(0, 2047, 10)
y = np.linspace(0, 2047, 10)
xv, yv = np.meshgrid(x, y)
names = np.array([-101] * len(xv.flatten())) # Try the integer ID (-101 for S200A1)

try:
    print(f"Checking {len(xv.flatten())} points...", flush=True)
    results = wcs(xv.flatten(), yv.flatten(), names)
    ra, dec, lam = results
    
    valid_mask = ~np.isnan(lam)
    valid_count = np.sum(valid_mask)
    print(f"Found {valid_count} valid wavelength points")
    if valid_count > 0:
        xv_valid = xv.flatten()[valid_mask]
        yv_valid = yv.flatten()[valid_mask]
        lam_valid = lam[valid_mask]
        print(f"X mean: {np.mean(xv_valid):.1f}, Y mean: {np.mean(yv_valid):.1f}")
        print(f"Wavelength range: {np.min(lam_valid):.4f} - {np.max(lam_valid):.4f} um")
        
        # Let's check another slit name
        names2 = np.array(['S1600A1'] * len(xv.flatten()))
        results2 = wcs(xv.flatten(), yv.flatten(), names2)
        valid_mask2 = ~np.isnan(results2[2])
        print(f"S1600A1 valid points: {np.sum(valid_mask2)}")
        
    else:
        print("No valid wavelengths anywhere on NRS2!")
except Exception as e:
    print(f"Centroid search failed: {e}")
