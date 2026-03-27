import os
import sys
import numpy as np
from jwst import datamodels

input_path = os.path.abspath('/Users/dcoe/NIRSpec/wavext/data/PID1492/jw01492001001_17101_00006_nrs2_wcs.fits')

print(f"Inspecting WCS in {input_path}")
model = datamodels.open(input_path)

# AssignWcsStep adds the WCS to the meta
if hasattr(model.meta, 'wcs'):
    wcs = model.meta.wcs
    print(f"WCS object: {wcs}")
    
    # Try to evaluate WCS at some points
    # For NRS2, G140H, F100LP
    # G140H: [0.97, 1.89] um nominal.
    # Our extended file should have [0.9, 3.3] um. (Actually from WAVELENGTH_RANGES: F100LP_G140M/H: [0.9, 3.3] µm)
    
    # Evaluate WCS for a few pixels to see if we get wavelengths
    # NIRSpec NRS2 detector is 2048x2048
    x = np.arange(0, 2048, 512)
    y = np.arange(0, 2048, 512)
    xv, yv = np.meshgrid(x, y)
    
    # Slit name for FS. From Extract2dStep output it was S200A1
    slit_name = 'S200A1'
    names = np.array([slit_name] * len(xv.flatten()))
    
    # For NIRSpec, evaluation is usually (x, y, name) -> (ra, dec, lam)
    try:
        results = wcs(xv.flatten(), yv.flatten(), names)
        # Check if results has 3 components
        ra, dec, lam = results
        print(f"WCS (x, y, {slit_name}) -> (ra, dec, lam) evaluation:")
        for i in range(len(xv.flatten())):
            print(f"  ({xv.flatten()[i]}, {yv.flatten()[i]}) -> RA: {ra[i]:.6f}, DEC: {dec[i]:.6f}, Lam: {lam[i]:.4f}")
    except Exception as e:
        print(f"WCS evaluation failed: {e}")

else:
    print("Model has no WCS attribute in meta.")

# Check the bounding box
try:
    print(f"Bounding box: {wcs.bounding_box}")
except:
    print("Bounding box not defined or failed to retrieve.")
