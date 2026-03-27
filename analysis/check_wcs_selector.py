import os
import sys
import numpy as np
from jwst import datamodels

input_path = os.path.abspath('/Users/dcoe/NIRSpec/wavext/data/PID1492/jw01492001001_17101_00006_nrs2_wcs.fits')

print(f"Inspecting WCS in {input_path}")
model = datamodels.open(input_path)
wcs = model.meta.wcs

bb = wcs.bounding_box
if hasattr(bb, 'bounding_boxes'):
    print(f"Bounding boxes dict: {bb.bounding_boxes}")
    for k in bb.bounding_boxes.keys():
        print(f"Key: {k}, type: {type(k)}")
        
if hasattr(bb, 'selector'):
    print(f"Selector type: {type(bb.selector)}")
    # Try different values for name to see which key matches
    candidates = ['S200A1', 'S200A2', 'S400A1', 'S1600A1', 'S200B1', '-1', '-2', '-3', '-4', '-5']
    for cand in candidates:
        try:
            matched_key = bb.selector(cand)
            print(f"Candidate '{cand}' matched key: {matched_key}")
        except:
            pass
else:
    print("Bounding box has no 'selector'")
