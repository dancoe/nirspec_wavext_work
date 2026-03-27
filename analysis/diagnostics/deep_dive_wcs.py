import os
import sys
import numpy as np
from jwst import datamodels

input_path = os.path.abspath('/Users/dcoe/NIRSpec/wavext/data/PID1492/jw01492001001_17101_00006_nrs2_wcs.fits')

print(f"Deep dive into WCS in {input_path}")
model = datamodels.open(input_path)

if hasattr(model.meta, 'wcs'):
    wcs = model.meta.wcs
    print(f"WCS inputs: {wcs.inputs}")
    print(f"WCS outputs: {wcs.outputs}")
    
    bb = wcs.bounding_box
    print(f"Bounding box type: {type(bb)}")
    
    if hasattr(bb, 'bounding_boxes'):
        print(f"Known keys: {bb.bounding_boxes.keys()}")
    
    if hasattr(bb, 'selector'):
        print(f"Selector: {bb.selector}")
        
    # Let's try to evaluate one point with name='S200A1' and check what the selector returns
    try:
        # Evaluate one point
        # x, y, name
        p = (1024, 1024, 'S200A1')
        print(f"Evaluating at {p}")
        # wcs(*p)
        # result = wcs(1024, 1024, 'S200A1') # This failed before
    except Exception as e:
        print(f"Evaluation failed: {e}")

    # Let's check the Slits in the model if it's a MultiSlitModel later
    # But this is still the input to extract_2d, so it's likely an ImageModel with wcs.
    
    # Check the WCS pipeline frames
    print("Pipeline frames:")
    for frame in wcs.available_frames:
        print(f"  - {frame}")

    # Inspect the 'detector' -> 'msa_frame' part
    try:
        det2msa = wcs.get_transform('detector', 'msa_frame')
        print(f"det2msa inputs: {det2msa.inputs}")
        # If it uses 'name', what does it do with it?
        # Usually it uses Mapping or selector to filter.
    except:
        print("Failed to get det2msa transform")

else:
    print("Model has no WCS attribute in meta.")
