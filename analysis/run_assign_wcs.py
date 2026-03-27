import sys
import os
from jwst.assign_wcs import AssignWcsStep
from jwst import datamodels

input_path = os.path.abspath('/Users/dcoe/NIRSpec/wavext/data/PID1492/jw01492001001_03108_00007_nrs2_rate.fits')
output_path = os.path.abspath('/Users/dcoe/NIRSpec/wavext/data/PID1492/jw01492001001_17101_00006_nrs2_wcs.fits')
# Point to actual PRISM NRS2 file I downloaded? No, use the nrs2_rate.fits
override_path = os.path.abspath('/Users/dcoe/NIRSpec/wavext/pipeline/wavelengthrange_extended.asdf')

print(f"Assigning WCS to {input_path}")
print(f"Override: {override_path}")

# Run natively in jwst_1.20.2 + asdf 4.5.0 environment. No silent catching.
model = datamodels.open(input_path)
print(f"Model {model.meta.filename} opened.")

step = AssignWcsStep()
step.override_wavelengthrange = override_path

# Execute step
result = step.run(model)
print("Step.run() finished.")

if result:
    result.save(output_path)
    print(f"SUCCESS: Saved to {output_path}")
else:
    print("FAILURE: result is None")
