import sys
import os
import logging
from jwst.assign_wcs import AssignWcsStep
from stdatamodels.jwst import datamodels

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
log = logging.getLogger("jwst.assign_wcs")
log.setLevel(logging.DEBUG)

# Setup inputs
if len(sys.argv) > 2:
    input_path = os.path.abspath(sys.argv[1])
    output_path = os.path.abspath(sys.argv[2])
else:
    input_path = os.path.abspath('/Users/dcoe/NIRSpec/wavext/data/PID1492/jw01492001001_0310a_00007_nrs2_rate.fits')
    output_path = os.path.abspath('/Users/dcoe/NIRSpec/wavext/data/PID1492/jw01492001001_0310a_00007_nrs2_g235h_wcs.fits')

override_path = os.path.abspath('/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/reference_files/wavelengthrange_extended.asdf')

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
