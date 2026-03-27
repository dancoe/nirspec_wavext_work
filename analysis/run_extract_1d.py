import os
import sys
import logging
from jwst.extract_1d import Extract1dStep
from stdatamodels.jwst import datamodels

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
log = logging.getLogger("jwst.extract_1d")
log.setLevel(logging.DEBUG)

input_path = os.path.abspath('/Users/dcoe/NIRSpec/wavext/data/PID1492/jw01492001001_0310a_00007_nrs2_g235h_extract_2d.fits')
output_path = os.path.abspath('/Users/dcoe/NIRSpec/wavext/data/PID1492/jw01492001001_0310a_00007_nrs2_g235h_extract_1d.fits')
# Override with the extended range file to ensure 1D extraction uses local wave_range
override_path = os.path.abspath('/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/reference_files/wavelengthrange_extended.asdf')

print(f"Running extract_1d for {input_path}")
model = datamodels.open(input_path)

step = Extract1dStep()
result = step.run(model)
print("Step.run() finished.")

if result:
    result.save(output_path)
    print(f"SUCCESS: Saved to {output_path}")
else:
    print("FAILURE: result is None")
