import os
import sys
import logging
from jwst.extract_1d import Extract1dStep
from stdatamodels.jwst import datamodels

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
log = logging.getLogger("jwst.extract_1d")
log.setLevel(logging.DEBUG)

input_path = os.path.abspath('/Users/dcoe/NIRSpec/wavext/data/PID1492/jw01492001001_17101_00006_nrs2_g395h_extract_2d.fits')
output_path = os.path.abspath('/Users/dcoe/NIRSpec/wavext/data/PID1492/jw01492001001_17101_00006_nrs2_g395h_extract_1d.fits')

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
