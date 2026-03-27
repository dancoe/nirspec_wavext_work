import os
import sys
import logging
from jwst.extract_2d import Extract2dStep
from stdatamodels.jwst import datamodels

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
log = logging.getLogger("jwst.extract_2d")
log.setLevel(logging.DEBUG)

input_path = os.path.abspath('/Users/dcoe/NIRSpec/wavext/data/PID1492/jw01492001001_0310a_00007_nrs2_g235h_wcs.fits')
output_path = os.path.abspath('/Users/dcoe/NIRSpec/wavext/data/PID1492/jw01492001001_0310a_00007_nrs2_g235h_extract_2d.fits')
# Use the override if needed, though wcs should already have it
override_path = os.path.abspath('/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/reference_files/wavelengthrange_extended.asdf')

print(f"Running extract_2d for {input_path}")

# Open the model
model = datamodels.open(input_path)
print(f"Model {model.meta.filename} opened.")

# Initialize the step with override
override_path = os.path.abspath('/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/reference_files/wavelengthrange_extended.asdf')
step = Extract2dStep()
step.override_wavelengthrange = override_path

# Execute step
result = step.run(model)
print("Step.run() finished.")

if result:
    result.save(output_path)
    print(f"SUCCESS: Saved to {output_path}")
    
    # Check the result
    if hasattr(result, 'slits'):
        print(f"Total slits: {len(result.slits)}")
        for i, slit in enumerate(result.slits):
            if slit.wavelength is not None:
                print(f"Slit {i} ({slit.name}) wavelength range: {slit.wavelength.min():.4f} - {slit.wavelength.max():.4f} um")
            else:
                print(f"Slit {i} ({slit.name}) wavelength is None!")
    else:
        print("Resulting model has no 'slits' attribute.")
else:
    print("FAILURE: result is None")
