import sys
from stdatamodels.jwst.datamodels import WavelengthrangeModel
import numpy as np

template = '/Users/dcoe/NIRSpec/wavext/pipeline/jwst_nirspec_wavext/jwst/assign_wcs/tests/data/waverange.asdf'
try:
    with WavelengthrangeModel(template) as model:
        print("Selector:", model.waverange_selector)
        for i, s in enumerate(model.waverange_selector):
            print(f"{i}: {s} -> {model.wavelengthrange[i]}")
except Exception as e:
    print(f"Error: {e}")
    sys.exit(1)
