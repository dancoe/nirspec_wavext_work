import sys
from stdatamodels.jwst.datamodels import WavelengthrangeModel
import numpy as np

template = '/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/reference_files/wavelengthrange_extended.asdf'
try:
    with WavelengthrangeModel(template) as model:
        print("Selector:", model.waverange_selector)
        for i, s in enumerate(model.waverange_selector):
            print(f"{i}: {s} -> {list(model.wavelengthrange[i])}")
except Exception as e:
    print(f"Error: {e}")
    sys.exit(1)
