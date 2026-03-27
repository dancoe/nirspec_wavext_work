import os
import sys
from stdatamodels.jwst.datamodels import WavelengthrangeModel

ref_path = os.path.abspath('/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/reference_files/wavelengthrange_extended.asdf')

print(f"Inspecting reference file: {ref_path}")
try:
    with WavelengthrangeModel(ref_path) as model:
        print(f"Selector: {model.waverange_selector}")
        for i, s in enumerate(model.waverange_selector):
            if 'G140H' in s:
                # Force load it
                wr = list(model.wavelengthrange[i])
                print(f"{s}: {wr}")
except Exception as e:
    print(f"Error: {e}")
