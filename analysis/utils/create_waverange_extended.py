from stdatamodels.jwst.datamodels import WavelengthrangeModel
import numpy as np

template = '/Users/dcoe/NIRSpec/wavext/pipeline/jwst_nirspec_wavext/jwst/assign_wcs/tests/data/waverange.asdf'
output = '/Users/dcoe/NIRSpec/wavext/pipeline/wavelengthrange_extended.asdf'

# Dictionary for updates (ranges in microns; will convert to meters)
updates = {
    'F070LP_G140M': [0.6, 3.3],
    'F070LP_G140H': [0.6, 3.3],
    'F100LP_G140M': [0.9, 3.3],
    'F100LP_G140H': [0.9, 3.3],
    'F170LP_G235M': [1.5, 5.3],
    'F170LP_G235H': [1.5, 5.3],
    'F290LP_G395M': [2.6, 5.6],
    'F290LP_G395H': [2.6, 5.6],
    'CLEAR_PRISM': [0.5, 5.6],
    'F110W_G395H': [0.7, 5.6]
}

with WavelengthrangeModel(template) as model:
    print(f"Opening template: {template}")
    
    # Iterate over the selector to find and update ranges
    for i, s in enumerate(model.waverange_selector):
        if s in updates:
            old_range = model.wavelengthrange[i]
            new_range = np.array(updates[s]) * 1e-6
            model.wavelengthrange[i] = new_range
            print(f"Updating {s}: {old_range} -> {new_range}")
        elif 'G235' in s or 'G140' in s or 'G395' in s or 'PRISM' in s:
            # Check for partial matches
            for key, val in updates.items():
                if key in s:
                    old_range = model.wavelengthrange[i]
                    new_range = np.array(val) * 1e-6
                    model.wavelengthrange[i] = new_range
                    print(f"Updating {s} (partial match): {old_range} -> {new_range}")
                    break

    model.save(output)
    print(f"Saved extended range file to: {output}")
