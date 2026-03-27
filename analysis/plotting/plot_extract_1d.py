import os
import matplotlib.pyplot as plt
from stdatamodels.jwst import datamodels
from stdatamodels.jwst.datamodels import WavelengthrangeModel
import numpy as np

# Use a non-interactive backend for server-side plotting
import matplotlib
matplotlib.use('Agg')

import sys

if len(sys.argv) > 1:
    input_file = os.path.abspath(sys.argv[1])
else:
    input_file = '/Users/dcoe/NIRSpec/wavext/data/PID1492/jw01492001001_17101_00006_nrs2_g395h_extract_1d.fits'

output_dir = '/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/plots'
nominal_ref = '/Users/dcoe/crds_cache/references/jwst/nirspec/jwst_nirspec_wavelengthrange_0006.asdf'
os.makedirs(output_dir, exist_ok=True)

print(f"Loading {input_file}...")
model = datamodels.open(input_file)

# Determine grating and filter
grating = model.meta.instrument.grating
filt = model.meta.instrument.filter
selector_key = f"{filt}_{grating}"
print(f"Grating: {grating}, Filter: {filt} (Selector: {selector_key})")

# Load nominal range
nominal_range = [None, None]
try:
    with WavelengthrangeModel(nominal_ref) as nm:
        if selector_key in nm.waverange_selector:
            idx = list(nm.waverange_selector).index(selector_key)
            nominal_range = [val * 1e6 for val in nm.wavelengthrange[idx]]  # Convert to microns
            print(f"Nominal range from CRDS: {nominal_range} um")
except Exception as e:
    print(f"Warning: Could not load nominal range from {nominal_ref}: {e}")

fig, ax = plt.subplots(figsize=(12, 6))

for spec in model.spec:
    wav = spec.spec_table['wavelength']
    flux = spec.spec_table['flux']
    
    # Filter out invalid wavelengths (<= 0)
    valid = (wav > 0) & (~np.isnan(flux))
    
    # Sort by wavelength for a clean line plot
    sort_idx = np.argsort(wav[valid])
    
    ax.plot(wav[valid][sort_idx], flux[valid][sort_idx], label=f"Slit {spec.name}")

# Highlight extended regions
if nominal_range[0] is not None:
    # Get total plot limits
    xlim = ax.get_xlim()
    
    # Shade below nominal_range[0] if relevant
    if xlim[0] < nominal_range[0]:
        ax.axvspan(xlim[0], nominal_range[0], color='lightyellow', alpha=0.5, label='Extended (Blue)')
    
    # Shade above nominal_range[1] if relevant
    if xlim[1] > nominal_range[1]:
        ax.axvspan(nominal_range[1], xlim[1], color='lightyellow', alpha=0.5, label='Extended (Red)')

ax.set_xlabel('Wavelength (um)')
ax.set_ylabel(f'Flux ({model.spec[0].spec_table.columns["flux"].unit})')
ax.set_title(f'NIRSpec 1D Extraction: {os.path.basename(input_file)}\n({grating}/{filt})')
ax.legend()
ax.grid(True, alpha=0.3)

# Save the plot
output_file = os.path.join(output_dir, f"{os.path.basename(input_file).replace('.fits', '.png')}")
plt.savefig(output_file, dpi=200)
print(f"SUCCESS: Plot saved to {output_file}")
plt.close()

