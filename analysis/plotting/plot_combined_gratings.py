import os
import matplotlib.pyplot as plt
import numpy as np
from stdatamodels.jwst import datamodels

def get_spectrum(filename):
    with datamodels.open(filename) as model:
        if len(model.spec) == 0:
            return None, None
        spec = model.spec[0]
        w = spec.spec_table['wavelength']
        f = spec.spec_table['flux']
        
        # Mask out NaNs and non-positive values for log plotting
        mask = np.isfinite(f) & (f > 1e-6)
        return w[mask], f[mask]

data_dir = '/Users/dcoe/NIRSpec/wavext/data/PID1492'
plots_dir = '/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/plots'

# Configuration
gratings = [
    {
        'name': 'G140M',
        'nominal_file': f'{data_dir}/jw01492-o001_t001-s000000001_nirspec_f100lp-g140m-s200a1-subs200a1_x1d.fits',
        'extended_file': f'{data_dir}/jw01492003001_03102_00005_nrs2_g140m_extract_1d.fits',
        'nom_color': 'blue',
        'ext_color': 'cornflowerblue'
    },
    {
        'name': 'G235M',
        'nominal_file': f'{data_dir}/jw01492-o001_t001-s000000001_nirspec_f170lp-g235m-s200a1-subs200a1_x1d.fits',
        'extended_file': f'{data_dir}/jw01492003001_03104_00004_nrs2_g235m_extract_1d.fits',
        'nom_color': 'darkgoldenrod',
        'ext_color': 'goldenrod'
    },
    {
        'name': 'G395M',
        'nominal_file': f'{data_dir}/jw01492-o001_t001-s000000001_nirspec_f290lp-g395m-s200a1-subs200a1_x1d.fits',
        'extended_file': f'{data_dir}/jw01492003001_03106_00002_nrs2_g395m_extract_1d.fits',
        'nom_color': 'red',
        'ext_color': 'lightcoral'
    }
]

prism_file = f'{data_dir}/jw01492-o001_t001-s000000001_nirspec_clear-prism-s200a1-subs200a1_x1d.fits'

plt.figure(figsize=(14, 8))

all_spectra = []
all_fluxes = []

# Collect data
# M-Gratings
for g in gratings:
    # Nominal
    nw, nf = get_spectrum(g['nominal_file'])
    if nw is not None and len(nw) > 0:
        all_spectra.append((nw, nf, g['nom_color'], f"{g['name']} Nominal", 0.6, 1.0, 5))
        all_fluxes.extend(nf)
    
    # Extended
    ew, ef = get_spectrum(g['extended_file'])
    if ew is not None and len(ew) > 0:
        all_spectra.append((ew, ef, g['ext_color'], f"{g['name']} Extended (S(λ))", 0.6, 0.8, 10))
        all_fluxes.extend(ef)

# PRISM baseline (Front)
pw, pf = get_spectrum(prism_file)
if pw is not None and len(pw) > 0:
    all_spectra.append((pw, pf, 'black', 'PRISM (Reference f(λ))', 0.6, 1.5, 20))
    all_fluxes.extend(pf)

# Plot everything in order of zorder
for w, f, c, l, a, lw, z in sorted(all_spectra, key=lambda x: x[-1]):
    plt.plot(w, f, color=c, label=l, alpha=a, lw=lw, zorder=z)

# Robust axis scaling
if all_fluxes:
    all_fluxes = np.array(all_fluxes)
    pos_fluxes = all_fluxes[all_fluxes > 0]
    if len(pos_fluxes) > 0:
        # Use 0.1% and 99.9% to avoid edge noise/spikes but capture real dynamic range
        ymin = np.percentile(pos_fluxes, 0.1) * 0.2
        ymax = np.percentile(pos_fluxes, 99.9) * 5.0
        plt.ylim(max(ymin, 1e-6), ymax)

plt.yscale('log')
plt.xlabel('Wavelength (µm)')
plt.ylabel('Flux (MJy/sr)')
plt.title('NIRSpec Multi-Resolution Spectral Overlap (PID 1492 FS)')
plt.legend(ncol=2, fontsize='small', loc='upper right', frameon=True)
plt.grid(False) # Removed interior grid lines
plt.xlim(0.6, 5.6)

output_plot = f'{plots_dir}/FS-1492_pre-cal.png'
plt.savefig(output_plot, dpi=200)
print(f"Combined plot saved to: {output_plot}")
