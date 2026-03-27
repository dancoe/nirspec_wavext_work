import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

# Paths to relevant reference files in crds_cache
crds_dir = "/Users/dcoe/crds_cache/references/jwst/nirspec"

flats_to_plot = [
    # MOS
    {"file": "jwst_nirspec_sflat_0202.fits", "mode": "MOS", "type": "S-flat", "label": "PRISM/CLEAR"},
    {"file": "jwst_nirspec_fflat_0152.fits", "mode": "MOS", "type": "F-flat", "label": "PRISM/CLEAR"},
    # FS
    {"file": "jwst_nirspec_sflat_0147.fits", "mode": "FS", "type": "S-flat", "label": "G140M/F100LP"},
    {"file": "jwst_nirspec_fflat_0154.fits", "mode": "FS", "type": "F-flat", "label": "G140M/F100LP"},
    # IFU
    {"file": "jwst_nirspec_sflat_0191.fits", "mode": "IFU", "type": "S-flat", "label": "G140M/F100LP"},
    {"file": "jwst_nirspec_fflat_0173.fits", "mode": "IFU", "type": "F-flat", "label": "G140M/F100LP"},
    # DFLAT (NRS1 for demonstration)
    {"file": "jwst_nirspec_dflat_0001.fits", "mode": "All", "type": "D-flat", "label": "NRS1 Detector Flat"}
]

def plot_single_flat_extension(ref_info, output_dir):
    fits_file = os.path.join(crds_dir, ref_info["file"])
    if not os.path.exists(fits_file):
        # Fallback to local if not in crds_cache
        fits_file = os.path.join("/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work", ref_info["file"])

    if not os.path.exists(fits_file):
        print(f"Skipping {ref_info['file']} - Not found")
        return

    with fits.open(fits_file) as hdul:
        table_name = 'FAST_VARIATION' if 'FAST_VARIATION' in hdul else 'FLAT_TABLE'
        
        if table_name in hdul:
            data = hdul[table_name].data
            nelem = data['nelem'][0] if 'nelem' in data.names else None
            
            if nelem:
                wl = data['wavelength'][0][:nelem]
                flat = data['data'][0][:nelem]
            else:
                wl = data['wavelength'][0]
                flat = data['data'][0]

            # Extension Strategy (Parlanti)
            last_wl = wl[-1]
            if last_wl < 1e-4: # meters
                ext_range = np.linspace(last_wl, 5.61e-6, 50)[1:]
                wl_plot = wl * 1e6
                ext_plot = ext_range * 1e6
            else: # microns
                ext_range = np.linspace(last_wl, 5.61, 50)[1:]
                wl_plot = wl
                ext_plot = ext_range
            
            ext_flat = np.full_like(ext_range, flat[-1])

            plt.figure(figsize=(9, 5))
            plt.step(wl_plot, flat, label=f'Nominal {ref_info["type"]}', color='blue', where='mid')
            plt.step(ext_plot, ext_flat, label='wavext Extended Region', color='red', linestyle='--', where='mid', alpha=0.8)
            
            plt.xlabel('Wavelength (µm)')
            plt.ylabel('Flat Field Correction')
            plt.title(f'NIRSpec {ref_info["mode"]} {ref_info["type"]} Extension ({ref_info["label"]})\n{ref_info["file"]}')
            
            # Highlight extended region
            plt.axvspan(ext_plot[0], 5.61, color='orange', alpha=0.1, label='Extrapolation (Reddest Value)')
            plt.grid(alpha=0.25)
            plt.legend()
            
            basename = f"{ref_info['mode'].lower()}_{ref_info['type'].lower().replace('-', '')}_{ref_info['label'].split('/')[0].lower()}.png"
            output_path = os.path.join(output_dir, basename)
            plt.savefig(output_path, dpi=300)
            print(f"Saved: {output_path}")
            plt.close()
        else:
            print(f"Skipping {ref_info['file']} - No fast variation table")

if __name__ == "__main__":
    out_dir = "/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/plots/flats"
    os.makedirs(out_dir, exist_ok=True)
    for info in flats_to_plot:
        plot_single_flat_extension(info, out_dir)
