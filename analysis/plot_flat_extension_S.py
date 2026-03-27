import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

# Example S-flat file found in CRDS cache
sflat_file = "/Users/dcoe/crds_cache/references/jwst/nirspec/jwst_nirspec_sflat_0202.fits"

def plot_flat_extension(fits_file):
    if not os.path.exists(fits_file):
        print(f"File not found: {fits_file}")
        return

    with fits.open(fits_file) as hdul:
        print(f"Opening {fits_file}")
        print(f"HDU names: {[h.name for h in hdul]}")
        
        # S-flats often have 'FLAT_TABLE' or 'FAST_VARIATION'
        table_name = 'FAST_VARIATION' if 'FAST_VARIATION' in hdul else 'FLAT_TABLE'
        if table_name in hdul:
            data = hdul[table_name].data
            print(f"Found {table_name} with {len(data)} rows")
            
            # For FAST_VARIATION, it might be a different structure
            if table_name == 'FAST_VARIATION':
                nelem = data['nelem'][0]
                wl = data['wavelength'][0][:nelem]
                flat = data['data'][0][:nelem]
                print(f"Slicing arrays at nelem={nelem}")
            else: # FLAT_TABLE
                wl = data['wavelength'][0]
                flat = data['data'][0]
            print(f"Wavelength range: {wl[0]} - {wl[-1]}")
            
            # Simulated extension (Parlanti strategy)
            # Extension range: last_wl to 5.6 microns
            last_wl = wl[-1]
            if last_wl < 1e-4: # meters
                ext_range = np.linspace(last_wl, 5.6e-6, 50)[1:]
            else: # microns?
                ext_range = np.linspace(last_wl, 5.6, 50)[1:]
            
            ext_flat = np.full_like(ext_range, flat[-1])
            
            plt.figure(figsize=(10, 6))
            plt.step(wl * 1e6 if last_wl < 1e-4 else wl, flat, label='Nominal Fast Variation', color='blue', where='mid')
            plt.step(ext_range * 1e6 if last_wl < 1e-4 else ext_range, ext_flat, label='wavext Extension (Reddest Value)', color='red', linestyle='--', where='mid', alpha=0.8)
            
            plt.xlabel('Wavelength (µm)')
            plt.ylabel('Flat Field Correction')
            plt.title(f'NIRSpec Flat Extension Strategy (Parlanti-style)\n{os.path.basename(fits_file)}')
            plt.grid(alpha=0.3)
            
            # Add yellow background for extended region
            xl = ext_range[0] * 1e6 if last_wl < 1e-4 else ext_range[0]
            plt.axvspan(xl, 5.6, color='yellow', alpha=0.1, label='Extended Region')
            plt.legend()
            
            output_plot = "/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/plots/flat_extension_strategy.png"
            plt.savefig(output_plot, dpi=300)
            print(f"Plot saved to: {output_plot}")
            plt.close()
        else:
            print("FLAT_TABLE not found in HDU list")

if __name__ == "__main__":
    plot_flat_extension(sflat_file)
