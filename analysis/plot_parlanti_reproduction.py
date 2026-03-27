import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

def create_parlanti_reproduction(ref_file, mode_label, output_path):
    fits_path = os.path.join("/Users/dcoe/crds_cache/references/jwst/nirspec", ref_file)
    if not os.path.exists(fits_path):
        print(f"File {ref_file} not found in cache.")
        return

    with fits.open(fits_path) as hdul:
        # 1. Fast Variation data
        fv_data = hdul['FAST_VARIATION'].data
        nelem = fv_data['nelem'][0]
        wl = fv_data['wavelength'][0][:nelem]
        flat = fv_data['data'][0][:nelem]
        
        # Extension
        last_wl = wl[-1]
        ext_wl = np.linspace(last_wl, 5.61e-6, 50)[1:]
        ext_flat = np.full_like(ext_wl, flat[-1])
        
        # 2. Slow Variation conceptual data
        # Nominal slow variation planes (if multiple)
        if 'WAVELENGTH' in hdul:
            nom_slow_wl = hdul['WAVELENGTH'].data['wavelength']
        else:
            # If 2D, nominal is the whole range
            nom_slow_wl = [wl[0], wl[-1]]
            
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
        
        # Top Panel: Fast Variation
        ax1.step(wl * 1e6, flat, color='blue', label='Nominal Fast Variation $f(\lambda)$')
        ax1.step(ext_wl * 1e6, ext_flat, color='red', linestyle='--', label='wavext extension (const)')
        ax1.axvspan(ext_wl[0]*1e6, 5.61, color='orange', alpha=0.1)
        ax1.set_ylabel('Fast Variation $f(\lambda)$')
        ax1.set_title(f'NIRSpec {mode_label} S-flat Extension Strategy (Parlanti Modeled)\n{ref_file}')
        ax1.legend(loc='lower left')
        ax1.grid(alpha=0.3)

        # Bottom Panel: Slow Variation
        # Show nominal slow variation as a "range" or lines
        ax2.hlines(1.0, wl[0]*1e6, last_wl*1e6, color='blue', linewidth=2, label='Nominal Slow Variation $S(x, y, \lambda)$')
        ax2.hlines(1.0, last_wl*1e6, 5.61, color='red', linestyle='--', linewidth=2, label='wavext extension $S(x, y, \lambda) = 1.0$')
        ax2.axvspan(ext_wl[0]*1e6, 5.61, color='orange', alpha=0.1)
        ax2.set_xlabel('Wavelength (µm)')
        ax2.set_ylabel('Slow Variation $S(x, y, \lambda)$')
        ax2.set_ylim(0.5, 1.5)
        ax2.legend(loc='lower left')
        ax2.grid(alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=300)
        print(f"Saved: {output_path}")
        plt.close()

if __name__ == "__main__":
    out_dir = "/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/plots/flats"
    os.makedirs(out_dir, exist_ok=True)
    
    # Generate for IFU (as in Parlanti A1)
    create_parlanti_reproduction("jwst_nirspec_sflat_0191.fits", "IFU G140M", os.path.join(out_dir, "parlanti_repro_ifu.png"))
    
    # Generate for MOS PRISM
    create_parlanti_reproduction("jwst_nirspec_sflat_0202.fits", "MOS PRISM", os.path.join(out_dir, "parlanti_repro_mos_prism.png"))
