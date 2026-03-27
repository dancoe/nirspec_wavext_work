import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

def plot_3d_sflat_planes(ref_file, mode_label, output_template):
    fits_path = os.path.join("/Users/dcoe/crds_cache/references/jwst/nirspec", ref_file)
    if not os.path.exists(fits_path):
        print(f"File {ref_file} not found.")
        return

    with fits.open(fits_path) as hdul:
        sci = hdul['SCI'].data
        if len(sci.shape) != 3:
            print(f"File {ref_file} is not 3D (shape {sci.shape}). Skipping multi-plane plot.")
            return
            
        n_planes = sci.shape[0]
        # Get wavelengths if available
        if 'WAVELENGTH' in hdul:
            wl = hdul['WAVELENGTH'].data['wavelength']
        else:
            wl = np.arange(n_planes)
            
        print(f"Plotting ALL {n_planes} planes for {ref_file}...")
        
        for idx in range(n_planes):
            plane_data = sci[idx]
            wl_um = wl[idx]
            # Handle meters to microns
            if wl_um < 1e-4: wl_um *= 1e6
            
            plt.figure(figsize=(10, 8))
            # Use 0.9 - 1.1 scale as requested
            plt.imshow(plane_data, origin='lower', cmap='viridis', vmin=0.9, vmax=1.1)
            plt.colorbar(label='S-flat correction')
            plt.title(f'{mode_label} S-flat Spatial Sensitivity\n{ref_file} - Plane {idx} ($\lambda$ = {wl_um:.3f} µm)')
            plt.xlabel('X Pixel')
            plt.ylabel('Y Pixel')
            
            basename = f"{mode_label.lower().replace(' ', '_')}_sflat_plane{idx:02d}.png"
            output_path = output_template.replace("TEMPLATE", basename)
            plt.savefig(output_path, dpi=300)
            print(f"Saved: {output_path}")
            plt.close()

if __name__ == "__main__":
    out_dir = "/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/plots/flats"
    os.makedirs(out_dir, exist_ok=True)
    
    # MOS PRISM NRS1 (0166)
    plot_3d_sflat_planes("jwst_nirspec_sflat_0166.fits", "MOS PRISM NRS1", os.path.join(out_dir, "TEMPLATE"))
    
    # MOS PRISM NRS2 (0202)
    plot_3d_sflat_planes("jwst_nirspec_sflat_0202.fits", "MOS PRISM NRS2", os.path.join(out_dir, "TEMPLATE"))
    
    # G395M NRS1 (0177)
    plot_3d_sflat_planes("jwst_nirspec_sflat_0177.fits", "MOS G395M NRS1", os.path.join(out_dir, "TEMPLATE"))
    
    # IFU NRS1 (0190)? check if 3D
    plot_3d_sflat_planes("jwst_nirspec_sflat_0190.fits", "IFU G140M NRS1", os.path.join(out_dir, "TEMPLATE"))
