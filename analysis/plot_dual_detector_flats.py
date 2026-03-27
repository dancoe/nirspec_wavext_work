import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

def plot_dual_chip_flats(ref_nrs1, ref_nrs2, flat_type, output_path):
    nrs1_path = os.path.join("/Users/dcoe/crds_cache/references/jwst/nirspec", ref_nrs1)
    nrs2_path = os.path.join("/Users/dcoe/crds_cache/references/jwst/nirspec", ref_nrs2)
    
    if not os.path.exists(nrs1_path) or not os.path.exists(nrs2_path):
        print(f"Partner files {ref_nrs1}/{ref_nrs2} not found.")
        return

    with fits.open(nrs1_path) as h1, fits.open(nrs2_path) as h2:
        # Get first SCI extension data
        d1 = h1['SCI'].data
        d2 = h2['SCI'].data
        
        # If 3D, take first plane
        if len(d1.shape) == 3:
            d1 = d1[0]
            d2 = d2[0]

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
        
        # D-flats have pixel variation, S-flats might be mostly 1.0 but can show edges
        cmap = 'gray' if 'D-flat' in flat_type else 'viridis'
        vmin, vmax = (0.9, 1.1) if 'D-flat' in flat_type else (0.98, 1.02)
        
        im1 = ax1.imshow(d1, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax)
        ax1.set_title(f'{flat_type} NRS1\n({ref_nrs1})')
        ax1.set_xticks([]); ax1.set_yticks([])
        
        im2 = ax2.imshow(d2, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax)
        ax2.set_title(f'{flat_type} NRS2\n({ref_nrs2})')
        ax2.set_xticks([]); ax2.set_yticks([])
        
        # Add colorbar for context
        fig.colorbar(im2, ax=ax2, label='Flat Correction')
        
        plt.suptitle(f'NIRSpec Dual-Detector {flat_type} (2D Spatial Variation)', fontsize=14)
        plt.tight_layout()
        plt.savefig(output_path, dpi=300)
        print(f"Saved: {output_path}")
        plt.close()

if __name__ == "__main__":
    out_dir = "/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/plots/flats"
    os.makedirs(out_dir, exist_ok=True)
    
    # 1. Dual Detector D-flat (Pixel Variation)
    plot_dual_chip_flats("jwst_nirspec_dflat_0001.fits", "jwst_nirspec_dflat_0002.fits", "D-flat", os.path.join(out_dir, "dual_detector_dflat.png"))
    
    # 2. Dual Detector S-flat (IFU G140M F100LP - available in cache)
    plot_dual_chip_flats("jwst_nirspec_sflat_0208.fits", "jwst_nirspec_sflat_0191.fits", "IFU G140M F100LP S-flat", os.path.join(out_dir, "dual_detector_ifu_sflat.png"))

    # 3. Dual Detector S-flat (MOS PRISM)
    plot_dual_chip_flats("jwst_nirspec_sflat_0166.fits", "jwst_nirspec_sflat_0202.fits", "MOS PRISM S-flat", os.path.join(out_dir, "dual_detector_mos_sflat.png"))

    # 4. Dual Detector S-flat (FS G140M)
    plot_dual_chip_flats("jwst_nirspec_sflat_0147.fits", "jwst_nirspec_sflat_0154.fits", "FS G140M S-flat", os.path.join(out_dir, "dual_detector_fs_sflat.png"))
