import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

def plot_dual_detector_image(ref_nrs1, ref_nrs2, mode_label, extension_name, output_path):
    nrs1_path = os.path.join("/Users/dcoe/crds_cache/references/jwst/nirspec", ref_nrs1)
    nrs2_path = os.path.join("/Users/dcoe/crds_cache/references/jwst/nirspec", ref_nrs2)
    
    if not os.path.exists(nrs1_path) or not os.path.exists(nrs2_path):
        print(f"Skipping {mode_label} - Missing pair: {ref_nrs1}/{ref_nrs2}")
        return

    with fits.open(nrs1_path) as hd1, fits.open(nrs2_path) as hd2:
        d1 = hd1[extension_name].data
        d2 = hd2[extension_name].data
        
        # Take first plane if 3D
        if len(d1.shape) == 3: d1 = d1[0]
        if len(d2.shape) == 3: d2 = d2[0]
            
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        # Scale 0.9 - 1.1 as requested
        im1 = ax1.imshow(d1, origin='lower', cmap='viridis', vmin=0.9, vmax=1.1)
        ax1.set_title(f'{mode_label} NRS1 ({ref_nrs1})')
        ax1.set_xticks([]); ax1.set_yticks([])
        
        im2 = ax2.imshow(d2, origin='lower', cmap='viridis', vmin=0.9, vmax=1.1)
        ax2.set_title(f'{mode_label} NRS2 ({ref_nrs2})')
        ax2.set_xticks([]); ax2.set_yticks([])
        
        fig.colorbar(im2, ax=ax2, label=f'{extension_name} Value')
        plt.suptitle(f'NIRSpec Dual-Detector {mode_label} (Scale 0.9 - 1.1)', fontsize=14)
        plt.tight_layout()
        plt.savefig(output_path, dpi=300)
        print(f"Saved: {output_path}")
        plt.close()

if __name__ == "__main__":
    out_dir = "/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/plots/flats"
    os.makedirs(out_dir, exist_ok=True)
    
    # D-flats (shows pixel-level noise/artifacts)
    plot_dual_detector_image("jwst_nirspec_dflat_0001.fits", "jwst_nirspec_dflat_0002.fits", "D-flat", "SCI", os.path.join(out_dir, "dual_detector_dflat_09_11.png"))
    
    # FS S-flats (shows slit footprints)
    plot_dual_detector_image("jwst_nirspec_sflat_0147.fits", "jwst_nirspec_sflat_0154.fits", "FS G140M S-flat", "SCI", os.path.join(out_dir, "dual_detector_fs_sflat_09_11.png"))
    
    # MOS PRISM S-flats (Checking if they show anything in this scale)
    plot_dual_detector_image("jwst_nirspec_sflat_0166.fits", "jwst_nirspec_sflat_0202.fits", "MOS PRISM S-flat", "SCI", os.path.join(out_dir, "dual_detector_mos_sflat_09_11.png"))
