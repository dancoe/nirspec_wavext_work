import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

def plot_two_chip_flat(ref_nrs1, ref_nrs2, mode_label, extension_name, wl_plane_idx, output_path):
    nrs1_path = os.path.join("/Users/dcoe/crds_cache/references/jwst/nirspec", ref_nrs1)
    nrs2_path = os.path.join("/Users/dcoe/crds_cache/references/jwst/nirspec", ref_nrs2)
    
    if not os.path.exists(nrs1_path) or not os.path.exists(nrs2_path):
        return

    with fits.open(nrs1_path) as h1:
        d1_full = h1[extension_name].data
        wl_table = h1['WAVELENGTH'].data['wavelength'] if 'WAVELENGTH' in h1 else None
        
    with fits.open(nrs2_path) as h2:
        d2_full = h2[extension_name].data
        
    # Get the specific plane
    if len(d1_full.shape) == 3:
        # Check if index is valid
        idx = min(wl_plane_idx, d1_full.shape[0]-1)
        d1 = d1_full[idx]
        d2 = d2_full[idx] if len(d2_full.shape) == 3 else d2_full
        wl_val = wl_table[idx] if wl_table is not None else idx
    else:
        d1 = d1_full
        d2 = d2_full
        wl_val = None

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Use seismic diverging colormap, centered on 1.0
    vmin, vmax = 0.9, 1.1
    im1 = ax1.imshow(d1, origin='lower', cmap='seismic', vmin=vmin, vmax=vmax)
    ax1.set_title(f'{mode_label} NRS1 ({ref_nrs1})')
    ax1.set_xticks([]); ax1.set_yticks([])
    
    im2 = ax2.imshow(d2, origin='lower', cmap='seismic', vmin=vmin, vmax=vmax)
    ax2.set_title(f'{mode_label} NRS2 ({ref_nrs2})')
    ax2.set_xticks([]); ax2.set_yticks([])
    
    title = f'NIRSpec {mode_label}'
    if wl_val is not None:
        title += f' (Plane {wl_plane_idx}, $\lambda$ = {wl_val:.5f} µm)'
    
    plt.suptitle(title, fontsize=14)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    print(f"Saved: {output_path}")
    plt.close()

if __name__ == "__main__":
    out_dir = "/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/plots/flats"
    os.makedirs(out_dir, exist_ok=True)
    
    # Pairs to plot
    pairs = [
        ("jwst_nirspec_dflat_0001.fits", "jwst_nirspec_dflat_0002.fits", "D-flat"),
        ("jwst_nirspec_sflat_0147.fits", "jwst_nirspec_sflat_0154.fits", "FS G140M S-flat"),
        ("jwst_nirspec_sflat_0166.fits", "jwst_nirspec_sflat_0202.fits", "MOS PRISM S-flat"),
        ("jwst_nirspec_sflat_0208.fits", "jwst_nirspec_sflat_0191.fits", "IFU G140M F100LP S-flat"),
    ]
    
    for r1, r2, label in pairs:
        # Check if 3D
        path = os.path.join("/Users/dcoe/crds_cache/references/jwst/nirspec", r1)
        if not os.path.exists(path): continue
        with fits.open(path) as h:
            shape = h['SCI'].data.shape
            wl_table = h['WAVELENGTH'].data['wavelength'] if 'WAVELENGTH' in h else None
            
        if len(shape) == 3:
            # Generate plots for every plane
            n_planes = shape[0]
            for i in range(n_planes):
                wl_str = f"_{wl_table[i]:.5f}um" if wl_table is not None else f"_plane{i}"
                fname = f"{label.lower().replace(' ', '_')}{wl_str}.png"
                plot_two_chip_flat(r1, r2, label, "SCI", i, os.path.join(out_dir, fname))
        else:
            fname = f"{label.lower().replace(' ', '_')}.png"
            plot_two_chip_flat(r1, r2, label, "SCI", 0, os.path.join(out_dir, fname))
