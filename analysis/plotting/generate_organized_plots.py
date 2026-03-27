import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

def plot_two_chip_flat_organized(ref_nrs1, ref_nrs2, mode_label, extension_name, wl_plane_idx, subdirectory):
    nrs1_path = os.path.join("/Users/dcoe/crds_cache/references/jwst/nirspec", ref_nrs1)
    nrs2_path = os.path.join("/Users/dcoe/crds_cache/references/jwst/nirspec", ref_nrs2)
    
    if not os.path.exists(nrs1_path) or not os.path.exists(nrs2_path):
        return

    out_root = "/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/plots/flats"
    out_dir = os.path.join(out_root, subdirectory)
    os.makedirs(out_dir, exist_ok=True)

    with fits.open(nrs1_path) as h1:
        d1_full = h1[extension_name].data
        if d1_full is None or d1_full.size == 0:
            print(f"Skipping {ref_nrs1} - Empty {extension_name}")
            return
        wl_table = h1['WAVELENGTH'].data['wavelength'] if 'WAVELENGTH' in h1 else None
        
    with fits.open(nrs2_path) as h2:
        d2_full = h2[extension_name].data
        if d2_full is None or d2_full.size == 0:
            d2_full = np.ones_like(d1_full) # Handle case where partner is empty
        
    # Get the specific plane
    if len(d1_full.shape) == 3:
        idx = min(wl_plane_idx, d1_full.shape[0]-1)
        d1 = d1_full[idx]
        d2 = d2_full[idx] if len(d2_full.shape) == 3 else d2_full
        wl_val = wl_table[idx] if wl_table is not None else None
    else:
        d1 = d1_full
        d2 = d2_full
        wl_val = None

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Use seismic diverging colormap, centered on 1.0
    vmin, vmax = 0.9, 1.1
    im1 = ax1.imshow(d1, origin='lower', cmap='seismic', vmin=vmin, vmax=vmax)
    ax1.set_title(f'{mode_label} NRS1 ({ref_nrs1})')
    ax1.set_xlabel('x-column, px'); ax1.set_ylabel('y-row, px')
    
    im2 = ax2.imshow(d2, origin='lower', cmap='seismic', vmin=vmin, vmax=vmax)
    ax2.set_title(f'{mode_label} NRS2 ({ref_nrs2})')
    ax2.set_xlabel('x-column, px'); ax2.set_ylabel('y-row, px')
    
    fig.colorbar(im2, ax=ax2, label=f'Correction ({extension_name})')
    
    title = f'NIRSpec {mode_label}'
    if wl_val is not None:
        title += f' ($\lambda$ = {wl_val:.5f} µm)'
    
    plt.suptitle(title, fontsize=14)
    plt.tight_layout()
    
    wl_str = f"_{wl_val:.5f}um" if wl_val is not None else ""
    basename = f"{mode_label.lower().replace(' ', '_')}{wl_str}.png"
    output_path = os.path.join(out_dir, basename)
    
    plt.savefig(output_path, dpi=200)
    print(f"Saved: {output_path}")
    plt.close()

if __name__ == "__main__":
    # Pairs for different modes
    sflat_pairs = [
        ("jwst_nirspec_sflat_0147.fits", "jwst_nirspec_sflat_0154.fits", "FS S-flat"),
        ("jwst_nirspec_sflat_0166.fits", "jwst_nirspec_sflat_0202.fits", "MOS PRISM S-flat"),
        ("jwst_nirspec_sflat_0208.fits", "jwst_nirspec_sflat_0191.fits", "IFU G140M F100LP S-flat"),
    ]
    
    dflat_pairs = [
        ("jwst_nirspec_dflat_0001.fits", "jwst_nirspec_dflat_0002.fits", "D-flat")
    ]
    
    fflat_pairs = [
        ("jwst_nirspec_fflat_0152.fits", "jwst_nirspec_fflat_0154.fits", "MOS PRISM F-flat"),
    ]

    for r1, r2, label in dflat_pairs:
        indices = [2, 11, 38] # ~0.8um, ~2.9um, ~5.3um
        for i in indices:
            plot_two_chip_flat_organized(r1, r2, label, "SCI", i, "d-flats")
            
    for r1, r2, label in sflat_pairs:
        path = os.path.join("/Users/dcoe/crds_cache/references/jwst/nirspec", r1)
        with fits.open(path) as h:
            shape = h['SCI'].data.shape
        if len(shape) == 3:
            for i in range(shape[0]):
                plot_two_chip_flat_organized(r1, r2, label, "SCI", i, "s-flats")
        else:
            plot_two_chip_flat_organized(r1, r2, label, "SCI", 0, "s-flats")
            
    for r1, r2, label in fflat_pairs:
        path = os.path.join("/Users/dcoe/crds_cache/references/jwst/nirspec", r1)
        with fits.open(path) as h:
            data = h['SCI'].data
            if data is not None: shape = data.shape
            else: shape = (0,)
        if len(shape) == 3 and shape[0] > 0:
            for i in range(shape[0]):
                plot_two_chip_flat_organized(r1, r2, label, "SCI", i, "f-flats")
        elif len(shape) >= 2:
            plot_two_chip_flat_organized(r1, r2, label, "SCI", 0, "f-flats")
        else:
            print(f"Skipping {label} F-flat - Empty image")
