import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

def plot_two_chip_flat_final(ref_nrs1, ref_nrs2, mode_label, wl_plane_idx, subdirectory):
    nrs1_path = os.path.join("/Users/dcoe/crds_cache/references/jwst/nirspec", ref_nrs1)
    nrs2_path = os.path.join("/Users/dcoe/crds_cache/references/jwst/nirspec", ref_nrs2)
    
    if not os.path.exists(nrs1_path) or not os.path.exists(nrs2_path):
        print(f"Skipping {mode_label} - Pair not available locally: {ref_nrs1}/{ref_nrs2}")
        return

    out_root = "/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/plots/flats"
    out_dir = os.path.join(out_root, subdirectory)
    os.makedirs(out_dir, exist_ok=True)

    with fits.open(nrs1_path) as h1:
        d1_full = h1['SCI'].data
        wl_table = h1['WAVELENGTH'].data['wavelength'] if 'WAVELENGTH' in h1 else None
        
    with fits.open(nrs2_path) as h2:
        d2_full = h2['SCI'].data
        
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
    
    # Diverging seismic colormap, 0.9 - 1.1 scale
    vmin, vmax = 0.9, 1.1
    im1 = ax1.imshow(d1, origin='lower', cmap='seismic', vmin=vmin, vmax=vmax)
    ax1.set_title(f'{mode_label} NRS1 ({ref_nrs1})')
    ax1.set_xlabel('x-column, px'); ax1.set_ylabel('y-row, px')
    
    im2 = ax2.imshow(d2, origin='lower', cmap='seismic', vmin=vmin, vmax=vmax)
    ax2.set_title(f'{mode_label} NRS2 ({ref_nrs2})')
    ax2.set_xlabel('x-column, px'); ax2.set_ylabel('y-row, px')
    
    fig.colorbar(im2, ax=ax2, label='Flat Correction (SCI)')
    
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
    # Updated pairs based on CRDS check
    # Many MOS pairs use 02XX series files
    pairs = [
        ("jwst_nirspec_dflat_0001.fits", "jwst_nirspec_dflat_0002.fits", "D-flat", "d-flats"),
        ("jwst_nirspec_sflat_0222.fits", "jwst_nirspec_sflat_0222.fits", "MOS PRISM S-flat", "s-flats"), # Temp partner check
        ("jwst_nirspec_sflat_0147.fits", "jwst_nirspec_sflat_0154.fits", "FS S-flat", "s-flats"),
        ("jwst_nirspec_sflat_0208.fits", "jwst_nirspec_sflat_0191.fits", "IFU G140M F100LP S-flat", "s-flats"),
    ]
    
    for r1, r2, label, subdir in pairs:
        path1 = os.path.join("/Users/dcoe/crds_cache/references/jwst/nirspec", r1)
        if not os.path.exists(path1): continue
        with fits.open(path1) as h:
            shape = h['SCI'].data.shape
        
        if len(shape) == 3:
            # For 3D files, plot a few planes
            indices = [min(0, shape[0]-1), shape[0]//2, shape[0]-1]
            for i in indices:
                plot_two_chip_flat_final(r1, r2, label, i, subdir)
        else:
            plot_two_chip_flat_final(r1, r2, label, 0, subdir)
