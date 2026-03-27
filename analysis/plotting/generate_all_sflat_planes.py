import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

def plot_two_chip_flat_final(ref_nrs1, ref_nrs2, mode_label, wl_plane_idx, subdirectory):
    nrs1_path = os.path.join("/Users/dcoe/crds_cache/references/jwst/nirspec", ref_nrs1)
    nrs2_path = os.path.join("/Users/dcoe/crds_cache/references/jwst/nirspec", ref_nrs2)
    
    if not os.path.exists(nrs1_path) and not os.path.exists(nrs2_path):
        return

    out_root = "/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/plots/flats"
    out_dir = os.path.join(out_root, subdirectory)
    os.makedirs(out_dir, exist_ok=True)

    d1 = None; d2 = None; wl_val = None

    if os.path.exists(nrs1_path):
        with fits.open(nrs1_path) as h1:
            d1_full = h1['SCI'].data
            wl_table = h1['WAVELENGTH'].data['wavelength'] if 'WAVELENGTH' in h1 else None
            if len(d1_full.shape) == 3:
                idx = min(wl_plane_idx, d1_full.shape[0]-1)
                d1 = d1_full[idx]
                wl_val = wl_table[idx] if wl_table is not None else None
            else:
                d1 = d1_full
    
    if os.path.exists(nrs2_path):
        with fits.open(nrs2_path) as h2:
            d2_full = h2['SCI'].data
            if len(d2_full.shape) == 3:
                idx = min(wl_plane_idx, d2_full.shape[0]-1)
                d2 = d2_full[idx]
            else:
                d2 = d2_full
                
    if d1 is None: d1 = np.ones_like(d2) if d2 is not None else np.zeros((2048,2048))
    if d2 is None: d2 = np.ones_like(d1)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    vmin, vmax = 0.9, 1.1
    im1 = ax1.imshow(d1, origin='lower', cmap='seismic', vmin=vmin, vmax=vmax)
    ax1.set_title(f'{mode_label} NRS1 ({ref_nrs1})')
    ax1.set_xlabel('x-column, px'); ax1.set_ylabel('y-row, px')
    
    im2 = ax2.imshow(d2, origin='lower', cmap='seismic', vmin=vmin, vmax=vmax)
    ax2.set_title(f'{mode_label} NRS2 ({ref_nrs2})')
    ax2.set_xlabel('x-column, px'); ax2.set_ylabel('y-row, px')
    
    fig.colorbar(im2, ax=ax2, label='Flat Correction (SCI)')
    title = f'NIRSpec {mode_label}'
    if wl_val is not None: title += f' ($\lambda$ = {wl_val:.5f} µm)'
    plt.suptitle(title, fontsize=14)
    plt.tight_layout()
    
    wl_str = f"_{wl_val:.5f}um" if wl_val is not None else ""
    basename = f"{mode_label.lower().replace(' ', '_')}{wl_str}.png"
    output_path = os.path.join(out_dir, basename)
    plt.savefig(output_path, dpi=200)
    plt.close()

if __name__ == "__main__":
    # S-flat pairs for MOS (02XX series)
    mos_pairs = [
        ("jwst_nirspec_sflat_0222.fits", "jwst_nirspec_sflat_0202.fits", "MOS PRISM S-flat"),
        ("jwst_nirspec_sflat_0221.fits", "jwst_nirspec_sflat_0228.fits", "MOS G140H S-flat"),
        ("jwst_nirspec_sflat_0231.fits", "jwst_nirspec_sflat_0226.fits", "MOS G140M S-flat"),
        ("jwst_nirspec_sflat_0229.fits", "jwst_nirspec_sflat_0224.fits", "MOS G235M S-flat"),
        ("jwst_nirspec_sflat_0232.fits", "jwst_nirspec_sflat_0225.fits", "MOS G395M S-flat"), # Swap check
    ]
    
    # FS and IFU
    other_pairs = [
        ("jwst_nirspec_sflat_0147.fits", "jwst_nirspec_sflat_0154.fits", "FS S-flat"),
        ("jwst_nirspec_sflat_0208.fits", "jwst_nirspec_sflat_0191.fits", "IFU G140M F100LP S-flat"),
    ]

    for r1, r2, label in mos_pairs:
        p1 = "/Users/dcoe/crds_cache/references/jwst/nirspec/" + r1
        if os.path.exists(p1):
            with fits.open(p1) as h:
                n_planes = h['SCI'].data.shape[0] if len(h['SCI'].data.shape) == 3 else 1
            for i in range(n_planes):
                # Put PRISM into subfolder if needed? User asked: "subdirectory for them plots/flats/f-flats/mos/"
                # I'll put PRISM into s-flats/mos for now but wait... user asked for f-flats/mos/
                subdir = "f-flats/mos" if "PRISM" in label else "s-flats"
                plot_two_chip_flat_final(r1, r2, label, i, subdir)

    for r1, r2, label in other_pairs:
        plot_two_chip_flat_final(r1, r2, label, 0, "s-flats")
