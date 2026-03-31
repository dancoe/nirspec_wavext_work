
import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Paths
BASE_DIR = '/Users/dcoe/NIRSpec/wavext/data/FS'
OUTPUT_REPORT_DIR = '/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/reports/v8/ifu_v8/extracted_ifu_v8'
os.makedirs(OUTPUT_REPORT_DIR, exist_ok=True)

TARGETS = [
    {
        'pid': '2186', 
        'name': 'UGC-5101', 
        'path': f'{BASE_DIR}/PID2186_UGC-5101/stage3_ext/g235m_f170lp_g235m-f170lp_s3d.fits',
        'x1d_path': f'{BASE_DIR}/PID2186_UGC-5101/stage3_ext/g235m_f170lp_g235m-f170lp_x1d.fits',
        'wavelengths': [2.0, 3.5, 5.0]
    },
    {
        'pid': '2654', 
        'name': 'SDSSJ0841', 
        'path': f'{BASE_DIR}/PID2654_SDSSJ0841/stage3_ext/g140m_f100lp_g140m-f100lp_s3d.fits',
        'x1d_path': f'{BASE_DIR}/PID2654_SDSSJ0841/stage3_ext/g140m_f100lp_g140m-f100lp_x1d.fits',
        'wavelengths': [1.0, 2.0, 3.0]
    }
]

def extract_circular(cube_path, radius_arcsec=0.5):
    with fits.open(cube_path) as hdul:
        sci = hdul['SCI'].data
        hdr = hdul['SCI'].header
        wl = np.arange(hdr['NAXIS3']) * hdr['CDELT3'] + hdr['CRVAL3']
        bunit = hdr.get('BUNIT', 'MJy/sr')
        pixar_sr = hdr.get('PIXAR_SR', 1.0)
        srctype = hdr.get('SRCTYPE', 'UNKNOWN')
        
        # Determine conversion to Jy
        conv = 1.0
        if bunit == 'MJy/sr':
            conv = 1e6 * pixar_sr
        
        # Mean image to find peak
        white_light = np.nanmean(sci, axis=0)
        peak_y, peak_x = np.unravel_index(np.nanargmax(white_light), white_light.shape)
        
        # Distance calculation
        cdelt1_arcsec = abs(hdr['CDELT1']) * 3600
        dist_arcsec = np.sqrt(((np.indices(white_light.shape)[1] - peak_x)*hdr['CDELT1'])**2 + \
                            ((np.indices(white_light.shape)[0] - peak_y)*hdr['CDELT2'])**2) * 3600
        
        mask = dist_arcsec <= radius_arcsec
        
        # Sum spectra in mask
        masked_sci = sci[:, mask]
        spectrum = np.nansum(masked_sci, axis=1) * conv
        
        return wl, spectrum, (peak_x, peak_y), srctype, sci, white_light, cdelt1_arcsec

def plot_vis(target):
    print(f"Visualizing {target['name']} (PID {target['pid']})...")
    wl, spec_v8, peak, srctype, sci, white_light, pix_scale = extract_circular(target['path'])
    
    # Load previous x1d for NPIXELS and metadata
    with fits.open(target['x1d_path']) as h:
        x1d_data = h['EXTRACT1D'].data
        wl_old = x1d_data['WAVELENGTH']
        spec_old = x1d_data['FLUX']
        npixels = x1d_data['NPIXELS']
        # Interpolate NPIXELS to find radius at specific wavelengths
        f_npix = interp1d(wl_old, npixels, bounds_error=False, fill_value='extrapolate')

    # 1. Plot slices
    fig_slices, axes = plt.subplots(1, len(target['wavelengths']), figsize=(15, 6))
    for ax, target_wl in zip(axes, target['wavelengths']):
        idx = np.argmin(np.abs(wl - target_wl))
        img = sci[idx]
        
        im = ax.imshow(img, origin='lower', cmap='viridis', interpolation='nearest', alpha=0.9)
        ax.set_title(f"{target_wl:.1f} um ({srctype})")
        
        # Draw v8 Aperture (Red)
        circ_v8 = plt.Circle(peak, 0.5/pix_scale, color='red', fill=False, lw=1.5, alpha=0.8, label='v8 (0.5")')
        ax.add_patch(circ_v8)
        
        # Draw Pipeline Aperture (Cyan dashed)
        if srctype == 'POINT':
            # Area = pi * r^2
            n_val = f_npix(target_wl)
            r_pix = np.sqrt(n_val / np.pi)
            circ_pipe = plt.Circle(peak, r_pix, color='cyan', linestyle='--', fill=False, lw=1.5, alpha=0.9, label='Pipeline')
            ax.add_patch(circ_pipe)
        elif srctype == 'EXTENDED':
            # Highlight full frame
            rect = plt.Rectangle((0, 0), img.shape[1]-1, img.shape[0]-1, color='cyan', linestyle='--', fill=False, lw=2, alpha=0.6, label='Pipeline (Whole Image)')
            ax.add_patch(rect)
        
        # Peak indicator
        ax.plot(peak[0], peak[1], 'r+', markersize=8, alpha=0.4)
        
        ax.set_xticks([])
        ax.set_yticks([])
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    fig_slices.suptitle(f"{target['name']} (PID {target['pid']}) - Spectral Slices\n[Red: v8 0.5\" r] [Cyan-Dashed: Default Pipeline ({srctype})]", fontsize=14)
    plt.tight_layout(rect=[0, 0.03, 1, 0.9])
    plt.savefig(f"{OUTPUT_REPORT_DIR}/slices_{target['pid']}.png", dpi=150)
    plt.close()

    # 2. Plot comparison from before
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]}, figsize=(10, 8))
    
    ax1.plot(wl_old, spec_old, label=f'Pipeline (SRCTYPE={srctype})', alpha=0.5, color='gray')
    ax1.plot(wl, spec_v8, label='v8 (0.5" radius)', color='blue')
    ax1.set_yscale('log')
    ax1.set_title(f"{target['name']} (PID {target['pid']}) - Extraction Comparison")
    ax1.set_ylabel("Flux (Jy)")
    ax1.legend()
    ax1.grid(True, which='both', alpha=0.2)
    
    # Ignore gaps
    mask_gap = ((wl >= 2.19) & (wl <= 2.23)) | \
               ((wl >= 3.65) & (wl <= 3.80))
    valid_v8 = (spec_v8 > 0) & np.isfinite(spec_v8) & (~mask_gap)
    
    if np.any(valid_v8):
        flux_vals = spec_v8[valid_v8]
        y_min = np.percentile(flux_vals, 1) * 0.5
        y_max = np.percentile(flux_vals, 99) * 2.0
        ax1.set_ylim(y_min, y_max)
    
    f_old = interp1d(wl_old, spec_old, bounds_error=False, fill_value=np.nan)
    spec_old_interp = f_old(wl)
    
    ratio = spec_v8 / spec_old_interp
    ax2.plot(wl, ratio, color='black', lw=1)
    ax2.axhline(1.0, color='red', linestyle='--', alpha=0.5)
    ax2.set_ylabel("Ratio (v8/old)")
    ax2.set_xlabel("Wavelength (um)")
    ax2.grid(True, which='both', alpha=0.2)
    
    valid_ratio = np.isfinite(ratio) & (~mask_gap)
    if np.any(valid_ratio):
        r_vals = ratio[valid_ratio]
        r_min = np.percentile(r_vals, 2) * 0.95
        r_max = np.percentile(r_vals, 98) * 1.05
        ax2.set_ylim(max(0, r_min), min(5, r_max))
    
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_REPORT_DIR}/extraction_{target['pid']}.png", dpi=150)
    plt.close()

def generate_report():
    report_path = f"{OUTPUT_REPORT_DIR}/REPORT_extracted_ifu_v8.md"
    with open(report_path, "w") as f:
        f.write("# IFU v8 Extracted Spectra Diagnostics\n\n")
        f.write("Comparison of the v8 0.5\" fixed circular aperture vs. the default pipeline extraction.\n\n")
        
        for target in TARGETS:
            f.write(f"## {target['name']} (PID {target['pid']})\n\n")
            f.write("### Spectral Slices and Extraction Regions\n")
            f.write("- **Red Solid Circle**: v8 extraction (r=0.5\")\n")
            f.write("- **Cyan Dashed**: Default Pipeline (Wavelength-dependent POINT or Whole-Image EXTENDED)\n\n")
            f.write(f"![{target['name']} Slices](slices_{target['pid']}.png)\n\n")
            f.write("### Spectrum and Ratio Comparison\n")
            f.write(f"![{target['name']} Extraction](extraction_{target['pid']}.png)\n\n")
            f.write("---\n\n")

if __name__ == "__main__":
    for target in TARGETS:
        plot_vis(target)
    generate_report()
    print("Done Visualization.")
