
import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from astropy.wcs import WCS

# Paths
BASE_FS = '/Users/dcoe/NIRSpec/wavext/data/FS'
BASE_IFU = '/Users/dcoe/NIRSpec/wavext/data/IFU'
OUTPUT_DIR = '/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/reports/v9/ifu_v9/extractions_all_ifu_v9'
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Reference File Constants (from jwst_nirspec_extract1d_0002.asdf)
PIPE_R = 0.45
PIPE_BKG_IN = 1.0
PIPE_BKG_OUT = 1.2

# Tweak for v9: r = 0.5" for all
# 2654 is EXCLUDED due to nearby contamination
TARGETS = [
    {'pid': '1536', 'name': 'J1743045', 'desc': 'Quasar', 's3d': f'{BASE_IFU}/PID1536_J1743045/stage3/f170lp_g235m-f170lp_s3d.fits'},
    {'pid': '1537', 'name': 'G191-B2B', 'desc': 'DA White Dwarf', 's3d': f'{BASE_IFU}/PID1537_G191-B2B/stage3/f170lp_g235m-f170lp_s3d.fits'},
    {'pid': '1538', 'name': 'P330E', 'desc': 'G-type Solar Analog', 's3d': f'{BASE_IFU}/PID1538_P330E/stage3/f170lp_g235m-f170lp_s3d.fits'},
    {'pid': '2186', 'name': 'UGC-5101', 'desc': 'ULIRG', 's3d': f'{BASE_FS}/PID2186_UGC-5101/stage3_ext/g235m_f170lp_g235m-f170lp_s3d.fits', 'is_science': True},
    {'pid': '6645', 'name': 'P330E-C3', 'desc': 'Solar Analog (Offset Visit)', 's3d': f'{BASE_IFU}/PID6645_P330E-C3/stage3/f170lp_g235m-f170lp_s3d.fits'}
]

def perform_extractions(target):
    s3d_path = target['s3d']
    x1d_path = s3d_path.replace('_s3d.fits', '_x1d.fits')
    
    if not os.path.exists(s3d_path):
        print(f"Skipping {target['name']} (PID {target['pid']}): s3d not found.")
        return None

    with fits.open(s3d_path) as hdul:
        sci = hdul['SCI'].data
        hdr = hdul['SCI'].header
        prim_hdr = hdul[0].header
        
        # WCS
        w = WCS(hdr, naxis=2)
        ra_ref = prim_hdr.get('RA_TARG') or prim_hdr.get('RA_REF') or hdr.get('CRVAL1')
        dec_ref = prim_hdr.get('DEC_TARG') or prim_hdr.get('DEC_REF') or hdr.get('CRVAL2')
        pipe_x, pipe_y = w.all_world2pix(ra_ref, dec_ref, 0)
        
        wl = np.arange(hdr['NAXIS3']) * hdr['CDELT3'] + hdr['CRVAL3']
        bunit = hdr.get('BUNIT', 'MJy/sr')
        pixar_sr = hdr.get('PIXAR_SR', 1.0)
        srctype = hdr.get('SRCTYPE', 'UNKNOWN')
        
        # Mean image to find peak (Project choice)
        white_light = np.nanmean(sci, axis=0)
        peak_y, peak_x = np.unravel_index(np.nanargmax(white_light), white_light.shape)
        
        # Coordinates in arcsec relative to peak (Project v9: r=0.5")
        yy, xx = np.indices(white_light.shape)
        dist_arcsec_v9 = np.sqrt(((xx - peak_x)*hdr['CDELT1'])**2 + \
                                ((yy - peak_y)*hdr['CDELT2'])**2) * 3600
        
        # Flux conversion to Jy
        conv = 1.0
        if bunit == 'MJy/sr':
            conv = 1e6 * pixar_sr

        # 1. extract1d (from file)
        spec_x1d = np.zeros_like(wl)
        wl_old = wl
        if os.path.exists(x1d_path):
            with fits.open(x1d_path) as h_x:
                x1d_data = h_x['EXTRACT1D'].data
                wl_old = x1d_data['WAVELENGTH']
                spec_x1d = x1d_data['FLUX']

        # 2. r = 0.5" (no bkg, peak-centered) — The v9 Primary Choice
        mask_05 = dist_arcsec_v9 <= 0.5
        spec_05 = np.nansum(sci[:, mask_05], axis=1) * conv

        # 3. Reference FOV (Entire spatial footprint)
        spec_fov = np.nansum(sci, axis=(1, 2)) * conv
        
        results = {
            'wl': wl,
            'wl_old': wl_old,
            'spec_x1d': spec_x1d,
            'spec_05': spec_05,
            'spec_fov': spec_fov,
            'srctype': srctype,
            'peak_xy': (peak_x, peak_y),
            'pipe_xy': (pipe_x, pipe_y),
            'white_light': white_light,
            'sci': sci
        }
        
        return results

def plot_target(target, results):
    if results is None: return
    
    # --- Main Spectral Plot ---
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]}, figsize=(12, 10))
    
    wl = results['wl']
    ax1.plot(results['wl_old'], results['spec_x1d'], label='Pipeline extract1d', alpha=0.6, lw=1.5, color='black', linestyle='--')
    ax1.plot(wl, results['spec_05'], label='v9 standard (r=0.5", no bkg)', alpha=0.9, color='blue', lw=2)
    ax1.plot(wl, results['spec_fov'], label='Entire FOV', alpha=0.4, color='gray')
    
    ax1.set_yscale('log')
    ax1.set_title(f"PID {target['pid']} – {target['name']} ({target['desc']})\nIFU Validation Suite v9 — Standard 0.5\" Aperture")
    ax1.set_ylabel("Flux (Jy)")
    ax1.legend(loc='upper right', fontsize='small')
    ax1.grid(True, which='both', alpha=0.2)
    
    # Scaling logic
    mask_gap = ((wl >= 2.19) & (wl <= 2.23)) | ((wl >= 3.65) & (wl <= 3.80))
    valid_all = (np.concatenate([results['spec_05'], results['spec_x1d']]) > 0) & \
                np.isfinite(np.concatenate([results['spec_05'], results['spec_x1d']])) & \
                (~np.concatenate([mask_gap, mask_gap]))
    
    if np.any(valid_all):
        all_fluxes = np.concatenate([results['spec_05'], results['spec_x1d']])
        y_min = np.percentile(all_fluxes[valid_all], 1) * 0.7
        y_max = np.percentile(all_fluxes[valid_all], 99) * 2.0
        ax1.set_ylim(y_min, y_max)
    
    if target['pid'] == '6645':
        ax1.set_ylim(1e-5, 0.1)
    
    f_old = interp1d(results['wl_old'], results['spec_x1d'], bounds_error=False, fill_value=np.nan)
    ref_flux = f_old(wl)
    
    ratio_v9 = results['spec_05'] / ref_flux
    ax2.plot(wl, ratio_v9, label='Ratio v9 / extract1d', color='blue', lw=2)
    ax2.axhline(1.0, color='red', linestyle='--', alpha=0.5)
    
    ax2.set_yscale('log')
    ax2.set_ylabel("Ratio to extract1d")
    ax2.set_xlabel("Wavelength (um)")
    if target['pid'] == '6645':
        ax2.set_ylim(0.5, 200.0)
    else:
        ax2.set_ylim(0.5, 2.0)
    from matplotlib.ticker import ScalarFormatter
    ax2.yaxis.set_major_formatter(ScalarFormatter())
    if target['pid'] == '6645':
        ax2.set_yticks([0.5, 1, 10, 100, 200])
    else:
        ax2.set_yticks([0.5, 0.7, 1.0, 1.5, 2.0])
    ax2.grid(True, which='both', alpha=0.1)
    ax2.legend(loc='lower left', fontsize='x-small')
    
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/extraction_all_{target['pid']}.png", dpi=150)
    plt.close()

    # --- Slices Plot ---
    fig_slice, axes = plt.subplots(1, 3, figsize=(15, 5))
    nwl = len(wl)
    indices = [nwl//4, nwl//2, 3*nwl//4]
    
    peak_x, peak_y = results['peak_xy']
    
    for ax, idx in zip(axes, indices):
        slice_img = results['sci'][idx]
        v_min, v_max = np.nanpercentile(slice_img, [5, 99])
        im = ax.imshow(slice_img, origin='lower', cmap='inferno', vmin=v_min, vmax=v_max)
        ax.set_title(f"λ = {wl[idx]:.2f} um")
        
        # Draw v9 aperture (Red)
        circ_v9 = plt.Circle((peak_x, peak_y), 0.5 / 0.1, color='red', fill=False, lw=2, label='v9 (0.5")')
        ax.add_patch(circ_v9)
        ax.plot(peak_x, peak_y, 'r+', markersize=5)
        
    axes[0].legend(loc='upper left', fontsize='xx-small')
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/slices_{target['pid']}.png", dpi=150)
    plt.close()

def run_all():
    summary_data = []
    
    for target in TARGETS:
        print(f"Processing {target['name']} (PID {target['pid']})...")
        results = perform_extractions(target)
        if results:
            plot_target(target, results)
            mask_gap = ((results['wl'] >= 2.19) & (results['wl'] >= 2.23)) | \
                       ((results['wl'] >= 3.65) & (results['wl'] <= 3.80))
            valid = np.isfinite(results['spec_05']) & (results['spec_05'] > 0) & (~mask_gap)
            med_v9 = np.nanmedian(results['spec_05'][valid]) if np.any(valid) else 0
            
            f_old = interp1d(results['wl_old'], results['spec_x1d'], bounds_error=False, fill_value=np.nan)
            ref_flux = f_old(results['wl'])
            ratio_v9 = np.nanmedian(results['spec_05'][valid] / ref_flux[valid]) if np.any(valid) else 0
            
            summary_data.append({
                'pid': target['pid'],
                'name': target['name'],
                'desc': target['desc'],
                'srctype': results['srctype'],
                'med_flux_v9': med_v9,
                'ratio_v9_to_x1d': ratio_v9
            })
    
    # Generate Report
    report_path = f"{OUTPUT_DIR}/REPORT_all_extractions_v9.md"
    with open(report_path, "w") as f:
        f.write("# Total IFU Extraction Baseline Comparison — v9 (Uniform 0.5\")\n\n")
        f.write("In v9, all extractions are standardized to a 0.5\" fixed-radius aperture centered on the brightest pixel.\n\n")
        
        f.write("## Summary Statistics\n\n")
        f.write("| PID | Name | SRCTYPE | median Flux v9 (Jy) | median Ratio (v9/x1d) |\n")
        f.write("| :--- | :--- | :--- | :--- | :--- |\n")
        for s in summary_data:
            f.write(f"| {s['pid']} | {s['name']} | {s['srctype']} | {s['med_flux_v9']:.6f} | {s['ratio_v9_to_x1d']:.3f} |\n")
        f.write("\n\n")
        
        for target in TARGETS:
            f.write(f"## PID {target['pid']} – {target['name']} ({target['desc']})\n\n")
            f.write(f"### Slices & Apertures\n")
            f.write(f"**Red Solid**: v9 extraction footprint (peak-centered, 0.5\")\n\n")
            f.write(f"![Slices](slices_{target['pid']}.png)\n\n")
            f.write(f"### Spectra\n")
            f.write(f"![Comparison](extraction_all_{target['pid']}.png)\n\n")
            f.write("---\n\n")

if __name__ == "__main__":
    run_all()
