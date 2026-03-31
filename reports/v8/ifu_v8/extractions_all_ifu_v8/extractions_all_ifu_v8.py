
import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Paths
BASE_FS = '/Users/dcoe/NIRSpec/wavext/data/FS'
BASE_IFU = '/Users/dcoe/NIRSpec/wavext/data/IFU'
OUTPUT_DIR = '/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/reports/v8/ifu_v8/extractions_all_ifu_v8'
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Reference File Constants (from jwst_nirspec_extract1d_0002.asdf)
PIPE_R = 0.45
PIPE_BKG_IN = 1.0
PIPE_BKG_OUT = 1.2

TARGETS = [
    {'pid': '1536', 'name': 'J1743045', 's3d': f'{BASE_IFU}/PID1536_J1743045/stage3/f170lp_g235m-f170lp_s3d.fits'},
    {'pid': '1537', 'name': 'G191-B2B', 's3d': f'{BASE_IFU}/PID1537_G191-B2B/stage3/f170lp_g235m-f170lp_s3d.fits'},
    {'pid': '1538', 'name': 'P330E', 's3d': f'{BASE_IFU}/PID1538_P330E/stage3/f170lp_g235m-f170lp_s3d.fits'},
    {'pid': '2186', 'name': 'UGC-5101', 's3d': f'{BASE_FS}/PID2186_UGC-5101/stage3_ext/g235m_f170lp_g235m-f170lp_s3d.fits'},
    {'pid': '2654', 'name': 'SDSSJ0841', 's3d': f'{BASE_FS}/PID2654_SDSSJ0841/stage3_ext/g140m_f100lp_g140m-f100lp_s3d.fits'},
    {'pid': '6645', 'name': 'P330E-C3', 's3d': f'{BASE_IFU}/PID6645_P330E-C3/stage3/f170lp_g235m-f170lp_s3d.fits'}
]

def perform_extractions(target):
    s3d_path = target['s3d']
    x1d_path = s3d_path.replace('_s3d.fits', '_x1d.fits')
    
    if not os.path.exists(s3d_path) or not os.path.exists(x1d_path):
        print(f"Skipping {target['name']} (PID {target['pid']}): missing files.")
        return None

    with fits.open(s3d_path) as hdul:
        sci = hdul['SCI'].data
        hdr = hdul['SCI'].header
        wl = np.arange(hdr['NAXIS3']) * hdr['CDELT3'] + hdr['CRVAL3']
        bunit = hdr.get('BUNIT', 'MJy/sr')
        pixar_sr = hdr.get('PIXAR_SR', 1.0)
        srctype = hdr.get('SRCTYPE', 'UNKNOWN')
        
        # Mean image to find peak
        white_light = np.nanmean(sci, axis=0)
        peak_y, peak_x = np.unravel_index(np.nanargmax(white_light), white_light.shape)
        
        # Coordinates in arcsec relative to peak
        yy, xx = np.indices(white_light.shape)
        dist_arcsec = np.sqrt(((xx - peak_x)*hdr['CDELT1'])**2 + \
                            ((yy - peak_y)*hdr['CDELT2'])**2) * 3600
        
        # Flux conversion to Jy
        conv = 1.0
        if bunit == 'MJy/sr':
            conv = 1e6 * pixar_sr

        # 1. extract1d (from file)
        with fits.open(x1d_path) as h_x:
            x1d_data = h_x['EXTRACT1D'].data
            wl_old = x1d_data['WAVELENGTH']
            spec_x1d = x1d_data['FLUX']

        # 2. r = 0.5" (no bkg)
        mask_05 = dist_arcsec <= 0.5
        spec_05 = np.nansum(sci[:, mask_05], axis=1) * conv

        # 3. r = 0.45" (with annulus if POINT)
        mask_045 = dist_arcsec <= 0.45
        spec_045 = np.nansum(sci[:, mask_045], axis=1) * conv
        
        if srctype == 'POINT' or srctype == 'UNKNOWN':
            mask_bkg = (dist_arcsec >= PIPE_BKG_IN) & (dist_arcsec <= PIPE_BKG_OUT)
            if np.any(mask_bkg):
                # Median background per pixel
                # sci shape (NWL, NY, NX) -> slice (NY, NX)
                # We need average background per wavelength slice
                masked_bkg = sci[:, mask_bkg]
                bkg_level = np.nanmedian(masked_bkg, axis=1) # (NWL)
                n_pix_aper = np.sum(mask_045)
                # Subtract bkg scaled to aperture area
                spec_045 -= bkg_level * n_pix_aper * conv
        
        # 4. Entire FOV
        spec_fov = np.nansum(sci, axis=(1, 2)) * conv
        
        return {
            'wl': wl,
            'wl_old': wl_old,
            'spec_x1d': spec_x1d,
            'spec_05': spec_05,
            'spec_045': spec_045,
            'spec_fov': spec_fov,
            'srctype': srctype
        }

def plot_target(target, results):
    if results is None: return
    
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]}, figsize=(12, 10))
    
    wl = results['wl']
    ax1.plot(results['wl_old'], results['spec_x1d'], label='Pipeline extract1d', alpha=0.7, lw=2, color='black')
    ax1.plot(wl, results['spec_05'], label='v8 (r=0.5", no bkg)', alpha=0.8, color='blue')
    ax1.plot(wl, results['spec_045'], label='r=0.45" (bkg sub if POINT)', alpha=0.8, color='cyan', linestyle='--')
    ax1.plot(wl, results['spec_fov'], label='Entire FOV', alpha=0.5, color='gray')
    
    ax1.set_yscale('log')
    ax1.set_title(f"{target['name']} (PID {target['pid']}) - Extraction Baseline Comparison ({results['srctype']})")
    ax1.set_ylabel("Flux (Jy)")
    ax1.legend(loc='upper right', fontsize='small')
    ax1.grid(True, which='both', alpha=0.2)
    
    # Range handling (ignore gaps)
    mask_gap = ((wl >= 2.19) & (wl <= 2.23)) | ((wl >= 3.65) & (wl <= 3.80))
    valid = (results['spec_05'] > 0) & np.isfinite(results['spec_05']) & (~mask_gap)
    if np.any(valid):
        y_min = np.percentile(results['spec_05'][valid], 1) * 0.5
        y_max = np.percentile(results['spec_05'][valid], 99) * 5.0
        ax1.set_ylim(y_min, y_max)
    
    # Ratio relative to extract1d
    f_old = interp1d(results['wl_old'], results['spec_x1d'], bounds_error=False, fill_value=np.nan)
    ref_flux = f_old(wl)
    
    ax2.plot(wl, results['spec_05'] / ref_flux, label='v8 / x1d', color='blue')
    ax2.plot(wl, results['spec_045'] / ref_flux, label='r=0.45 / x1d', color='cyan', linestyle='--')
    ax2.axhline(1.0, color='red', linestyle='--', alpha=0.5)
    
    ax2.set_ylabel("Ratio to extract1d")
    ax2.set_xlabel("Wavelength (um)")
    ax2.set_ylim(0.5, 3.0)
    ax2.grid(True, which='both', alpha=0.2)
    ax2.legend(loc='upper right', ncol=2, fontsize='x-small')
    
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/extraction_all_{target['pid']}.png", dpi=150)
    plt.close()

def run_all():
    summary_data = []
    
    for target in TARGETS:
        print(f"Processing {target['name']} (PID {target['pid']})...")
        results = perform_extractions(target)
        if results:
            plot_target(target, results)
            
            # Simple stats for report
            mask_gap = ((results['wl'] >= 2.19) & (results['wl'] >= 2.23)) | \
                       ((results['wl'] >= 3.65) & (results['wl'] <= 3.80))
            valid = np.isfinite(results['spec_05']) & (results['spec_05'] > 0) & (~mask_gap)
            med_05 = np.nanmedian(results['spec_05'][valid]) if np.any(valid) else 0
            
            f_old = interp1d(results['wl_old'], results['spec_x1d'], bounds_error=False, fill_value=np.nan)
            ref_flux = f_old(results['wl'])
            ratio_05 = np.nanmedian(results['spec_05'][valid] / ref_flux[valid]) if np.any(valid) else 0
            
            summary_data.append({
                'pid': target['pid'],
                'name': target['name'],
                'srctype': results['srctype'],
                'med_flux_v8': med_05,
                'ratio_v8_to_x1d': ratio_05
            })
    
    # Generate Report
    report_path = f"{OUTPUT_DIR}/REPORT_all_extractions.md"
    with open(report_path, "w") as f:
        f.write("# Total IFU Extraction Baseline Comparison\n\n")
        f.write("Comparison of 4 extraction methods across all v8 catalog datasets:\n")
        f.write("1. **extract1d**: Official pipeline output (includes aperture correction for POINT).\n")
        f.write("2. **v8 (r=0.5\")**: Project custom fixed circular aperture, no background subtraction, no aperture correction.\n")
        f.write("3. **r=0.45\"**: Reference fixed aperture with background annulus subtraction (1.0-1.2\") if SRCTYPE=POINT.\n")
        f.write("4. **Entire FOV**: Summation of the entire spatial footprint of the cube.\n\n")
        
        f.write("## Summary Statistics\n\n")
        f.write("| PID | Name | SRCTYPE | median Flux v8 (Jy) | median Ratio (v8/x1d) |\n")
        f.write("| :--- | :--- | :--- | :--- | :--- |\n")
        for s in summary_data:
            f.write(f"| {s['pid']} | {s['name']} | {s['srctype']} | {s['med_flux_v8']:.6f} | {s['ratio_v8_to_x1d']:.3f} |\n")
        f.write("\n\n")
        
        for target in TARGETS:
            f.write(f"## {target['name']} (PID {target['pid']})\n\n")
            f.write(f"![Comparison](extraction_all_{target['pid']}.png)\n\n")
            f.write("---\n\n")

if __name__ == "__main__":
    run_all()
