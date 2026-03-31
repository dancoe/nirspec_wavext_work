
import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

# Paths
BASE_DIR = '/Users/dcoe/NIRSpec/wavext/data/FS'
OUTPUT_REPORT_DIR = '/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/reports/v8/ifu_v8'
os.makedirs(OUTPUT_REPORT_DIR, exist_ok=True)

TARGETS = [
    {
        'pid': '2186', 
        'name': 'UGC-5101', 
        'path': f'{BASE_DIR}/PID2186_UGC-5101/stage3_ext/g235m_f170lp_g235m-f170lp_s3d.fits',
        'x1d_path': f'{BASE_DIR}/PID2186_UGC-5101/stage3_ext/g235m_f170lp_g235m-f170lp_x1d.fits'
    },
    {
        'pid': '2654', 
        'name': 'SDSSJ0841', 
        'path': f'{BASE_DIR}/PID2654_SDSSJ0841/stage3_ext/g140m_f100lp_g140m-f100lp_s3d.fits',
        'x1d_path': f'{BASE_DIR}/PID2654_SDSSJ0841/stage3_ext/g140m_f100lp_g140m-f100lp_x1d.fits'
    }
]

def extract_circular(cube_path, radius_arcsec=0.5):
    with fits.open(cube_path) as hdul:
        sci = hdul['SCI'].data
        hdr = hdul['SCI'].header
        wl = np.arange(hdr['NAXIS3']) * hdr['CDELT3'] + hdr['CRVAL3']
        bunit = hdr.get('BUNIT', 'MJy/sr')
        pixar_sr = hdr.get('PIXAR_SR', 1.0)
        
        # Determine conversion to Jy
        conv = 1.0
        if bunit == 'MJy/sr':
            conv = 1e6 * pixar_sr
        
        # Mean image to find peak
        white_light = np.nanmean(sci, axis=0)
        peak_y, peak_x = np.unravel_index(np.nanargmax(white_light), white_light.shape)
        
        # Grid of pixel distances
        yy, xx = np.indices(white_light.shape)
        
        # Distance calculation
        # Simplified since cdelt1/2 are roughly equal
        cdelt1_arcsec = abs(hdr['CDELT1']) * 3600
        dist_arcsec = np.sqrt(((xx - peak_x)*hdr['CDELT1'])**2 + ((yy - peak_y)*hdr['CDELT2'])**2) * 3600
        
        mask = dist_arcsec <= radius_arcsec
        
        # Sum spectra in mask
        masked_sci = sci[:, mask]
        spectrum = np.nansum(masked_sci, axis=1) * conv
        
        return wl, spectrum, (peak_x, peak_y), mask

def run_extraction():
    results = {}
    for target in TARGETS:
        print(f"Processing {target['name']} (PID {target['pid']})...")
        wl, spec_v8, peak, mask = extract_circular(target['path'])
        
        # Load previous x1d
        with fits.open(target['x1d_path']) as h:
            x1d_data = h['EXTRACT1D'].data
            wl_old = x1d_data['WAVELENGTH']
            spec_old = x1d_data['FLUX']
            
        results[target['pid']] = {
            'wl': wl,
            'spec_v8': spec_v8,
            'wl_old': wl_old,
            'spec_old': spec_old,
            'peak': peak,
            'name': target['name']
        }
        
    return results

def plot_comparison(results):
    for pid, data in results.items():
        plt.figure(figsize=(10, 6))
        plt.plot(data['wl_old'], data['spec_old'], label='Previous extract1d', alpha=0.5, color='gray')
        plt.plot(data['wl'], data['spec_v8'], label='v8 (0.5" radius)', color='blue')
        plt.yscale('log')
        plt.title(f"{data['name']} (PID {pid}) - Extraction Comparison")
        plt.xlabel("Wavelength (um)")
        plt.ylabel("Flux (Jy)")
        plt.legend()
        plt.grid(True, which='both', alpha=0.2)
        plt.savefig(f"{OUTPUT_REPORT_DIR}/extraction_comparison_{pid}.png", dpi=150)
        plt.close()

def generate_report(results):
    report_path = f"{OUTPUT_REPORT_DIR}/EXTRACTED_ifu_v8.md"
    with open(report_path, "w") as f:
        f.write("# IFU v8 Extracted Spectra Comparison\n\n")
        f.write("This report compares the v8 0.5\" radius circular aperture extractions (summed in Jy) to the previous standard `extract1d` results.\n\n")
        
        for pid, data in results.items():
            f.write(f"## {data['name']} (PID {pid})\n\n")
            f.write(f"**Peak Spaxel (pixel coords):** {data['peak']}\n\n")
            f.write(f"![{data['name']} Comparison](extraction_comparison_{pid}.png)\n\n")
            
            # Simple summary
            # Exclude nans for comparison
            mask_v8 = np.isfinite(data['spec_v8']) & (data['spec_v8'] > 0)
            mask_old = np.isfinite(data['spec_old']) & (data['spec_old'] > 0)
            med_v8 = np.nanmedian(data['spec_v8'][mask_v8]) if np.any(mask_v8) else 0
            med_old = np.nanmedian(data['spec_old'][mask_old]) if np.any(mask_old) else 1e-10
            
            ratio = med_v8 / med_old
            f.write(f"- Median flux level (v8): {med_v8:.6f} Jy\n")
            f.write(f"- Median flux level (previous): {med_old:.6f} Jy\n")
            f.write(f"- Median flux ratio (v8 / previous): {ratio:.3f}\n")
            f.write(f"- Note: The v8 extraction uses a fixed circular aperture of 0.5\" radius centered on the brightest pixel, with proper `MJy/sr` to `Jy` conversion via `PIXAR_SR`.\n\n")

if __name__ == "__main__":
    results = run_extraction()
    plot_comparison(results)
    generate_report(results)
    print("Done.")
