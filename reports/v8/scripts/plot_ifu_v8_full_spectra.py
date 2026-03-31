
import os
import glob
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

# Paths
IFU_DIR = '/Users/dcoe/NIRSpec/wavext/data/IFU'
FS_DIR = '/Users/dcoe/NIRSpec/wavext/data/FS'
OUTPUT_DIR = '/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/reports/v8/ifu_v8'
os.makedirs(OUTPUT_DIR, exist_ok=True)

TARGETS = [
    {'pid': '1536', 'name': 'J1743045', 'dir': f'{IFU_DIR}/PID1536_J1743045/stage3_ext'},
    {'pid': '1537', 'name': 'G191-B2B', 'dir': f'{IFU_DIR}/PID1537_G191-B2B/stage3_ext'},
    {'pid': '1538', 'name': 'P330E', 'dir': f'{IFU_DIR}/PID1538_P330E/stage3_ext'},
    {'pid': '2186', 'name': 'UGC-5101', 'dir': f'{FS_DIR}/PID2186_UGC-5101/stage3_ext', 'science': True},
    {'pid': '2654', 'name': 'SDSSJ0841', 'dir': f'{FS_DIR}/PID2654_SDSSJ0841/stage3_ext', 'science': True},
    {'pid': '6645', 'name': 'P330E-C3', 'dir': f'{IFU_DIR}/PID6645_P330E-C3/stage3_ext'},
]

def extract_05_radius(s3d_file):
    try:
        with fits.open(s3d_file) as hdul:
            sci = hdul['SCI'].data
            hdr = hdul['SCI'].header
            wl = np.arange(hdr['NAXIS3']) * hdr['CDELT3'] + hdr['CRVAL3']
            pixar_sr = hdr.get('PIXAR_SR', 1.0)
            conv = 1e6 * pixar_sr if hdr.get('BUNIT') == 'MJy/sr' else 1.0
            
            white_light = np.nanmean(sci, axis=0)
            peak_y, peak_x = np.unravel_index(np.nanargmax(white_light), white_light.shape)
            yy, xx = np.indices(white_light.shape)
            dist_arcsec = np.sqrt(((xx - peak_x)*hdr['CDELT1'])**2 + ((yy - peak_y)*hdr['CDELT2'])**2) * 3600
            mask = dist_arcsec <= 0.5
            
            spec = np.nansum(sci[:, mask], axis=1) * conv
            return wl, spec
    except Exception as e:
        print(f"Error extracting {s3d_file}: {e}")
        return None, None

def plot_target(target):
    pid = target['pid']
    print(f"Plotting PID {pid}...")
    
    # Gratings: G140M (NRS1/2 gap 2.17-2.28), G235M (NRS1/2 gap 3.66-3.82)
    plt.figure(figsize=(12, 6))
    
    # Standard Gratings
    gratings = ['g140m', 'g235m', 'g395m']
    colors = {'g140m': 'blue', 'g235m': 'green', 'g395m': 'red'}
    
    for grating in gratings:
        # Look for s3d or x1d files
        files = glob.glob(f"{target['dir']}/*{grating}*x1d.fits")
        if not files:
            continue
            
        file = files[0]
        if target.get('science'):
            # Use 0.5" radius for science targets from s3d cube
            s3d_files = glob.glob(f"{target['dir']}/*{grating}*s3d.fits")
            if s3d_files:
                wl, spec = extract_05_radius(s3d_files[0])
            else:
                # Fallback to x1d
                with fits.open(file) as h:
                    wl, spec = h[1].data['WAVELENGTH'], h[1].data['FLUX']
        else:
            # Use standard x1d for standard stars
            with fits.open(file) as h:
                wl, spec = h[1].data['WAVELENGTH'], h[1].data['FLUX']

        if wl is None: continue
        
        # Gap aware plotting
        if grating == 'g140m':
            # Gap at ~2.17-2.28
            mask1 = wl < 2.17
            mask2 = wl > 2.28
            plt.plot(wl[mask1], spec[mask1], color=colors[grating], label=f'{grating.upper()} (v8)')
            plt.plot(wl[mask2], spec[mask2], color=colors[grating])
        elif grating == 'g235m':
            # Gap at ~3.66-3.82
            mask1 = wl < 3.66
            mask2 = wl > 3.82
            plt.plot(wl[mask1], spec[mask1], color=colors[grating], label=f'{grating.upper()} (v8)')
            plt.plot(wl[mask2], spec[mask2], color=colors[grating])
        else:
            plt.plot(wl, spec, color=colors[grating], label=f'{grating.upper()} (v8)')
            
    plt.yscale('log')
    plt.xlabel("Wavelength (um)")
    plt.ylabel("Flux (Jy)")
    plt.title(f"{target['name']} – PID {pid} – IFU wavext v8")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.xlim(0.6, 5.6)
    
    label = "Full Spectrum"
    if target.get('science'):
        label += " (0.5\" radius extraction)"
    plt.text(0.7, 0.1, label, transform=plt.gca().transAxes, fontsize=10, bbox=dict(facecolor='white', alpha=0.5))
    
    plt.savefig(f"{OUTPUT_DIR}/full_spectrum_v8_{pid}.png", dpi=150)
    plt.close()

if __name__ == "__main__":
    for target in TARGETS:
        try:
            plot_target(target)
        except Exception as e:
            print(f"Failed PID {target['pid']}: {e}")
    print("All plots generated.")
