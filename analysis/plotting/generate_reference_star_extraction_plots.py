import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits
from jwst import datamodels
from jwst.datamodels import SlitModel
import sys

# Ensure project code is in path
sys.path.append('/Users/dcoe/NIRSpec/wavext/jwst_nirspec_wavext')
from analysis.util_mos import mos

def find_files():
    """
    Scans data directory to find rate files and their corresponding 2D extracts.
    Groups them by (Program ID / Star) and Grating.
    """
    data_dir = '../data'
    # List of known star/program folders
    folders = {
        'PID 1492': 'PID1492',
        'PID 1537': 'PID1537_G191-B2B',
        'PID 1538': 'PID1538_P330E'
    }
    
    star_gratings = {}
    
    for label, folder in folders.items():
        folder_path = os.path.join(data_dir, folder)
        if not os.path.exists(folder_path):
            continue
            
        star_gratings[label] = {}
        
        # Find all extract_2d files
        extract_files = [f for f in os.listdir(folder_path) if f.endswith('extract_2d.fits')]
        extract_files.sort()
        
        # Find all rate files
        rate_files = [f for f in os.listdir(folder_path) if f.endswith('rate.fits')]
        
        # Group by Grating
        for extract_f in extract_files:
            # For PID 1492, prioritize 01001 as requested
            if label == 'PID 1492' and 'jw01492001001' not in extract_f:
                continue
            
            det = 'nrs1' if 'nrs1' in extract_f else 'nrs2'
            extract_path = os.path.join(folder_path, extract_f)
            
            try:
                with fits.open(extract_path) as h:
                    grating = h[0].header.get('GRATING', 'Unknown').upper()
            except:
                continue
                
            if grating not in star_gratings[label]:
                star_gratings[label][grating] = {
                    'nrs1': None, 'nrs2': None,
                    'extract_nrs1': None, 'extract_nrs2': None
                }
            
            # Find matching rate file
            base = extract_f.split('_' + det)[0]
            rate_f = base + '_' + det + '_rate.fits'
            rate_path = os.path.join(folder_path, rate_f)
            
            if os.path.exists(rate_path):
                star_gratings[label][grating][det] = rate_path
                star_gratings[label][grating]['extract_' + det] = extract_path
                
        # Now try to find any MISSING pairs (e.g. rate file exists but no extract file yet)
        for rate_f in rate_files:
            if label == 'PID 1492' and 'jw01492001001' not in rate_f:
                continue
                
            det = 'nrs1' if 'nrs1' in rate_f else 'nrs2'
            rate_path = os.path.join(folder_path, rate_f)
            
            try:
                with fits.open(rate_path) as h:
                    grating = h[0].header.get('GRATING', 'Unknown').upper()
            except:
                continue
                
            if grating not in star_gratings[label]:
                star_gratings[label][grating] = {
                    'nrs1': None, 'nrs2': None,
                    'extract_nrs1': None, 'extract_nrs2': None
                }
            
            if star_gratings[label][grating][det] is None:
                star_gratings[label][grating][det] = rate_path
                # Check if extract exists but wasn't found (e.g. name mismatch)
                ext_f = rate_f.replace('_rate.fits', '_extract_2d.fits')
                ext_path = os.path.join(folder_path, ext_f)
                if os.path.exists(ext_path):
                    star_gratings[label][grating]['extract_' + det] = ext_path
                    
    return star_gratings

def plot_all():
    star_gratings = find_files()
    if not star_gratings:
        print("No files discovered!")
        return

    os.makedirs('plots/extraction', exist_ok=True)

    for label, gratings in star_gratings.items():
        print(f"Plotting for {label}...")
        
        # Sort gratings: G140, G235, G395, PRISM
        sorted_gratings = sorted(gratings.keys(), key=lambda x: (0 if '140' in x else (1 if '235' in x else (2 if '395' in x else 3))))
        
        rate_files = []
        slit_models = []
        titles = []
        
        for grating in sorted_gratings:
            info = gratings[grating]
            nrs1 = info['nrs1']
            nrs2 = info['nrs2']
            ext1_path = info['extract_nrs1']
            ext2_path = info['extract_nrs2']
            
            # Setup for 2 columns: NRS1, NRS2
            rate_files.append(nrs1)
            rate_files.append(nrs2)
            
            # Load extract models for WCS
            ext1 = datamodels.open(ext1_path) if ext1_path else None
            ext2 = datamodels.open(ext2_path) if ext2_path else None
            slit_models.append(ext1)
            slit_models.append(ext2)
            
            # Titles
            t1 = f"{grating} NRS1 | {os.path.basename(nrs1) if nrs1 else '(Missing)'}"
            t2 = f"{grating} NRS2 | {os.path.basename(nrs2) if nrs2 else '(Missing)'}"
            titles.append(t1)
            titles.append(t2)

        if not rate_files:
            continue

        # Shared scale estimation
        all_data = []
        for rf in rate_files:
            if rf and os.path.exists(rf):
                with fits.open(rf) as h:
                    data = h['SCI'].data
                    if data.ndim == 3: data = data[0]
                    # Sample some areas
                    all_data.append(data.ravel()[::100])
        
        vmin, vmax = None, None
        if all_data:
            combined = np.concatenate(all_data)
            combined = combined[np.isfinite(combined)]
            if len(combined) > 0:
                vmin = np.percentile(combined, 20)
                vmax = np.percentile(combined, 99.5)
                print(f"Using shared scale: vmin={vmin:.4f}, vmax={vmax:.4f}")

        plot_path = f"plots/extraction/{label.replace(' ', '_')}_All_Gratings_NRS1_NRS2_v6.png"
        nrows = len(sorted_gratings)
        ncols = 2
        print(f"  Generating {nrows}x{ncols} plot...")
        
        mos.show_MOS_rate_files(rate_files, slit_models, titles=titles,
                                save_plot=plot_path, close_plot=True,
                                autocrop=True, vmin=vmin, vmax=vmax,
                                show_labels=False, # Suppress red labels
                                nrows=nrows, ncols=ncols)
        print(f"SAVING {plot_path}")
        
    # Cleanup open models
    for sm in slit_models:
        if sm: sm.close()

if __name__ == '__main__':
    plot_all()
