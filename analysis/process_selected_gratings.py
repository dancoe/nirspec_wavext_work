import os
import subprocess
import glob
import logging
from stdatamodels.jwst import datamodels
import numpy as np
import matplotlib.pyplot as plt

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger("wavext_workflow")

DATA_DIR = '/Users/dcoe/NIRSpec/wavext/data/PID1492'
PLOT_DIR = '/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/plots'
REF_FILE = '/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/reference_files/wavelengthrange_extended.asdf'
ANALYSIS_DIR = '/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/analysis'

# Configuration for processing
targets = [
    {
        'name': 'G140M',
        'rate': 'jw01492003001_03102_00005_nrs2_rate.fits',
        'nominal_x1d': 'jw01492-o001_t001-s000000001_nirspec_f100lp-g140m-s200a1-subs200a1_x1d.fits',
        'overlap_range': [1.88, 3.3]
    },
    {
        'name': 'G235M',
        'rate': 'jw01492003001_03104_00004_nrs2_rate.fits',
        'nominal_x1d': 'jw01492-o001_t001-s000000001_nirspec_f170lp-g235m-s200a1-subs200a1_x1d.fits',
        'overlap_range': [3.1, 5.3]
    },
    {
        'name': 'G395M',
        'rate': 'jw01492003001_03106_00002_nrs2_rate.fits',
        'nominal_x1d': 'jw01492-o001_t001-s000000001_nirspec_f290lp-g395m-s200a1-subs200a1_x1d.fits',
        'overlap_range': [5.1, 5.6]
    },
    {
        'name': 'PRISM',
        'rate': 'jw01492005001_17101_00006_nrs2_rate.fits',
        'nominal_x1d': 'jw01492001001_0310e_00002_nrs1_x1d.fits',
        'overlap_range': [1.2, 5.6]
    }
]

def run_step(script, input_file, output_file):
    cmd = ['python3', os.path.join(ANALYSIS_DIR, script), input_file, output_file]
    logger.info(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

def get_flux_data(file_path):
    if not os.path.exists(file_path):
        return None, None
    try:
        model = datamodels.open(file_path)
        all_wav, all_flux = [], []
        specs = getattr(model, 'spec', [])
        if not specs and hasattr(model, 'spec_table'):
            specs = [model]
            
        for s in specs:
            wav = s.spec_table['wavelength']
            flx = s.spec_table['flux']
            msk = (wav > 0) & (flx > 0) & (~np.isnan(flx)) 
            all_wav.append(wav[msk])
            all_flux.append(flx[msk])
        if len(all_wav) > 0:
            return np.concatenate(all_wav), np.concatenate(all_flux)
    except Exception as e:
        logger.error(f"Error reading {file_path}: {e}")
    return None, None

def main():
    os.makedirs(PLOT_DIR, exist_ok=True)
    
    # First, ensure we have the PRISM baseline
    prism_nom_path = os.path.join(DATA_DIR, targets[3]['nominal_x1d'])
    pw, pf = get_flux_data(prism_nom_path)
    
    for t in targets:
        logger.info(f"Processing {t['name']}...")
        input_rate = os.path.join(DATA_DIR, t['rate'])
        wcs_out = input_rate.replace('_rate.fits', f'_{t["name"].lower()}_wcs.fits')
        e2d_out = input_rate.replace('_rate.fits', f'_{t["name"].lower()}_extract_2d.fits')
        e1d_out = input_rate.replace('_rate.fits', f'_{t["name"].lower()}_extract_1d.fits')
        
        # 1. Assign WCS
        run_step('run_assign_wcs.py', input_rate, wcs_out)
        
        # 2. Extract 2D
        run_step('run_extract_2d.py', wcs_out, e2d_out)
        
        # 3. Extract 1D
        run_step('run_extract_1d.py', e2d_out, e1d_out)
        
        # 4. Plot
        nw, nf = get_flux_data(os.path.join(DATA_DIR, t['nominal_x1d']))
        ew, ef = get_flux_data(e1d_out)
        
        plt.figure(figsize=(10, 6))
        if nw is not None:
             idx = np.argsort(nw)
             plt.plot(nw[idx], nf[idx], label='Nominal', color='blue', alpha=0.5, linewidth=1)
        if ew is not None:
             idx = np.argsort(ew)
             plt.plot(ew[idx], ef[idx], label='Extended', color='red', alpha=0.7, linewidth=0.8)
        if pw is not None:
             idx = np.argsort(pw)
             plt.plot(pw[idx], pf[idx], label='PRISM Nominal', color='black', alpha=0.3, linewidth=0.5)
        
        if t['overlap_range']:
            plt.axvspan(t['overlap_range'][0], t['overlap_range'][1], color='yellow', alpha=0.1, label='Extended/Overlap')
            
        plt.yscale('log')
        plt.xlabel('Wavelength (um)')
        plt.ylabel('Flux (Jy)')
        plt.title(f'Comparison: {t["name"]} (PID 1492)')
        plt.legend()
        plt.grid(alpha=0.2)
        
        plot_path = os.path.join(PLOT_DIR, f'{t["name"].lower()}_comparison.png')
        plt.savefig(plot_path)
        logger.info(f"Saved plot: {plot_path}")
        plt.close()

if __name__ == '__main__':
    main()
