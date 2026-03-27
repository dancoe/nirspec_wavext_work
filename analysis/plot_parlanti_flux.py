import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from stdatamodels.jwst import datamodels

def get_level3_flux(file_path):
    # Extracts all valid fluxes from a MAST x1d/x1dints file
    try:
        model = datamodels.open(file_path)
        all_wav, all_flux = [], []
        # Support both 'spec' list (MultiSpecModel) or direct .spec (SpecModel)
        specs = getattr(model, 'spec', [])
        if not specs and hasattr(model, 'spec_table'):
            specs = [model]
            
        for s in specs:
            wav = s.spec_table['wavelength']
            flx = s.spec_table['flux']
            mask = (wav > 0) & (~np.isnan(flx)) & (flx != 0.0)
            all_wav.append(wav[mask])
            all_flux.append(flx[mask])
        if len(all_wav) > 0:
            return np.concatenate(all_wav), np.concatenate(all_flux), model.meta.instrument.grating
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
    return None, None, None

def main():
    data_dir = '/Users/dcoe/NIRSpec/wavext/data/PID1492'
    
    # 1. PRISM (serves as f(lambda))
    prism_files = glob.glob(f'{data_dir}/*0310e*x1d.fits')
    prism_wav, prism_flux, _ = get_level3_flux(prism_files[0]) if prism_files else (None, None, None)
    
    # 2. G235H (Nominal, from x1d downloaded from MAST)
    g235h_files = glob.glob(f'{data_dir}/*0310a*x1d.fits')
    g235n_wav, g235n_flux, _ = get_level3_flux(g235h_files[0]) if g235h_files else (None, None, None)

    # 3. G140M (Nominal)
    g140m_files = glob.glob(f'{data_dir}/*03102*x1d.fits')
    g140mn_wav, g140mn_flux, _ = get_level3_flux(g140m_files[0]) if g140m_files else (None, None, None)
    
    # Let's plot PRISM and G140M/G235H nominal to see what they look like
    plt.figure(figsize=(12, 7))
    if prism_wav is not None:
        idx = np.argsort(prism_wav)
        plt.plot(prism_wav[idx], prism_flux[idx], label='PRISM (Reference f(λ))', color='black', alpha=0.3, linewidth=1)
        
    if g235n_wav is not None:
        idx = np.argsort(g235n_wav)
        plt.plot(g235n_wav[idx], g235n_flux[idx], label='G235H (Nominal)', color='blue', alpha=0.7)
        
    if g140mn_wav is not None:
        idx = np.argsort(g140mn_wav)
        plt.plot(g140mn_wav[idx], g140mn_flux[idx], label='G140M (Nominal)', color='green', alpha=0.7)

    # 4. Our custom extended extractions
    custom_g235h = f'{data_dir}/jw01492001001_0310a_00007_nrs2_g235h_extract_1d.fits'
    cw, cf, _ = get_level3_flux(custom_g235h)
    if cw is not None:
        idx = np.argsort(cw)
        plt.plot(cw[idx], cf[idx], label='G235H (Extended)', color='blue', linestyle='--', linewidth=2)
        
    custom_g140m = f'{data_dir}/jw01492003001_03102_00005_nrs2_extract_1d.fits'
    cw2, cf2, _ = get_level3_flux(custom_g140m)
    if cw2 is not None:
        idx = np.argsort(cw2)
        plt.plot(cw2[idx], cf2[idx], label='G140M (Extended S(λ))', color='red', linestyle='--', linewidth=2)

    plt.xlabel(r'Wavelength ($\mu$m)')
    plt.ylabel('Flux (Jy)')
    plt.yscale('log')
    plt.xlim(0.6, 5.6)
    plt.grid(True, which='both', linestyle='--', alpha=0.5)
    plt.legend(loc='lower right', frameon=True, facecolor='white', framealpha=0.9)
    plt.title('NIRSpec Multi-Order Overlap Analysis (PID 1492 FS)')
    
    # Shade the 2nd order region for G140M (approx > 1.9 um)
    plt.axvspan(1.88, 3.3, alpha=0.1, color='red', label='2nd Order Contamination Region')

    
    plot_path = '/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/plots/parlanti_flux_demo.png'
    os.makedirs(os.path.dirname(plot_path), exist_ok=True)
    plt.savefig(plot_path)
    print(f"Saved plot to {plot_path}")

if __name__ == '__main__':
    main()
