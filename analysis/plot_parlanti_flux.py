import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from stdatamodels.jwst import datamodels

def get_level3_flux(file_path):
    try:
        model = datamodels.open(file_path)
        all_wav, all_flux = [], []
        specs = getattr(model, 'spec', [])
        if not specs and hasattr(model, 'spec_table'):
            specs = [model]
            
        for s in specs:
            wav = s.spec_table['wavelength']
            flx = s.spec_table['flux']
            # Clip values nearly zero or negative to nan for plotting
            msk = (wav > 0) & (flx > 1e-12) & (~np.isnan(flx)) 
            all_wav.append(wav[msk])
            all_flux.append(flx[msk])
        if len(all_wav) > 0:
            return np.concatenate(all_wav), np.concatenate(all_flux), model.meta.instrument.grating
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
    return None, None, None

def plot_grating(grating_name, nominal_wav, nominal_flux, extended_wav, extended_flux, prism_wav, prism_flux, plot_filename, order_shade_range=None):
    plt.figure(figsize=(10, 6))
    
    # PRISM baseline
    if prism_wav is not None:
        idx = np.argsort(prism_wav)
        plt.plot(prism_wav[idx], prism_flux[idx], label='PRISM (Reference f(λ))', color='black', alpha=0.3, linewidth=1)
        
    # Nominal
    if nominal_wav is not None:
        idx = np.argsort(nominal_wav)
        color = 'red' if '140M' in grating_name else 'blue'
        plt.plot(nominal_wav[idx], nominal_flux[idx], label=f'{grating_name} (Nominal)', color=color, alpha=1.0, linewidth=1)
        
    # Extended
    if extended_wav is not None:
        idx = np.argsort(extended_wav)
        color = 'lightcoral' if '140M' in grating_name else 'cornflowerblue'
        plt.plot(extended_wav[idx], extended_flux[idx], label=f'{grating_name} (Extended S(λ))', color=color, alpha=0.6, linewidth=0.5)

    if order_shade_range:
        plt.axvspan(order_shade_range[0], order_shade_range[1], alpha=0.1, color='gray', label='Overlap Region')

    plt.xlabel(r'Wavelength ($\mu$m)')
    plt.ylabel('Flux (Jy)')
    plt.yscale('log')
    plt.xlim(0.6, 5.6)
    # Autoscale Y based on data
    # Remove interior grid lines
    plt.grid(False) 
    plt.legend(loc='lower right')
    plt.title(f'NIRSpec NIRSpec/FS {grating_name} Analysis (PID 1492)')

    
    os.makedirs(os.path.dirname(plot_filename), exist_ok=True)
    plt.savefig(plot_filename)
    print(f"Saved plot to {plot_filename}")
    plt.close()

def main():
    data_dir = '/Users/dcoe/NIRSpec/wavext/data/PID1492'
    plot_dir = '/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/plots'

    # PRISM
    prism_files = glob.glob(f'{data_dir}/*0310e*x1d.fits')
    pw, pf, _ = get_level3_flux(prism_files[0]) if prism_files else (None, None, None)
    
    # G140M
    g140mn_files = glob.glob(f'{data_dir}/*03102*x1d.fits')
    g140mn_w, g140mn_f, _ = get_level3_flux(g140mn_files[0]) if g140mn_files else (None, None, None)
    g140me_file = f'{data_dir}/jw01492003001_03102_00005_nrs2_extract_1d.fits'
    g140me_w, g140me_f, _ = get_level3_flux(g140me_file) if os.path.exists(g140me_file) else (None, None, None)
    
    plot_grating('G140M', g140mn_w, g140mn_f, g140me_w, g140me_f, pw, pf, 
                 f'{plot_dir}/g140m_parlanti_overlap.png', order_shade_range=[1.88, 3.3])

    # G235H
    g235hn_files = glob.glob(f'{data_dir}/*0310a*x1d.fits')
    g235hn_w, g235hn_f, _ = get_level3_flux(g235hn_files[0]) if g235hn_files else (None, None, None)
    g235he_file = f'{data_dir}/jw01492001001_0310a_00007_nrs2_g235h_extract_1d.fits'
    g235he_w, g235he_f, _ = get_level3_flux(g235he_file) if os.path.exists(g235he_file) else (None, None, None)

    plot_grating('G235H', g235hn_w, g235hn_f, g235he_w, g235he_f, pw, pf, 
                 f'{plot_dir}/g235h_parlanti_overlap.png', order_shade_range=[3.1, 5.3])

if __name__ == '__main__':
    main()
