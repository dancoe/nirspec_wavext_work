import os, sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Paths
COEFF_DIR = '/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/plots/Parlanti/cal/153678_v3'

def plot_unified_coeffs():
    c140 = pd.read_csv(f'{COEFF_DIR}/coeffs_G140M.csv')
    c235 = pd.read_csv(f'{COEFF_DIR}/coeffs_G235M.csv')
    
    plt.figure(figsize=(12, 10))
    
    plt.subplot(2,1,1)
    plt.plot(c140['wav'], c140['k'], 'b-', lw=2, label='k(λ) 1st Order')
    plt.plot(c140['wav'], c140['alpha'], 'c-', label='α̃(λ) 2nd')
    plt.plot(c140['wav'], c140['beta'], 'm-', label='β̃(λ) 3rd')
    plt.axhline(1.0, color='gray', linestyle='--', alpha=0.5)
    plt.title('G140M Final Extension Coefficients (V3)')
    plt.ylabel('Efficiency / Coefficient'); plt.legend(); plt.grid(alpha=0.1)
    
    plt.subplot(2,1,2)
    plt.plot(c235['wav'], c235['k'], color='#B8860B', lw=2, label='k(λ) 1st Order')
    plt.plot(c235['wav'], c235['alpha'], color='orange', label='α̃(λ) 2nd')
    plt.plot(c235['wav'], c235['beta'], color='gold', label='β̃(λ) 3rd')
    plt.axhline(1.0, color='gray', linestyle='--', alpha=0.5)
    plt.title('G235M Final Extension Coefficients (V3)')
    plt.ylabel('Efficiency / Coefficient'); plt.legend(); plt.grid(alpha=0.1)
    plt.xlabel('Wavelength [µm]')
    
    plt.tight_layout()
    plt.savefig(f'{COEFF_DIR}/FINAL_COEFFS_V3.png')
    plt.close()

plot_unified_coeffs()
