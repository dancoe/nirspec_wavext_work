import os, sys
import numpy as np
import matplotlib.pyplot as plt
from stdatamodels.jwst import datamodels
from scipy.interpolate import interp1d

# Paths
DATA_DIR = '/Users/dcoe/NIRSpec/wavext/data'
WAVEXT_DIR = '/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work'

from astropy.io import fits

def load_spec(path):
    print(f"Loading {os.path.basename(path)}...")
    with fits.open(path) as hdul:
        if 'EXTRACT1D' in hdul:
            data = hdul['EXTRACT1D'].data
            w = data['wavelength']
            f = data['flux']
            msk = (w > 0) & np.isfinite(f) & (f > 0)
            return w[msk], f[msk]
    return np.array([]), np.array([])

def solve_grating(grating, wav_grid, boundary_wav):
    print(f"\n--- Solving for {grating} ---")
    
    # Sources configuration
    sources = {
        'IRAS-05248': {
            'prism': f'{DATA_DIR}/PID1492/jw01492-o001_t001-s000000001_nirspec_clear-prism-s200a1-subs200a1_x1d.fits',
            'ext_n1': f'{DATA_DIR}/PID1492/jw01492-o001_t001-s000000001_nirspec_f100lp-g140m-s200a1-subs200a1_x1d.fits' if grating == 'G140M' else f'{DATA_DIR}/PID1492/jw01492-o001_t001-s000000001_nirspec_f170lp-g235m-s200a1-subs200a1_x1d.fits',
            'ext_n2': f'{DATA_DIR}/PID1492/jw01492003001_03102_00005_nrs2_g140m_extract_1d.fits' if grating == 'G140M' else f'{DATA_DIR}/PID1492/jw01492003001_03104_00004_nrs2_g235m_extract_1d.fits'
        },
        'G191-B2B': {
            'prism': f'{DATA_DIR}/PID1537_G191-B2B/jw01537-o009_t001-s000000001_nirspec_clear-prism-s1600a1-sub512_x1d.fits',
            'ext_n1': f'{DATA_DIR}/PID1537_G191-B2B/jw01537-o003_t001-s000000001_nirspec_f100lp-g140m-s1600a1-sub2048_x1d.fits' if grating == 'G140M' else f'{DATA_DIR}/PID1537_G191-B2B/jw01537-o007_t001-s000000001_nirspec_f170lp-g235m-s1600a1-sub2048_x1d.fits',
            'ext_n2': f'{DATA_DIR}/PID1537_G191-B2B/jw01537007001_07101_00003_nrs2_extract_1d.fits' if grating == 'G140M' else f'{DATA_DIR}/PID1537_G191-B2B/jw01537007001_09101_00003_nrs2_extract_1d.fits'
        },
        'P330E': {
            'prism': f'{DATA_DIR}/PID1538_P330E/jw01538-o107_t002-s000000001_nirspec_clear-prism-s1600a1-sub512_x1d.fits',
            'ext_n1': f'{DATA_DIR}/PID1538_P330E/jw01538-o130_t002-s000000001_nirspec_f100lp-g140m-s1600a1-sub2048_x1d.fits' if grating == 'G140M' else f'{DATA_DIR}/PID1538_P330E/jw01538-o160_t002-s000000001_nirspec_f170lp-g235m-s1600a1-sub2048_x1d.fits',
            'ext_n2': f'{DATA_DIR}/PID1538_P330E/jw01538160001_06101_00001_nrs2_extract_1d.fits' if grating == 'G140M' else f'{DATA_DIR}/PID1538_P330E/jw01538160001_08101_00005_nrs2_extract_1d.fits'
        }
    }

    # Build Matrices
    all_Cs = [] 
    all_As = [] 

    for name, paths in sources.items():
        w_prism, f_prism = load_spec(paths['prism'])
        w_ext, c_ext = load_spec(paths['ext_n2'])
        
        # Scaling
        mask_prism = (w_prism > boundary_wav - 0.1) & (w_prism < boundary_wav + 0.1)
        mask_ext = (w_ext > boundary_wav) & (w_ext < boundary_wav + 0.2)
        
        scale = np.median(f_prism[mask_prism]) / np.median(c_ext[mask_ext])
        print(f"Source {name} scale factor: {scale:.2e}")
        c_ext_scaled = c_ext * scale
        
        interp_prism = interp1d(w_prism, f_prism, bounds_error=False, fill_value=0)
        interp_ext = interp1d(w_ext, c_ext_scaled, bounds_error=False, fill_value=0)
        
        source_A, source_C = [], []
        for L in wav_grid:
            source_A.append([interp_prism(L), interp_prism(L/2), interp_prism(L/3)])
            source_C.append(interp_ext(L))
        
        all_As.append(source_A)
        all_Cs.append(source_C)

    # Solve 3x3 at each wavelength
    ks, alphas, betas = [], [], []
    from scipy.optimize import nnls
    for i in range(len(wav_grid)):
        A_mat = np.array([all_As[0][i], all_As[1][i], all_As[2][i]])
        C_vec = np.array([all_Cs[0][i], all_Cs[1][i], all_Cs[2][i]])
        try:
            sol, _ = nnls(A_mat, C_vec)
            ks.append(sol[0]); alphas.append(sol[1]); betas.append(sol[2])
        except:
            ks.append(np.nan); alphas.append(np.nan); betas.append(np.nan)

    # Smoothing
    def smooth(y, box_pts):
        box = np.ones(box_pts)/box_pts
        return np.convolve(y, box, mode='same')

    ks, alphas, betas = smooth(np.array(ks), 11), smooth(np.array(alphas), 11), smooth(np.array(betas), 11)

    # Plot results
    plt.figure(figsize=(10,6))
    plt.plot(wav_grid, ks, label='k (1st)')
    plt.plot(wav_grid, alphas, label='alpha (2nd)')
    plt.plot(wav_grid, betas, label='beta (3rd)')
    plt.ylim(-0.1, 1.2) # Added from original plotting section
    plt.title(f'Direct 3-Source Solver Results ({grating})') # Adjusted title
    plt.xlabel('Wavelength (um)'); plt.ylabel('Derived Coefficients'); plt.legend(); plt.grid() # Adjusted ylabel
    
    output_plot = os.path.join(WAVEXT_DIR, f'plots/Parlanti/cal/direct_3source_solve_{grating}.png')
    os.makedirs(os.path.dirname(output_plot), exist_ok=True)
    plt.savefig(output_plot, dpi=200); plt.close()
    print(f"SUCCESS: Plot saved to {output_plot}")

    # Save coefficients
    coeff_file = os.path.join(WAVEXT_DIR, f'plots/Parlanti/cal/direct_coefficients_{grating}.txt')
    with open(coeff_file, 'w') as f:
        f.write("# wav_um k alpha_tilde beta_tilde\n")
        for w, k, a_val, b_val in zip(wav_grid, ks, alphas, betas):
            f.write(f"{w:.4f} {k:.6f} {a_val:.6f} {b_val:.6f}\n")
    print(f"SUCCESS: Saved {grating} coefficients and plot.")

if __name__ == "__main__":
    solve_grating('G235M', np.linspace(3.4, 5.2, 250), 3.3)
    solve_grating('G140M', np.linspace(2.0, 3.2, 250), 1.9)
