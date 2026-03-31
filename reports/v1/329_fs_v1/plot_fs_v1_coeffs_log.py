import os
import numpy as np
import matplotlib.pyplot as plt

# Paths
BASE        = '/Users/dcoe/NIRSpec/wavext'
COEFFS_DIR  = f'{BASE}/nirspec_wavext_work/plots/Parlanti/cal/fs_v1'
OUTPUT_DIR  = f'{BASE}/nirspec_wavext_work/reports/329_fs_v1'
os.makedirs(OUTPUT_DIR, exist_ok=True)

def plot_log_coeffs():
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
    fig.suptitle('FS Calibration v1 (329) — k/α/β coefficients', fontsize=14)

    colors = {'k': 'violet', 'alpha': 'darkorange', 'beta': '#3498db'} # Using Parlanti-like colors

    for ax, grating in zip(axes, ['G140M', 'G235M']):
        csv_path = f'{COEFFS_DIR}/coeffs_fs_v1_{grating}.csv'
        if not os.path.exists(csv_path):
            print(f'Warning: {csv_path} not found')
            continue
            
        data = np.loadtxt(csv_path, delimiter=',', skiprows=1)
        wl = data[:, 0]
        k = data[:, 1]
        a = data[:, 2]
        b = data[:, 3]

        ax.plot(wl, k, color=colors['k'], lw=2, label=r'$k(\lambda)$')
        ax.plot(wl, a, color=colors['alpha'], lw=1.5, label=r'$\tilde{\alpha}(\lambda)$')
        ax.plot(wl, b, color=colors['beta'], lw=1.5, label=r'$\tilde{\beta}(\lambda)$')

        ax.set_yscale('log')
        ax.set_ylim(1e-4, 2.0)
        ax.set_xlabel('Wavelength [µm]', fontsize=12)
        if ax == axes[0]:
            ax.set_ylabel('Value', fontsize=12)
        
        ax.set_title(f'{grating}', fontsize=12)
        ax.grid(True, which='both', alpha=0.2)
        ax.legend(fontsize=10)

    plt.tight_layout()
    outpath = f'{OUTPUT_DIR}/fs_v1_coeffs_log.png'
    plt.savefig(outpath, dpi=200, bbox_inches='tight')
    plt.close()
    print(f'Saved: {outpath}')

if __name__ == '__main__':
    plot_log_coeffs()
