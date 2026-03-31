"""Plot FS v2 k/alpha/beta coefficients (log scale), mirroring Parlanti Fig 3 style."""
import os
import numpy as np
import matplotlib.pyplot as plt

BASE       = '/Users/dcoe/NIRSpec/wavext'
COEFFS_DIR = f'{BASE}/nirspec_wavext_work/plots/Parlanti/cal/fs_v2'
OUTPUT_DIR = f'{BASE}/nirspec_wavext_work/reports/329_fs_v2'
os.makedirs(OUTPUT_DIR, exist_ok=True)

fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=False)
fig.suptitle('FS Calibration v2 (329) — k/α/β coefficients\n'
             'k: median R from P330E + J1743045 | α, β: Parlanti published',
             fontsize=12)

colors = {'k': 'k', 'alpha': 'darkorange', 'beta': '#3498db'}

for ax, grating in zip(axes, ['G140M', 'G235M']):
    csv_path = f'{COEFFS_DIR}/coeffs_fs_v2_{grating}.csv'
    if not os.path.exists(csv_path):
        ax.text(0.5, 0.5, f'{csv_path} not found', transform=ax.transAxes, ha='center')
        continue
    data = np.loadtxt(csv_path, delimiter=',', skiprows=1)
    wl, k, a, b = data[:, 0], data[:, 1], data[:, 2], data[:, 3]

    # Replace zero/nan alpha,beta with a small floor for log-scale plotting
    a_plot = np.where(a > 0, a, 1e-5)
    b_plot = np.where(b > 0, b, 1e-5)

    ax.plot(wl, k,      color=colors['k'],     lw=2.0, label=r'$k(\lambda)$ [FS v2]')
    ax.plot(wl, a_plot, color=colors['alpha'], lw=1.5, label=r'$\tilde{\alpha}(\lambda)$ [Parlanti]')
    ax.plot(wl, b_plot, color=colors['beta'],  lw=1.5, label=r'$\tilde{\beta}(\lambda)$ [Parlanti]')
    ax.axhline(1.0, color='r', lw=0.8, ls='--', alpha=0.5, label='k=1')

    ax.set_yscale('log')
    ax.set_ylim(1e-4, 3.0)
    ax.set_xlabel('Wavelength [µm]', fontsize=11)
    ax.set_ylabel('Value', fontsize=11)
    ax.set_title(f'{grating}', fontsize=12)
    ax.grid(True, which='both', alpha=0.2)
    ax.legend(fontsize=9)

plt.tight_layout()
outpath = f'{OUTPUT_DIR}/fs_v2_coeffs_log.png'
plt.savefig(outpath, dpi=200, bbox_inches='tight')
plt.close()
print(f'SAVING: {outpath}')
