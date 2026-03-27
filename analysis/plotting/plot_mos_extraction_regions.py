
import os
import sys
import matplotlib.pyplot as plt

# Mock tqdm and astroquery to avoid ModuleNotFoundError
import sys
from unittest.mock import MagicMock

sys.modules['tqdm'] = MagicMock()
sys.modules['astroquery'] = MagicMock()
sys.modules['astroquery.mast'] = MagicMock()

from jwst import datamodels

# Add utilities and mos to path
sys.path.append('/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/analysis/util_mos')
import mos

def main():
    # Data directory
    data_dir = '/Users/dcoe/NIRSpec/wavext/data/PID1492'
    plot_dir = '/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/plots/extraction'
    os.makedirs(plot_dir, exist_ok=True)

    # Selected files for plotting side-by-side
    cases = [
        {
            'name': 'G140M',
            'rate': os.path.join(data_dir, 'jw01492003001_03102_00005_nrs2_rate.fits'),
            'e2d': os.path.join(data_dir, 'jw01492003001_03102_00005_nrs2_g140m_extract_2d.fits')
        },
        {
            'name': 'G235M',
            'rate': os.path.join(data_dir, 'jw01492003001_03104_00004_nrs2_rate.fits'),
            'e2d': os.path.join(data_dir, 'jw01492003001_03104_00004_nrs2_g235m_extract_2d.fits')
        },
        {
            'name': 'G395M',
            'rate': os.path.join(data_dir, 'jw01492003001_03106_00002_nrs2_rate.fits'),
            'e2d': os.path.join(data_dir, 'jw01492003001_03106_00002_nrs2_g395m_extract_2d.fits')
        }
    ]

    rate_files = []
    slit_models = []

    for case in cases:
        if os.path.exists(case['rate']) and os.path.exists(case['e2d']):
            rate_files.append(case['rate'])
            # Load the MultiSlitModel to get extraction regions
            model = datamodels.open(case['e2d'])
            slit_models.append(model)
        else:
            print(f"Skipping {case['name']} - missing files")

    if not rate_files:
        print("No files found to plot.")
        return

    # Call the mos utility function
    # We'll use a manually specified output path
    output_png = os.path.join(plot_dir, 'MOS_rate_extractions_NRS2.png')
    
    print(f"Generating side-by-side plot for {len(rate_files)} gratings...")
    
    # show_MOS_rate_files(rate_files, slit_models=[], save_plot=False, close_plot=False, integration=None, vmin=None, vmax=None, show_colorbar=True)
    mos.show_MOS_rate_files(
        rate_files, 
        slit_models=slit_models, 
        save_plot=output_png, 
        vmin=-0.001, 
        vmax=0.01, # Adjusted for visibility of faint light
        show_colorbar=True
    )
    
    print(f"Plot saved to: {output_png}")

if __name__ == '__main__':
    main()
