#!/usr/bin/env python
"""Quick inspection of spectra wavelength ranges and flux levels."""
import os, sys
import numpy as np
sys.path.insert(0, '/Users/dcoe/NIRSpec/wavext/jwst_nirspec_wavext')
os.environ.setdefault('CRDS_PATH', os.path.expanduser('~/crds_cache'))
os.environ.setdefault('CRDS_SERVER_URL', 'https://jwst-crds.stsci.edu')

from stdatamodels.jwst import datamodels

DATA = '/Users/dcoe/NIRSpec/wavext/data/PID1492'

files = {
    'PRISM':    f'{DATA}/jw01492-o001_t001-s000000001_nirspec_clear-prism-s200a1-subs200a1_x1d.fits',
    'G140M_nom':f'{DATA}/jw01492-o001_t001-s000000001_nirspec_f100lp-g140m-s200a1-subs200a1_x1d.fits',
    'G235M_nom':f'{DATA}/jw01492-o001_t001-s000000001_nirspec_f170lp-g235m-s200a1-subs200a1_x1d.fits',
    'G395M_nom':f'{DATA}/jw01492-o001_t001-s000000001_nirspec_f290lp-g395m-s200a1-subs200a1_x1d.fits',
    'G140M_ext':f'{DATA}/jw01492003001_03102_00005_nrs2_g140m_extract_1d.fits',
    'G235M_ext':f'{DATA}/jw01492003001_03104_00004_nrs2_g235m_extract_1d.fits',
    'G395M_ext':f'{DATA}/jw01492003001_03106_00002_nrs2_g395m_extract_1d.fits',
}

def read_spec(path):
    m = datamodels.open(path)
    specs = getattr(m, 'spec', [])
    all_w, all_f = [], []
    for s in specs:
        w = s.spec_table['wavelength']
        f = s.spec_table['flux']
        msk = (w > 0) & np.isfinite(f) & (f > 0)
        all_w.append(np.array(w[msk]))
        all_f.append(np.array(f[msk]))
    m.close()
    if all_w:
        w = np.concatenate(all_w)
        f = np.concatenate(all_f)
        i = np.argsort(w)
        return w[i], f[i]
    return np.array([]), np.array([])

for name, path in files.items():
    try:
        w, f = read_spec(path)
        print(f'{name:15s}: {len(w):5d} pts  λ=[{w[0]:.3f},{w[-1]:.3f}] µm  '
              f'flux=[{f.min():.3e},{f.median() if hasattr(f,"median") else np.median(f):.3e},{f.max():.3e}] Jy')
    except Exception as e:
        print(f'{name}: ERROR - {e}')
