import os
import numpy as np
from astropy.io import fits

cache_dir = "/Users/dcoe/crds_cache/references/jwst/nirspec/"
files = [f for f in os.listdir(cache_dir) if f.startswith("jwst_nirspec_sflat_")]

for f in files:
    path = os.path.join(cache_dir, f)
    with fits.open(path) as hdul:
        if 'SCI' in hdul:
            data = hdul['SCI'].data
            if data.size > 0:
                std = np.nanstd(data)
                min_val = np.nanmin(data)
                max_val = np.nanmax(data)
                print(f"{f}: std={std:.6f} range=[{min_val:.3f}, {max_val:.3f}] size={data.shape}")
            else:
                print(f"{f}: EMPTY SCI")
