import os
import numpy as np
from astropy.io import fits

cache_dir = "/Users/dcoe/crds_cache/references/jwst/nirspec/"
all_files = [f for f in os.listdir(cache_dir) if f.endswith(".fits")]

for f in all_files:
    try:
        with fits.open(os.path.join(cache_dir, f)) as hdul:
            if 'SCI' in hdul:
                data = hdul['SCI'].data
                if data.size > 0:
                    std = np.nanstd(data)
                    if std > 1e-4:
                        print(f"{f}: std={std:.6f} range=[{np.nanmin(data):.3f}, {np.nanmax(data):.3f}] shape={data.shape}")
    except:
        pass
