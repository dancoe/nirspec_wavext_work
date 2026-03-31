import astropy.io.fits as f
import numpy as np
import os
RESULTS_DIR = "/Users/dcoe/NIRSpec/wavext/results/v5"
for g in ["g140m_f100lp", "g235m_f170lp"]:
    path = os.path.join(RESULTS_DIR, f"calib_v5_{g}.fits")
    if not os.path.exists(path):
        print(f"Skipping {g}, file not found.")
        continue
    with f.open(path) as hdul:
        d = hdul[1].data
        ok = np.isfinite(d["K"]) & (d["K"] > 0)
        k_med = np.nanmedian(d["K"][ok])
        a_med = np.nanmedian(d["ALPHA"][ok])
        b_med = np.nanmedian(d["BETA"][ok])
        print(f"{g}: k={k_med:.3f}, alpha={a_med:.4f}, beta={b_med:.4f}")
