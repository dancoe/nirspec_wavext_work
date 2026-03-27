from astroquery.mast import Observations
import os
import numpy as np

prog_id = '1492'
obs_table = Observations.query_criteria(
    proposal_id=prog_id,
    instrument_name='NIRSPEC/SLIT',
    intentType='science'
)

output_dir = '/Users/dcoe/NIRSpec/wavext/data/PID1492'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Gratings found in DATA.md: G140H, G235H, G395H, PRISM, G140M, G235M, G395M
gratings = ['G140H', 'G235H', 'G395H', 'PRISM', 'G140M', 'G235M', 'G395M']
selected_obs = []

for g in gratings:
    # Use a safe list comprehension for filtering columns
    g_mask = [g in str(f) for f in obs_table['filters']]
    g_obs = obs_table[g_mask]
    if len(g_obs) > 0:
        # Pick the one with longest exp time for SN
        g_obs.sort('t_exptime', reverse=True)
        selected_obs.append(g_obs[0])

print(f"Selected {len(selected_obs)} representative observations for download.")

for obs in selected_obs:
    print(f"Processing obs_id: {obs['obs_id']} ({obs['filters']})")
    try:
         # Use ROW directly to get product list (worked in -c test)
         products = Observations.get_product_list(obs)
         
         # 1. DOWNLOAD the x1d.fits if calib_level is 3
         if obs['calib_level'] == 3:
              x1d_p = [p for p in products if p['productFilename'].endswith('x1d.fits')]
              if x1d_p:
                   Observations.download_products(x1d_p[0:1], download_dir=output_dir, flat=True)
         
         # 2. ALSO DOWNLOAD the rate.fits for corresponding exposures
         # This might be in a different obs entry since obs_table rows vary by level
         # Let's just find the rate.fits in the current product list (if any)
         rate_p = [p for p in products if p['productFilename'].endswith('rate.fits')]
         if rate_p:
              # Download only first 2 rate files (NRS1, NRS2)
               print(f"  Downloading Level 2 rate.fits: {[rp['productFilename'] for rp in rate_p[:2]]}")
               Observations.download_products(rate_p[:2], download_dir=output_dir, flat=True)
              
    except Exception as e:
         print(f"  Error handling {obs['obs_id']}: {e}")
