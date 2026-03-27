from astroquery.mast import Observations
import os

prog_id = '1492'
obs_table = Observations.query_criteria(
    proposal_id=prog_id,
    instrument_name='NIRSPEC/SLIT',
    intentType='science'
)

output_dir = '/Users/dcoe/NIRSpec/wavext/data/PID1492'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

gratings = ['G140M', 'G235M', 'G395M']
selected_obs = []

for g in gratings:
    g_mask = [g in str(f) for f in obs_table['filters']]
    g_obs = obs_table[g_mask]
    if len(g_obs) > 0:
        g_obs.sort('t_exptime', reverse=True)
        # pick the Level 3 (combined) or longest to get its products
        # actually, just get the products for the longest level 3 exposure
        l3_obs = g_obs[g_obs['calib_level'] == 3]
        if len(l3_obs) > 0:
            selected_obs.append(l3_obs[0])

print(f"Selected {len(selected_obs)} representative observations for download.")

for obs in selected_obs:
    try:
        products = Observations.get_product_list(obs['obsid'])
        # download level 3 x1d.fits
        import numpy as np
        
        # Only download for G140M or G235M rate/x1d files
        x1d_mask = np.char.endswith(products['productFilename'].astype(str), 'x1d.fits')
        x1d_p = products[x1d_mask]
        if len(x1d_p) > 0:
            print(f"Downloading x1d.fits...")
            Observations.download_products(x1d_p, download_dir=output_dir, flat=True)
            
        rate_mask = np.char.endswith(products['productFilename'].astype(str), 'rate.fits')
        rate_p = products[rate_mask]
        if len(rate_p) > 0:
            print(f"Downloading rate.fits files...")
            Observations.download_products(rate_p, download_dir=output_dir, flat=True)  # get all rate files for this obs
    except Exception as e:

        print(f"Error handling {obs['obs_id']}: {e}")
