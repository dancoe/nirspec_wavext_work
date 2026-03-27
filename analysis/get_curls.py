from astroquery.mast import Observations
import os

prog_id = '1492'
obs_table = Observations.query_criteria(
    proposal_id=prog_id,
    instrument_name='NIRSPEC/SLIT',
    intentType='science'
)

# Filters for each grating
gratings = ['G140H', 'G235H', 'G395H', 'PRISM']
output_dir = '/Users/dcoe/NIRSpec/wavext/data/PID1492'

for g in gratings:
    mask = [g in str(f) for f in obs_table['filters']]
    g_obs = obs_table[mask & (obs_table['calib_level'] == 3)]
    if len(g_obs) > 0:
        g_obs.sort('t_exptime', reverse=True)
        # Get product list for ONE obs
        products = Observations.get_product_list(g_obs[0])
        # Find x1d and rate
        x1d_p = [p for p in products if p['productFilename'].endswith('x1d.fits')]
        rate_p = [p for p in products if p['productFilename'].endswith('rate.fits')]
        
        # Download x1d if found
        if x1d_p:
            url = f"https://mast.stsci.edu/api/v0.1/Download/file?uri={x1d_p[0]['dataURI']}"
            opath = os.path.join(output_dir, x1d_p[0]['productFilename'])
            print(f"curl -o {opath} \"{url}\"")
        
        # Download one NRS2 rate file if found
        nrs2_rate = [p for p in rate_p if 'nrs2' in p['productFilename']]
        if nrs2_rate:
            url = f"https://mast.stsci.edu/api/v0.1/Download/file?uri={nrs2_rate[0]['dataURI']}"
            opath = os.path.join(output_dir, nrs2_rate[0]['productFilename'])
            print(f"curl -o {opath} \"{url}\"")
