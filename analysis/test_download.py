from astroquery.mast import Observations
import os

# Test for one observation in PID 1537 (G191-B2B)
prog_id = '1537'
target = 'G191-B2B'
obs_id = 'jw01537009001_03102_00003_nrs2'

print(f"Searching for {obs_id} in PID {prog_id}...")
obs_table = Observations.query_criteria(
    proposal_id=prog_id,
    target_name=target,
    instrument_name='NIRSPEC/SLIT',
    intentType='science'
)

# Find matching row
mask = obs_table['obs_id'] == obs_id
if any(mask):
    obs = obs_table[mask]
    print(f"Found observation row.")
    products = Observations.get_product_list(obs)
    print(f"Found {len(products)} products.")

    if len(products) > 0:
        rate_mask = [str(f).endswith('rate.fits') for f in products['productFilename']]
        rate_products = products[rate_mask]
        print(f"Found {len(rate_products)} rate products.")
        for p in rate_products:
            print(f"  {p['productFilename']}")

        if len(rate_products) > 0:
            output_dir = '/Users/dcoe/NIRSpec/wavext/data/test_download'
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            
            # Use only first 2 for test
            to_download = rate_products[:2]
            print(f"Downloading {len(to_download)} products...")
            # Use clear=True to avoid the TESS FFI check if possible? 
            # Actually download_products is a method of Observations
            Observations.download_products(to_download, download_dir=output_dir, flat=True)
            print(f"Download complete.")
    else:
        print(f"No products found.")
else:
    print(f"Observation ID {obsid} not found in search results.")
