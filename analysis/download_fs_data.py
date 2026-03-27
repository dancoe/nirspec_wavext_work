from astroquery.mast import Observations
import os

prog_id = '1492'
obs_table = Observations.query_criteria(
    proposal_id=prog_id,
    instrument_name='NIRSPEC/SLIT',
    intentType='science'
)

# Subset table for first observation
obs_subset = obs_table[0:1]
products = Observations.get_product_list(obs_subset)
print(f"Products: {len(products)}")

rate_products = [p for p in products if p['productFilename'].endswith('rate.fits')]
print(f"Found {len(rate_products)} rate products.")

output_dir = '/Users/dcoe/NIRSpec/wavext/data/PID1492'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

if rate_products:
     print(f"Downloading {rate_products[0]['productFilename']}...")
     manifest = Observations.download_products(rate_products[0:1], download_dir=output_dir)
     print(f"Manifest: {manifest}")
