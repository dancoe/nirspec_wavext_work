from astroquery.mast import Observations
import os

prog_id = '1492'
obs_table = Observations.query_criteria(
    proposal_id=prog_id,
    instrument_name='NIRSPEC/SLIT',
    intentType='science'
)

obs_subset = obs_table[0:1]
products = Observations.get_product_list(obs_subset)
print(f"Total products found: {len(products)}")
for p in products:
    print(f"  {p['productFilename']}")
