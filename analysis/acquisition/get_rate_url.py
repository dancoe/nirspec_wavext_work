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

rate_products = [p for p in products if p['productFilename'].endswith('rate.fits')]

if rate_products:
    url = f"https://mast.stsci.edu/api/v0.1/Download/file?uri={rate_products[0]['dataURI']}"
    print(f"URL: {url}")
    print(f"Filename: {rate_products[0]['productFilename']}")
else:
    print("No rate products found.")
