"""Debug MAST download - check what rate files look like."""
import os
os.environ['CRDS_PATH'] = os.path.expanduser('~/crds_cache')
os.environ['CRDS_SERVER_URL'] = 'https://jwst-crds.stsci.edu'
from astroquery.mast import Observations
import numpy as np

pid, target = '1538', 'P330E'

obs_table = Observations.query_criteria(
    instrument_name=['NIRSPEC/IFU'],
    provenance_name=['CALJWST'],
    obs_id=[f'*{pid}*'],
    filters=['F100LP;G140M', 'F170LP;G235M'],
)
print(f'Obs found: {len(obs_table)}')

row = obs_table[0]
print(f'obs_id: {row["obs_id"]}')
products = Observations.get_product_list(row)
print(f'Total products: {len(products)}')
print('Columns:', products.colnames)

rate_files = Observations.filter_products(
    products,
    calib_level=[2],
    productSubGroupDescription='RATE',
)
print(f'\nRate files (calib_level=2): {len(rate_files)}')
if len(rate_files):
    for f in rate_files:
        print(f'  {f["productFilename"]}  size={f["size"]}  calib={f["calib_level"]}')

# Also try calib_level 1
rate1 = Observations.filter_products(
    products,
    calib_level=[1],
    productSubGroupDescription='RATE',
)
print(f'\nRate files (calib_level=1): {len(rate1)}')
if len(rate1):
    for f in rate1[:5]:
        print(f'  {f["productFilename"]}  size={f["size"]}  calib={f["calib_level"]}')

# Show all product types
print('\nAll product types:')
from collections import Counter
types = Counter(zip(products['calib_level'], products['productSubGroupDescription']))
for k, v in sorted(types.items()):
    print(f'  calib={k[0]}  type={k[1]}  count={v}')
