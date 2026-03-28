"""Query MAST for IFU rate files for PIDs 1537 and 1538."""
import os
os.environ['CRDS_PATH'] = os.path.expanduser('~/crds_cache')
os.environ['CRDS_SERVER_URL'] = 'https://jwst-crds.stsci.edu'
from astroquery.mast import Observations
import numpy as np

target_filters = ['F100LP;G140M', 'F170LP;G235M']

for pid, target in [('1537', 'G191-B2B'), ('1538', 'P330E')]:
    print(f'\n=== PID {pid} ({target}) ===')
    obs_table = Observations.query_criteria(
        instrument_name=['NIRSPEC/IFU'],
        provenance_name=['CALJWST'],
        obs_id=[f'*{pid}*'],
        filters=target_filters,
    )
    print(f'  G140M/G235M IFU obs: {len(obs_table)}')

    all_uris = []
    for row in obs_table:
        print(f'  obs: {row["obs_id"]}')
        products = Observations.get_product_list(row)
        rate_files = Observations.filter_products(
            products,
            calib_level=[2],
            productSubGroupDescription='RATE',
        )
        print(f'    rate files: {len(rate_files)}')
        for f in rate_files[:5]:
            print(f'      {f["productFilename"]}  size={f["size"]/1e6:.1f}MB')
        all_uris.extend(rate_files['dataURI'].tolist())
    print(f'  Total rate URIs: {len(all_uris)}')
