from astroquery.mast import Observations
import os
import numpy as np

# We focus on targets that have PRISM, G140M, G235M, G395M
targets_to_fetch = {
    '1536': 'J1743045',
    '1537': 'G191-B2B',
    '1538': 'P330E'
}

gratings = ['PRISM', 'G140M', 'G235M', 'G395M']
# Maps for specific filters to ensure we get the right ones
grating_filters = {
    'PRISM': 'CLEAR',
    'G140M': 'F100LP',
    'G235M': 'F170LP',
    'G395M': 'F290LP'
}

base_data_dir = '/Users/dcoe/NIRSpec/wavext/data'

for prog_id, target in targets_to_fetch.items():
    print(f"\n=== Processing PID {prog_id} Target {target} ===")
    
    output_dir = os.path.join(base_data_dir, f"PID{prog_id}_{target}")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    obs_table = Observations.query_criteria(
        proposal_id=prog_id,
        target_name=target,
        instrument_name='NIRSPEC/SLIT',
        intentType='science'
    )
    
    print(f"Found {len(obs_table)} observations for {target}.")
    
    selected_indices = []
    for g, f in grating_filters.items():
        # Mask for grating and filter
        mask = np.array([g in str(row['filters']) and (f == 'ANY' or f in str(row['filters'])) for row in obs_table])
        if not any(mask):
            continue
            
        g_obs = obs_table[mask]
        # Sort by calib_level to get Level 3 first, then by exposure time
        # We need to find the indices in the original table
        l3_mask = g_obs['calib_level'] == 3
        if any(l3_mask):
             l3_selected = g_obs[l3_mask][0]
             # find original index by obs_id
             orig_idx = np.where(obs_table['obs_id'] == l3_selected['obs_id'])[0][0]
             selected_indices.append(orig_idx)
             print(f"  Selected Level 3 for {g}/{f}: {l3_selected['obs_id']}")
             
        l2_mask = g_obs['calib_level'] == 2
        if any(l2_mask):
             l2_selected = g_obs[l2_mask][0]
             orig_idx = np.where(obs_table['obs_id'] == l2_selected['obs_id'])[0][0]
             selected_indices.append(orig_idx)
             print(f"  Selected Level 2 for {g}/{f}: {l2_selected['obs_id']} (and siblings)")

    for idx in selected_indices:
        # Wrap the row in a 1-row Table to avoid MAST bigint error
        obs_row = obs_table[idx:idx+1]
        print(f"  Downloading products for {obs_row['obs_id'][0]}...")
        try:
            products = Observations.get_product_list(obs_row)
            
            if obs_row['calib_level'][0] == 3:
                # x1d.fits
                x1d_mask = [str(f).endswith('x1d.fits') for f in products['productFilename']]
                Observations.download_products(products[x1d_mask], download_dir=output_dir, flat=True)
            else:
                # rate.fits
                rate_mask = [str(f).endswith('rate.fits') for f in products['productFilename']]
                Observations.download_products(products[rate_mask], download_dir=output_dir, flat=True)
            print(f"    Done.")
        except Exception as e:
            print(f"    Error: {e}")
