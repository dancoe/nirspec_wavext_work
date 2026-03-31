from astroquery.mast import Observations
import os
import numpy as np

# v5 Programs
programs = [
    {'pid': '6644', 'targets': ['NGC2506G31'], 'mode': 'NIRSPEC/SLIT'},
    {'pid': '2654', 'targets': ['SDSSJ0749', 'SDSSJ0841'], 'mode': 'NIRSPEC/IFU'},
    {'pid': '2186', 'targets': ['UGC5101'], 'mode': 'NIRSPEC/IFU'},
]

# We want specific configurations
# Format: (grating, filter)
configs = {
    'G140M': 'F100LP',
    'G235M': 'F170LP',
    'G395M': 'F290LP',
    'PRISM': 'CLEAR'
}

base_dir = '/Users/dcoe/NIRSpec/wavext/data'

def download_program(pid, targets, mode):
    print(f"\n=== Querying PID {pid} [{mode}] ===")
    
    obs = Observations.query_criteria(
        proposal_id=pid,
        instrument_name=mode,
        intentType='science'
    )
    
    if len(obs) == 0:
        print(f"  No science observations found for PID {pid}.")
        return

    for target in targets:
        print(f"  Filtering for target {target}...")
        # fuzzy match target
        t_mask = [target.lower() in str(t).lower().replace('-', '').replace(' ', '') for t in obs['target_name']]
        t_obs = obs[t_mask]
        
        if len(t_obs) == 0:
            print(f"    Target {target} not found in this program.")
            continue
            
        print(f"    Found {len(t_obs)} total observations for {target}.")
        
        output_dir = os.path.join(base_dir, f"PID{pid}_{target}")
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # Get all products for these observations once
        print(f"    Listing products...")
        products = Observations.get_product_list(t_obs)
        
        # Filter for rate files across the configurations of interest
        # We need rate files that match the (grating/filter) pairs AND are science (not TA)
        rate_mask = [
            str(f).endswith('rate.fits') and 
            'target-acquisition' not in str(f).lower() and
            'nrs_taconfirm' not in str(f).lower()
            for f in products['productFilename']
        ]
        
        all_rate_products = products[rate_mask]
        print(f"    Found {len(all_rate_products)} total science rate files.")
        
        if len(all_rate_products) == 0:
            continue
            
        # Download in batches
        print(f"    Starting batch download to {output_dir}...")
        for i in range(0, len(all_rate_products), 10):
            chunk = all_rate_products[i:i+10]
            Observations.download_products(chunk, download_dir=output_dir, flat=True)

for p in programs:
    download_program(p['pid'], p['targets'], p['mode'])
