from astroquery.mast import Observations
import os
import numpy as np

programs = ['1536', '1537', '1538']
gratings = ['PRISM', 'G140M', 'G235M', 'G395M']

for prog_id in programs:
    print(f"\n--- Searching for Program ID {prog_id} NIRSpec FS data ---")
    try:
        obs_table = Observations.query_criteria(
            proposal_id=prog_id,
            instrument_name='NIRSPEC/SLIT',
            intentType='science'
        )
        print(f"Found {len(obs_table)} observations.")
        if len(obs_table) > 0:
            for obs in obs_table:
                filters = str(obs['filters'])
                target = obs['target_name']
                obs_id = obs['obs_id']
                level = obs['calib_level']
                if any(g in filters for g in gratings):
                    print(f"ObsID: {obs_id}, Target: {target}, Filters: {filters}, Level: {level}")
    except Exception as e:
        print(f"Error searching program {prog_id}: {e}")
