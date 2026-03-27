from astroquery.mast import Observations
import numpy as np

pids = ['1534', '1535', '1536', '1537', '1538']
for pid in pids:
    print(f"\n--- Searching for Program ID {pid} ---")
    obs_table = Observations.query_criteria(proposal_id=pid, instrument_name='NIRSPEC/SLIT')
    targets = set(obs_table['target_name'])
    for t in targets:
        rows = obs_table[obs_table['target_name'] == t]
        filters = " ".join(rows['filters'])
        if all(g in filters for g in ['G140M', 'G235M', 'G395M', 'PRISM']):
            print(f"MATCH: {t}")
        else:
            # Maybe check for just M gratings
            if all(g in filters for g in ['G140M', 'G235M']):
                print(f"Partial Match (M-gratings only): {t}")
