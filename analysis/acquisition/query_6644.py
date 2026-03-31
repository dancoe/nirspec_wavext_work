import os
from astroquery.mast import Observations
import pandas as pd

pids = ['6644']

print(f"=== Science Configurations in PID {pids[0]} ===")

obs_table = Observations.query_criteria(
    proposal_id=pids[0],
    instrument_name='NIRSPEC/SLIT',
    intentType='science'
)

if len(obs_table) > 0:
    df = obs_table.to_pandas()
    # List unique obs ids and their EXP_TYPE if available
    print(df[['obs_id', 'target_name', 'filters', 't_exptime']].head(20))
