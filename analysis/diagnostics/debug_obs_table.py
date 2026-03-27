from astroquery.mast import Observations
import os

prog_id = '1492'
print(f"Searching for Program ID {prog_id} NIRSPEC/SLIT data...")

obs_table = Observations.query_criteria(
    proposal_id=prog_id,
    instrument_name='NIRSPEC/SLIT',
    intentType='science'
)

print(f"Columns: {obs_table.colnames}")
if 'obsid' not in obs_table.colnames and 'obs_id' in obs_table.colnames:
    # Some versions might use different keys or the table has different column names
    print("No 'obsid' column, checking for 'obs_id' or others...")

# Re-run search for first few to minimize
obs_table = obs_table[:2]

products = Observations.get_product_list(obs_table)
# ... filter ...
