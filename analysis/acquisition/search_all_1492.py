from astroquery.mast import Observations

prog_id = '1492'
print(f"Searching for Program ID {prog_id}...")

obs_table = Observations.query_criteria(
    proposal_id=prog_id,
    intentType='science'
)

print(f"Found {len(obs_table)} observations.")
for obs in obs_table:
    print(f"ObsID: {obs['obs_id']}, Instrument: {obs['instrument_name']}, Project: {obs['project']}")
