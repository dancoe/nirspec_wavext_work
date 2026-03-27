from astroquery.mast import Observations

prog_id = '1492'
print(f"Searching for Program ID {prog_id} NIRSpec data...")

obs_table = Observations.query_criteria(
    proposal_id=prog_id,
    instrument_name='NIRSPEC',
    intentType='science'
)

print(f"Found {len(obs_table)} observations.")
if len(obs_table) > 0:
    for obs in obs_table:
        print(f"ObsID: {obs['obs_id']}, Target: {obs['target_name']}, Type: {obs['obs_collection']}")

    products = Observations.get_product_list(obs_table)
    filtered_products = Observations.filter_products(products, productType="SCIENCE", extension="fits")
    print(f"Products to download ({len(filtered_products)}):")
    for p in filtered_products[:10]: # Print first 10
         print(f"Product: {p['productFilename']}")
