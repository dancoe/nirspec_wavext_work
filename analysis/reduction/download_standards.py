import os
from astroquery.mast import Observations

PIDS = ['1536', '1537', '1538']
BASE_DIR = '/Users/dcoe/NIRSpec/wavext/data'

def download_pid_rate(pid):
    print(f"\n--- Downloading PID {pid} ---")
    # Search for the program
    obs_table = Observations.query_criteria(proposal_id=pid, instrument_name='NIRSPEC')
    
    # Filter for FS and Fixed Slit
    # S1600A1 is usually the central slit
    mask = (obs_table['intentType'] == 'science') 
    science_obs = obs_table[mask]
    
    if len(science_obs) == 0:
        print(f"No science observations found for PID {pid}")
        return

    # Filter for rate files
    data_products = Observations.get_product_list(science_obs)
    # We want rate.fits
    mask_rate = (data_products['productSubGroupDescription'] == 'RATE') & \
                (data_products['extension'] == 'fits')
    
    rate_products = data_products[mask_rate]
    print(f"  Found {len(rate_products)} rate files.")
    
    if len(rate_products) == 0:
        return
        
    # Download to data/PIDxxxx_TargetName
    # We find the target name from the first obs
    target_name = science_obs['target_name'][0]
    out_dir = os.path.join(BASE_DIR, f'PID{pid}_{target_name}')
    os.makedirs(out_dir, exist_ok=True)
    
    # We only want nrs1 and nrs2 rate files for actual slits
    manifest = Observations.download_products(rate_products, 
                                             download_dir=out_dir, 
                                             product_type="SCIENCE",
                                             extension="fits")
    
    # Move files from mastDownload subdirs to out_dir
    import shutil
    for row in manifest:
        src = row['Local Path']
        dest = os.path.join(out_dir, os.path.basename(src))
        if os.path.exists(src):
            shutil.move(src, dest)
            print(f"  Moved {os.path.basename(dest)}")
            
    # Clean up mastDownload
    mast_dir = os.path.join(out_dir, 'mastDownload')
    if os.path.exists(mast_dir):
        shutil.rmtree(mast_dir)

if __name__ == '__main__':
    for pid in PIDS:
        download_pid_rate(pid)
