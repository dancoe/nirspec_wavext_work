"""
Audit local data downloads against MAST archive expectations.
Checks counts of *_rate.fits files for science and calibration PIDs.
"""
from astroquery.mast import Observations
import os
import glob

# PIDs and targets from v5 and general calibration
PROGRAMS = [
    {'pid': '2654', 'targets': ['SDSSJ0749', 'SDSSJ0841'], 'mode': 'NIRSPEC/IFU'},
    {'pid': '2186', 'targets': ['UGC-5101'], 'mode': 'NIRSPEC/IFU'},
    {'pid': '1536', 'targets': ['J1743045'], 'mode': 'NIRSPEC/IFU'},
    {'pid': '1537', 'targets': ['G191-B2B'], 'mode': 'NIRSPEC/IFU'},
    {'pid': '1538', 'targets': ['P330E'], 'mode': 'NIRSPEC/IFU'},
    {'pid': '6645', 'targets': ['P330E-C3'], 'mode': 'NIRSPEC/IFU'},
    {'pid': '6644', 'targets': ['NGC2506G31'], 'mode': 'NIRSPEC/SLIT'},
]

BASE_DIR = '/Users/dcoe/NIRSpec/wavext/data'

def get_expected_count(pid, target_name, mode):
    """Query MAST for expected science rate file count."""
    try:
        obs = Observations.query_criteria(
            proposal_id=pid,
            target_name=f'*{target_name}*',
            instrument_name=mode,
            intentType='science'
        )
        if len(obs) == 0:
            return 0
        products = Observations.get_product_list(obs)
        rate_mask = [
            str(f).endswith('rate.fits') and 
            'target-acquisition' not in str(f).lower() and
            'nrs_taconfirm' not in str(f).lower()
            for f in products['productFilename']
        ]
        return len(products[rate_mask])
    except Exception as e:
        print(f"  Error querying MAST for PID {pid}: {e}")
        return -1

def get_local_count(pid, target_name):
    """Count unique rate files on disk across potential directory variations."""
    # Handle dir name variations (hyphens, underscores)
    search_patterns = [
        os.path.join(BASE_DIR, f"PID{pid}_{target_name}", "**", "*rate.fits"),
        os.path.join(BASE_DIR, f"PID{pid}_{target_name.replace('-', '')}", "**", "*rate.fits"),
        os.path.join(BASE_DIR, "IFU", f"PID{pid}_{target_name}", "**", "*rate.fits"),
        os.path.join(BASE_DIR, "IFU", f"PID{pid}_{target_name.replace('-', '')}", "**", "*rate.fits"),
    ]
    
    all_files = []
    for pattern in search_patterns:
        all_files.extend(glob.glob(pattern, recursive=True))
    
    # Return count of unique basenames
    unique_files = set([os.path.basename(f) for f in all_files])
    return len(unique_files)

if __name__ == "__main__":
    print(f"{'PID':<6} | {'Target':<15} | {'MAST':<5} | {'Local':<5} | {'Status'}")
    print("-" * 55)
    
    for p in PROGRAMS:
        pid = p['pid']
        for target in p['targets']:
            mast = get_expected_count(pid, target, p['mode'])
            local = get_local_count(pid, target)
            
            status = "✅ OK" if local >= mast and mast > 0 else "⏳ Partial"
            if mast == 0: status = "❓ No Obs"
            if mast == -1: status = "❌ Error"
            
            print(f"{pid:<6} | {target:<15} | {mast:<5} | {local:<5} | {status}")
