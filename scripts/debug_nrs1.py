import os
import sys

def debug():
    data_root = '../data'
    all_nrs1 = []
    for root, dirs, files in os.walk(data_root):
        for f in files:
            if f.endswith('rate.fits') and 'nrs1' in f: all_nrs1.append(os.path.join(root, f))
    
    print(f"Total NRS1 rate files found: {len(all_nrs1)}")
    
    nrs2_rate_path = '../data/PID1492/jw01492001001_03108_00007_nrs2_rate.fits'
    parts = os.path.basename(nrs2_rate_path).split('_')
    print(f"Parts for {os.path.basename(nrs2_rate_path)}: {parts}")
    
    full_obs_id = parts[0] + '_' + parts[1] + '_' + parts[2]
    print(f"Full Obs ID: {full_obs_id}")
    
    matches = [p for p in all_nrs1 if full_obs_id in os.path.basename(p)]
    print(f"Matches for {full_obs_id}: {matches}")

if __name__ == '__main__':
    debug()
