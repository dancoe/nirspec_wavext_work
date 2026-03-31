import os
import glob
import json
from jwst.pipeline import Spec3Pipeline
from jwst.associations import asn_from_list as afl
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base
from astropy.io import fits

def run_spec3_1492():
    cwd = os.getcwd()
    cal_dir = os.path.join(cwd, 'data/FS/PID1492/spec2_full_cal')
    out_dir = os.path.join(cwd, 'data/FS/PID1492/spec3_full_cal')
    os.makedirs(out_dir, exist_ok=True)
    
    cal_files = glob.glob(os.path.join(cal_dir, '*_cal.fits'))
    print(f"Found {len(cal_files)} calibrated files.")
    
    # Group files by (DETECTOR, GRATING, FILTER, FXD_SLIT)
    groups = {}
    for f in cal_files:
        with fits.open(f) as hd:
            hdr = hd[0].header
            det = hdr.get('DETECTOR')
            grat = hdr.get('GRATING')
            filt = hdr.get('FILTER')
            slit = hdr.get('FXD_SLIT')
            
            key = (det, grat, filt, slit)
            if key not in groups:
                groups[key] = []
            groups[key].append(f)
            
    print(f"Found {len(groups)} unique configurations.")
    
    # Process each group
    for (det, grat, filt, slit), files in groups.items():
        prod_name = f"jw01492_{grat}_{filt}_{slit}_{det}".lower()
        asn_path = os.path.join(out_dir, f"{prod_name}_asn.json")
        
        print(f"Processing {prod_name} with {len(files)} files...")
        
        # Create ASN
        asn = afl.asn_from_list(files, rule=DMS_Level3_Base, product_name=prod_name)
        with open(asn_path, 'w') as f:
            f.write(asn.dump()[1])
            
        # Run Spec3
        try:
            # Configure Step
            spec3 = Spec3Pipeline()
            spec3.output_dir = out_dir
            spec3.save_results = True
            
            # Skip outlier detection if few files or just for speed now? 
            # No, let's run it properly.
            
            result = spec3.run(asn_path)
            print(f"  Finished {prod_name}")
        except Exception as e:
            print(f"  Error processing {prod_name}: {e}")

if __name__ == '__main__':
    run_spec3_1492()
