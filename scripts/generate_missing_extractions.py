import os
import glob
from jwst.assign_wcs import AssignWcsStep
from jwst.extract_2d import Extract2dStep
from jwst import datamodels

def process_nrs1(pid_dir, override_file):
    print(f"Processing NRS1 in {pid_dir}...")
    # Find all nrs1_rate.fits
    rate_files = glob.glob(os.path.join(pid_dir, "*nrs1*rate.fits"))
    for rf in rate_files:
        ext_file = rf.replace("_rate.fits", "_extract_2d.fits")
        if os.path.exists(ext_file):
            print(f"  {ext_file} already exists, skipping.")
            continue
        
        print(f"  Processing {rf}...")
        try:
            # Step 1: Assign WCS with extended override
            wcs_out = AssignWcsStep.call(rf, override_wavelengthrange=override_file)
            
            # Step 2: Extract 2D
            ext_out = Extract2dStep.call(wcs_out)
            
            # Save the result
            ext_out.save(ext_file)
            print(f"  Saved {ext_file}")
            
        except Exception as e:
            print(f"  Error processing {rf}: {e}")

if __name__ == "__main__":
    override = "/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/reference_files/wavelengthrange_extended.asdf"
    
    # Process PIDs with known nrs1_rate files
    # PID 1492
    process_nrs1("../data/PID1492", override)
    
    # Check others
    process_nrs1("../data/PID1537_G191-B2B", override)
    process_nrs1("../data/PID1538_P330E", override)
