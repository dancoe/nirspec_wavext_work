import os
import sys
from jwst.assign_wcs import AssignWcsStep
from jwst.extract_2d import Extract2dStep
from jwst import datamodels

# Set environment
os.environ['CRDS_SERVER_URL'] = 'https://jwst-crds.stsci.edu'
os.environ['CRDS_PATH'] = os.path.expanduser('~/crds_cache')

def run():
    data_dir = '../data/PID1492'
    ref_file = '/Users/dcoe/NIRSpec/wavext/jwst_nirspec_wavext/reference_files/wavelengthrange_extended.asdf'
    
    # Process both NRS1 and NRS2 for 01001
    files = sorted([f for f in os.listdir(data_dir) if f.startswith('jw01492001001') and f.endswith('rate.fits')])
    
    # Initialize steps
    aw_step = AssignWcsStep()
    ex_step = Extract2dStep()
    ex_step.wavelengthrange_url = ref_file
    
    for f in files:
        rate_path = os.path.join(data_dir, f)
        ext_path = rate_path.replace('_rate.fits', '_extract_2d.fits')
        
        # We always regenerate for 01001 to ensure WCS is present
        print(f"Processing {f}...")
        try:
            # 1. Assign WCS
            print(f"  Running assign_wcs...")
            model_wcs = aw_step.run(rate_path)
            
            # 2. Extract 2D
            print(f"  Running extract_2d...")
            result = ex_step.run(model_wcs)
            
            # Save
            result.save(ext_path)
            print(f"  Saved {ext_path}")
            
            model_wcs.close()
            result.close()
        except Exception as e:
            print(f"  Error processing {f}: {e}")

if __name__ == '__main__':
    run()
