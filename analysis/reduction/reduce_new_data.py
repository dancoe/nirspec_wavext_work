import os
import glob
from jwst.assign_wcs import AssignWcsStep
from jwst.extract_2d import Extract2dStep
from jwst.extract_1d import Extract1dStep
from stdatamodels.jwst import datamodels

# Paths
base_data_dir = '/Users/dcoe/NIRSpec/wavext/data'
override_path = os.path.abspath('/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/reference_files/wavelengthrange_extended.asdf')

# Targets downloaded
targets = ['PID1537_G191-B2B', 'PID1538_P330E']

for target_dir in targets:
    tpath = os.path.join(base_data_dir, target_dir)
    if not os.path.exists(tpath):
        print(f"Target directory {tpath} not found. Skipping.")
        continue
    
    print(f"\n--- Reducing data in {target_dir} ---")
    
    # We only need to process NRS2 for the extension, but let's do both for completeness if they are there
    rate_files = glob.glob(os.path.join(tpath, '*_rate.fits'))
    
    for rate_file in sorted(rate_files):
        print(f"\nProcessing {os.path.basename(rate_file)}")
        try:
            # 1. Assign WCS
            wcs_output = rate_file.replace('_rate.fits', '_wcs.fits')
            if not os.path.exists(wcs_output):
                print(f"  Step 1: AssignWcsStep...")
                model = datamodels.open(rate_file)
                step_wcs = AssignWcsStep()
                step_wcs.override_wavelengthrange = override_path
                res_wcs = step_wcs.run(model)
                if res_wcs:
                    res_wcs.save(wcs_output)
                    print(f"  Saved WCS to {os.path.basename(wcs_output)}")
                model.close()
            else:
                print(f"  WCS already exists: {os.path.basename(wcs_output)}")

            # 2. Extract 2D
            e2d_output = wcs_output.replace('_wcs.fits', '_extract_2d.fits')
            if not os.path.exists(e2d_output):
                print(f"  Step 2: Extract2dStep...")
                model_wcs = datamodels.open(wcs_output)
                step_e2d = Extract2dStep()
                step_e2d.override_wavelengthrange = override_path
                res_e2d = step_e2d.run(model_wcs)
                if res_e2d:
                    res_e2d.save(e2d_output)
                    print(f"  Saved 2D to {os.path.basename(e2d_output)}")
                    # Quick check of wavelengths
                    if hasattr(res_e2d, 'slits'):
                        for slit in res_e2d.slits:
                             if slit.wavelength is not None:
                                  print(f"    Slit {slit.name} wave range: {slit.wavelength.min():.4f}-{slit.wavelength.max():.4f} um")
                model_wcs.close()
            else:
                print(f"  2D already exists: {os.path.basename(e2d_output)}")

            # 3. Extract 1D
            x1d_output = e2d_output.replace('_extract_2d.fits', '_extract_1d.fits')
            if not os.path.exists(x1d_output):
                print(f"  Step 3: Extract1dStep...")
                model_e2d = datamodels.open(e2d_output)
                step_x1d = Extract1dStep()
                # step_x1d.override_wavelengthrange = override_path # Not strictly needed as 2D has it
                res_x1d = step_x1d.run(model_e2d)
                if res_x1d:
                    res_x1d.save(x1d_output)
                    print(f"  Saved 1D to {os.path.basename(x1d_output)}")
                model_e2d.close()
            else:
                print(f"  1D already exists: {os.path.basename(x1d_output)}")
        except Exception as e:
            print(f"  ERROR processing {rate_file}: {e}")

print("\nReduction complete.")
