import os
import numpy as np
from stdatamodels.jwst import datamodels
from jwst.assign_wcs.assign_wcs_step import AssignWcsStep

files = [
    '/Users/dcoe/NIRSpec/wavext/data/PID1492/jw01492001001_03102_00001_nrs2_rate.fits', # G140M
    '/Users/dcoe/NIRSpec/wavext/data/PID1492/jw01492001001_03104_00001_nrs2_rate.fits', # G235M
    '/Users/dcoe/NIRSpec/wavext/data/PID1492/jw01492001001_03106_00001_nrs2_rate.fits', # G395M
    '/Users/dcoe/NIRSpec/wavext/data/PID1492/jw01492001001_03108_00007_nrs2_rate.fits', # G140H
    '/Users/dcoe/NIRSpec/wavext/data/PID1492/jw01492001001_0310a_00007_nrs2_rate.fits', # G235H
]
override = '/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/reference_files/wavelengthrange_extended.asdf'

for f in files:
    try:
        # We don't have all these downloaded in our data directory yet, except maybe G140H and G235H.
        # Oh wait, are they downloaded? Let's check if they exist.
        if not os.path.exists(f):
            continue
        model = datamodels.open(f)
        
        # assign wcs
        step = AssignWcsStep()
        step.override_wavelengthrange = override
        res = step.process(model)
        
        slits_iter = getattr(res, 'slits', [res])
        
        # Check S200A1
        for slit in slits_iter:
            slit_name = getattr(slit, 'name', 'UNKNOWN')
            if slit_name in ['S200A1', 'S1600A1']:
                wcs = slit.meta.wcs
                bbox = wcs.bounding_box
                if bbox is not None:
                    try:
                        # CompoundBoundingBox might have a dict structure
                        bbox_x = bbox.bounding_boxes[slit.name].intervals['x']
                        nx = int(bbox_x.upper)
                    except:
                        # Single BoundingBox
                        nx = int(bbox[0][1])
                else:
                    nx = 2047
                
                # Check center of slit
                x = np.arange(nx + 1)
                y = np.ones(nx + 1) * 28 # roughly the center vertically
                _, _, wav = wcs(x, y)
                wav = wav[wav > 0]
                if len(wav) > 0:
                    print(f'{os.path.basename(f)} ({slit.name}): {slit.meta.instrument.grating} {slit.meta.instrument.filter} NRS2 Max X = {nx}, Max Wavelength = {np.nanmax(wav):.3f} um')

    except Exception as e:
        print(f'{f} error: {e}')
