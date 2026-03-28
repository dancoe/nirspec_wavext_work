import sys
sys.path.append('analysis/util_mos')
from astropy.io import fits
import os
import mos
from jwst import datamodels
from unittest.mock import MagicMock
sys.modules['tqdm'] = MagicMock()
sys.modules['astroquery'] = MagicMock()
sys.modules['astroquery.mast'] = MagicMock()

# Import the actual logic to find files
import analysis.plotting.generate_reference_star_extraction_plots as gen

def test():
    files = gen.find_files()
    pid = 'PID 1492'
    if pid in files:
        gratings = files[pid]
        print(f"Gratings for {pid}: {list(gratings.keys())}")
        for g, d in gratings.items():
            print(f"  {g}: NRS1: {d['nrs1'] != None}, NRS2: {d['nrs2'] != None}")
            print(f"      EXT1: {d['extract_nrs1'] != None}, EXT2: {d['extract_nrs2'] != None}")

if __name__ == '__main__':
    test()
