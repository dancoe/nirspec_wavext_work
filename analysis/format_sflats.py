import glob
from astropy.io import fits
import os

files = glob.glob('/Users/dcoe/crds_cache/references/jwst/nirspec/jwst_nirspec_sflat_*.fits')

mos_files = []
fs_files = []
ifu_files = []

for f in sorted(files):
    fname = os.path.basename(f)
    try:
        with fits.open(f) as h:
            exp = h[0].header.get('EXP_TYPE', '')
            det = h[0].header.get('DETECTOR', '')
            grat = h[0].header.get('GRATING', '')
            filt = h[0].header.get('FILTER', '')
            info = f"- `{fname}`: {det} | {grat} | {filt}"
            if 'MOS' in exp or 'MSASPEC' in exp:
                mos_files.append(info)
            elif 'FIXED' in exp:
                fs_files.append(info)
            elif 'IFU' in exp:
                ifu_files.append(info)
            else:
                pass
    except:
        pass

with open('/tmp/flat_list.md', 'w') as out:
    out.write("### MOS S-flats\n")
    for m in mos_files: out.write(m + "\n")
    out.write("\n### FS S-flats\n")
    for m in fs_files: out.write(m + "\n")
    out.write("\n### IFU S-flats\n")
    for m in ifu_files: out.write(m + "\n")
