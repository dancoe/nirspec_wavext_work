import crds
from collections import defaultdict

import os

# Set CRDS paths and server
os.environ["CRDS_SERVER_URL"] = "https://jwst-crds.stsci.edu"
os.environ["CRDS_PATH"] = "/Users/dcoe/crds_cache"

ctx = crds.get_default_context()
rmap = crds.get_cached_mapping(ctx)
imap = rmap.get_imap('nirspec')
sflat_rmap = imap.get_rmap('sflat')

all_refs = sflat_rmap.reference_names()

def get_parkeys(ref_name):
    # Retrieve the match keys for a given reference file name in the rmap.
    for match, ref in sflat_rmap.get_matched_files():
        if ref == ref_name:
            return match
    return None

mos_refs = []
fs_refs = []
ifu_refs = []

# Collect everything
for match, ref in sflat_rmap.get_matched_files():
    # match is a tuple like (EXP_TYPE, DETECTOR, GRATING, FILTER, ...)
    # The order of parkeys for nirspec sflat:
    # parkeys = sflat_rmap.parkey
    meta = dict(zip(sflat_rmap.parkey, match))
    
    exp = meta.get('META.EXPOSURE.TYPE', '').upper()
    det = meta.get('META.INSTRUMENT.DETECTOR', '')
    grat = meta.get('META.INSTRUMENT.GRATING', '')
    filt = meta.get('META.INSTRUMENT.FILTER', '')
    
    # Exclude ambiguous "ANY", focus on actual configurations
    if "ANY" in det or "ANY" in grat or "ANY" in filt:
        continue
    
    info = {"ref": ref, "det": det, "grat": grat, "filt": filt, "exp": exp}
    
    if "MOS" in exp or "MSASPEC" in exp:
        mos_refs.append(info)
    elif "FIXED" in exp:
        fs_refs.append(info)
    elif "IFU" in exp:
        ifu_refs.append(info)

def generate_table(title, refs_list):
    # Group by Grating, Filter, then get NRS1 and NRS2
    table = {}
    for r in refs_list:
        if "PRISM" in r["grat"]:
            key = (r["grat"], "CLEAR")
        else:
            key = (r["grat"], r["filt"])
            
        if key not in table:
            table[key] = {"NRS1": "Missing", "NRS2": "Missing"}
            
        # extract just the number (e.g. jwst_nirspec_sflat_0222.fits -> 0222)
        ref_num = r["ref"].split("_sflat_")[-1].replace(".fits", "")
        
        table[key][r["det"]] = ref_num

    print(f"\n### {title}")
    # Sort
    sorted_keys = sorted(table.keys(), key=lambda x: (x[0], x[1]))
    print("Grating | Filter | NRS1 | NRS2")
    print("------- | ------ | ---- | ----")
    for grat, filt in sorted_keys:
        nrs1 = table[(grat, filt)]["NRS1"]
        nrs2 = table[(grat, filt)]["NRS2"]
        nrs1_fmt = f"`{nrs1}`" if nrs1 != "Missing" else ""
        nrs2_fmt = f"`{nrs2}`" if nrs2 != "Missing" else ""
        print(f"{grat} | {filt} | {nrs1_fmt} | {nrs2_fmt}")
        
    return table

mos_table = generate_table("MOS S-flats", mos_refs)
generate_table("FS S-flats", fs_refs)
generate_table("IFU S-flats", ifu_refs)

