import os
os.environ["CRDS_SERVER_URL"] = "https://jwst-crds.stsci.edu"
os.environ["CRDS_PATH"] = "/Users/dcoe/crds_cache"

from crds.client.api import get_best_references

combinations = [
    ("PRISM", "CLEAR"),
    ("G140M", "F070LP"),
    ("G140M", "F100LP"),
    ("G140H", "F070LP"),
    ("G140H", "F100LP"),
    ("G235M", "F170LP"),
    ("G235H", "F170LP"),
    ("G395M", "F290LP"),
    ("G395H", "F290LP"),
]

modes = [
    ("MOS", "NRS_MSASPEC"),
    ("FS", "NRS_FIXEDSLIT"),
    ("IFU", "NRS_IFU"),
]

def format_ref(ref_str):
    if not ref_str: return ""
    return f"`{ref_str.split('_sflat_')[-1].replace('.fits','')}`"

with open("/tmp/new_tables.md", "w") as f:
    for mode_name, exp_type in modes:
        f.write(f"### {mode_name} S-flats\n\n")
        f.write("`jwst_nirspec_sflat_0###.fits`\n\n")
        f.write("Grating | Filter | NRS1 | NRS2\n")
        f.write("------- | ------ | ---- | ----\n")
        for grat, filt in combinations:
            res = {}
            for det in ["NRS1", "NRS2"]:
                cfg = {
                    "INSTRUME": "NIRSPEC",
                    "META.INSTRUMENT.NAME": "NIRSPEC",
                    "EXP_TYPE": exp_type,
                    "DETECTOR": det,
                    "GRATING": grat,
                    "FILTER": filt,
                    "DATE-OBS": "2026-03-27",
                    "TIME-OBS": "00:00:00",
                    "META.OBSERVATION.DATE": "2026-03-27",
                    "META.OBSERVATION.TIME": "00:00:00"
                }
                try:
                    refs = get_best_references("jwst-latest", cfg, reftypes=["sflat"])
                    res[det] = refs.get("sflat", "NOT_FOUND")
                except Exception as e:
                    res[det] = "ERROR"
            
            f.write(f"{grat} | {filt} | {format_ref(res['NRS1'])} | {format_ref(res['NRS2'])}\n")
        f.write("\n")
