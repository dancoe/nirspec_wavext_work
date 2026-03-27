import crds
import os

os.environ["CRDS_SERVER_URL"] = "https://jwst-crds.stsci.edu"
os.environ["CRDS_PATH"] = "/Users/dcoe/crds_cache"

# use crds bestref to find the exact file name
missing_configs = [
    {"EXP_TYPE": "NRS_MSASPEC", "DETECTOR": "NRS1", "GRATING": "G140M", "FILTER": "F070LP", "DATE-OBS": "2026-03-27", "TIME-OBS": "00:00:00", "META.OBSERVATION.DATE": "2026-03-27", "META.OBSERVATION.TIME": "00:00:00"},
    {"EXP_TYPE": "NRS_MSASPEC", "DETECTOR": "NRS1", "GRATING": "G140H", "FILTER": "F070LP", "DATE-OBS": "2026-03-27", "TIME-OBS": "00:00:00", "META.OBSERVATION.DATE": "2026-03-27", "META.OBSERVATION.TIME": "00:00:00"},
    {"EXP_TYPE": "NRS_MSASPEC", "DETECTOR": "NRS1", "GRATING": "G235H", "FILTER": "F170LP", "DATE-OBS": "2026-03-27", "TIME-OBS": "00:00:00", "META.OBSERVATION.DATE": "2026-03-27", "META.OBSERVATION.TIME": "00:00:00"},
]

from crds.client.api import get_best_references

for cfg in missing_configs:
    try:
        refs = get_best_references("jwst-latest", cfg, reftypes=["sflat"])
        print(f"Missing config {cfg['GRATING']} {cfg['FILTER']} {cfg['DETECTOR']} -> {refs['sflat']}")
    except Exception as e:
        print(f"Error for {cfg['GRATING']} {cfg['FILTER']} {cfg['DETECTOR']}: {e}")
