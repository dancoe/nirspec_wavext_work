import crds
import sys

# We can query crds context directly via command line wrapped in python
import subprocess

out = subprocess.check_output(["crds", "list", "--latest-context", "--file-properties", "EXP_TYPE,DETECTOR,GRATING,FILTER"]).decode()
with open("/tmp/crds_all_files.txt", "w") as f:
    f.write(out)

# Parse only nirspec sflats
sflats = []
for line in out.splitlines():
    if "jwst_nirspec_sflat" in line:
        sflats.append(line)

with open("/tmp/crds_sflats.txt", "w") as f:
    for s in sflats:
        f.write(s + "\n")
print(f"Found {len(sflats)} sflats.")
