import os
import glob
import shutil

root = "/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/plots/flats"

# Map of mode to subdirectory
mappings = {
    "mos_prism": "s-flats/mos/prism",
    "mos_g140h": "s-flats/mos/g140h",
    "mos_g140m": "s-flats/mos/g140m",
    "mos_g235m": "s-flats/mos/g235m",
    "mos_g395m": "s-flats/mos/g395m",
    "fs_s-flat": "s-flats/fs",
    "ifu_s-flat": "s-flats/ifu",
    "d-flat": "d-flats"
}

# Recursively find all PNG files
for src in glob.glob(root + "/**/*.png", recursive=True):
    fname = os.path.basename(src).lower()
    moved = False
    for key, sub in mappings.items():
        if key in fname:
            target_dir = os.path.join(root, sub)
            os.makedirs(target_dir, exist_ok=True)
            dst = os.path.join(target_dir, os.path.basename(src))
            if src != dst:
                shutil.move(src, dst)
            moved = True
            break
    
# Clean up empty directories
for dirpath, dirnames, filenames in os.walk(root, topdown=False):
    if not dirnames and not filenames:
        os.rmdir(dirpath)
