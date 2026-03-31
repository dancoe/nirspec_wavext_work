import os
import shutil
import glob

def reorganize_pid_dir(pid_dir):
    print(f"Reorganizing {pid_dir}")
    if not os.path.isdir(pid_dir):
        return

    # Target directories
    dirs = {
        'rate': os.path.join(pid_dir, 'rate'),
        'spec2_nrs1_cal': os.path.join(pid_dir, 'spec2_nrs1_cal'),
        'spec2_nrs2_cal': os.path.join(pid_dir, 'spec2_nrs2_cal'),
        'spec3_nrs1_nom': os.path.join(pid_dir, 'spec3_nrs1_nom'),
        'spec3_nrs2_ext': os.path.join(pid_dir, 'spec3_nrs2_ext'),
    }

    # Ensure target directories exist (they will be created as needed)
    
    # Mapping old to new names for existing subdirectories
    mapping = {
        'nrs1_spec2_cal': 'spec2_nrs1_cal',
        'nrs2_spec2_cal': 'spec2_nrs2_cal',
        'nrs1_spec3_nom': 'spec3_nrs1_nom',
        'nrs2_spec3_ext': 'spec3_nrs2_ext',
    }

    for old_sub, new_sub in mapping.items():
        old_path = os.path.join(pid_dir, old_sub)
        new_path = os.path.join(pid_dir, new_sub)
        if os.path.isdir(old_path):
            if not os.path.exists(new_path):
                os.makedirs(new_path)
            for item in os.listdir(old_path):
                s = os.path.join(old_path, item)
                d = os.path.join(new_path, item)
                shutil.move(s, d)
            os.rmdir(old_path)

    # Now move files from the root of the PID directory
    files = glob.glob(os.path.join(pid_dir, "*"))
    for f in files:
        if os.path.isdir(f):
            continue
        
        name = os.path.basename(f)
        
        # 1. Level 3 (Spec3)
        # Association-level products usually have o001-t001 or similar
        # and contain the grating/filter name in the filename
        if ("-o" in name and "-t" in name) or "f100lp" in name.lower() or "f170lp" in name.lower() or "f290lp" in name.lower() or "clear-prism" in name.lower():
             # We need to distinguish nrs1 vs nrs2 if possible, but spec3 often merges them.
             # However, the user asked for spec3_nrs1_nom and spec3_nrs2_ext.
             # Usually, the ones in the root are the standard MAST L3.
             # Let's put them in spec3_nrs1_nom as "nominal" if they don't specify nrs2.
             dest = dirs['spec3_nrs1_nom']
             if "nrs2" in name.lower():
                 dest = dirs['spec3_nrs2_ext']
             
             if not os.path.exists(dest): os.makedirs(dest)
             shutil.move(f, os.path.join(dest, name))
             continue

        # 2. rate files
        if "_rate.fits" in name:
            dest = dirs['rate']
            if not os.path.exists(dest): os.makedirs(dest)
            shutil.move(f, os.path.join(dest, name))
            continue

        # 3. Individual nrs1/nrs2 files
        if "_nrs1_" in name:
            dest = dirs['spec2_nrs1_cal']
            if not os.path.exists(dest): os.makedirs(dest)
            shutil.move(f, os.path.join(dest, name))
        elif "_nrs2_" in name:
            dest = dirs['spec2_nrs2_cal']
            if not os.path.exists(dest): os.makedirs(dest)
            shutil.move(f, os.path.join(dest, name))

def main():
    root_data = "/Users/dcoe/NIRSpec/wavext/data"
    pids = ["PID1492", "PID1536_J1743045", "PID1537_G191-B2B", "PID1538_P330E", "PID2186_UGC-5101", "PID2186_UGC5101", "PID2654_SDSSJ0749", "PID2654_SDSSJ0841", "PID6644", "PID6644_NGC2506G31"]
    
    for pid in pids:
        pid_dir = os.path.join(root_data, pid)
        if os.path.isdir(pid_dir):
            reorganize_pid_dir(pid_dir)

if __name__ == "__main__":
    main()
