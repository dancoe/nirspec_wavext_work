import os
import shutil
import glob

# --- Paths ---
LOCAL_DATA = "/Users/dcoe/NIRSpec/wavext/data"
REMOTE_DATA = "/Volumes/wit4/nirspec/dcoe/wavext/data"

FS_TARGETS = [
    "PID1492", "PID1536_J1743045", "PID1537_G191-B2B", "PID1538_P330E",
    "PID6644", "PID6644_NGC2506G31"
]

IFU_TARGETS = [
    "PID1536_J1743045", "PID1537_G191-B2B", "PID1538_P330E", 
    "PID6645_P330E-C3", "PID2186_UGC-5101", "PID2186_UGC5101",
    "PID2654_SDSSJ0749", "PID2654_SDSSJ0841"
]

def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)

def offload_fs():
    ensure_dir(os.path.join(REMOTE_DATA, "FS"))
    for pid in FS_TARGETS:
        src_dir = os.path.join(LOCAL_DATA, pid)
        if not os.path.isdir(src_dir): continue
        
        dest_dir = os.path.join(REMOTE_DATA, "FS", pid)
        print(f"Offloading FS {pid} to central store...")
        
        # 1. Copy everything to remote
        if os.path.exists(dest_dir):
            shutil.rmtree(dest_dir)
        shutil.copytree(src_dir, dest_dir)
        
        # 2. Cleanup local: Keep spec3_nrs1_nom and spec3_nrs2_ext
        # We also need rate files if we want to re-run reductions, 
        # but the user said "keep Level 3 on my laptop". 
        # For now, let's delete rate/ and spec2_ folders locally to save space.
        for sub in ['rate', 'spec2_nrs1_cal', 'spec2_nrs2_cal']:
            local_sub = os.path.join(src_dir, sub)
            if os.path.exists(local_sub):
                print(f"  Removing local {sub}...")
                shutil.rmtree(local_sub)

def offload_ifu():
    ensure_dir(os.path.join(REMOTE_DATA, "IFU"))
    
    # Check the data/IFU subdirectory structure
    local_ifu_root = os.path.join(LOCAL_DATA, "IFU")
    if not os.path.exists(local_ifu_root): return

    for pid in os.listdir(local_ifu_root):
        src_dir = os.path.join(local_ifu_root, pid)
        if not os.path.isdir(src_dir): continue
        
        dest_dir = os.path.join(REMOTE_DATA, "IFU", pid)
        print(f"Offloading IFU {pid} to central store...")
        
        # 1. Copy everything to remote
        if os.path.exists(dest_dir):
            shutil.rmtree(dest_dir)
        shutil.copytree(src_dir, dest_dir)
        
        # 2. Cleanup local: Keep stage3 and stage3_ext
        # Delete stage1, stage2, stage2_ext, etc.
        for sub in ['stage1', 'stage2', 'stage2_ext', 'asn', 'asn_ext']:
            local_sub = os.path.join(src_dir, sub)
            if os.path.exists(local_sub):
                print(f"  Removing local {sub}...")
                shutil.rmtree(local_sub)

if __name__ == "__main__":
    if not os.path.exists("/Volumes/wit4"):
        print("ERROR: /Volumes/wit4 not mounted.")
    else:
        offload_fs()
        offload_ifu()
        print("Offload complete.")
