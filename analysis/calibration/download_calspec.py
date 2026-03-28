"""
Download CALSPEC stellar spectral energy distribution models for the 4 NIRSpec
IFU calibration standards used by Parlanti et al. 2025.

Sources:
    G191-B2B  — hot white dwarf (PID 1537)
    P330E     — gray A-type star (PIDs 1538, 6645)
    J1743-045 — solar analog (PID 1536)

STScI CALSPEC archive:
    https://archive.stsci.edu/hlsps/reference-atlases/cdbs/current_calspec/

Output: data/CALSPEC/*.fits
"""
import os
import urllib.request
import urllib.error
import astropy.io.fits as fits
import numpy as np

CALSPEC_BASE = (
    "https://archive.stsci.edu/hlsps/reference-atlases/cdbs/current_calspec/"
)
# Fallback URL for older/archived files
CALSPEC_BASE_FULL = (
    "https://archive.stsci.edu/hlsps/reference-atlases/cdbs/calspec/"
)

OUTPUT_DIR = "/Users/dcoe/NIRSpec/wavext/data/CALSPEC"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Map from our target labels to CALSPEC filenames.
# We list multiple candidate filenames (newest first) and download the first
# one that exists.
TARGETS = {
    "G191-B2B": {
        # From CALSPEC table 1a: prefix g191b2b, mod_012, stiswfcnic_004
        "candidates": [
            "g191b2b_mod_012.fits",
            "g191b2b_stiswfcnic_004.fits",
            "g191b2b_mod_010.fits",
            "g191b2b_stisnic_003.fits",
        ],
        "pids": ["1537"],
    },
    "P330E": {
        # From CALSPEC table 1a: prefix p330e, mod_008, stiswfcnic_007
        "candidates": [
            "p330e_mod_008.fits",
            "p330e_stiswfcnic_007.fits",
            "p330e_mod_003.fits",
            "p330e_stisnic_010.fits",
        ],
        "pids": ["1538", "6645"],
    },
    "J1743045": {
        # From CALSPEC table 1a: prefix 1743045, mod_007, stisnic_009
        "candidates": [
            "1743045_mod_007.fits",
            "1743045_stisnic_009.fits",
            "j1743045_mod_007.fits",
            "j1743045_stisnic_009.fits",
        ],
        "pids": ["1536"],
    },
}


def try_download(url, outfile):
    """Attempt to download url → outfile. Returns True on success."""
    try:
        print(f"    Trying {url.split('/')[-1]} ... ", end="", flush=True)
        urllib.request.urlretrieve(url, outfile)
        size = os.path.getsize(outfile)
        if size < 1000:
            # May be an HTML error page
            os.remove(outfile)
            print("FAILED (too small)")
            return False
        print(f"OK ({size/1024:.0f} KB)")
        return True
    except urllib.error.HTTPError as e:
        print(f"HTTP {e.code}")
        return False
    except Exception as e:
        print(f"ERROR: {e}")
        return False


def download_calspec():
    results = {}
    for target, info in TARGETS.items():
        print(f"\n{target} (PIDs: {', '.join(info['pids'])})")
        # Check if already downloaded
        existing = [f for c in info["candidates"]
                    if os.path.exists(os.path.join(OUTPUT_DIR, c))
                    for f in [os.path.join(OUTPUT_DIR, c)]]
        if existing:
            print(f"  Already present: {os.path.basename(existing[0])}")
            results[target] = existing[0]
            continue
        # Try each candidate filename
        success = False
        for candidate in info["candidates"]:
            outfile = os.path.join(OUTPUT_DIR, candidate)
            # Try current_calspec first, then full archive
            for base in [CALSPEC_BASE, CALSPEC_BASE_FULL]:
                url = base + candidate
                if try_download(url, outfile):
                    results[target] = outfile
                    success = True
                    break
            if success:
                break
        if not success:
            print(f"  WARNING: Could not download CALSPEC for {target}")
            results[target] = None

    return results


def verify_calspec(results):
    """Print a summary of the downloaded CALSPEC files."""
    print("\n=== CALSPEC Summary ===")
    for target, path in results.items():
        if path is None:
            print(f"  {target}: MISSING")
            continue
        try:
            h = fits.open(path)
            d = h[1].data
            wl = d['WAVELENGTH']
            flux = d['FLUX']
            print(f"  {target}: {os.path.basename(path)}")
            print(f"    Wavelength: {wl.min()/1e4:.3f} - {wl.max()/1e4:.3f} µm "
                  f"(n={len(wl)})")
            print(f"    Flux [flam]: {flux.min():.3e} - {flux.max():.3e}")
            h.close()
        except Exception as e:
            print(f"  {target}: {path} — ERROR: {e}")


if __name__ == "__main__":
    print(f"Downloading CALSPEC models to {OUTPUT_DIR}")
    print(f"Base URL: {CALSPEC_BASE}\n")
    results = download_calspec()
    verify_calspec(results)
    print("\nDone.")
