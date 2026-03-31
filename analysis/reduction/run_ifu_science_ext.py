"""
Extended-range IFU pipeline for science targets (PIDs 2654, 2186).

Applies Parlanti et al. (2025) reference files to recover extended wavelengths:
  G140M/F100LP:  ~1.88 – 3.57 µm  (nominal NRS2 coverage begins here)
  G235M/F170LP:  ~3.15 – 5.27 µm

Targets:
  PID 2654: SDSSJ0749, SDSSJ0841 (AGN; z~2; G140M validation — H-alpha recovery)
  PID 2186: UGC-5101 (ULIRG; G235M extended continuum validation)

Usage:
    micromamba run -n jwst_1.20.2 python run_ifu_science_ext.py --pid 2654 --target SDSSJ0841
    micromamba run -n jwst_1.20.2 python run_ifu_science_ext.py --pid 2186 --target UGC5101
    micromamba run -n jwst_1.20.2 python run_ifu_science_ext.py --all

Output directories (under data/PID{pid}_{target}/):
    asn_ext/        — Level 2 and Level 3 ASN JSON
    stage2_ext/     — cal files (Parlanti overrides applied)
    stage3_ext/     — s3d cubes and x1d spectra

Reference:
    Parlanti et al. 2025, arXiv:2512.14844
"""
import os, sys, glob, json, time, argparse, warnings
import numpy as np
from astropy.io import fits

warnings.simplefilter("ignore", RuntimeWarning)

# ── Environment ────────────────────────────────────────────────────────────────
os.environ.setdefault('CRDS_PATH', os.path.expanduser('~/crds_cache'))
os.environ.setdefault('CRDS_SERVER_URL', 'https://jwst-crds.stsci.edu')
PYTHONPATH_EXT = '/Users/dcoe/NIRSpec/wavext/jwst_nirspec_wavext'
if PYTHONPATH_EXT not in sys.path:
    sys.path.insert(0, PYTHONPATH_EXT)

from jwst.pipeline import Spec2Pipeline, Spec3Pipeline
from jwst.associations import asn_from_list as afl
from jwst.associations.lib.rules_level2_base import DMSLevel2bBase
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base
import jwst
print(f"JWST pipeline: {jwst.__version__}")

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR      = '/Users/dcoe/NIRSpec/wavext/data'
PARLANTI_DIR  = f'{BASE_DIR}/parlanti_repo/CRDS_1364'
WAVELENGTHRANGE = f'{PARLANTI_DIR}/jwst_nirspec_wavelengthrange_0008.asdf'
CUBEPAR         = f'{PARLANTI_DIR}/jwst_nirspec_cubepar_0009.fits'

SFLAT = {
    ('G140M', 'F100LP', 'NRS1'): f'{PARLANTI_DIR}/jwst_nirspec_sflat_0208.fits',
    ('G140M', 'F100LP', 'NRS2'): f'{PARLANTI_DIR}/jwst_nirspec_sflat_0191.fits',
    ('G235M', 'F170LP', 'NRS1'): f'{PARLANTI_DIR}/jwst_nirspec_sflat_0211.fits',
    ('G235M', 'F170LP', 'NRS2'): f'{PARLANTI_DIR}/jwst_nirspec_sflat_0192.fits',
}
FFLAT = {
    'G140M': f'{PARLANTI_DIR}/jwst_nirspec_fflat_0105.fits',
    'G235M': f'{PARLANTI_DIR}/jwst_nirspec_fflat_0091.fits',
}
TARGET_CONFIGS = {('G140M', 'F100LP'), ('G235M', 'F170LP')}
EXP_TYPE = 'NRS_IFU'

# ── Science Target Registry ────────────────────────────────────────────────────
# Maps (pid, target_key) → list of rate file directories
# For each target, we list ALL directories that may contain its rate files.
SCIENCE_TARGETS = {
    ('2654', 'SDSSJ0749'): [f'{BASE_DIR}/PID2654_SDSSJ0749'],
    ('2654', 'SDSSJ0841'): [f'{BASE_DIR}/PID2654_SDSSJ0841'],
    # For 2186, rate files live in the rate/ subdirectory; base_dir is still
    # PID2186_UGC-5101 so that stage2_ext/stage3_ext land in the right place.
    ('2186', 'UGC5101'):   [f'{BASE_DIR}/PID2186_UGC-5101',
                             f'{BASE_DIR}/PID2186_UGC-5101/rate',
                             f'{BASE_DIR}/PID2186_UGC5101',
                             f'{BASE_DIR}/PID2186_UGC5101/rate'],
}

# ── Helpers ────────────────────────────────────────────────────────────────────

def get_header(f, key=None):
    h = fits.getheader(f)
    return h if key is None else h.get(key)


def check_reffiles():
    missing = []
    for p in [WAVELENGTHRANGE, CUBEPAR] + list(SFLAT.values()) + list(FFLAT.values()):
        if not os.path.exists(p) or os.path.getsize(p) == 0:
            missing.append(os.path.basename(p))
    if missing:
        raise FileNotFoundError(f"Missing Parlanti refs: {missing}\n  Expected in: {PARLANTI_DIR}")
    print(f"  All Parlanti reference files present.")


def collect_rate_files(dirs, grating=None, filt=None):
    """Collect IFU rate files from given directories, optionally filtering by config."""
    files = []
    for d in dirs:
        files += sorted(glob.glob(os.path.join(d, '*_rate.fits')))
    # Filter to IFU science exposures only
    sel = []
    for f in files:
        try:
            h = get_header(f)
            if h.get('EXP_TYPE') != EXP_TYPE:
                continue
            if grating and h.get('GRATING') != grating:
                continue
            if filt and h.get('FILTER') != filt:
                continue
            sel.append(f)
        except Exception:
            pass
    return sorted(sel)


def writel3asn(cal_files, asn_dir, product_name):
    """Create a Level-3 ASN combining all cal files for one grating/filter combo."""
    asnf = os.path.join(asn_dir, f'{product_name}_l3asn.json')
    asn = afl.asn_from_list(cal_files, rule=DMS_Level3_Base, product_name=product_name)
    _, s = asn.dump()
    with open(asnf, 'w') as fh:
        fh.write(s)
    return asnf

# ── Stage 2 ───────────────────────────────────────────────────────────────────

def run_stage2(pid, target, grating, filt, rate_files, out_dir, asn_dir):
    """Run Spec2 with Parlanti overrides on the given rate files."""
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(asn_dir, exist_ok=True)

    cal_files = []
    for rate_file in rate_files:
        base = os.path.basename(rate_file).replace('_rate.fits', '')
        det  = get_header(rate_file, 'DETECTOR')
        out_cal = os.path.join(out_dir, f'{base}_cal.fits')

        if os.path.exists(out_cal):
            print(f"    EXISTS: {base}_cal.fits")
            cal_files.append(out_cal)
            continue

        sflat_path = SFLAT.get((grating, filt, det))
        fflat_path = FFLAT.get(grating)
        if not sflat_path:
            print(f"    SKIP {base}: no sflat for {grating}/{filt}/{det}")
            continue

        # Build a minimal Level-2 ASN for this single exposure
        asnf = os.path.join(asn_dir, f'{base}_l2asn.json')
        asn  = afl.asn_from_list([rate_file], rule=DMSLevel2bBase, product_name=base)
        asn.data['program'] = str(pid)
        _, s = asn.dump()
        with open(asnf, 'w') as fh:
            fh.write(s)

        spec2_steps = {
            'assign_wcs': {'override_wavelengthrange': WAVELENGTHRANGE},
            'flat_field':  {
                'override_sflat': sflat_path,
                'override_fflat': fflat_path,
            },
            'bkg_subtract': {'skip': True},
            'cube_build':   {'skip': True},
            'extract_1d':   {'skip': True},
        }

        print(f"    Spec2: {base} [{det}]")
        try:
            Spec2Pipeline.call(asnf, save_results=True,
                               steps=spec2_steps, output_dir=out_dir)
            if os.path.exists(out_cal):
                cal_files.append(out_cal)
                print(f"      → {base}_cal.fits")
            else:
                # Try finding the output with different suffix
                produced = sorted(glob.glob(os.path.join(out_dir, f'*{base}*_cal.fits')))
                cal_files.extend(produced)
                if produced:
                    print(f"      → {[os.path.basename(p) for p in produced]}")
                else:
                    print(f"      WARNING: no cal file produced")
        except Exception as e:
            print(f"      ERROR: {e}")

    return cal_files

# ── Stage 3 ───────────────────────────────────────────────────────────────────

def run_stage3(pid, target, grating, filt, cal_files, out_dir, asn_dir):
    """Run Spec3 cube build + extraction on the cal files."""
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(asn_dir, exist_ok=True)

    grat_lower = grating.lower()
    filt_lower = filt.lower()
    product_name = f'{grat_lower}_{filt_lower}'
    out_s3d = os.path.join(out_dir, f'{product_name}_s3d.fits')
    out_x1d = os.path.join(out_dir, f'{product_name}_x1d.fits')

    if os.path.exists(out_x1d):
        print(f"    EXISTS: {product_name}_x1d.fits")
        return out_s3d, out_x1d

    # Filter cal files to this grating/filter
    sel_cal = []
    for f in cal_files:
        try:
            g = get_header(f, 'GRATING')
            fi = get_header(f, 'FILTER')
            if g == grating and fi == filt:
                sel_cal.append(f)
        except Exception:
            pass

    if not sel_cal:
        print(f"    No cal files for {grating}/{filt}")
        return None, None

    n_nrs1 = sum(1 for f in sel_cal if 'nrs1' in f.lower())
    n_nrs2 = sum(1 for f in sel_cal if 'nrs2' in f.lower())
    print(f"    Stage3 {grating}/{filt}: {len(sel_cal)} cal ({n_nrs1} NRS1, {n_nrs2} NRS2)")
    print(f"    cubepar: {os.path.basename(CUBEPAR)}")

    asnf = writel3asn(sel_cal, asn_dir, product_name)
    spec3_steps = {
        'master_background': {'skip': True},
        'outlier_detection': {'kernel_size': '3 3'},
        'cube_build': {'override_cubepar': CUBEPAR},
        'extract_1d': {},
    }
    try:
        Spec3Pipeline.call(asnf, save_results=True,
                           steps=spec3_steps, output_dir=out_dir)
    except Exception as e:
        print(f"    ERROR stage3: {e}")
        return None, None

    s3d = out_s3d if os.path.exists(out_s3d) else None
    x1d = out_x1d if os.path.exists(out_x1d) else None
    # Also check alternate naming
    if x1d is None:
        x1d_alts = glob.glob(os.path.join(out_dir, f'*{grat_lower}*x1d.fits'))
        if x1d_alts:
            x1d = sorted(x1d_alts)[0]
    print(f"    → s3d={os.path.basename(s3d) if s3d else 'None'}, "
          f"x1d={os.path.basename(x1d) if x1d else 'None'}")
    return s3d, x1d


# ── Main ──────────────────────────────────────────────────────────────────────

def process_target(pid, target, gratings=None):
    """Full Stage 2+3 extended pipeline for one science target."""
    key = (pid, target)
    if key not in SCIENCE_TARGETS:
        print(f"Unknown target: PID {pid}, {target}")
        print(f"Available: {sorted(SCIENCE_TARGETS.keys())}")
        return

    rate_dirs = SCIENCE_TARGETS[key]
    # Output base — use the first rate directory
    base_dir = rate_dirs[0]

    out_stage2 = os.path.join(base_dir, 'stage2_ext')
    out_stage3 = os.path.join(base_dir, 'stage3_ext')
    asn_dir    = os.path.join(base_dir, 'asn_ext')

    # Determine which configs to process
    if gratings is None:
        # Use all TARGET_CONFIGS present in the rate files
        available = set()
        for d in rate_dirs:
            for f in glob.glob(os.path.join(d, '*nrs1_rate.fits')):
                try:
                    h = get_header(f)
                    if h.get('EXP_TYPE') == EXP_TYPE:
                        available.add((h['GRATING'], h['FILTER']))
                except Exception:
                    pass
        configs = [(g, fi) for (g, fi) in available if (g, fi) in TARGET_CONFIGS]
    else:
        configs = [(g, fi) for (g, fi) in TARGET_CONFIGS
                   if g in gratings or fi in gratings]

    print(f"\n{'='*60}")
    print(f"PID {pid} — {target}")
    print(f"Configs: {configs}")
    print(f"Rate dirs: {[os.path.basename(d) for d in rate_dirs]}")
    print(f"{'='*60}")

    all_cal = []
    for grating, filt in sorted(configs):
        rate_files = collect_rate_files(rate_dirs, grating=grating, filt=filt)
        print(f"\n  [{grating}/{filt}] {len(rate_files)} rate files")
        if not rate_files:
            print(f"  No rate files for {grating}/{filt} in {rate_dirs}")
            continue

        cal = run_stage2(pid, target, grating, filt, rate_files, out_stage2, asn_dir)
        all_cal.extend(cal)
        print(f"  [{grating}/{filt}] Stage 2 done: {len(cal)} cal files")

    print(f"\n  [{pid}/{target}] Running Stage 3 ...")
    results = {}
    for grating, filt in sorted(configs):
        s3d, x1d = run_stage3(pid, target, grating, filt, all_cal, out_stage3, asn_dir)
        results[(grating, filt)] = (s3d, x1d)

    return results


def parse_args():
    p = argparse.ArgumentParser(description='Science IFU extended pipeline (PIDs 2654, 2186)')
    p.add_argument('--pid',    default=None, help='PID number (2654 or 2186)')
    p.add_argument('--target', default=None, help='Target name (SDSSJ0749, SDSSJ0841, UGC5101)')
    p.add_argument('--all',    action='store_true', help='Process all science targets')
    p.add_argument('--stage',  default='23', help='Stage(s): 2, 3, or 23')
    p.add_argument('--grating', default=None,
                   help='Grating filter: G140M or G235M (default: all)')
    return p.parse_args()


if __name__ == '__main__':
    args = parse_args()

    print("Checking Parlanti reference files ...")
    check_reffiles()

    t0 = time.perf_counter()

    if args.all:
        targets_to_run = list(SCIENCE_TARGETS.keys())
    elif args.pid and args.target:
        targets_to_run = [(args.pid, args.target)]
    elif args.pid == '2654':
        targets_to_run = [('2654', 'SDSSJ0749'), ('2654', 'SDSSJ0841')]
    elif args.pid == '2186':
        targets_to_run = [('2186', 'UGC5101')]
    else:
        print("Specify --pid 2654 --target SDSSJ0841, or --all")
        sys.exit(1)

    gratings = None
    if args.grating:
        gratings = [args.grating]

    for pid, target in targets_to_run:
        process_target(pid, target, gratings=gratings)

    elapsed = time.perf_counter() - t0
    print(f"\nTotal time: {elapsed/60:.1f} min")
