"""Quick check of v6 data files: wavelengths, flux units."""
import os
from astropy.io import fits
import numpy as np

BASE = '/Users/dcoe/NIRSpec/wavext'

files_fs_standards = [
    ('FS J1743 G140M NRS2 ext', f'{BASE}/data/PID1536_J1743045/nrs2_spec3_ext/nrs2_l3_1536_g140m_s1600a1_s000000001_x1d.fits'),
    ('FS J1743 G235M NRS2 ext', f'{BASE}/data/PID1536_J1743045/nrs2_spec3_ext/nrs2_l3_1536_g235m_s1600a1_s000000001_x1d.fits'),
    ('FS P330E G140M NRS2 ext', f'{BASE}/data/PID1538_P330E/nrs2_spec3_ext/nrs2_l3_1538_g140m_s1600a1_s000000001_x1d.fits'),
    ('FS G191 G140M NRS2 ext',  f'{BASE}/data/PID1537_G191-B2B/nrs2_spec3_ext/nrs2_l3_1537_g140m_s1600a1_s000000001_x1d.fits'),
    ('FS 6644 G140M NRS2 ext',  f'{BASE}/data/PID6644_NGC2506G31/nrs2_spec3_ext/nrs2_l3_6644_g140m_s1600a1_s000000001_x1d.fits'),
    ('FS 6644 G235M NRS2 ext',  f'{BASE}/data/PID6644_NGC2506G31/nrs2_spec3_ext/nrs2_l3_6644_g235m_s1600a1_s000000001_x1d.fits'),
]

files = [
    ('IFU J1743 G140M ext', f'{BASE}/data/IFU/PID1536_J1743045/stage3_ext/f100lp_g140m-f100lp_x1d.fits'),
    ('IFU J1743 G235M ext', f'{BASE}/data/IFU/PID1536_J1743045/stage3_ext/f170lp_g235m-f170lp_x1d.fits'),
    ('IFU P330E G140M ext', f'{BASE}/data/IFU/PID1538_P330E/stage3_ext/f100lp_g140m-f100lp_x1d.fits'),
    ('IFU P330E G235M ext', f'{BASE}/data/IFU/PID1538_P330E/stage3_ext/f170lp_g235m-f170lp_x1d.fits'),
    ('IFU G191 G140M ext',  f'{BASE}/data/IFU/PID1537_G191-B2B/stage3_ext/f100lp_g140m-f100lp_x1d.fits'),
    ('IFU G191 G235M ext',  f'{BASE}/data/IFU/PID1537_G191-B2B/stage3_ext/f170lp_g235m-f170lp_x1d.fits'),
    ('IFU UGC5101 G235M ext', f'{BASE}/data/PID2186_UGC-5101/stage3_ext/g235m_f170lp_g235m-f170lp_x1d.fits'),
    ('IFU UGC5101 G395M nom', f'{BASE}/data/PID2186_UGC5101/stage3_nom/g395m_f290lp_g395m-f290lp_x1d.fits'),
    ('IFU SDSSJ0841 G140M ext', f'{BASE}/data/PID2654_SDSSJ0841/stage3_ext/g140m_f100lp_g140m-f100lp_x1d.fits'),
    ('FS PID1492 G140M MAST',  f'{BASE}/data/PID1492/jw01492-o001_t001-s000000001_nirspec_f100lp-g140m-s200a1-subs200a1_x1d.fits'),
    ('FS PID1492 G235M MAST',  f'{BASE}/data/PID1492/jw01492-o001_t001-s000000001_nirspec_f170lp-g235m-s200a1-subs200a1_x1d.fits'),
    ('FS PID1492 G395M MAST',  f'{BASE}/data/PID1492/jw01492-o001_t001-s000000001_nirspec_f290lp-g395m-s200a1-subs200a1_x1d.fits'),
    ('FS PID1492 G140M NRS2 ext', f'{BASE}/data/PID1492/jw01492003001_03102_00005_nrs2_g140m_extract_1d.fits'),
    ('FS PID1492 G235M NRS2 ext', f'{BASE}/data/PID1492/jw01492003001_03104_00004_nrs2_g235m_extract_1d.fits'),
]

all_files = files_fs_standards + files
for label, path in all_files:
    if not os.path.exists(path):
        print(f'MISSING: {label}')
        continue
    try:
        with fits.open(path) as h:
            d = h[1].data
            wl = d['WAVELENGTH'].astype(float)
            fl = d['FLUX'].astype(float)
            ok = (wl > 0) & np.isfinite(fl) & (fl > 0)
            u_wl = h[1].header.get('TUNIT1', '?')
            u_fl = h[1].header.get('TUNIT2', '?')
            med = np.nanmedian(fl[ok]) if ok.any() else float('nan')
            print(f'{label}: {wl[ok].min():.3f}-{wl[ok].max():.3f} {u_wl}, '
                  f'unit={u_fl}, median={med:.3e}, n={ok.sum()}')
    except Exception as e:
        print(f'ERROR {label}: {e}')
