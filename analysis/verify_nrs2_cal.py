"""Verify new NRS2 x1d files have valid Jy fluxes."""
import numpy as np
from astropy.io import fits
import glob

calspec = {
    'g191b2b': '/Users/dcoe/NIRSpec/wavext/data/CALSPEC/g191b2b_mod_012.fits',
}

def load_calspec_jy(fname):
    C = 2.99792458e18
    with fits.open(fname) as h:
        d = h[1].data
        wl = d['WAVELENGTH'].astype(float)
        fl = d['FLUX'].astype(float)
    fnu = fl * wl**2 / C * 1e23
    return wl / 1e4, fnu  # um, Jy

wl_cs, fj_cs = load_calspec_jy(calspec['g191b2b'])

base = '/Users/dcoe/NIRSpec/wavext/data/PID1537_G191-B2B/nrs2_spec2_cal/'
for f in sorted(glob.glob(base + '*_x1d.fits')):
    with fits.open(f) as h:
        wl = h[1].data['WAVELENGTH'].astype(float)
        fl = h[1].data['FLUX'].astype(float)
        unit = str(h[1].columns['FLUX'].unit)
        ok = np.isfinite(fl) & (fl != 0)
        n_ok = ok.sum()
        if n_ok > 0:
            med = float(np.nanmedian(fl[ok]))
            # Compare with CALSPEC at overlap wavelengths
            w_lo, w_hi = wl[ok].min(), wl[ok].max()
            cs_mask = (wl_cs > w_lo) & (wl_cs < w_hi)
            med_cs = float(np.nanmedian(fj_cs[cs_mask])) if cs_mask.sum() > 0 else float('nan')
            ratio = med / med_cs if med_cs > 0 else float('nan')
            print(f'{f.split("/")[-1]}: unit={unit}, wl=[{w_lo:.3f},{w_hi:.3f}], '
                  f'n_ok={n_ok}/2048, med={med:.4g} Jy, CALSPEC_med={med_cs:.4g} Jy, ratio={ratio:.3f}')
        else:
            print(f'{f.split("/")[-1]}: unit={unit}, NO valid flux')
