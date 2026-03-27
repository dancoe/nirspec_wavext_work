from util import *

import numpy as np
import os

#import jwst
from jwst import datamodels
from jwst.associations import asn_from_list   # create association file

import matplotlib  # as mpl
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import matplotlib.patheffects as pe

from astropy.io import fits
import astropy.wcs as wcs
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import Table, vstack, unique
from astropy.visualization import simple_norm, ImageNormalize, AsinhStretch
from astropy.visualization import LogStretch, LinearStretch, ManualInterval
from astropy.stats import sigma_clip

import inspect
from copy import deepcopy

# Color image
import PIL  # Python Image Library
from PIL import Image, ImageEnhance
PIL.Image.MAX_IMAGE_PIXELS = 933120000  # allow it to load large image

# https://docs.google.com/spreadsheets/d/1ZItzi28T-njSIn4fpmXi1g9BKqUCSGtDwc4glgC-25g
line_list = '''
Ly$\\alpha$    1215.67
[OII]  3727.0635
[OII]  3729.8472
H$\\beta$  4862.6739
[OIII] 4960.2785
[OIII] 5008.2236
[NeIII] 3868.76
[NeIII] 3967.47
#[NII] 5754.59
#[NII] 6548.05
#[NII] 6583.46
H$\\alpha$ 6564.6237
'''.split('\n')[1:-1]

def show_MOS_rate(rate_file, slits_model=None,
                  save_plot=False, close_plot=False, integration=None,
                  cmap='viridis', bad_color=(1, 0.7, 0.7),
                  vmin=-0.003, vmax=0.022, show_colorbar=True, units='DN/s',
                  title_prefix=None, title_path=False):
    """
    Parameters
    ----------
    rate_file: str
        Path to the rate/rateints FITS file (count rate data).
    slit_model: MultiSlitModel or None
        Slit models for all the sources.
        Can come from the CAL or S2D product datamodel (from Spec2Pipeline).
    save_plot: bool or str
        If True, save the plot as a PNG file.
        If a string is provided, it is used as the filename. Default is False.
    close_plot: bool
        If True, close the plot after displaying.
    integration: str or int or None
        If the data is 3D (rateints), define an integration to plot.
        If 'min', extracts the minimum slice of the cube for plotting.
        Default is None (assumes 2D data).
    cmap: str
        Colormap for the plot.
    bad_color: str or tuple
        Color of "bad" pixels (e.g., NaNs) in the colormap.
        Can be specified as a color name or RGB tuple (values 0 to 1).
    vmin/vmax: float
        Minimum/Maximum value for colormap range.
    show_colorbar: bool
        Show the figure colorbar? Default is True.
    title_prefix: str or None
        Prefix for the plot title.
    title_path: bool
        If True, include the full path of the FITS file in the plot title.
    """

    # Open the rate FITS file
    with fits.open(rate_file) as hdu_list:
        data = hdu_list['SCI'].data  # count rate data

        # Is the data 3D (_rateints.fits)?
        if integration == 'min':
            data = data.min(axis=0)  # min slice of cube
        elif integration is not None:
            data = data[integration]

    # Setup the figure and colorbar
    cmap = matplotlib.colormaps[cmap]
    cmap.set_bad(bad_color, 1.)  # color/opacity of bad pixels

    fig, ax = plt.subplots(1, 1, figsize=(12, 12))
    norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=AsinhStretch())
    im = ax.imshow(data, origin='lower', cmap=cmap,
                   norm=norm, interpolation='nearest')
    ax.set_ylabel('Pixel Row')
    ax.set_xlabel('Pixel Column')

    if slits_model is not None:  # Draw slits and label source ids

        # add white outline to text below
        # path_effects=[pe.withStroke(linewidth=3, foreground="w", alpha=0.9)]
        path_effects = []  # no outline
        fontsize = 7
        color = 'w'

        # For each slitlet, draw the slit patch
        slit_patches = []
        for slit in slits_model.slits:
            slit_patch = Rectangle((slit.xstart, slit.ystart),
                                   slit.xsize, slit.ysize)
            slit_patches.append(slit_patch)

            y = slit.ystart + slit.ysize/2  # slit center location
            va = 'center'  # vertically centered text

            # Label the spectra on the left hand side for NRS1:
            if 'nrs1' in rate_file:
                x = slit.xstart
                ha = 'right'  # align test to the right
            else:  # Label the spectra on the right hand side for NRS2:
                x = slit.xstart + slit.xsize
                ha = 'left'  # align test to the left

            plt.text(x, y, slit.source_id, color=color, ha=ha, va=va,
                     fontsize=fontsize, path_effects=path_effects,
                     weight='bold')

        # plot the slit patches
        ax.add_collection(PatchCollection(slit_patches, ec='r', fc='None'))

    # plot title
    title = title_prefix + '  ' if title_prefix else ' '
    title += rate_file if title_path else os.path.basename(rate_file)

    if integration is not None:
        title = title.replace('rateints', 'rateints[%s]' % integration)

    plt.title(title)
    print(title)

    # https://stackoverflow.com/questions/18195758/set-matplotlib-colorbar-size-to-match-graph
    if show_colorbar:
        plt.subplots_adjust(left=0.05, right=0.85)
        #units = 'DN/s'
        cbar_dx = 0.02
        cbar_ax = fig.add_axes([ax.get_position().x1+cbar_dx,
                                ax.get_position().y0, cbar_dx,
                                ax.get_position().height])
        cbar = fig.colorbar(im, label=units, cax=cbar_ax)
        cbar_ticks = cbar_ax.get_yticks()
        cbar_ticks = np.concatenate([cbar_ticks, [vmin, vmax]])
        cbar_ticks = np.compress(between(vmin, cbar_ticks, vmax),
                                 cbar_ticks)
        cbar_ticks = np.sort(cbar_ticks)
        cbar.set_ticks(cbar_ticks)
        # print(cbar_ticks)

    # Save the plot?
    if save_plot:
        if type(save_plot) != type('a.fits'):
            save_plot = rate_file.replace('fits', 'png')
            if integration != None:
                save_plot = save_plot.replace('.png', '%s.png' % integration)
        plt.savefig(save_plot, dpi=200)
    if close_plot:
        plt.close()
        
    
def show_MOS_rate_files(rate_files, slit_models=[], save_plot=False, close_plot=False, integration=None,
                        vmin=None, vmax=None, show_colorbar=True):
                  #vmin=-0.003, vmax=0.022, 

    fig, axs = plt.subplots(1, len(rate_files), figsize=(15,8), sharey=True)
    plt.subplots_adjust(wspace=0.02)

    cmap = 'viridis'
    bad_color = 1, 0.7, 0.7
    cmap = matplotlib.colormaps[cmap]
    cmap.set_bad(bad_color, 1.)
    
    for ifile, rate_file in enumerate(rate_files):        
        if not os.path.exists(rate_file):
            continue
            
        with fits.open(rate_file) as hdu_list:
            data = hdu_list['SCI'].data
            if integration == 'min':
                data = data.min(axis=0)
            elif integration != None:
                data = data[integration]

        ax = axs[ifile]
        if vmax:
            norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=AsinhStretch())
        else:
            norm = simple_norm(data, 'asinh', min_percent=20, max_percent=98)
        #norm = ImageNormalize(vmin=-0.003, vmax=0.022, stretch=AsinhStretch())
        #print(norm.vmin, norm.vmax)
        # Turn off interpolation! Or else bad pixels will appear to grow in the plot
        im  = ax.imshow(data, origin='lower', cmap=cmap, norm=norm, interpolation='nearest')
        
        if len(slit_models):
            # Labels 2D extraction regions on MSA detector rate image
            slit_model = slit_models[ifile]
            slit_patches = []
            path_effects=[pe.withStroke(linewidth=3, foreground="w", alpha=0.9)]
            fontsize = 8
            for slit in slit_model.slits:
                slit_patch = Rectangle((slit.xstart, slit.ystart), slit.xsize, slit.ysize)
                slit_patches.append(slit_patch)
                #plt.text(slit.xstart + slit.xsize, slit.ystart + slit.ysize/2, slit.source_id, color='w', va='center', fontsize=12, weight='bold')
                #print('nrs1' in rate_file)
                if 'nrs1' in rate_file: # Label the spectra on the left hand side for NRS1:
                    ax.text(slit.xstart, slit.ystart + slit.ysize/2, slit.source_id, color='r', ha='right', va='center', fontsize=fontsize, path_effects=path_effects)
                else:  # 'nrs2' # Label the spectra on the right hand side for NRS2:
                    ax.text(slit.xstart + slit.xsize, slit.ystart + slit.ysize/2, slit.source_id, color='r', va='center', fontsize=fontsize, path_effects=path_effects)

            ax.add_collection(PatchCollection(slit_patches, ec='r', fc='None'))
        
        title = os.path.basename(rate_file)
        if integration != None:
            title = title.replace('rateints', 'rateints[%s]' % integration)
    
        ax.set_title(title)
        print(title)

    # https://stackoverflow.com/questions/18195758/set-matplotlib-colorbar-size-to-match-graph
    if show_colorbar:
        plt.subplots_adjust(left=0.05, right=0.85)
        units = 'DN/s'
        cbar_dx = 0.02
        cbar_ax = fig.add_axes([ax.get_position().x1+cbar_dx,
                                ax.get_position().y0, cbar_dx,
                                ax.get_position().height])
        cbar = fig.colorbar(im, label=units, cax=cbar_ax)
        
        # Add min and max to ticks on colorbar
        cbar_ticks = cbar_ax.get_yticks()
        cbar_ticks = np.concatenate([cbar_ticks, [norm.vmin, norm.vmax]])
        cbar_ticks = np.compress(between(norm.vmin, cbar_ticks, norm.vmax), cbar_ticks)
        cbar_ticks = np.sort(cbar_ticks)
        cbar.set_ticks(cbar_ticks)
        # print(cbar_ticks)

    if save_plot:
        if type(save_plot) != type('a.fits'):
            save_plot = rate_file.replace('fits', 'png')
            save_plot = save_plot.replace('_nrs1', '')
            save_plot = save_plot.replace('_nrs2', '')
            if integration != None:
                save_plot = save_plot.replace('.png', '%s.png' % integration)
        print('SAVING', save_plot)
        plt.savefig(save_plot, dpi=200)

    if close_plot:
        plt.close()
        

def create_spec2_asn_files(rate_files, output_dir, nods=3):
    """
    Generate association (ASN) files for sets of 3 nodded rate files.
    Each rate file will have an associated ASN file where it serves
    as the science exposure not background.

    Parameters
    ----------
    rate_files: list of str
        A list of strings representing the paths to the rate files.
    output_dir: str
        Path to the directory where the generated ASN files will be saved.

    Returns
    -------
    asn_files: list str
        A list of paths to the newly generated ASN files.
    """
    asn_files = []
    rate_files = [os.path.basename(rate_file) for rate_file in rate_files]
    
    # Blocks of 3 nods
    rate_file_groups = [rate_files[i:i+nods] for i in range(0, len(rate_files), nods)]
    for rate_file_group in rate_file_groups:
        # Loop through the rate file list
        for science_index, science_rate_file in enumerate(rate_file_group):
            # Create empty ASN file
            product_name = os.path.basename(science_rate_file)
            asn_data = asn_from_list.asn_from_list(rate_file_group,
                                                   product_name=product_name)
    
            # Add science and background exposures to the ASN file
            for i, member in enumerate(asn_data['products'][0]['members']):
                if i != science_index:
                    member['exptype'] = 'background'
    
            # Name and save the ASN files
            asn_file = science_rate_file.replace('_rate.fits', '_asn.json')
            asn_file = os.path.join(output_dir, asn_file)
            asn_files.append(asn_file)
            print(asn_file)
            with open(asn_file, 'w') as outfile:
                name, serialized = asn_data.dump(format='json')
                outfile.write(serialized)

    return asn_files


def load_NIRCam_image(color_image_file, image_file=None):
    # color_rgb.fits
    # color.png, image.fits
    # image.fits
    if '.fits' in color_image_file:  # or could be black & white, still works
        print('Loading', color_image_file)
        color_image_hdulist = fits.open(color_image_file)
        image_wcs = wcs.WCS(color_image_hdulist[0].header, color_image_hdulist)
        NIRCam_image = np.stack([hdu.data for hdu in color_image_hdulist[1:]], axis=-1)
    else:
        print('Loading', color_image_file)
        im = Image.open(color_image_file)
        im = im.transpose(method=Image.FLIP_TOP_BOTTOM)
        NIRCam_image = np.asarray(im)  # (24576, 16384, 3)
        
        print('Loading', image_file)
        image_hdulist = fits.open(image_file)
        image_wcs = wcs.WCS(image_hdulist[0].header, image_hdulist)

    # Make these variables available globally within the module
    # (or from outside via mos.NIRCam_image, mos.image_wcs)
    for variable in 'NIRCam_image image_wcs'.split():
        globals()[variable] = locals()[variable]
        

def spectra_similarity(spectrum1, spectrum2):
    # Remove NaN values
    mask = ~np.isnan(spectrum1) & ~np.isnan(spectrum2)
    spectrum1 = spectrum1[mask]
    spectrum2 = spectrum2[mask]
    
    if len(spectrum1) == 0:  # no overlap
        return 0, 0

    # Correlation coefficient
    correlation = np.corrcoef(spectrum1, spectrum2)[0, 1]

    # Multiplicative factor from areas under the curve
    factor_from_integrals = np.trapz(spectrum1) / np.trapz(spectrum2)

    return correlation, factor_from_integrals


open_slit_x_size = 0.20  # open slit width  in arcseconds (dispersion direction)
open_slit_y_size = 0.46  # open slit height in arcseconds (cross-dispersion direction)
slit_bar_width   = 0.07  # bar width=height in arcseconds (between the slits)

open_slit_aspect = open_slit_y_size / open_slit_x_size

full_slit_x_size = open_slit_x_size + slit_bar_width
full_slit_y_size = open_slit_y_size + slit_bar_width
full_slit_aspect = full_slit_y_size / full_slit_x_size
#full_slit_x_size, full_slit_y_size

x_scale_open_to_full = full_slit_x_size / open_slit_x_size
y_scale_open_to_full = full_slit_y_size / open_slit_y_size
#x_scale_open_to_full, y_scale_open_to_full

# Slit corners in coordinates that range from (0,0) to (1,1)
open_slit_x_corners = 0, 0, 1, 1 
open_slit_y_corners = 0, 1, 1, 0

# Convert coordinates to slit centroid at (0,0)
open_slit_x_corners = np.array(open_slit_x_corners) - 0.5
open_slit_y_corners = np.array(open_slit_y_corners) - 0.5

# Convert to full slit (not just open area)
full_slit_x_corners = x_scale_open_to_full * open_slit_x_corners
full_slit_y_corners = y_scale_open_to_full * open_slit_y_corners


class SlitSpectrum:
    
    """
    Plot the 2D (S2D) and corresponding 1D (X1D) MOS spectra.

    Parameters:
    ----------
    s2d_model: MultiSlitModel
        Slit model(s) from an S2D product datamodel
        either for all the sources (from Spec2Pipeline)
        or from a single source (from Spec3Pipeline).
    x1d_model: MultiSlitModel
        Slit model(s) from an X1D product datamodel
        either for all the sources (from Spec2Pipeline)
        or from a single source (from Spec3Pipeline).
    source_id: int or None
        Source to inspect (may be left as None for 
        a single source from Spec3Pipeline).
    cmap: str
        Colormap for the plot.
    bad_color: str or tuple
        Color of "bad" pixels (e.g., NaNs) in the colormap.
        Can be specified as a color name or RGB tuple (values 0 to 1).
    expand_wavelength_gap: bool
        Expand the wavelength gap to fill missing wavelength values?
        Default is True.
    """
    
    def __init__(self, s2d_files):
        self.z = None
        self.msa_metadata_id = 1  # default; can change
        self.s2d_files = s2d_files
        #self.load_data()
        
    def assign(self, *args, set_value=False):
        frame = inspect.currentframe().f_back
        local_vars = frame.f_locals

        for value in args:
            for var_name, var_value in local_vars.items():
                if var_value is value and not var_name.startswith('_'):
                    if set_value is not False:
                        value = set_value
                    setattr(self, var_name, value)
                    break

    def assign_none(self, variables):
        frame = inspect.currentframe().f_back
        local_vars = frame.f_locals

        # Set each variable to None
        for var_name in variables:
            setattr(self, var_name, None)

    #def extract(self, *vars):
    #    frame = inspect.currentframe().f_back
    #    for var in vars:
    #        frame.f_locals[var] = getattr(self, var.__name__)

    def extract(self, var_names):
        frame = inspect.currentframe().f_back
        for var_name in var_names:
            frame.f_locals[var_name] = getattr(self, var_name)

    def load_data(self, source_id=None):
        s2d_file = self.s2d_files
        if source_id:
            try:
                #s2d_file = select_files(self.s2d_files, ['_s%09d' % source_id])[0]
                s2d_file = select_files(self.s2d_files, ['%09d' % source_id])[0]
            except:
                print('Source', source_id, 'not found in S2D files')
                pass
            
        s2d_filename = os.path.basename(s2d_file)
        s2d_dirname  = os.path.dirname(s2d_file)
        
        print('Loading %s...' % s2d_file)
        s2d_model = datamodels.open(s2d_file)
        x1d_file  = s2d_file.replace('s2d', 'x1d')
        print('Loading %s...' % x1d_file)
        x1d_model = datamodels.open(x1d_file)

        self.assign(s2d_file, s2d_model, s2d_dirname, s2d_filename,
                    x1d_file, x1d_model,
                    source_id)

    def load_source(self, source_id=None):
        #variables = 's2d_model x1d_model source_id'.split()
        #for variable in variables: 
        #    locals()[variable] = None
        #self.extract(variables)
    
        # 2D spectrum
        # s2d has all the objects; extract the one with source_id
        if 'slits' in list(self.s2d_model):
            source_ids = [slit.source_id for slit in self.s2d_model.slits]
            i_slit = source_ids.index(source_id)
            slit_model = self.s2d_model.slits[i_slit]
        else:  # s2d only has one object
            slit_model = deepcopy(self.s2d_model)
            # if background, then source_id will be assigned slitlet_id
            # so may need to override that with input source_id
            source_id = source_id or slit_model.source_id
            i_slit = 0

        s2d_data = slit_model.data + 0  # load and make copy
        # Replace zeros with nan where there is no data
        s2d_data = np.where(slit_model.err, s2d_data, np.nan)

        # 1D spectrum
        x1d_wave    = self.x1d_model.spec[i_slit].spec_table.WAVELENGTH
        x1d_flux    = self.x1d_model.spec[i_slit].spec_table.FLUX
        x1d_fluxerr = self.x1d_model.spec[i_slit].spec_table.FLUX_ERROR
        # fluxerr all nan?? pipeline bug
        if np.sum(np.isnan(x1d_fluxerr)) == len(x1d_fluxerr):
            # Replace zeros with nan where there is no data
            x1d_flux = np.where(x1d_flux, x1d_flux, np.nan)
        else:
            # Replace zeros with nan where there is no data
            x1d_flux = np.where(np.isnan(x1d_fluxerr), np.nan, x1d_flux)
    
        num_waves = len(x1d_wave)
        
        self.assign(source_id, slit_model, s2d_data, x1d_wave, x1d_flux, x1d_fluxerr, num_waves)
        
    # Expand the wavelength array for the gap?
    def expand_wavelength_gap(self):
        # calculate differences between consecutive wavelengths and
        # find if there is a gap to fill
        dx1d_wave = self.x1d_wave[1:] - self.x1d_wave[:-1]
        igap = np.argmax(dx1d_wave)
        dx1d_max = np.max(dx1d_wave)
        dx_replace = (dx1d_wave[igap-1] + dx1d_wave[igap+1]) / 2.
        num_fill = int(np.round(dx1d_max / dx_replace))
        print("Expanding wavelength gap %.2f -- %.2f microns"
              % (self.x1d_wave[igap], self.x1d_wave[igap+1]))

        #print(self.s2d_data.shape, self.x1d_flux.shape, num_fill)
        
        if num_fill > 1:  # There is a gap to fill
            wave_fill = np.mgrid[self.x1d_wave[igap]:self.x1d_wave[igap+1]:(num_fill+1)*1j]
            x1d_wave = np.concatenate([self.x1d_wave[:igap+1],wave_fill[1:-1],self.x1d_wave[igap+1:]])

            num_rows, num_waves = self.s2d_data.shape
            s2d_fill = np.zeros(shape=(num_rows, num_fill-1)) * np.nan
            self.s2d_data = np.concatenate([self.s2d_data[:, :igap+1], s2d_fill, 
                                            self.s2d_data[:, igap+1:]], axis=1)

            x1d_fill = np.zeros(shape=(num_fill-1)) * np.nan
            self.x1d_flux = np.concatenate([self.x1d_flux[:igap+1], x1d_fill, 
                                            self.x1d_flux[igap+1:]])
            self.x1d_fluxerr = np.concatenate([self.x1d_fluxerr[:igap+1], x1d_fill, 
                                                self.x1d_fluxerr[igap+1:]])
   
            # Fill the wavelength gap linearly (good enough since there's no data)
            wave_fill = np.linspace(self.x1d_wave[igap], self.x1d_wave[igap+1], num_fill+2)[1:-1]
            #print(self.x1d_wave[igap], self.x1d_wave[igap+1])
            self.x1d_wave = np.concatenate([self.x1d_wave[:igap+1], wave_fill,
                                            self.x1d_wave[igap+1:]])
        
        
        self.num_waves = len(self.x1d_flux)  # update to new length
        #print(self.s2d_data.shape, self.x1d_flux.shape, num_fill)

    def load_msa_metafile(self, output_dir=None):
        msa_metafile = fits.getval(self.s2d_file, 'MSAMETFL')
        output_dir = output_dir or self.s2d_dirname
        msa_metafile = os.path.join(output_dir, msa_metafile)
        msa_hdu_list = fits.open(msa_metafile)
        #msa_hdu_list.info()
        source_table  = Table(msa_hdu_list['SOURCE_INFO'].data)
        shutter_table = Table(msa_hdu_list['SHUTTER_INFO'].data)
        msa_metadata_ids = sorted(list(set(shutter_table['msa_metadata_id'])))
        slitlet_ids = np.sort(list(set(shutter_table['slitlet_id'])))        

        self.assign(msa_metafile, source_table, shutter_table, msa_metadata_ids, slitlet_ids)
        
    def make_figure(self, figsize=None, height_ratios=None):
        print('figsize', figsize)
        print('height_ratios', height_ratios)
        fig, (ax_2d, ax_1d) = plt.subplots(2, 1, figsize=figsize, height_ratios=height_ratios,
                                         squeeze=True, sharex=True)
        self.assign(fig, ax_1d, ax_2d)

    def set_wavelength_range(self, wave_min=None, wave_max=None):
        self.wave_min = wave_min or self.x1d_wave[0]
        self.wave_max = wave_max or self.x1d_wave[-1]
        #xtick_pos = np.interp(xticks, self.x1d_wave, np.arange(len(self.x1d_wave)))
        self.xmin = np.interp(self.wave_min, self.x1d_wave, np.arange(len(self.x1d_wave)))
        self.xmax = np.interp(self.wave_max, self.x1d_wave, np.arange(len(self.x1d_wave)))

    def plot_spec2d(self, cmap='RdBu', bad_color='w', sigma=5, maxiters=3, symmetric_range=False, title=None):
        cmap = matplotlib.colormaps[cmap]
        cmap.set_bad(bad_color, 1.)
        
        # ymax = 2 * np.nanpercentile(s2d_data, 90)
        # ymin = -0.2 * ymax

        sigma_clipped_data = sigma_clip(self.s2d_data, sigma=sigma, maxiters=maxiters)
        ymin = np.min(sigma_clipped_data)
        # ymin = 0.2 * ymin
        ymax = np.max(sigma_clipped_data)
        if symmetric_range:
            #center range about zero; but there may be some bias positive or negative
            ymin = -ymax

        # Plot the rectified 2D spectrum
        # norm = simple_norm(s2dsci, 'linear', min_percent=2, max_percent=99.9)
        # norm = simple_norm(s2d_data, 'asinh',
        # min_percent=min_percent, max_percent=max_percent)
        #norm = simple_norm(s2d_data, 'linear', min_cut=ymin, max_cut=ymax)
        norm = ImageNormalize(vmin=ymin, vmax=ymax, stretch=LinearStretch())
        print('s2d range:', norm.vmin, norm.vmax)

        self.ax_2d.imshow(self.s2d_data, origin='lower', cmap=cmap,
                    aspect='auto', norm=norm, interpolation='nearest')
        num_cross, num_dispersion = self.s2d_data.shape
        ny, nx = self.s2d_data.shape
        self.ax_2d.yaxis.set_ticks_position('both')
        # ax_2d.axhline((ny-1)/2., c='w', lw=0.5, alpha=0.5, ls='--')
        #self.ax_2d.set_ylabel("Pixel Row")

        #if 'extraction_ystart' in slit_spec.x1d_model.spec[0].instance.keys():
        if hasattr(self.x1d_model.spec[0], 'extraction_ystart'):
            # Extraction region: 1 pixel higher than input
            ystart = self.x1d_model.spec[0].extraction_ystart - 1 - 0.5
            ystop  = self.x1d_model.spec[0].extraction_ystop  - 1
            self.ax_2d.axhline(ystart-0.5, c='m', lw=0.5)
            self.ax_2d.axhline(ystop+0.5,  c='m', lw=0.5)
            self.ax_2d.set_yticks([0, ystart, ystop, ny-1])
            self.ax_2d.set_yticklabels([f'{y:.1f}' for y in self.ax_2d.get_yticks()])
        
        self.ax_2d.tick_params(labelbottom=False)  # Hide x-axis labels on the 2D plot
        self.ax_2d.yaxis.set_ticks_position('right')
        
        self.ax_2d.set_xlim(self.xmin, self.xmax)
        
        if title is None:
            title = self.s2d_filename.replace('_s2d.fits', '')
            title = title.replace('_', '  ')
            title += '  slit %d'   % self.slit_model.slitlet_id
            title += '  source %d' % self.slit_model.source_id

        self.ax_2d.set_title(title)


    def label_spec_lines(self, line_color=(1, 0.85, 0.7), missing_line_color=(0.9, 0.9, 0.9), line_label_color='k', detection_snr=3):
        label_boxes = []  # List to store label bounding boxes
        for line in line_list[::-1]:
            line_name, line_lam = line.split()
            if line_name[0] == '#':
                continue
            line_lam_rest = float(line_lam) * u.AA
            line_lam_obs = (line_lam_rest * (1+self.z)).to(u.um).value
            line_x = np.interp(line_lam_obs, self.x1d_wave, np.arange(len(self.x1d_wave)))
            
            if between(self.xmin+1, line_x, self.xmax-1):
                # Check for detection
                line_index = int(line_x)
                line_flux = self.x1d_flux[line_index]
                line_error = self.x1d_fluxerr[line_index]
                line_snr = line_flux / line_error
                #print(line_name, line_lam_obs, line_flux, line_error, line_snr)
                is_detected = line_snr > detection_snr

                color = line_color if is_detected else missing_line_color
                print(line_name, '%5.1f' % line_snr, color)
                self.ax_1d.axvline(line_x, c=color, zorder=-10)

                line_label = ' ' + line_name
                lam_str = '%d' % roundint(float(line_lam)) + '${\\rm \AA}$\n'
                line_label += '\n ' + lam_str
                
                # Start with the default y position
                y_pos = 0.95
                
                # Create a text object
                text = self.ax_1d.text(line_x, y_pos, line_label, color=line_label_color,
                                    va='top', transform=self.ax_1d.get_xaxis_transform(),
                                    fontsize=10)
                
                # Get the bounding box of the text in display coordinates
                bbox = text.get_window_extent(renderer=self.fig.canvas.get_renderer())
                
                # Convert the bbox to data coordinates
                bbox_data = bbox.transformed(self.ax_1d.transData.inverted())
                
                # Check for overlaps and adjust position if necessary
                while any(bbox_data.overlaps(box) for box in label_boxes):
                    y_pos -= 0.01
                    text.set_position((line_x, y_pos))
                    bbox = text.get_window_extent(renderer=self.fig.canvas.get_renderer())
                    bbox_data = bbox.transformed(self.ax_1d.transData.inverted())
                
                # Store the final bounding box
                label_boxes.append(bbox_data)

    def plot_spec1d(self, wave_tick_interval=0.2, wave_tick_fmt='%.1f',
                    ymargin_1d_pos = 0.1,  # margin above max
                    ymargin_1d_neg = 0.07, # margin below zero
                    sigma_1d_pos=1000, maxiters_1d_pos=1,  # 1D spectrum clipping
                    #sigma_1d_neg=10,   maxiters_1d_neg=3,  # 1D spectrum clipping
                    sigma_1d_neg=3,    maxiters_1d_neg=3,  # 1D spectrum clipping
                    wave_min=None, wave_max=None,
                    ymin=None, ymax=None):
        
        # Plot the 1D extraction x1d vs. indices, same as s2d array above
        self.ax_1d.axhline(0, c='0.50', lw=0.5, alpha=0.66, ls='-')
        self.ax_1d.step(np.arange(self.num_waves), self.x1d_fluxerr, lw=0.5, c='r', alpha=0.66)
        self.ax_1d.step(np.arange(self.num_waves), self.x1d_flux, lw=1)
        self.ax_1d.set_xlim(0, self.num_waves)
        self.ax_1d.yaxis.set_ticks_position('both')
        
        #def set_tick_labels(self, major_tick_interval=0.2):
        # Tick Labels
        #wave_min, wave_max = self.x1d_wave[0], self.x1d_wave[-1]
        #self.wave_min = wave_min or self.x1d_wave[0]
        #self.wave_max = wave_max or self.x1d_wave[-1]
        #wave_min, wave_max = self.x1d_wave[0], self.x1d_wave[-1]
        eps = 1e-7
        #major_tick_interval = 0.2  # 0.5
        xtick_min = (np.ceil((self.wave_min - eps) / wave_tick_interval) * wave_tick_interval)
        xticks = np.arange(xtick_min, self.wave_max, wave_tick_interval)
        # xticks = np.append(xticks, [5.2])
        #xtick_pos = np.interp(xticks, self.x1d_wave, np.arange(self.num_waves))
        xtick_pos = np.interp(xticks, self.x1d_wave, np.arange(len(self.x1d_wave)))
        xtick_labels = [wave_tick_fmt % xtick for xtick in xticks]
        
        if sigma_1d_pos:
            sigma_clipped_data = sigma_clip(self.x1d_flux, sigma=sigma_1d_pos, maxiters=maxiters_1d_pos)
            ymax_1d = np.max(sigma_clipped_data) * (1 + ymargin_1d_pos)

            sigma_clipped_data = sigma_clip(self.x1d_flux, sigma=sigma_1d_neg, maxiters=maxiters_1d_neg)
            ymin_1d = np.min(sigma_clipped_data)

            if ymin_1d > 0:
                ymin_1d = -ymargin_1d_neg * ymax_1d
            
            print('1D limits:', ymin_1d, ymax_1d)
            self.ax_1d.set_ylim(ymin_1d, ymax_1d)

        self.ax_1d.set_xlim(self.xmin, self.xmax)
        if ymin or ymax:
            self.ax_1d.set_ylim(ymin, ymax)
    
        #self.ax_1d.set_title('X1D Flux (Jy)')
        plt.xticks(xtick_pos, xtick_labels)
        plt.xlabel('wavelength (microns)')    
        self.ax_1d.yaxis.set_ticks_position('right')
        
        # Remove overlapping tick labels
        self.fig.canvas.draw()
        renderer = self.fig.canvas.get_renderer()
        tick_labels = self.ax_1d.get_xticklabels()
        prev_bbox = None
        pruned_xtick_labels = []
        for label in tick_labels:
            bbox = label.get_window_extent(renderer)
            if prev_bbox is None or not bbox.overlaps(prev_bbox):
                pruned_xtick_labels.append(label)
                prev_bbox = bbox
            else:
                label.set_visible(False)
    
    
    def master_background(self, background_wave, background_flux):
        self.background = np.interp(slit_spec.x1d_wave, background_wave, background_flux)    
        self.x1d_flux = self.x1d_flux - self.background
        #self.s2d_data = self.s2d_data - self.background[np.newaxis,:]    
        
        
    def save_plot(self, outfile=None, save_subdir='plots', save_dir=None, just_get_name=False):
        save_dir = save_dir or os.path.join(self.s2d_dirname, save_subdir)
        os.makedirs(save_dir, exist_ok=True)

        if not outfile:
            outfile = self.s2d_filename.replace('_s2d.fits', '.png')

        outfile = os.path.join(save_dir, outfile)
        
        if just_get_name:
            return outfile
        else:
            print('SAVING', outfile)
            plt.savefig(outfile, dpi=200, bbox_inches='tight')

        
    def extract_slit_from_table(self): # shutter_table, source_table, source_id, slit_to_sky):
        # WCS transformations
        slit_wcs    = self.slit_model.meta.wcs
        det_to_sky  = slit_wcs.get_transform('detector', 'world') # coordinate transform from detector pixels to sky 
        slit_to_sky = slit_wcs.get_transform('slit_frame', 'world')
        sky_to_slit = slit_wcs.get_transform('world', 'slit_frame')
        self.assign(slit_wcs, det_to_sky, slit_to_sky, sky_to_slit)
                
        #num_dithers = np.max(self.shutter_table['dither_point_index'])

        source_shutter_table = filter_table(self.shutter_table, source_id=self.source_id, 
                                            dither_point_index=1) #, msa_metadata_id=1)
                                            
        if not source_shutter_table:
            variables = 'i_primary, dy_columns, dx_obs, dy_obs, estimated_source_in_shutter_x, estimated_source_in_shutter_y'.split(', ')
            self.assign_none(variables)  # set them all to None
            self.dx_obs = self.dy_obs = 0
            print('Source %d not found in MSA metafile shutter table' % self.source_id)
            return

        #msa_metadata_ids = sorted(list(set(source_shutter_table['msa_metadata_id'])))
        #msa_metadata_id = int(msa_metadata_ids[1])
        #source_shutter_table = filter_table(source_shutter_table, msa_metadata_id=self.msa_metadata_id)
        source_shutter_table = filter_table(source_shutter_table)
        source_metadata_ids = list(set(source_shutter_table['msa_metadata_id']))
        if len(source_metadata_ids) > 1:
            print('Source %d has %d metadata IDs:' % (self.source_id, len(source_metadata_ids)))
            print(source_metadata_ids)
            print('Using the first one: %d' % source_metadata_ids[0])
            source_shutter_table = filter_table(source_shutter_table, msa_metadata_id=source_metadata_ids[0])

        num_dithers = np.max(source_shutter_table['dither_point_index'])

        i_primary = 0
        primary_source_list = list(source_shutter_table['primary_source'])
        if len(primary_source_list):
                if 'Y' in primary_source_list:
                    i_primary = primary_source_list.index('Y')

        select_source_table = self.source_table[self.source_table['source_id'] == self.source_id]
        source_ra  = select_source_table['ra'][0]
        source_dec = select_source_table['dec'][0]

        estimated_source_in_shutter_x = source_shutter_table['estimated_source_in_shutter_x'][i_primary]
        estimated_source_in_shutter_y = source_shutter_table['estimated_source_in_shutter_y'][i_primary]

        # Shift coordinate centroid to (0,0)
        estimated_source_in_shutter_x -= 0.5
        estimated_source_in_shutter_y -= 0.5

        # Coordinates are actually for full slit (not just open area)
        estimated_source_in_shutter_x *= x_scale_open_to_full
        estimated_source_in_shutter_y *= y_scale_open_to_full

        # Transform to sky (RA,Dec) using S2D WCS transformation
        estimated_source_ra, estimated_source_dec, zero = self.slit_to_sky(
            estimated_source_in_shutter_x, estimated_source_in_shutter_y, 0)

        # Transform to image pixels (x,y) using image WCS
        estimated_source_coordinates = SkyCoord(ra=estimated_source_ra*u.deg, dec=estimated_source_dec*u.deg)
        estimated_source_x, estimated_source_y = image_wcs.world_to_pixel(estimated_source_coordinates)

        # Calculate offset between S2D WCS transformation and actual coordinates in input catalog
        # The input catalog coordinates (RA,Dec) are more accurate than the pipeline (RA,Dec)
        source_coordinates = SkyCoord(ra=source_ra*u.deg, dec=source_dec*u.deg)
        source_xy = source_x, source_y = image_wcs.world_to_pixel(source_coordinates)

        # We'll correct for this offset below
        dx_obs = estimated_source_x - source_x
        dy_obs = estimated_source_y - source_y

        # Indices to iterate along slit with multiple shutters
        # (e.g., dy_columns = -1,0,1 for 3-shutter slitlet)
        dx_rows    = source_shutter_table['shutter_row']    - source_shutter_table['shutter_row'][i_primary]
        dy_columns = source_shutter_table['shutter_column'] - source_shutter_table['shutter_column'][i_primary]

        print(dy_columns, 'dy_columns')

        # Scale to full slit (not just open area)
        # No, do this later
        #dx_rows    = np.array(dx_rows)    * x_scale_open_to_full
        #dy_columns = np.array(dy_columns) * y_scale_open_to_full

        # Note the cross-dispersion direction is defined as columns in the MSA metafile
        # even though we normally show them as rows in the MSA

        self.assign(i_primary, dy_columns, dx_obs, dy_obs, estimated_source_in_shutter_x, estimated_source_in_shutter_y, num_dithers)
        #return i_primary, dy_columns, dx_obs, dy_obs, estimated_source_in_shutter_x, estimated_source_in_shutter_y

        
    def plot_image(self, figsize=(12,6), height_ratios=[1, 3], x_slit_hi=1.2, slit_color = 'r', 
                   full_slit_color = 'gray', source_color = (0, 0.8, 0), slit_bkg_color = (0.88, 0.7, 0.7),
                   y_slit_lo_min=3, y_slit_hi_max=-2):

        # Extract slit from tables
        #slit_list = 
        #self.extract_slit_from_table() #self.shutter_table, self.source_table, self.source_id, slit_to_sky)
        #if slit_list and plot_image:
        #i_primary, dy_columns, dx_obs, dy_obs, \
        #estimated_source_in_shutter_x, estimated_source_in_shutter_y = slit_list        

        # Make FIGURE
        fig = plt.figure(figsize=figsize)
        fig_width, fig_height = figsize
        grid = plt.GridSpec(2, 2, width_ratios=[1, 12 * 1.2/x_slit_hi * 6/fig_height], 
                        height_ratios=height_ratios, wspace=0.02, hspace=0.12)

        ax_image = fig.add_subplot(grid[:, 0])
        ax_2d    = fig.add_subplot(grid[0, 1])
        ax_1d    = fig.add_subplot(grid[1, 1], sharex=ax_2d)

        # IMAGE
        ax_image.axis('off')
    
        y_s2d = np.arange(self.s2d_data.shape[0]) # grid of pixel y indices: spatial cross-dispersion
        x_s2d = y_s2d * 0  # dispersion direction irrelevant for RA, Dec
        ra_s2d, dec_s2d, s2d_waves = self.det_to_sky(x_s2d, y_s2d) # RA, Dec, wavelength (microns) for each pixel 
        slit_x, slit_y, zero = self.sky_to_slit(ra_s2d, dec_s2d, 0)  # slit coordinates

        # Define the image y limits and make sure it's tall enough
        y_slit_lo = max(slit_y[0], y_slit_lo_min)
        y_slit_hi = min(slit_y[-1], y_slit_hi_max)
        print('y_slit:', y_slit_lo, y_slit_hi)

        # x_slit_hi = 1.2  # 1.0 is just the slit + half bar; add extra for padding
        x_slit_lo = - (x_slit_hi - 1)
        x_slit_lo, x_slit_hi = x_slit_hi, x_slit_lo  # need to transpose?

        # Create a high-resolution coordinate grid in the slit plane
        nx_slit_image = 100
        ny_slit_image = nx_slit_image * full_slit_aspect * (y_slit_hi - y_slit_lo)
        y_slit, x_slit = np.mgrid[y_slit_lo:y_slit_hi:ny_slit_image*1j, x_slit_lo:x_slit_hi:nx_slit_image*1j] - 0.5
        x_slit *= x_scale_open_to_full
        y_slit *= y_scale_open_to_full
        ny_slit, nx_slit = y_slit.shape

        # Transform these to the image plane
        ra_slit, dec_slit, zero = self.slit_to_sky(x_slit, y_slit, 0)
        coords_slit = SkyCoord(ra=ra_slit*u.deg, dec=dec_slit*u.deg)
        x_image, y_image = image_wcs.world_to_pixel(coords_slit)

        # Extract the image values (colors) at each coordinate
        slit_stamp = NIRCam_image[np.round(y_image-self.dy_obs).astype(int), 
                                       np.round(x_image-self.dx_obs).astype(int)]
        if np.all(slit_stamp == slit_stamp[0]):  # all the same value probably means no data
            slit_stamp = 128 + 0 * slit_stamp  # make gray for no data
            
        slit_extent = x_slit[0,0], x_slit[0,-1], y_slit[0,0], y_slit[-1,0]
        #print('ymin, xmin', np.min(y_image-dy_obs), np.min(x_image-dx_obs))
        #print(source_id, i_primary, dy_columns, dx_obs, dy_obs)

            
        # Plot image
        ax_image.imshow(slit_stamp, origin='lower', aspect=open_slit_aspect, extent=slit_extent)
    
        if self.num_dithers == 3:
            dy_columns_extended = list(self.dy_columns) + [self.dy_columns[0] - 1] +  [self.dy_columns[-1] + 1]
        else:
            dy_columns_extended = self.dy_columns

        print(self.num_dithers, 'dithers')
        for dy in dy_columns_extended:
            #print('dy', dy)
            if between(self.dy_columns[0], dy, self.dy_columns[-1]):
                color = slit_color
            else:
                color = slit_bkg_color
            xy_corners = np.array([open_slit_x_corners, open_slit_y_corners + dy * y_scale_open_to_full]).T
            patch = matplotlib.patches.Polygon(xy_corners, fc='None', ec=color, alpha=1, zorder=100)
            ax_image.add_patch(patch)

        ax_image.plot(self.estimated_source_in_shutter_x, 
                      self.estimated_source_in_shutter_y, 'o', mec=source_color, mfc='None')

        self.assign(fig, ax_image, ax_2d, ax_1d)
        
    def show_MOS_spectrum(self, source_id=None, show_image=True,
            symmetric_range=True, wave_min=None, wave_max=None, wave_tick_interval=0.2, wave_tick_fmt='%.1f',
            save=False, save_subdir='plots', save_dir=None, output_dir=None, figsize=(12,6), height_ratios=[1,3], 
            detection_snr=3, x_slit_hi=1.2, title=None, line_color=(1, 0.85, 0.7), line_label_color='k'):
        self.load_source(source_id)
        self.expand_wavelength_gap()
        output_dir = output_dir or self.s2d_dirname
        self.load_msa_metafile(output_dir)
        self.extract_slit_from_table()
        if show_image and self.dy_obs != None:
            # color image, 2D & 1D spectra
            self.plot_image(figsize=figsize, height_ratios=height_ratios, x_slit_hi=x_slit_hi)
        else:
            # no color image, just 2D & 1D spectra
            print('Slits not defined in MSA metafile')
            self.make_figure(figsize=figsize, height_ratios=height_ratios)
        
        self.set_wavelength_range(wave_min=wave_min, wave_max=wave_max)
        self.plot_spec2d(symmetric_range=symmetric_range, title=title)
        self.plot_spec1d(wave_tick_interval=wave_tick_interval, wave_tick_fmt=wave_tick_fmt)
        if self.z != None:
            self.label_spec_lines(line_color=line_color, line_label_color=line_label_color, detection_snr=detection_snr)
        if save:
            outfile = None
            if type(save) == str:
                outfile = save
            self.save_plot(outfile=outfile, save_subdir=save_subdir, save_dir=save_dir)

    def show_MOS_CAL(self, ext='SCI', source_id=None, vmin=None, vmax=None):
        fig, axs = plt.subplots(3, 1, figsize=(15,7.5), squeeze=True, sharex=True)
        units = 'MJy/sr'
        
        for iplot, cal_model in enumerate(self.cal_models):
            #cal_model = cal_models[iplot]

            if 'slits' in list(cal_model):  # s2d has all the objects; extract the one with source_id
                source_ids = [slit.source_id for slit in cal_model.slits]
                i_slit = source_ids.index(source_id)
                slit_model = cal_model.slits[i_slit]
            else:  # s2d only has one object
                slit_model = cal_model
                i_slit = 0
            
            if ext == 'SCI':
                data = slit_model.data + 0
                cmap = 'RdBu'
                bad_color = 'w'
            elif ext == 'ERR':
                data = slit_model.err + 0
                cmap = 'viridis'
                bad_color = 1, 0.7, 0.7
            elif ext == 'DQ':
                data = slit_model.dq + 0
                data = np.log2(data)  # 'white' is good: no flags set: log(0) = -inf
                vmin, vmax = 0, 32
                cmap = 'turbo'
                bad_color = 'w'
                units = 'flag bit'
                
            if vmin == vmax == None:
                sigma_clipped_data = sigma_clip(data, sigma=5, maxiters=3)
                vmin = np.min(sigma_clipped_data)
                vmax = np.max(sigma_clipped_data)

            cmap = matplotlib.colormaps[cmap]
            cmap.set_bad(bad_color, 1.)
                    
            left   = slit_model.xstart
            right  = slit_model.xsize + left
            bottom = slit_model.ystart
            top    = slit_model.ysize + bottom
            extent = left, right, bottom, top
                            
            ax = axs[iplot]
            im = ax.imshow(data, extent=extent, vmin=vmin, vmax=vmax, origin='lower', cmap=cmap, aspect='auto', interpolation='nearest') # , norm=norm
            plt.colorbar(im, ax=ax, label=units, aspect=10, pad=0.02)
            ax.set_aspect(5)
            
            #ax.set_aspect(1)
            #ax.set_xlim(900, 1150)

            filename = cal_model.meta.filename
            filename_ext = filename.split('_')[-1].split('.')[0]
            title = filename_ext.upper()
            if ext != 'SCI':
                title += ' ' + ext
            title += '  ' + filename
            title += '  slit %d'   % slit_model.slitlet_id
            title += '  source %d' % slit_model.source_id
            ax.set_title(title)