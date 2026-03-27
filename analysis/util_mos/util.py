import numpy as np
import os

#import jwst
from jwst import datamodels

import matplotlib  # as mpl
import matplotlib.pyplot as plt

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

import astroquery # To retrieve data from MAST
from astroquery.mast import Observations  # MAST
from astropy.table import Table, vstack, unique
#print('astroquery version', astroquery.__version__)

import requests  # download large file
from tqdm import tqdm  # shows progress bar while downloading file
import gzip  # decompress .gz file
import shutil

tqdm_params = {
    #'desc': url,
    #'total': total,
    'miniters': 1,
    'unit': 'B',
    'unit_scale': True,
    'unit_divisor': 1024,
    'bar_format': '{l_bar}{bar:20}{r_bar}{bar:-20b}',  # progress bar length 20 pixels
}

def show_progress(*args, **kwargs):
    return tqdm(*args, **{**tqdm_params, **kwargs})

# each label will appear only once in legend; this function checks
def label_once_in_legend(label, ax):
    ax.legend()  # generate it to update
    legend = ax.get_legend()
    label_is_in_legend = False
    if legend is not None:
        label_is_in_legend = label in [text.get_text() for text in legend.get_texts()]
        #print(label, label_is_in_legend, legend.get_texts())

    if label_is_in_legend:
        return None  # label is in legend already
    else:
        return label
        

def print_dict(d, indent=0, print_key=False):
    """
    Another way to view dictionary contents.

    Parameters:
    ----------
    d: dict
        Dictionary
    indent: int
        Number of indents (default = 0)

    Return:
    ------
    None
    """

    for key in d.keys():
        val = d[key]
        try:
            # keys = val.keys()
            if print_key:
                print(' ' * indent, key)
            print_dict(val, indent+4)
        except Exception:
            if isinstance(val, list):
                for i, item in enumerate(val):
                    print(' ' * indent, '[%d]' % i)
                    print_dict(item, indent+4)
            else:
                print(' ' * indent, key, ':', d[key])

def all_same(lst):
    if all(x == lst[0] for x in lst):
        return lst[0]
    else:
        return False
                
def between(lo, x, hi):
    """
    Check if a value falls within a specified range.

    Parameters:
    ----------
    lo: float or int
        The lower bound of the range.
    x: float or int
        The value to check if it falls within the specified range.
    hi: float or int
        The higher bound of the range.

    Returns:
    -------
    int: 1 if the value x is between lo and hi (inclusive), 0 otherwise.
    """
    return (lo <= x) * (x <= hi)
    
def roundint(x):
    return int(round(x))

# Helper functions to select subset of files from list
def all_in(elements, list_or_string):
    """
    Check if all elements are present in a list or string.

    Parameters:
    ----------
    elements: str or list
        Elements to check.
    list_or_string: list or str
        List or string to search for elements.

    Returns:
    -------
    bool: True if all elements are present, False otherwise.
    """
    if isinstance(elements, str):
        elements = elements.split()
    for element in elements:
        if element not in list_or_string:
            return False
    return True


def any_in(elements, list_or_string):
    """
    Check if any element is present in a list or string.

    Parameters:
    ----------
    elements: str or list
        Elements to check.
    list_or_string: list or str
        List or string to search for elements.

    Returns:
    -------
        bool: True if any element is present, False otherwise.
    """
    if isinstance(elements, str):
        elements = elements.split()
    for element in elements:
        if element in list_or_string:
            return True
    return False


# Select subset of files containing all search strings
def select_files(all_files, search_strings=[], exclude_strings=[]):
    """
    Select a subset of files containing all search strings.

    Parameters:
    ----------
    all_files: list
        List of all files.
    search_strings: list
        List of search strings to match.

    Returns:
    -------
    chosen_files: list or str
        List of chosen files or a single chosen file if only one matches.
    """
    chosen_files = [file for file in all_files if all_in(search_strings, file)]
    chosen_files = [file for file in chosen_files if not any_in(exclude_strings, file)]
    #if len(chosen_files) == 1:
    #    chosen_files = chosen_files[0]
    return chosen_files
    
def fits_header_table(fits_files, 
    header_keywords = 'FILENAME OBSERVTN VISIT EFFEXPTM EXP_TYPE NOD_TYPE NUMDTHPT PATT_NUM NINTS NGROUPS READPATT'.split()):
    # PROGRAM EXP_TYPE

    extracted_data = {keyword: [] for keyword in header_keywords}

    # Loop through each FITS file
    for fits_file in fits_files:
        # Open the FITS file
        with fits.open(fits_file) as hdu_list:
            # Extract the header from the primary HDU
            header = hdu_list[0].header

            # Get the values of the desired keywords and append to the respective lists
            for keyword in header_keywords:
                value = header.get(keyword, 'N/A')
                extracted_data[keyword].append(value)

    fits_table = Table(extracted_data)
    return fits_table

def single_value(x):  # True = one number; False = multiple numbers (list / tuple / array / set)
    return isinstance(x, (int, float))

def filter_table1(full_table, **kwargs):
    """
    Filters an Astropy Table based an arbitrary number of input column-value pairs.
    Each value can be either a single value or a list (or tuple, array, or set).
    Example:
    select_shutter_table = filter_table(shutter_table, msa_metadata_id=1, dither_point_index=1, source_id=[6355,5144])
    """
    filtered_table = full_table
    for column, value in kwargs.items():
        if single_value(value):
            filtered_table = filtered_table[filtered_table[column] == value]
        else: # list
            filtered_table = filtered_table[[(item in value) for item in filtered_table[column]]]
    return filtered_table

def filter_table(full_table, **kwargs):
    filtered_table = full_table
    for column, value in kwargs.items():
        if isinstance(value, (list, tuple, np.ndarray)):
            filtered_table = filtered_table[np.isin(filtered_table[column], value)]
        else:
            filtered_table = filtered_table[filtered_table[column] == value]
    return filtered_table
    
#def fits_header_filter(fits_files, **kwargs):
def select_fits_headers(fits_files, **kwargs):
    fits_table = fits_header_table(fits_files)
    chosen_fits_files = list(filter_table(fits_table, **kwargs)['FILENAME'])
    fits_dir = os.path.dirname(fits_files[0])  # assume all files are in the same directory?
    chosen_fits_files = [os.path.join(fits_dir, fits_file) for fits_file in chosen_fits_files]
    return chosen_fits_files


# Helper function to download JWST files from MAST
def download_jwst_files(filenames, download_dir, mast_dir='mast:jwst/product'):
    """
    Helper function to download JWST files from MAST.

    Parameters:
    ----------
    filenames: list of str
        List of filenames to download.
    download_dir: str
        Directory where the files will be downloaded.
    mast_dir: str
        MAST directory containing JWST products.

    Returns:
    -------
    downloaded_files: list of str
        List of downloaded file paths.
    """
    # Download data
    downloaded_files = []
    os.makedirs(download_dir, exist_ok=True)
    for filename in filenames:
        filename = os.path.basename(filename)
        mast_path = os.path.join(mast_dir, filename)
        local_path = os.path.join(download_dir, filename)
        if os.path.exists(local_path):
            print(local_path, 'EXISTS')
        else:
            # Can let this command check if local file exists
            # However, it will delete it if it's there
            # and the wrong size (e.g., reprocessed)
            Observations.download_file(mast_path,   local_path=local_path)
        downloaded_files.append(local_path)

    return downloaded_files
    

def decompress_file(filename, decompressed_file=''):
    decompressed_file = decompressed_file or filename[:-3]
    print('Decompressing to:', decompressed_file, '...')
    
    with gzip.open(filename, 'rb') as f_in:
        with open(decompressed_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    # Clean up: remove compressed file, now that you've decompressed it
    if os.path.exists(filename):
        if os.path.getsize(filename):
            os.remove(filename)

def download_large_file(url, filename='', decompress=True):
    if filename:
        if os.path.isdir(filename):
            filename = os.path.join(filename, os.path.basename(url))
        # else just use the filename as is
    else:
        filename = os.path.basename(url)  # if left blank, just save the filename as in the URL

    decompressed_file = filename
    if filename[-3:] == '.gz':
        if decompress:
            decompressed_file = filename[:-3]
    else:
        if url[-3:] == '.gz':
            filename += '.gz'
    
    if os.path.exists(decompressed_file):
        if os.path.getsize(decompressed_file):
            print(decompressed_file, 'EXISTS')
            return decompressed_file
    
    print('Downloading', url)
    print('to:', filename)
    #print('(d)', decompressed_file)

    with open(filename, 'wb') as f:
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            total = int(r.headers.get('content-length', 0))
            tqdm_params['total'] = total
            with tqdm(**tqdm_params) as pb:
                for chunk in r.iter_content(chunk_size=8192):
                    pb.update(len(chunk))
                    f.write(chunk)
                      
    if decompress:
        if filename[-3:] == '.gz':
            decompress_file(filename, decompressed_file)        
            filename = decompressed_file
            
    return filename
    

mast_url_path = 'https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:jwst/product/'

# Download large file from MAST with progress bar
def download_MAST_file(filename, download_dir):
    filename = os.path.basename(filename)
    mast_file_path = os.path.join(mast_url_path, filename)
    local_path = os.path.join(download_dir, filename)
    download_large_file(mast_file_path, local_path)
    
# Create links in output_dir directory to source files
def link_to_files(source_files, output_dir, overwrite=False):
    links = []
    for source_file in source_files:
        link = os.path.join(output_dir, os.path.basename(source_file))
        links.append(link)
        print(link, '->', source_file)
        if os.path.exists(link) and overwrite:
            print('Removing', link)
            os.remove(link)
        if not os.path.exists(link):
            os.symlink(os.path.abspath(source_file), link)
    return links
            
def download_file(url, download_dir=''):  # simple for small files
    filename = os.path.basename(url)
    if download_dir:
        if os.path.isdir(download_dir):
            filename = os.path.join(download_dir, filename)
    if os.path.exists(filename):
        print(filename, 'EXISTS')
    else:
        print('DOWNLOADING', filename)
        r = requests.get(url, allow_redirects=True)
        open(filename, 'wb').write(r.content)
    return filename

def table_with_units(fits_table):
    # Create a new Astropy Table
    table = Table(fits_table.data)
    
    # Apply units to each column
    for col in table.colnames:
        unit_str = fits_table.columns[col].unit
        if unit_str:
            try:
                unit = u.Unit(unit_str)
                table[col].unit = unit
                print(f"Column {col} has unit {unit}")
            except ValueError:
                print(f"Warning: Unrecognized unit '{unit_str}' for {col}. Leaving as dimensionless.")
        else:
            print(f"Warning: No unit specified for {col}. Leaving as dimensionless.")
    
    return table