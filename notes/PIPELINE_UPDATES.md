# Notes and Plans for Pipeline Updates

Set up a branch of the pipeline where we remove the lines of code that throw errors for M-grating NRS2 processing and see what happens if we do the following:
* run the pipeline on default MAST associations (we need to make sure this wouldn't cause new errors for MAST/DMS)
* run the pipeline on NRS2 data without giving it an extended wavelength range (to see if that would lead to confusing errors for users)
* run the pipeline on NRS2 data and give it an extended wavelength range with an override to the wavelengthrange referenfe file, and make sure that we NRS2 gets processed through spec2 so that we will be able to do calibration work.

For IFU data there is a check in `assign_wcs` that will explicitly throw a `NoDataOnDetectorError` when processing NRS2 data for M-gratings and F070LP/G140H:
*  https://github.com/spacetelescope/jwst/blob/eca9c5711265352834cd46a06130c593608a1281/jwst/assign_wcs/nirspec.py#L284-L292

We would like to remove the hard-coded error, so that calibration data can be custom processed on NRS2.  I believe that associations are only made for NRS1 in cases such as these, but we may still want a data driven safeguard to not continue processing NRS2 data if the processed wavelength range doesn't fall on NRS2 (perhaps similar to checking if there are no pixels with valid wavelengths as is done for MOS/FS slits in extract_2d).

Currently the pipeline does support processing MOS and FS data (on both NRS1 and NRS2) when given extended wavelength ranges (though of course the calibrations aren't available).

IFU is the only mode that throws an error when processing NRS2 for M gratings (or F070LP G140H).  I don't know what would be the best way to remove that barrier for cal work at the moment, that might need to be investigated (we might be able to simply remove the check that throws an error, but we'd need to make sure this wouldn't cause unexpected issues for DMS in MAST processing)
