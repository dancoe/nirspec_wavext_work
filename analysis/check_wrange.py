from stdatamodels.jwst.datamodels import WavelengthrangeModel
m = WavelengthrangeModel('/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/reference_files/wavelengthrange_extended.asdf')
print(m.wavelengthrange[0])
print(type(m.wavelengthrange[0]))

