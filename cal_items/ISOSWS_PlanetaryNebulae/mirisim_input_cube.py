import os
import numpy as np
from astropy.io import fits

# set these
in_spec = 'spec_in_mJy.txt'
in_image = 'N20080702S0123.fits'
pixel_scale = 0.003  # in arcsec

out_filename = 'mirisim_input_cube.fits'
if os.path.isfile(out_filename): os.remove(out_filename)

# sort out the image
in_hdulist = fits.open(in_image)

# normalise the image
image_total = np.sum(np.sum(in_hdulist[0].data, axis=0), axis=0)
norm_image = in_hdulist[0].data / image_total
image_dim = norm_image.shape

# sort out the spectrum
sed_data = np.loadtxt(in_spec)
wave_length = len(sed_data[:,0])

# create the datacube
datacube = np.zeros([wave_length,image_dim[0],image_dim[1]])
for n,row in enumerate(sed_data):
    frame = ((row[1] * norm_image) * 1.e3) * (1/pixel_scale**2)
    datacube[n,:,:] = frame

# add header
# determine the wavelength scale
wave_scale = sed_data[1,0] - sed_data[0,0]

prihdr = fits.Header()
prihdr['CRVAL1'] = 0
prihdr['CRPIX1'] = 516
prihdr['CDELT1'] = pixel_scale
prihdr['CTYPE1'] = 'RA---TAN'
prihdr['CUNIT1'] = 'arcsecond'
prihdr['CRVAL2'] = 0
prihdr['CRPIX2'] = 664
prihdr['CDELT2'] = pixel_scale
prihdr['CTYPE2'] = 'DEC--TAN'
prihdr['CUNIT2'] = 'arcsecond'
prihdr['CRVAL3'] = sed_data[0,0]
prihdr['CRPIX3'] = 0
prihdr['CDELT3'] = wave_scale
prihdr['CTYPE3'] = 'WAVE'
prihdr['CUNIT3'] = 'micron'
prihdr['UNITS'] = 'uJy arcsec-2'

# save fits
out_hdu = fits.PrimaryHDU(datacube,header=prihdr)
hdulist = fits.HDUList([out_hdu])
hdulist.writeto(out_filename)
hdulist.close()
