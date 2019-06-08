import numpy as np
from astropy.io import fits
from astropy.constants import c
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

# MRS spectral bands
allbands = ['1A','1B','1C','2A','2B','2C','3A','3B','3C','4A','4B','4C']

MRS_bands = {'1A':[4.83,5.82],
    '1B':[5.62,6.73],
    '1C':[6.46,7.76],
    '2A':[7.44,8.90],
    '2B':[8.61,10.28],
    '2C':[9.94,11.87],
    '3A':[11.47,13.67],
    '3B':[13.25,15.80],
    '3C':[15.30,18.24],
    '4A':[17.54,21.10],
    '4B':[20.44,24.72],
    '4C':[23.84,28.82]} # microns

# define paths
workDir = '/Users/ioannisa/Desktop/python/mirisim/stellar_spectra/'
target = 'HD50083'
input_file_LOW = '{}_LOW.txt'.format(target)
input_file_SH = '{}_SH_final.fits'.format(target)
input_file_LH = '{}_LH_final.fits'.format(target)

# load data
wavs_LOW,fluxes_LOW = np.genfromtxt(workDir+input_file_LOW,usecols=(0,1),unpack=True) # wavs in micron,fluxes in mJy
hdulist_SH,hdulist_LH = fits.open(workDir+input_file_SH),fits.open(workDir+input_file_LH)
wavs_SH,fluxes_SH = hdulist_SH[1].data['WAVE'],hdulist_SH[1].data['SPEC'] # wavs in micron, fluxes in Jansky
wavs_LH,fluxes_LH = hdulist_LH[1].data['WAVE'],hdulist_LH[1].data['SPEC'] # wavs in micron, fluxes in Jansky
hdulist_SH.close(),hdulist_LH.close()

# plt.figure()
# plt.plot(wavs_LOW,fluxes_LOW/1000.)
# plt.plot(wavs_SH[~np.isnan(fluxes_SH)],fluxes_SH[~np.isnan(fluxes_SH)])
# plt.plot(wavs_LH[~np.isnan(fluxes_LH)],fluxes_LH[~np.isnan(fluxes_LH)])
# plt.show()

for band in allbands:
    if band in ['1A','1B','1C','2A','2B']:
        input_file = input_file_LOW.replace('_LOW.txt','')
        sel = (wavs_LOW>=MRS_bands[band][0]) & (wavs_LOW<=MRS_bands[band][-1])
        wavs = wavs_LOW[sel]
        fluxes = fluxes_LOW[sel]/1000. # convert mJy back to Jy
    elif band == '2C':
        input_file = input_file_SH.replace('_SH_final.fits','')
        sel = (wavs_SH>=wavs_SH[0]) & (wavs_SH<=MRS_bands['2C'][-1])
        wavs = wavs_SH[sel]
        fluxes = fluxes_SH[sel]
    elif band in ['3A','3B','3C','4A']:
        input_file = input_file_SH.replace('_SH_final.fits','')
        sel = (wavs_SH>=MRS_bands[band][0]) & (wavs_SH<=MRS_bands[band][-1])
        wavs = wavs_SH[sel]
        fluxes = fluxes_SH[sel]
    elif band in ['4B','4C']:
        input_file = input_file_LH.replace('_LH_final.fits','')
        sel = (wavs_LH>=MRS_bands[band][0]) & (wavs_LH<=MRS_bands[band][-1])
        wavs = wavs_LH[sel]
        fluxes = fluxes_LH[sel]

    # convert Jansky to microJansky
    fluxes *= 1e6

    # omit zero wavelengths and NaNs
    sel = np.nonzero(wavs)
    wavs,fluxes = wavs[sel],fluxes[sel]
    sel = ~np.isnan(fluxes)
    wavs,fluxes = wavs[sel],fluxes[sel]

    step = np.diff(wavs)[0]
    wrange = np.arange(wavs[0],wavs[-1],step)
    ifluxes = interp1d(wavs,fluxes)(wrange)

    f = open(workDir+input_file+'_band'+band+'.txt', "w")

    f.write("wavelength flux\n")
    for i in range(len(wrange)):
        f.write("{}\t{}\n".format(wrange[i],ifluxes[i]))
    f.close()
