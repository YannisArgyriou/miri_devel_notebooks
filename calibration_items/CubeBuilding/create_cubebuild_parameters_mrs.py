#  Program determine the Cube Defaults reference file
#  THis program sets the defaults for plate scale and wavelength resolution
import numpy as np
from astropy.io import fits
import datetime

def create_cubebuild_parameters_mrs(workDir=None):

    if workDir is None:
        workDir = '/Users/ioannisa/Desktop/python/miri_devel/notebooks/CubeBuilding/'
    
    instrument = 'MIRI'

    type = 'CUBEPAR'
    modelname = 'FM'
    version = '07.04.00'

    out_file = instrument+'_'+modelname+'_'+type+'_'+version+'.fits'
    out_filepath='/Users/ioannisa/Desktop/'+out_file
    print  'Output filename: {}'.format(out_filepath)

    hdu0 = fits.PrimaryHDU()

    print datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    # add the standard header values

    hdu0.header["DATE"] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    hdu0.header['REFTYPE'] = type
    hdu0.header['DESCRIP'] = 'Default IFU Cube Sampling and ROI size'
    hdu0.header['PEDIGREE'] = 'GROUND'
    hdu0.header['TELESCOP'] = 'JWST'
    hdu0.header['INSTRUME'] = instrument

    hdu0.header['MODELNAM'] = modelname

    hdu0.header['DETECTOR'] = 'N/A'

    hdu0.header['EXP_TYPE'] = 'MIR_MRS'
    hdu0.header['BAND'] = 'N/A'
    hdu0.header['CHANNEL'] = 'N/A'

    hdu0.header['FILENAME'] = out_file
    hdu0.header['USEAFTER'] = '2000-01-01T00:00:00'

    hdu0.header['VERSION'] = version
    hdu0.header['AUTHOR'] = 'Ioannis Argyriou'
    hdu0.header['ORIGIN'] = 'KU Leuven'
    hdu0.header['FILENAME'] = out_file
    hdu0.header['HISTORY'] = 'IFU Cube defaults'
    hdu0.header['HISTORY'] = 'DOCUMENT: TBD'
    hdu0.header['HISTORY'] = 'SOFTWARE: IDL J Morrison create_cubebuild_parameters.pro'
    hdu0.header['HISTORY'] = 'DATA USED: Simulated Data created by Ben Sargent'
    hdu0.header['HISTORY'] = 'DIFFERENCES: The format of the CUBEPARS file was changed'
    hdu0.header['HISTORY'] = 'DIFFERENCES: per advice from David Law.  The new format'
    hdu0.header['HISTORY'] = 'DIFFERENCES: enables greater control over parameters '
    hdu0.header['HISTORY'] = 'DIFFERENCES: that define the cube-building process.'
    hdu0.header['HISTORY'] = 'DIFFERENCES: It also contains parameters that vary '
    hdu0.header['HISTORY'] = 'DIFFERENCES: smoothly as a function of wavelength to be '
    hdu0.header['HISTORY'] = 'DIFFERENCES: used when building multi-band and multi-'
    hdu0.header['HISTORY'] = 'DIFFERENCES: channel cubes.'

    # ********************************************************************************
    #  Make a binary table to hold the values

    #  For MIRI the data has 12 bands defined by CH# & band
    #  For NIRSPEC the data as 11 bands defined by GRISM and FILTER

    eps=0.01# 1e-2

    chan=np.array([1,1,1,2,2,2,3,3,3,4,4,4])
    bnd=np.array([1,2,3,1,2,3,1,2,3,1,2,3]) # ['A','B','C','A','B','C','A','B','C','A','B','C']
    wmin=np.array([4.89,5.65,6.52,7.49,8.65,9.99,11.53,13.37,15.44,17.66,20.54,23.95])
    wmax=np.array([5.75,6.64,7.66,8.78,10.14,11.71,13.48,15.63,18.05,20.92,24.40,28.45])
    pscl=np.array([0.13,0.13,0.13,0.17,0.17,0.17,0.20,0.20,0.20,0.40,0.40,0.40])
    roispat=np.array([0.10,0.10,0.10,0.15,0.15,0.15,0.20,0.20,0.20,0.40,0.40,0.40])
    wsamp=np.array([0.0025,0.0025,0.0030,0.0040,0.0045,0.0050,0.0060,0.0060,0.0080,0.0140,0.0140,0.0140])
    roispec=np.array([0.0025,0.0025,0.0030,0.0040,0.0045,0.0050,0.0060,0.0060,0.0080,0.0140,0.0140,0.0140])
    power=np.array([2,2,2,2,2,2,2,2,2,2,2,2])
    softrad=np.zeros(len(chan))+0.01 # eps**(1/power)

    col1 = fits.Column(name='CHANNEL', format='J', array=chan)
    col2 = fits.Column(name='BAND', format='J', array=bnd)
    col3 = fits.Column(name='WAVEMIN', format='D', unit='microns', array=wmin)
    col4 = fits.Column(name='WAVEMAX', format='D', unit='microns', array=wmax)
    col5 = fits.Column(name='SPAXELSIZE', format='D', unit='arcseconds', array=pscl)
    col6 = fits.Column(name='SPECTRALSTEP', format='D', unit='microns', array=wsamp)
    col7 = fits.Column(name='ROISPATIAL', format='D', unit='arcseconds', array=roispat)
    col8 = fits.Column(name='ROISPECTRAL', format='D', unit='microns', array=roispec)
    col9 = fits.Column(name='POWER', format='D', unit='unitless', array=power)
    col10 = fits.Column(name='SOFTRAD', format='D', unit='unitless', array=softrad)

    coldefs = fits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8, col9, col10])
    hdu1 = fits.BinTableHDU().from_columns(coldefs,name='CUBEPAR_MSM')

    optpsfile=workDir+'optimization/opt_platescale_revised.txt'
    lamlaw,optps = np.genfromtxt(optpsfile,unpack=True,usecols=(0,1))
    optroispatfile=workDir+'optimization/opt_roispatial_revised.txt'
    lamlaw,optroispat = np.genfromtxt(optroispatfile,unpack=True,usecols=(0,1))
    optwsampfile=workDir+'optimization/opt_wavelengthsampling_revised.txt'
    lamlaw,optwsamp = np.genfromtxt(optwsampfile,unpack=True,usecols=(0,1))
    optroispecfile=workDir+'optimization/opt_roispectral_revised.txt'
    lamlaw,optroispec = np.genfromtxt(optroispecfile,unpack=True,usecols=(0,1))

    powercontinuous=np.zeros(len(lamlaw))+2.
    softrad2ext=np.zeros(len(lamlaw))+0.01

    col1 = fits.Column(name='WAVELENGTH',format='D',array=lamlaw,unit='microns')
    col2 = fits.Column(name='ROISPATIAL',format='D',array=optroispat,unit='arcseconds')
    col3 = fits.Column(name='ROISPECTRAL',format='D',array=optroispec,unit='microns')
    col4 = fits.Column(name='POWER',format='D',array=powercontinuous,unit='unitless')
    col5 = fits.Column(name='SOFTRAD',format='D',array=softrad2ext,unit='unitless')

    coldefs = fits.ColDefs([col1, col2, col3, col4, col5])
    hdu2 = fits.BinTableHDU().from_columns(coldefs,name='MULTICHANNEL_MSM')

    hdulist = fits.HDUList([hdu0,hdu1,hdu2])
    hdulist.writeto(out_filepath,overwrite=True)

create_cubebuild_parameters_mrs()
