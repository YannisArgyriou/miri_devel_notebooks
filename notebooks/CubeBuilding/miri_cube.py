"""
NAME:
miri_cube

PURPOSE:
Builds a data cube from CV ground test data or mirisim simulated
data. Assumes that WCS keywords have already been added to the data, either using
mmrs_cv_preprocess.pro (or Jane's equivalent routine) or mmrs_mirisim_preprocess

Can be run on Lvl2b slope data

Can be run to only produce a single image slice using /imonly
keyword.

Can be stopped at a particular x,y,z location in the cube for
debugging purposes by specifying slice, stopx, and stopy

Note that the actual heavy lifting code is mmrs_cube.py;
miri_cube is really the calling script.

CALLING SEQUENCE:
miri_cube

INPUTS:

OPTIONAL INPUTS:

OUTPUT:
Data cubes and slices thereof

OPTIONAL OUTPUT:

COMMENTS:
Works with CV2, CV3, and mirisim data instead of requiring
seperate codes for each one.

/cvint is Lvl2 data from CV testing that is neither
ramp nor Lvl2b format, and has no SCI extension

EXAMPLES:

BUGS:

PROCEDURES CALLED:

INTERNAL SUPPORT ROUTINES:

REVISION HISTORY:
Early 2016    Written by David Law (dlaw@stsci.edu)
Oct 2016:     Changed sign on RA WCS, was wrong
Dec 2016:     Update for cube testing with Jane
07-Mar-2017   Ported to miri_cube from cv3cube and simcube.  Update to 1-index
            distortion mapping.
07-Jun-2017:  Update to use x,y in 0-indexed detector frame to
           match python code

------------------------------------------------------------------------------"""

import tel_tools

import os
import numpy as np
from astropy.io import fits
from mmrs_cube import mmrs_cube
from mmrs_xytoabl import mmrs_xytoabl
from mmrs_abtov2v3 import mmrs_abtov2v3
from scipy.optimize import curve_fit

import matplotlib.pyplot as plt

def twoD_Gaussian_ravel( xy, amplitude, xo, yo, sigma_x, sigma_y, theta, offset ):
    x,y = xy
    xo = float(xo)
    yo = float(yo)
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo)
                            + c*((y-yo)**2)))
    return g.ravel()

def twoD_Gaussian_noravel( xy, amplitude, xo, yo, sigma_x, sigma_y, theta, offset ):
    x,y = xy
    xo = float(xo)
    yo = float(yo)
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo)
                            + c*((y-yo)**2)))
    return g

def miri_cube(file,type,band,imonly=None,cvint=None,slice=None,stopx=None,stopy=None):

    channel=int(band[0])
    if (channel == 1): det_name='MIRIFUSHORT'
    elif (channel == 2): det_name='MIRIFUSHORT'
    elif (channel == 3): det_name='MIRIFULONG'
    elif (channel == 4): det_name='MIRIFULONG'
    subband = band[1] # A,B, or C
    if (subband == 'A'): subband_name='SHORT'
    elif (subband == 'B'): subband_name='MEDIUM'
    elif (subband == 'C'): subband_name='LONG'

    indir = os.path.dirname(os.path.realpath(file))

    outdir=indir+'/stack/'

    if os.path.exists(outdir) is False:
        os.system('mkdir {}'.format(outdir))
    outcube=outdir+'cube.fits'
    outslice=outdir+'slice.fits'
    outcollapse=outdir+'collapse.fits'

    hdr=fits.open(file)[0].header
    thisdet=hdr['DETECTOR'] # MIRIFUSHORT or MIRIFULONG

    # CV data use DGAA_POS and DGAB_POS
    if (type == 'cv'):
        thisband=hdr['DGAA_POS'] # Assume no crossed-setups.  SHORT,MEDIUM,LONG
    # mirisim data use BAND
    if (type == 'mirisim'):
        thisband=hdr['BAND']

    # Loop over inputs reading dither positions
    hdu0=fits.open(file)[0].header
    hdu1=fits.open(file)[1].header

    ny = hdu1['NAXIS2']

    raref= hdu0['RA_REF']
    decref= hdu0['DEC_REF']
    v2ref= hdu0['V2_REF']
    v3ref= hdu0['V3_REF']
    rollref= hdu0['ROLL_REF']

    # Band-specific information about pixel/slice size
    # and cube-building parameters
    if band[0] == '1':
        pwidth=0.196 # pixel size along alpha in arcsec
        swidth=0.176 # slice width in arcsec
        xmin = 8 # Minimum x pixel
        xmax = 509 # Maximum x pixel

        # Output cube parameters
        rlim_arcsec=0.1 # in arcseconds
        rlimz_mic=0.0025
        ps_x=0.13 # arcsec
        ps_y=0.13 # arcsec
        ps_z=0.0025 # micron

    elif band[0] == '2':
        pwidth=0.196 # pixel size along alpha in arcsec
        swidth=0.277 # slice width in arcsec
        xmin=520 # Minimum x pixel
        xmax=1020 # Maximum x pixel

        # Output cube parameters
        rlim_arcsec=0.1 # in arcseconds
        rlimz_mic=0.004
        ps_x=0.13 # arcsec
        ps_y=0.13 # arcsec
        ps_z=0.004 # micron

    elif band[0] == '3':
        pwidth=0.245 # pixel size along alpha in arcsec
        swidth=0.387 # slice width in arcsec
        xmin=510 # Minimum x pixel
        xmax=1025 # Maximum x pixel

        # Output cube parameters
        rlim_arcsec=0.1 # in arcseconds
        rlimz_mic=0.004
        ps_x=0.1 # arcsec
        ps_y=0.1 # arcsec
        ps_z=0.002 # micron

    elif band[0] == '4':
        pwidth=0.273 # pixel size along alpha in arcsec
        swidth=0.645 # slice width in arcsec
        xmin=14 # Minimum x pixel
        xmax=480 # Maximum x pixel

        # Output cube parameters
        rlim_arcsec=0.4 # in arcseconds
        rlimz_mic=0.004
        ps_x=0.2 # arcsec
        ps_y=0.2 # arcsec
        ps_z=0.002 # micron

    # Ramps data are not pixel area corrected, while Lvl2b
    # data are
    parea=1.0 # pixel area

    # Define 0-indexed base x and y pixel number
    basex,basey = np.meshgrid(np.arange(1032),np.arange(1024))
    # Convert to base alpha,beta,lambda
    basealpha,basebeta,baselambda,slicenum = mmrs_xytoabl(band)

    basex = basex[:,xmin:xmax]
    basey = basey[:,xmin:xmax]
    basealpha = basealpha[:,xmin:xmax]
    basebeta = basebeta[:,xmin:xmax]
    baselambda = baselambda[:,xmin:xmax]
    slicenum = slicenum[:,xmin:xmax]

    # Crop to only pixels on a real slice
    index0=np.where(slicenum > 0)
    slicenum=slicenum[index0]
    basex=basex[index0]
    basey=basey[index0]
    basebeta=basebeta[index0]
    basealpha=basealpha[index0]
    baselambda=baselambda[index0]
    npix=len(slicenum)

    # Convert all alpha,beta base locations to v2,v3 base locations
    basev2,basev3,xan,yan = mmrs_abtov2v3(basealpha,basebeta,band)

    # Create a master vector of fluxes and v2,v3 locations
    nfiles = 1
    master_flux=np.zeros(npix*nfiles)
    master_ra=np.zeros(npix*nfiles)
    master_dec=np.zeros(npix*nfiles)
    master_lam=np.zeros(npix*nfiles)
    master_expnum=np.zeros(npix*nfiles)
    master_dq=np.zeros(npix*nfiles)
    # Extra vectors for debugging
    master_detx=np.zeros(npix*nfiles) # 0-indexed
    master_dety=np.zeros(npix*nfiles) # 0-indexed
    master_v2=np.zeros(npix*nfiles)
    master_v3=np.zeros(npix*nfiles)

    # Loop over input files reading them into master vectors
    hdulist = fits.open(file)

    thisimg=hdulist[0].data[0,:,:]
    thisdq=hdulist[0].data[2,:,:]

    # If dimensionality is not 2, something went wrong
    if (len(thisimg.shape) != 2):
        raise KeyError('Error: wrong input file dimensions, is this ramp data?')

    # Crop to correct 1/2 of detector
    thisflux=thisimg[:,xmin:xmax]
    thisdq=thisdq[:,xmin:xmax]
    # Crop to only pixels on real slices
    thisflux=thisflux[index0]
    thisdq=thisdq[index0]
    i=0
    master_flux[i*npix:(i+1)*npix]=thisflux
    master_dq[i*npix:(i+1)*npix]=thisdq
    # Coordinate transform
    ra,dec,newroll = tel_tools.jwst_v2v3toradec(basev2,basev3,hdr=hdu0)

    master_ra[i*npix:(i+1)*npix]=ra
    master_dec[i*npix:(i+1)*npix]=dec
    master_lam[i*npix:(i+1)*npix]=baselambda
    master_expnum[i*npix:(i+1)*npix]=i
    master_detx[i*npix:(i+1)*npix]=basex
    master_dety[i*npix:(i+1)*npix]=basey
    master_v2[i*npix:(i+1)*npix]=basev2
    master_v3[i*npix:(i+1)*npix]=basev3

    # hack
    # newcen0=0.;0.005 *15. / 3600.
    # newcen1=0.;0.14 / 3600.
    # rad=sqrt((master_ra-newcen0)**2+(master_dec-newcen1)**2)*3600.
    # newflux=1./(rad)
    # master_flux=3.

    # Safety case deal with 0-360 range to ensure no problems
    # around ra=0 with coordinate wraparound
    medra=np.median(master_ra)
    wrapind=np.where(abs(master_ra - medra) > 180.)
    nwrap = len(wrapind[0])
    if ((nwrap != 0) & (medra < 180.)):
        master_ra[wrapind] = master_ra[wrapind]-360.
    if ((nwrap != 0) & (medra >= 180.)):
        master_ra[wrapind] = master_ra[wrapind]+360.

    # Declare maxima/minima of the cube range *before* doing any QA cuts for specific exposures
    lmin=min(master_lam)
    lmax=max(master_lam)
    ra_min=min(master_ra)
    ra_max=max(master_ra)
    dec_min=min(master_dec)
    dec_max=max(master_dec)
    dec_ave=(dec_min+dec_max)/2.
    ra_ave=(ra_min+ra_max)/2.

    # Crop any pixels with bad DQ flags
    good=np.where(master_dq == 0)
    master_flux=master_flux[good]
    master_ra=master_ra[good]
    master_dec=master_dec[good]
    master_lam=master_lam[good]
    master_expnum=master_expnum[good]
    master_dq=master_dq[good]
    master_detx=master_detx[good]
    master_dety=master_dety[good]
    master_v2=master_v2[good]
    master_v3=master_v3[good]

    # Trim to eliminate any nan fluxes
    index1 = ~np.isnan(master_flux)
    master_flux=master_flux[index1]
    master_ra=master_ra[index1]
    master_dec=master_dec[index1]
    master_lam=master_lam[index1]
    master_expnum=master_expnum[index1]
    master_dq=master_dq[index1]
    master_detx=master_detx[index1]
    master_dety=master_dety[index1]
    master_v2=master_v2[index1]
    master_v3=master_v3[index1]

    # Tangent plane projection to xi/eta
    xi_min=3600.*(ra_min-ra_ave)*np.cos(dec_ave*np.pi/180.)
    xi_max=3600.*(ra_max-ra_ave)*np.cos(dec_ave*np.pi/180.)
    eta_min=3600.*(dec_min-dec_ave)
    eta_max=3600.*(dec_max-dec_ave)

    # Define cube sizes
    n1a=np.ceil(abs(xi_min)/ps_x)
    n1b=np.ceil(abs(xi_max)/ps_x)
    n2a=np.ceil(abs(eta_min)/ps_y)
    n2b=np.ceil(abs(eta_max)/ps_y)
    cube_xsize=int(n1a+n1b)
    cube_ysize=int(n2a+n2b)

    # Redefine xi/eta minima/maxima to exactly
    # match pixel boundaries
    xi_min = -n1a*ps_x - ps_x/2.
    xi_max = n1b*ps_x + ps_x/2.
    eta_min = -n2a*ps_y - ps_y/2.
    eta_max = n2b*ps_y + ps_y/2.

    xi=3600.*(master_ra-ra_ave)*np.cos(dec_ave*np.pi/180.)
    eta=3600.*(master_dec-dec_ave)
    cube_x=(xi-xi_min-ps_x/2.)/ps_x
    cube_y=(eta-eta_min-ps_y/2.)/ps_y

    racen=ra_ave
    decen=dec_ave
    xcen=n1a
    ycen=n2a

    zrange=lmax-lmin
    cube_zsize=int(np.ceil(zrange/ps_z))
    lamcen=(lmax+lmin)/2.
    lamstart=lamcen-(cube_zsize/2.)*ps_z
    lamstop=lamstart+cube_zsize*ps_z
    cube_z=(master_lam-lamstart)/ps_z # Z output cube location in pixels
    wavevec=np.arange(cube_zsize)*ps_z+min(master_lam)

    # squash factors
    xpsf_arcsec=0.3
    ypsf_arcsec=0.24
    zpsf_micron=1 # NOT DONE
    xpsf=1. # 0.3 # xpsf_arcsec/ps_x
    ypsf=1. # 0.24 # ypsf_arcsec/ps_y
    # Note- arbitrarily squashing x and y too much effectively
    # downranks the importance of the z distancedecreasing spectral
    # resolution to improve spatial resolution.

    # roi
    rlimx=rlim_arcsec/ps_x # in pixels
    rlimy=rlimx # in pixels
    rlimz=rlimz_mic/ps_z
    # (Gives about 1-2 spec elements at each spatial element)
    rlim=[1.3*rlimx,1.3*rlimy,rlimz] # in pixels


    # Scale correction factor is the ratio of area between an input pixel
    # (in arcsec^2) and the output pixel size in arcsec^2
    # The result means that the cube will be in calibrated units/pixel
    # scale=ps_x*ps_y/(parea)
    scale=1.0 # Output is in flux/solid angle

    if ((slice is not None) & (stopx is not None) & (stopy is not None)):
        print('Will stop at {}, {}, {} microns'.format(stopx,stopy,wavevec[slice]))

    # Make images at specified slice
    if slice is None: slice=None
    if stopx is None: stopx=-1
    if stopy is None: stopy=-1
    dim_out = [cube_xsize,cube_ysize,cube_zsize]
    print(dim_out)
    print(rlim)
    print(master_flux)
    im = mmrs_cube(cube_x,cube_y,cube_z,master_flux,master_expnum,dim_out,rlim,xsquash=xpsf,ysquash=ypsf,scale=scale,slice=slice,wtype=2,detx=master_detx,dety=master_dety,detlam=master_lam,stopx=stopx,stopy=stopy)

    if slice is not None:
        slice = im.shape[2]/2
        # Recover gaussian FWHM
        x = np.arange(im.shape[0])
        y = np.arange(im.shape[1])
        x, y = np.meshgrid(x, y)

        initial_guess = (3,im.shape[0]/2,im.shape[1]/2,1,1,0,0)

        popt, pcov = curve_fit(twoD_Gaussian_ravel, (x, y), im[:,:,0].ravel(), p0=initial_guess)
        data_fitted = twoD_Gaussian_noravel((x, y), *popt)

        # plt.figure()
        # plt.imshow(im[:,:,0])
        # plt.contour(x, y, data_fitted, 4, colors='w')
        # plt.show()

        # plt.figure()
        # plt.plot(im[:,:,0][15,:])
        # plt.plot(data_fitted[15,:])
        # plt.show()

        fwhmx=np.round(popt[3]*2.355*ps_x*1e3)/1e3
        fwhmy=np.round(popt[4]*2.355*ps_y*1e3)/1e3
        print('Wavelength (micron): {}'.format(slice*ps_z+min(baselambda)))
        print('X FWHM (arcsec): {}'.format(fwhmx))
        print('Y FWHM (arcsec): {}'.format(fwhmy))

    # Write file
    hdu0 = fits.PrimaryHDU()
    hdu1 = fits.ImageHDU(im.transpose(2,1,0))
    cdarr=np.zeros((2,2))
    cdarr[0,0]=2.77778e-4*ps_x
    cdarr[1,1]=2.77778e-4*ps_y
    # Make astrometry solution, don't forget to put xcen/ycen in
    # 1-indexed FITS convention
    hdu1.header['CD1_1'] = cdarr[0,0]
    hdu1.header['CD1_2'] = 0
    hdu1.header['CD1_3'] = 0
    hdu1.header['CD2_1'] = 0
    hdu1.header['CD2_2'] = 0
    hdu1.header['CD2_2'] = cdarr[1,1]
    hdu1.header['CD3_1'] = 0
    hdu1.header['CD3_2'] = 0
    hdu1.header['CD3_3'] = ps_z
    hdu1.header['CRPIX1'] = xcen+1
    hdu1.header['CRPIX2'] = ycen+1
    hdu1.header['CRPIX3'] = 1
    hdu1.header['CRVAL1'] = racen
    hdu1.header['CRVAL2'] = decen
    hdu1.header['CRVAL3'] = lamstart
    hdu1.header['CDELT3'] = ps_z
    hdu1.header['CUNIT3'] = 'um'
    hdu1.header['CTYPE3'] = 'WAVE'

    hdulist = fits.HDUList([hdu0,hdu1])
    if im.shape[2] == 1:
        hdulist.writeto(outslice,overwrite=True)
    elif im.shape[2] > 1:
        hdulist.writeto(outcube,overwrite=True)

    print('Successful completion of MIRI_CUBE')

miri_cube('/Users/ioannisa/Desktop/python/miri_devel/CV2_data/LVL2/processed/MIRM0363-P0-LONG-4230100846_1_495_SE_2014-08-18T10h28m39_LVL2.fits','cv','1C',imonly=None,cvint=None,slice=200,stopx=None,stopy=None)
