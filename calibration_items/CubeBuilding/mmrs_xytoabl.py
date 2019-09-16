"""
+
NAME:
  mmrs_xytoabl

PURPOSE:
  Convert MRS detector coordinates to MRS local alpha,beta coordinates.
  Convention is that x,y pixel locations follow the JWST pipeline
  convention where the detector has 1032x1024 pixels and (0,0) is
  the middle of the lower left detector pixel.

CALLING SEQUENCE:
  mmrs_xytoabl,x,y,a,b,l,channel,[slicenum=,slicename=,refdir=]

INPUTS:
  x      - X coordinate in 0-indexed pixels
  y      - Y coordinate in 0-indexed pixels
  channel - channel name (e.g, '1A')

OPTIONAL INPUTS:
  refdir - Root directory for distortion files
  /trim  - Return only valid pixels in a slice

OUTPUT:
  a       - Alpha coordinate in arcsec
  b       - Beta coordinate in arcsec
  l       - Lambda coordinate in microns

OPTIONAL OUTPUT:
  slicenum - Slice number (e.g., 11)
  slicename - Slice name (e.g., 211A for ch2, slice 11, sub-band A)

COMMENTS:
  Works with CDP5 delivery files.  Inverse function is mmrs_abtoxy.pro
  Not all input x,y can actually map to alpha,beta because some pixels
  fall between slices.  alpha,beta,lambda for these are set to -999.

  CDP transforms provided by A. Glauser assume a 1-indexed detector
  frame convention, whereas this code uses the JWST pipeline
  convention that (0,0) is the middle of the lower-left detector pixel.
  Therefore this code also does the 0- to 1-indexed transform prior
  to calling the Glauser distortions.

EXAMPLES:

BUGS:

PROCEDURES CALLED:

INTERNAL SUPPORT ROUTINES:

REVISION HISTORY:
  30-July-2015  Written by David Law (dlaw@stsci.edu)
  30-Sep-2015   Include lambda axis (D. Law)
  24-Jan-2016   Update extension names for CDP5 (D. Law)
  07-Jun-2017  Update to use x,y in 0-indexed detector frame to
               match python code
  13-Dec-2017  Update directory path for new STScI-MIRI workspace
  10-Oct-2018  Update directory path for new miricoord structure
-
------------------------------------------------------------------------------"""

import numpy as np
from astropy.io import fits

def mmrs_xytoabl(band,cdpDir=None,slice_transmission='80pc',fileversion = "7B.05.01"):

    if cdpDir is None:
        cdpDir='/Users/ioannisa/Desktop/python/miri_devel/cdp_data/'

    # all MRS distortion files
    distcdp = {}
    distcdp["3C"] = "MIRI_FM_MIRIFULONG_34LONG_DISTORTION_%s.fits" %fileversion
    distcdp["3B"] = "MIRI_FM_MIRIFULONG_34MEDIUM_DISTORTION_%s.fits" %fileversion
    distcdp["3A"] = "MIRI_FM_MIRIFULONG_34SHORT_DISTORTION_%s.fits" %fileversion

    distcdp["1C"] = "MIRI_FM_MIRIFUSHORT_12LONG_DISTORTION_%s.fits" %fileversion
    distcdp["1B"] = "MIRI_FM_MIRIFUSHORT_12MEDIUM_DISTORTION_%s.fits" %fileversion
    distcdp["1A"] = "MIRI_FM_MIRIFUSHORT_12SHORT_DISTORTION_%s.fits" %fileversion

    distcdp["4C"] = distcdp["3C"]
    distcdp["4B"] = distcdp["3B"]
    distcdp["4A"] = distcdp["3A"]

    distcdp["2C"] = distcdp["1C"]
    distcdp["2B"] = distcdp["1B"]
    distcdp["2A"] = distcdp["1A"]


    # import parameters needed for d2c mapping
    if fileversion[:5] == '7B.05':
        from astropy.io import fits
        dist = fits.open(cdpDir+distcdp[band])
        alphaPoly = dist['Alpha_CH{}'.format(band[0])].data
        lambdaPoly = dist['Lambda_CH{}'.format(band[0])].data
        bdel = dist[0].header['B_DEL{}'.format(band[0])]
        bzero = dist[0].header['B_ZERO{}'.format(band[0])]

        slice_idx = int(slice_transmission[0])-1
        sliceMap = dist['Slice_Number'].data[slice_idx,:,:]

    else:
        from miri.datamodels.cdp import MiriMrsDistortionModel12, MiriMrsDistortionModel34
        if band in ["1A", "1B", "1C"]:
            dist = MiriMrsDistortionModel12(cdpDir + distcdp[band])
            alphaPoly = dist.alpha_ch1
            lambdaPoly = dist.lambda_ch1
            sliceMap = dist.slicenumber
            bdel = dist.meta.instrument.bdel1
            bzero = dist.meta.instrument.bzero1

        if band in ["2A", "2B", "2C"]:
            dist = MiriMrsDistortionModel12(cdpDir + distcdp[band])
            alphaPoly = dist.alpha_ch2
            lambdaPoly = dist.lambda_ch2
            sliceMap = dist.slicenumber
            bdel = dist.meta.instrument.bdel2
            bzero = dist.meta.instrument.bzero2

        if band in ["3A", "3B", "3C"]:
            dist = MiriMrsDistortionModel34(cdpDir + distcdp[band])
            alphaPoly = dist.alpha_ch3
            lambdaPoly = dist.lambda_ch3
            sliceMap = dist.slicenumber
            bdel = dist.meta.instrument.bdel3
            bzero = dist.meta.instrument.bzero3

        if band in ["4A", "4B", "4C"]:
            dist = MiriMrsDistortionModel34(cdpDir + distcdp[band])
            alphaPoly = dist.alpha_ch4
            lambdaPoly = dist.lambda_ch4
            sliceMap = dist.slicenumber
            bdel = dist.meta.instrument.bdel4
            bzero = dist.meta.instrument.bzero4

    # create maps with wavelengths, alpha and beta coordinates and pixel size

    channel = int(band[0])
    #> slice numbers in the slice map of the distortion CDP for this band
    sliceInventory = np.unique(sliceMap)
    slicesInBand = sliceInventory[np.where( (sliceInventory >= 100*channel ) & (sliceInventory <100*(channel+1)))]

    #> initialise the maps with wavelengths, alpha and beta coordinates of corners for every pixel
    lambdaMap  = np.zeros(sliceMap.shape)
    alphaMap   = np.zeros(sliceMap.shape)
    betaMap    = np.zeros(sliceMap.shape)

    for ss in slicesInBand:
        s = int(ss - 100*channel)
        #> construct a list of y,x coordinates of detector pixels belonging to slices of this band
        pixels = np.where(sliceMap == ss)
        #> for all pixels within the band, construct arrays with center y,x coordinates, and y,z
        # coordinates of the corners of the pixel
        pixelCtry = pixels[0]
        pixelCtrx = pixels[1]

        # Calculate wavelengths for center of the pixels, following (Eq 3) in MIRI-TN-00001-ETH
        # lambda(x,y) = SUM_i(SUM_j ( (K_lam(i,j)*(x-xs)**j * y**i)))
        lambdas   = np.zeros(len(pixelCtry))
        lp = lambdaPoly[s-1]
        xs = lp[0]
        for i in range(5):
            for j in range(5):
                cIndex = 1 + i*5 + j
                lambdas   = lambdas   + lp[cIndex]*(pixelCtrx-xs)**j * pixelCtry**i
        lambdaMap[pixels]   = lambdas

        #> Calculate alpha coordinate for the corners of the pixels, following (Eq 2) in
        # MIRI-TN-00001-ETH
        # alpha(x,y) = SUM_i(SUM_j ( (K_alpha(i,j)*(x-xs)**j * y**i)))
        alphas   = np.zeros(len(pixelCtry))

        ap = alphaPoly[s-1]
        xs = lp[0]
        for i in range(5):
            for j in range(5):
                cIndex = 1 + i*5 + j
                alphas   = alphas   + ap[cIndex]*(pixelCtrx-xs)**j * pixelCtry**i
        alphaMap[pixels]   = alphas

        #> Calculate beta coordinate for the corners of the pixels, following (Eq 4) in
        # MIRI-TN-00001-ETH
        # Beta(s) = Beta_zero + (s-1) * Delta_Beta
        betas   = bzero + (s-1)*bdel
        betaMap[pixels]   = betas

    return alphaMap,betaMap,lambdaMap,sliceMap-int(band[0])*100
