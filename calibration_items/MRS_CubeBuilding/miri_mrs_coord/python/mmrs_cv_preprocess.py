"""
+
NAME:
  mmrs_cv_preprocess

PURPOSE:
  Add WCS keywords to CV ground test data.  Note that while early
  versions of this code originally changed the FORMAT of the data too,
  applying a very simple minded ramps-to-slopes module.  It no
  longer does so.

  Note also that the WCS keyword addition is unnecessary if
  Jane's code, to do the same, thing has already been run.

CALLING SEQUENCE:
  mmrs_cv_preprocess

INPUTS:
  directory - cv input directory

OPTIONAL INPUTS:

OUTPUT:
  Outputs are the input files cloned into the new directory with
  modified headers.

OPTIONAL OUTPUT:

COMMENTS:
  Works with CV2 and CV3 data.  Used to be CV3-specific.

EXAMPLES:

BUGS:

PROCEDURES CALLED:

INTERNAL SUPPORT ROUTINES:

REVISION HISTORY:
  Early 2016   Written by David Law (dlaw@stsci.edu)
  07-Mar-2017  Overhaul to only modify header information
-
------------------------------------------------------------------------------"""

import os
import numpy as np
from astropy.io import fits

import matplotlib
matplotlib.use('TkAgg')

from tkinter import filedialog
from tkinter import *

def mmrs_cv_preprocess():

    root = Tk()
    root.withdraw()
    filenames = filedialog.askopenfilenames(initialdir = "/", title = "Read Files to Process",filetypes = (("FITS files","*.fits"),("all files","*.*")) )
    nfiles = len(filenames)

    # Pretend we're pointing at RA=45 degrees, dec=0 degrees
    # with local roll 0
    razp=45.
    deczp=0.
    ROLLREF=0.

    # Define V2REF, V3REF at middle of 1A field (alpha=beta=0.0)
    V2REF = -8.3942412*60. # In arcsec
    V3REF = -5.3123744*60. # In arcsec

    for i in range(nfiles):
        filepath = filenames[i]
        # Load fits file
        hdulist=fits.open(filepath)
        hdr=hdulist[0].header

        # Derive the offsets from the XACTPOS, YACTPOS in the headers
        osim_xact=hdr['XACTPOS'] # arcminutes
        osim_yact=hdr['YACTPOS'] # arcminutes

        # Convert to V2,V3 positions of the OSIM point source
        v2=osim_xact*60. # arcseconds
        v3=-(osim_yact+7.8)*60. # arcseconds

        # Calculate the deltav2 and deltav3 relative to the reference point
        dv2 = (v2-V2REF)/3600. * np.cos(V3REF/3600. * np.pi/180.)# Arc offset in degrees
        dv3 = (v3-V3REF)/3600. # Arc offset in degrees
        dra = dv2 * np.cos(ROLLREF*np.pi/180.) + dv3 * np.sin(ROLLREF*np.pi/180.)
        ddec = -dv2 * np.sin(ROLLREF*np.pi/180.) + dv3 * np.cos(ROLLREF*np.pi/180.)
        # Figure out corresponding RA,DEC location
        # (subtract because we're pretending that the telescope moved
        # instead of the point source?)
        RAREF=razp-dra/np.cos(deczp*np.pi/180.)
        DECREF=deczp-ddec

        # Add dither offset keywords
        hdulist[0].header['V2_REF'] = V2REF
        hdulist[0].header['V3_REF'] = V3REF
        hdulist[0].header['ROLL_REF'] = ROLLREF
        hdulist[0].header['RA_REF'] = RAREF
        hdulist[0].header['DEC_REF'] = DECREF

        directory = os.path.dirname(os.path.realpath(filepath))
        filename = filepath.split(directory)[1][1:]

        if os.path.exists(directory+'/processed/') is False:
            os.system('mkdir {}/processed'.format(directory))
        hdulist.writeto(directory+'/processed/'+filename,overwrite=True)

    return

mmrs_cv_preprocess()
