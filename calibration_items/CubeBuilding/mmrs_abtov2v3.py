"""
+
NAME:
  mmrs_abtov2v3

PURPOSE:
  Convert MRS local alpha,beta coordinates to JWST v2,v3 coordinates

CALLING SEQUENCE:
  mmrs_abtov2v3,a,b,band,[refdir=]

INPUTS:
  a       - Alpha coordinate in arcsec
  b       - Beta coordinate in arcsec
  band    - band name (e.g, '1A')

OPTIONAL INPUTS:
  refdir - Root directory for distortion files

OUTPUT:
  v2      - V2 coordinate in arcsec
  v3      - V3 coordinate in arcsec

OPTIONAL OUTPUT:
  XAN     - XAN coordinate in arcsec
  YAN     - YAN coordinate in arcsec

COMMENTS:
  Works with CDP5 delivery files.  Inverse function is mmrs_v2v3toab.pro

EXAMPLES:

BUGS:

PROCEDURES CALLED:

INTERNAL SUPPORT ROUTINES:

REVISION HISTORY:
  30-Jul-2015  Written by David Law (dlaw@stsci.edu)
  27-Oct-2015  Add conversion to REAL V2,V3 (D. Law)
  16-Nov-2015  Add conversion to XAN,YAN (D. Law)
  24-Jan-2016  Update reference files to CDP5 (D. Law)
  17-Oct-2016  Input/output v2/v3 in arcsec (D. Law)
  13-Dec-2017  Update directory path for new STScI-MIRI workspace
  10-Oct-2018  Update directory path for new miricoord structure
-
------------------------------------------------------------------------------"""

import numpy as np
from astropy.io import fits

def mmrs_abtov2v3(a,b,band,cdpDir=None,fileversion = "06.04.00"):

    if cdpDir is None:
        cdpDir='/Users/ioannisa/Desktop/python/miri_devel/cdp_data/'

    # Ensure we're not using integer inputs
    adbl=a.astype('float')
    bdbl=b.astype('float')

    # Determine input reference FITS file
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

    # Read alpha,beta -> XAN,YAN table
    convtable=fits.open(cdpDir+distcdp[band])['albe_to_XANYAN'].data
    # Determine which rows we need
    v2index=np.where(convtable.field(0)=='T_CH{}_XAN'.format(band))[0][0]
    v3index=np.where(convtable.field(0)=='T_CH{}_YAN'.format(band))[0][0]
    # Trim to relevant v2, v3 rows for this band
    conv_v2=list(convtable[v2index])
    conv_v3=list(convtable[v3index])

    # WARNING index 0 is just text, the coefficients start from index 1
    xan=conv_v2[1]+conv_v2[2]*adbl + conv_v2[3]*bdbl + conv_v2[4]*adbl*bdbl
    yan=conv_v3[1]+conv_v3[2]*adbl + conv_v3[3]*bdbl + conv_v3[4]*adbl*bdbl

    # convert to arcsec from arcmin
    v2=xan*60.
    v3=-(yan+7.8)*60.
    xan=xan*60.
    yan=yan*60.

    return v2,v3,xan,yan
