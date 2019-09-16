"""
+
NAME:
  jwst_v2v3toradec

PURPOSE:
  Convert v2,v3 coordinates in a JWST frame to RA,DEC coordinates
  given a JWST attitude matrix (or relevant attitude keywords)
  describing the relative orientation of the V3,V3 and RA,DEC
  reference frames.  These can be derived from JWST file FITS headers.

  This constructs the attitude matrix using the keywords V2_REF,
  V3_REF, RA_REF, DEC_REF, and ROLL_REF where the first four
  associate a fixed reference location in V2,V3 with a location in RA,DEC
  and the ROLL_REF specifies the local roll (defined as the position
  angle of the V3 axis measured from N towards E) of the V2,V3 coordinate
  system at the reference location.

  Note that all v2,v3 locations are given in arcseconds while all
  RA,DEC information is given in degrees

  In normal operation this function computes and uses the full JWST
  attitude matrixit can also be run in a /local approximation
  that neglects the full matrix formalism for a local approximation
  with simpler math.

  The full attitude matrix calculations are based on section 6 of
  technical report JWST-STScI-001550 'Description and Use of
  the JWST Science Instrument Aperture File', author C. Cox.

CALLING SEQUENCE:
  jwst_v2v3toradec,v2,v3,ra,dec,hdr=hdr,V2REF=V2REF,V3REF=V3REF,$
     RAREF=RAREF, DECREF=DECREF, ROLLREF=ROLLREF,NEWROLL=NEWROLL,/local

INPUTS:
  v2       - v2 location in arcseconds
  v3       - v3 location in arcseconds

OPTIONAL INPUTS:
  (note that either hdr or the 5 other keywords MUST be provided)
  hdr    - JWST type FITS header containing the V2_REF, V3_REF,
           RA_REF, DEC_REF, and ROLL_REF keywords
  V2_REF - V2_REF location in arcseconds if hdr not given
  V3_REF - V3_REF location in arcseconds if hdr not given
  RA_REF - RA_REF location in degrees (0 -> 360) if hdr not given
  DEC_REF - DEC_REF location in degrees (-90 -> +90) if hdr not given
  ROLL_REF - ROLL_REF location in degrees (0 -> 360) if hdr not given

  /local  - Use the local approximation algorithms

OUTPUT:
  ra      - Right ascension corresponding to v2 in degrees (0 -> 360)
  dec     - Declination corresponding to v3 in degrees (-90 -> +90)

OPTIONAL OUTPUT:
  NEWROLL - Local roll in degrees (0 -> 360) at the location v2,v3

COMMENTS:
  The input v2,v3 locations may be provided as vectors of values, in
  which case the output ra,dec,NEWROLL will also be vectors.  The
  current loop to do this in the full matrix transforms is a little slow.

EXAMPLES:

BUGS:
  This has not been exhaustively coded to catch all edge cases
  (e.g., exactly at the coordinate system poles, exactly at
  RA 360 degrees, out of bound inputs, etc.).

  Note that not specifying input values makes them default to zero

PROCEDURES CALLED:

INTERNAL SUPPORT ROUTINES:
  jwst_att1()
  jwst_att2()
  jwst_attmatrix()
  jwst_localroll()

REVISION HISTORY:
  12-Apr-2016  Written by David Law (dlaw@stsci.edu)
  17-Oct-2016  Deal with zero inputs, v2/v3 in arcsec (D. Law)
-
------------------------------------------------------------------------------"""

##########################################################
import numpy as np

# JWST M1 attitude matrix (V2 and V3 rotations)
# V2REF and V3REF should be in radians
def jwst_att1(V2REF,V3REF):
    # M1=  a00  a01  a02
    #      a10  a11  a12
    #      a20  a21  a22
    thematrix=np.zeros(3,3)
    thematrix[0,0]=np.cos(V2REF)*np.cos(V3REF)
    thematrix[0,1]=np.sin(V2REF)*np.cos(V3REF)
    thematrix[0,2]=np.sin(V3REF)
    thematrix[1,0]=-np.sin(V2REF)
    thematrix[1,1]=np.cos(V2REF)
    thematrix[1,2]=0.
    thematrix[2,0]=-np.cos(V2REF)*np.sin(V3REF)
    thematrix[2,1]=-np.sin(V2REF)*np.sin(V3REF)
    thematrix[2,2]=np.cos(V3REF)
    return thematrix

##########################################################

# JWST M2 attitude matrix (RA,DEC,ROLL rotations)
# RAREF, DECREF, ROLLREF should be in radians
def jwst_att2(RAREF,DECREF,ROLLREF):
    # M2=  a00  a01  a02
    #      a10  a11  a12
    #      a20  a21  a22
    thematrix=np.zeros(3,3)
    thematrix[0,0]=np.cos(RAREF)*np.cos(DECREF)
    thematrix[0,1]=-np.sin(RAREF)*np.cos(ROLLREF)+np.cos(RAREF)*np.sin(DECREF)*np.sin(ROLLREF)
    thematrix[0,2]=-np.sin(RAREF)*np.sin(ROLLREF)-np.cos(RAREF)*np.sin(DECREF)*np.cos(ROLLREF)
    thematrix[1,0]=np.sin(RAREF)*np.cos(DECREF)
    thematrix[1,1]=np.cos(RAREF)*np.cos(ROLLREF)+np.sin(RAREF)*np.sin(DECREF)*np.sin(ROLLREF)
    thematrix[1,2]=np.cos(RAREF)*np.sin(ROLLREF)-np.sin(RAREF)*np.sin(DECREF)*np.cos(ROLLREF)
    thematrix[2,0]=np.sin(DECREF)
    thematrix[2,1]=-np.cos(DECREF)*np.sin(ROLLREF)
    thematrix[2,2]=np.cos(DECREF)*np.cos(ROLLREF)
return thematrix

##########################################################

# JWST M = (M2 # M1) attitude matrix
# V2REF, V3REF, RAREF, DECREF, ROLLREF should be in radians
def jwst_attmatrix(V2REF,V3REF,RAREF,DECREF,ROLLREF):

    thematrix = jwst_att2(RAREF,DECREF,ROLLREF) # jwst_att1(V2REF,V3REF)

return thematrix

##########################################################

# Compute the local roll (the position angle measured from N
# towards E of the V3 axis) at any V2,V3 given an attitude matrix
# V2, V3 must be in radians, result is in radians
def jwst_localroll(V2,V3,ATTMATRIX):

    X=-(ATTMATRIX[2,0]*np.cos(V2)+ATTMATRIX[2,1]*np.sin(V2))*np.sin(V3)+ATTMATRIX[2,2]*np.cos(V3)
    Y=(ATTMATRIX[0,0]*ATTMATRIX[1,2]-ATTMATRIX[1,0]*ATTMATRIX[0,2])*np.cos(V2)+(ATTMATRIX[0,1]*ATTMATRIX[1,2]-ATTMATRIX[1,1]*ATTMATRIX[0,2])*np.sin(V2)

return np.arctan(Y/X)

##########################################################

# Everything is provided and returned in degrees
# Note that v2 has range 0 -> 360
# v3 has range -90 -> +90
def jwst_v2v3toradec(v2,v3,hdr=None,V2REF=None,V3REF=None,RAREF=None,DECREF=None,ROLLREF=None,NEWROLL=None,local=None)

    # Read in attitude keywords and convert from arcseconds to radians
    # If a header was provided start by using those attitude keywords
    if (hdr is not None):
        thisV2REF=hdr['V2_REF']/3600.*np.pi/180.
        thisV3REF=hdr['V3_REF']/3600.*np.pi/180.
        thisRAREF=hdr['RA_REF']*np.pi/180.
        thisDECREF=hdr['DEC_REF']*np.pi/180.
        thisROLLREF=hdr['ROLL_REF']*np.pi/180.
    # Case where no keywords set
    elif ((V2REF is None) and (V3REF is None) and (RAREF is None) and (DECREF is None) and (ROLLREF is None)):
        # Fail out
        raise KeyError, 'No attitude information provided!'
    elif:
        # If individual keywords not set (or zero), default to zero
        if (V2REF is None): V2REF=0.
        if (V3REF is None): V3REF=0.
        if (RAREF is None): RAREF=0.
        if (DECREF is None): DECREF=0.
        if (ROLLREF is None): ROLLREF=0.
        # If attitude keywords were provided, use them
        thisV2REF=V2REF/3600.*np.pi/180.
        thisV3REF=V3REF/3600.*np.pi/180.
        thisRAREF=RAREF*np.pi/180.
        thisDECREF=DECREF*np.pi/180.
        thisROLLREF=ROLLREF*np.pi/180.

    # If running in /local mode, use the local approximate transform
    if (local is not None):
        dv2=(v2/3600.*np.pi/180.-thisV2REF)*np.cos(thisV3REF) # Offset from V2REF in radians
        dv3=v3/3600.*np.pi/180.-thisV3REF # Offset from V3REF in radians
        dra=dv2*np.cos(thisROLLREF)+dv3*np.sin(thisROLLREF) # Offset from RAREF in radians
        ddec=-dv2*np.sin(thisROLLREF)+dv3*np.cos(thisROLLREF) # Offset from DECREF in radians
        ra=(thisRAREF+dra/np.cos(thisDECREF))*180./np.pi # New RA in degrees
        dec=(thisDECREF+ddec)*180./np.pi # New DEC in degrees
        NEWROLL=thisROLLREF*180./np.pi # New roll (identical to old) in degrees
    # If running in normal mode, use the full attitude matrix transform
    else:
        # Compute the JWST attitude matrix from the 5 attitude keywords
        attmat=jwst_attmatrix(thisV2REF,thisV3REF,thisRAREF,thisDECREF,thisROLLREF)

        # Make empty vectors to hold the output ra,dec,NEWROLL
        ra=np.zeros(len(v2))
        dec=np.zeros(len(v2))
        NEWROLL=np.zeros(len(v2))

        # If the input was a vector, loop over elements in the simplest way
        for i in range(len(v2)):
            # Compute the vector describing the input location
            invector=[np.cos(v2[i]/3600.*np.pi/180.)*np.cos(v3[i]/3600.*np.pi/180.),np.sin(v2[i]/3600.*np.pi/180.)*np.cos(v3[i]/3600.*np.pi/180.),np.sin(v3[i]/3600.*np.pi/180.)]

            # Compute the output vector (np.cos(RA)np.cos(dec),np.sin(RA)np.cos(dec),np.sin(dec))
            # by applying the attitude matrix
            outvector=attmat # invector

            # Split the output vector into RA and DEC components and convert
            # back to degrees
            ra[i]=np.arctan(outvector[1]/outvector[0])*180./np.pi

            # Ensure 0-360 degrees
            if (ra[i] < 0.): ra[i]=ra[i]+360.
            dec[i]=anp.sin(outvector[2])*180./np.pi

            # Compute the local roll at this location and convert
            # back to degrees
            NEWROLL[i]=jwst_localroll(v2[i]/3600.*np.pi/180.,v3[i]/3600.*np.pi/180.,attmat)*180./np.pi
            # Ensure 0-360 degrees
            if (NEWROLL[i] < 0.): NEWROLL[i]=NEWROLL[i]+360.

    return ra,dec
