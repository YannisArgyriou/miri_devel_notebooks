"""
+
NAME:
  mmrs_mirisim_preprocess

PURPOSE:
  Convert mirisim output files to a format that can be further processed
  by DRL cube building code or the JWST pipeline.

  Originally, this meant digging through the logs and mirisim code
  structure to put the WCS reference keywords in the headers.
  Later this information was already in the headers, but was
  incorrectly defined so had to be tweaked.

  Puts the results in a preproc/ directory.

  Note that this code originally changed the FORMAT of the data too,
  applying a very simple minded ramps-to-slopes module.  It no
  longer does so.

CALLING SEQUENCE:
  mmrs_mirisim_preprocess

INPUTS:
  none- prompts for files

OPTIONAL INPUTS:

OUTPUT:
  Outputs are the input files cloned into the new directory with
  modified headers.

OPTIONAL OUTPUT:

COMMENTS:

EXAMPLES:

BUGS:
  The dither-table dependent method hasn't been checked recently.

PROCEDURES CALLED:

INTERNAL SUPPORT ROUTINES:
  mmrs_mirisim_dloc()

REVISION HISTORY:
  Early 2016  Written by David Law (dlaw@stsci.edu)
  17-Oct-2016  Use coordinate keywords in header
  10-Nov-2016  Remove format changes, use arcsec for v2,v3, general overhaul
-
------------------------------------------------------------------------------"""

import os
import numpy as np

# Convert a nominal dither offset in alpha,beta
# to offset in ra, dec
def mmrs_mirisim_dloc(aoff,boff,ROLLREF):

    # Convert alpha,beta offsets to v2, v3
    v2ref,v3ref = mmrs_abtov2v3(0.,0.,'1A')
    v2off,v3off = mmrs_abtov2v3(aoff,boff,'1A')
    v2off=-(v2off-v2ref)
    v3off=-(v3off-v3ref)

    raoff=v2off*np.cos(ROLLREF*np.pi/180.)+v3off*np.sin(ROLLREF*np.pi/180.)
    decoff=-v2off*np.sin(ROLLREF*np.pi/180.)+v3off*np.cos(ROLLREF*np.pi/180.)

return raoff,decoff

#------------------------------------------------------------------------------

def mmrs_mirisim_preprocess(filepath,oldwcs=None):

    # Look in the first file to see if we have header keywords
    hdulist=fits.open(filepath)
    hdu0 = hdulist[0].header
    hdu1 = hdulist[1].header
    try:
        # Translate from Oct 2016 mirisim header keywords
        # to my conventional keywords.  Mirisim currently uses
        # an odd frame with reference point nearly but not quite XAN,YAN
        # that is flipped in the DEC direction.
        # Hard code this for now in a way that will probably fail for
        # ROLL ne 0
        #temp1=-1.138#fxpar(hdr,'V2_REF')*60.
        #temp2=-468.363#fxpar(hdr,'V3_REF')*60.-7.8*60.
        temp1=hdu1['V2_REF']*60.
        temp2=hdu1['V3_REF']*60.-7.8*60.
        temp3=hdu1['RA_REF']
        temp4=-hdu1['DEC_REF']
        temp5=hdu1['ROLL_REF']

        V2REF = -8.3942412*60. # In arcsec
        V3REF = -5.3123744*60. # In arcsec
        ra,dec=jwst_v2v3toradec(V2REF,V3REF,V2REF=temp1,V3REF=temp2,RAREF=temp3,DECREF=temp4,ROLLREF=temp5)

        RAREF=ra[0]
        DECREF=dec[0]
        ROLLREF=temp5
    except KeyError:
        # If not, read in reference dither table
        ditherfile=os.environ['MIRISIM_DIR']+'obssim/data/mrs_recommended_dither.dat'
        # Table of possible dither positions (2 x N)
        dithertable=np.genfromtxt(ditherfile,delimiter=',',unpack=True,usecols=(0,1))

        # Read the log file to pull out key information
        directory = os.path.dirname(os.path.realpath(filename))
        filename = filepath.split(directory)[1]
        logfile=directory+'mirisim.log'
        with open(filename) as f:
            lines = f.readlines()

        idx = []
        for i,line in enumerate(lines):
            if 'Starting index of executed' in line:
                idx.append(i)

        # Look for a line of info about the dither pattern start index (1-indexed)
        dstart=int(lines[idx[0]].split('\n')[0])-1
        # Pretend we're pointing at RA=45 degrees, dec=0 degrees
        # with local roll 0
        razp=45.
        deczp=0.
        ROLLREF=0.
        # Define V2REF, V3REF at middle of 1A field (alpha=beta=0.0)
        V2REF = -8.3942412*60 # In arcsec
        V3REF = -5.3123744*60 # In arcsec

        ditherno=dstart+hdu0['PNTG_SEQ']-1# Add starting point, convert to 0 indexed
        # Add dither offset keywords.
        aoffset=dithertable[ditherno,0]# In alpha
        boffset=dithertable[ditherno,1]# In beta
        # Convert alpha,beta offsets to ra,dec
        dra,ddec = mmrs_mirisim_dloc(aoffset,boffset,ROLLREF)
        # Figure out corresponding RA,DEC location
        # (subtract because we're pretending that the telescope moved
        # instead of the point source)
        RAREF=razp-dra/np.cos(deczp*np.pi/180.)/3600.
        DECREF=deczp-ddec/3600.
    else:   # If there are header keywords, process them

        # Add (or overwrite) these keywords into header
        hdu0.header['V2_REF'] = V2REF
        hdu0.header['V3_REF'] = V3REF
        hdu0.header['ROLL_REF'] = ROLLREF
        hdu0.header['RA_REF'] = RAREF
        hdu0.header['DEC_REF'] = DECREF

        # Modify 0th extension headers
        try:
            os.system('mkdir {}/processed'.format(directory))
        except:
            pass
        hdulist.writeto(directory+'processed/'+filename,overwrite=True)


    return
