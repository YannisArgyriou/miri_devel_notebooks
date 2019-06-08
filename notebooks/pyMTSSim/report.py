"""
NAME:
       REPORT

PURPOSE:

       To store a hard-copy of the spectra of the target and the background
       sources at the input plane of MIRI (Irradiances).

CALLING SEQUENCE:

       REPORT, VDATA, CONFIG, STT

INPUTS:

       VDATA  - Structure with the spectral information (See initialize.pro
                for more information).
       CONFIG - structure with the configuration of the System (MTS [+ MIRI]).
       STT    - "internal" structure with a bunch of parameters which
                characterize the instruments and the code.

OUTPUTS:

       A fits file with the spectra in columns is produced, its name being:

              - [config.outpath]/[config.spcfile]

OPTIONAL INPUT KEYWORDS:

      None.
"""

def REPORT(Vdata, config, stt):
    """
    Dump text file with:
    wavelength
    target spectrum
    target uncertainty
    BKG
    BKG_PNH
    BKG_EXT
    ** BKG sources
    [target in ADUs]
    [target uncer. in ADUs]
    [BKG in ADUs]
    [aperture-target]
    [uncert. in aperture-target]
    [aperture-target in ADUs]
    [uncert. in aperture-target in ADUs]
    """
    from _mkhdr import _MKHDR
    from astropy.io import fits

    file = config['outpath'] + config['spcfile']

    hdr0 = _MKHDR(Vdata, config, stt)

    tags = ['wave','SPC','E_SPC','SPC_APER','E_SPC_APER','BKG','BKG_PNH','BKG_EXT',
            'BKG_BBFRAME','BKG_BBSTRAP','BKG_VAS','BKG_TRGT','BKG_POM','BKG_POM_M',
            'BKG_FM1','BKG_FM2','BKG_FM3','BKG_FM4','BKG_OBA','SPC_ADUs','E_SPC_ADUs',
            'BKG_ADUs','SPC_APER_ADUs','E_SPC_APER_ADUs']

    cols = []
    for tag in tags:
        cols.append(fits.Column(name=tag,format='D', array=Vdata[tag]))
    hdr1 = fits.BinTableHDU.from_columns(fits.ColDefs(cols))
    hdulist = fits.HDUList([hdr0,hdr1])

    hdulist.writeto(file, overwrite=True)
