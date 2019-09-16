"""
NAME:

    INITIALIZE

PURPOSE:

    To set-up the main data structure that holds all the information about
    the radiometric output of the MTS.

EXPLANATION:

    The data structure "Vdata" goes through different modules named after
    the optical component they model (BB, VAS, FW, etc.), and the spectra
    of the MTS Calibration Source and the background is modified according
    to the design and characteristics of such elements.

CALLING SEQUENCE:

    VDATA = INITIALIZE(CONFIG, STT)

INPUTS:

    CONFIG - structure with the configuration of the System (MTS [+ MIRI]).
    STT    - "internal" structure with a bunch of parameters which
             characterize the instruments and the code (internal processing
             variables, physical measures of components, names of files, etc.).

OUTPUTS:

    VDATA  - Structure with following parameters, set to "start"
             values. Spectra are set to a "NaN" value. These are
             the values:

                 - wave : array with wavenlengths (um).
                 - dw : spectral sampling (um).
                 - nw : number of spectral bins.
                 - SPC : spectrum of the Calibration Source
                            (BB or the PSS).
                 - E_SPC : uncertainty in SPC.
                 - spc_id : physical quantity that SPC represents.
                          E : Irradiance, L : Radiance,
                          M : Emittance, P : spectral power.
                 - [SPC_APER]: If the point source is used, flux integrated
                   over an aperture that encloses "all" the flux.
                 - [E_SPEC_APER]: uncertainty in SPC_APER.
                 - units : units of the spectrum (Irradiance)
                          MKS : W m-2 um-1
                          cgs : erg s-1 cm-2 um-1
                          ph : ph s-1 cm-2 um-1
                          Jy : 1.E-26 W m-2 Hz-1
                 - BKGtags
                       List of Background sources.
                 - BKG_id
                       Analogue to spc_id, but for the Background spectra.

                 - BKG : Total Background.
                 - BKG_PNH : Total Background of sources before target
                             pinhole (point-like).
                 - BKG_EXT : Extended Background.
                 - BKG_BBFRAME : Background contribution of the BB-Frame.
                 - BKG_BBSTRAP : BKG. of the BB-Strap.
                 - BKG_VAS : BKG. of the VAS.
                 - BKG_TRGT : BKG. of the Source Scanning Subsystem.
                 - BKG_POM : BKG. of the Pick-Off-Mirror frame.
                 - BKG_POM_m : BKG. of the Pick-Off-Mirror.
                 - BKG_FM[1-4] : BKG. of the Folding Mirror #n
                 - BKG_OBA : BKG. of the Optical Bench Assembly.
                 - SPC_ADUs : spectrum of the "target" in ADU s-1 pix-1.
                        Depends on the MIRI configuration.
                 - E_SPC_ADUs: uncertainty in SPC_ADUs.
                 - BKG_ADUs : "BKG" given in ADU s-1 pix-1.
                 - [FWHM] : size of the "spot" at the MIRI Input plane if the point source
                            is used, in meters.
                 - [FWHM_PIX] : size of the "spot" at the MIRI Input plane if the
                            point source is used, in pixels. It is a 2xnw matrix, to allow
                            for differences in the size, in pixels, of the PSF in "x" and "y"  directions when the MRS is used.
                 - [PX]   : pixel size, in meters PROJECTED at the MIRI Input plane.
                            It is a 2 x nw matrix, to allow for the rectangular
                            pixel size of MRS (the column with the largest dimension is
                            the projection of the slit, and the other is the projection
                            of the detector pixels). If the imager is used, then
                            both columns have the same value.
"""
def INITIALIZE(config, stt):
    import numpy as np
    # If the etalons are "on" the spectral resolution is increased
    # to sample the lines properly.

    if config['EtalonON']: dw = stt['dw_hires']
    else: dw = stt['dw_lores']

    wmin = config['wmin']
    wmax = config['wmax']           # wavelength limits of the spectra

    nw = int(np.floor((wmax-wmin) / dw))+1
                                 # Number of spectral bins

    spc = np.full(nw,0.)        # spectrum of BB
    e_spc = np.full(nw,0.)             # uncertainty of spc
    spc_id = 'None'              # physical quantity that spectrum represents.
                                 # It changes along the algorightms: BEWARE!
    BKG_id = 'None'              # Same as spc_id, but for background sources.
    units = 'MKS'                # units of E = W m-2 (A-1)


    wave = wmin + np.arange(nw) * dw # wavelength vector

    # List of Background contributors:

    BKGtags = ['BKG','BKG_PNH','BKG_EXT','BKG_BBFRAME','BKG_BBSTRAP',\
    'BKG_VAS','BKG_TRGT','BKG_POM','BKG_POM_M','BKG_FM1','BKG_FM2',\
    'BKG_FM3','BKG_FM4','BKG_OBA']

    fkarr = np.full(nw,np.nan)

    # Some parameters must be given values at this early stage.

    px1d = (stt['FOV_Lineal'] / stt['num_px']) * 1.E-3 # meters!
    px = np.full((2,nw),0.)

    if config['MIRI'] == 'MRS' :
        chans = stt['mrs_chans']

        for ichan in stt['mrs_chans']:
            pb = stt['mrs_pbs'][ichan]

            if ichan == chans[0]: pb = [np.min(wave),pb[1]]
            if ichan == chans[-1]: pb = [pb[0],np.max(wave)]

            ixw = np.where((wave >= pb[0]) & (wave <= pb[1]))
            if len(ixw[0]) != 0 :
                px[0,ixw] = stt['mrs_pixrel'][ichan][0] * px1d
                px[1,ixw] = stt['mrs_pixrel'][ichan][1] * px1d
    else:
        px[0,:] = px1d
        px[1,:] = px1d # squared pixels by default (Imager)

    Vdata = {'wave':wave, 'dw':dw, 'nw':nw, 'SPC':spc, 'E_SPC':e_spc,\
    'SPC_APER': fkarr, 'E_SPC_APER': fkarr,\
    'spc_id':spc_id, 'units':units,\
    'BKGtags':BKGtags, 'BKG_id':BKG_id,\
    'BKG':fkarr,'BKG_PNH':fkarr,'BKG_EXT':fkarr,\
    'BKG_BBFRAME':fkarr,\
    'BKG_BBSTRAP':fkarr, 'BKG_VAS':fkarr,\
    'BKG_TRGT':fkarr, 'BKG_POM':fkarr,\
    'BKG_POM_M':fkarr, 'BKG_FM1':fkarr,\
    'BKG_FM2':fkarr, 'BKG_FM3':fkarr,\
    'BKG_FM4':fkarr, 'BKG_OBA':fkarr,\
    'SPC_ADUs':fkarr,\
    'E_SPC_ADUs':fkarr, 'BKG_ADUs':fkarr, \
    'SPC_APER_ADUs':fkarr, 'E_SPC_APER_ADUs':fkarr,\
    'FWHM':fkarr,'FWHM_PIX':[np.transpose(fkarr),np.transpose(fkarr)],\
    'PX':px}

    return Vdata
