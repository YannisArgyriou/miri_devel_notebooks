"""
NAME:

    VAS

PURPOSE:

    Include the effects of the "Variable Aperture System" (VAS)
    on the beam.

EXPLANATION:

    This function multiplies the input spectrum by the attenuation
    corresponding to the aperture of the VAS. It also provides, optionally,
    the background contribution of the VAS.

CALLING SEQUENCE:

    VDATA = VAS(VDATA, CONFIG, STT)

INPUTS:

    VDATA  - Structure with the spectral information (See initialize.pro
             for more information).
    CONFIG - structure with the configuration of the System (MTS [+ MIRI]).
    STT    - "internal" structure with a bunch of parameters which
             characterize the instruments and the code.

OUTPUTS:

    VDATA  - (after modifications).

OPTIONAL INPUT KEYWORDS:

    IF config.VASEMIT EQ True, the contribution to the background of
    the VAS is also computed. See set_up.pro.


    IF config.corrVAS EQ True, then a correction for the difference
    between commanded and actual aperture of the VAS, as learnt from
    VM2 data, is applied.
"""
def VAS(Vdata, config, stt):
    from numpy import genfromtxt,pi
    from scipy.interpolate import interp1d
    from planck import planck

    VAScomm = config['VASap']
    config['VASCOMM'] = VAScomm

    if VAScomm <= 0:
         VASap = 0.
         config['VASap'] = 0.

    else:
        VASap = VAScomm

    if config['corrVAS']:
        VAScorrf = stt['Opath'] + stt['VAScorrf']
        COMM,F560W = genfromtxt(VAScorrf,skip_header=3,usecols=(0,1),unpack=True)

        config['VASap'] = interp1d(COMM,F560W)(VASap)

    spcP = Vdata['SPC']        # Power

    assert Vdata['spc_id'] == 'P','Expecting units of "P", but got {} instead.'.format(Vdata['spc_id'])

    #spcP0 = spcP
    spcP *= VASap / 100.
    Vdata['SPC'] = spcP
    Vdata['spc_id'] = 'P'

    # BACKGROUND

    # Attenuation of BKG sources

    Vdata['BKG_BBFRAME'] *= VASap / 100.
    Vdata['BKG_BBSTRAP'] *= VASap / 100.

    # VAS BKG

    if config['VASemit']:

        wave = Vdata['wave']

        Vis_VAS = 1. - VASap / 100.

        spcM_VAS = planck(wave * 1.E4, stt['T_VAS']) * stt['epsilon_VAS']

	    #spcM_VAS = PLANCK(wave * 1.E4, 50.) * stt.epsilon_VAS

                                            # erg s-1 cm-2 A-1
        spcM_VAS *= 1.E1                    # W m-2 um-1

        spcP_VAS = spcM_VAS * pi * (stt['ap_collimator']/2. * 1.E-3)**2.
                                            # area of VAS ~ area of collimator
        spcP_VAS *= Vis_VAS

        spcP_VAS *= stt['VASFctr']             # "Fine" tunning parameter.


    	# BEWARE!!
    	# The beam coming from the collimator is collimated, but the VAS
    	# contribution to the background is not, and thus not all of its
    	# emission goes into the IS (as assumed), only a fraction. But as the
    	# VAS emission is quite uncertain, by now we won't add any
    	# geometric factor for the angle subtended by the VAS from the
    	# entrance port of the IS, and let that be adjusted with the
    	# "fiddling" factor "VASFctr".

        Vdata['BKG_VAS'] = spcP_VAS
        Vdata['BKG_id'] = 'P'

    return Vdata
