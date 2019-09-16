"""
NAME:

    IS

PURPOSE:

    Include the effects of the Integrating Sphere on the beam:
            diffusion, attenuation.

EXPLANATION:

    The integration sphere provides a spatially flat, extended
    light source, by multiple reflections in its flat golden
    inner surface.

CALLING SEQUENCE:

    VDATA = IS(VDATA, CONFIG, STT)

INPUTS:

    VDATA  - Structure with the spectral information (See initialize.pro
             for more information).
    CONFIG - structure with the configuration of the System (MTS [+ MIRI]).
    STT    - "internal" structure with a bunch of parameters which
             characterize the instruments and the code.

OUTPUTS:

    VDATA  - (after modifications)

OPTIONAL INPUT KEYWORDS:

    None
"""
def IS(Vdata, config, stt):
    from numpy import pi,tan
    from loadtrans import LoadTrans

    spcP = Vdata['SPC']
    wave = Vdata['wave']

    assert Vdata['spc_id'] == 'P','Expecting units of "P", but got {} instead.'.format(Vdata['spc_id'])

    beam_divergence = config['beam_divergence']

    Area_sphere = 4. * pi * (stt['R_sphere']*1.E-3)**2. # m2
    IS_aperture = stt['ap_collimator'] + tan(beam_divergence) * stt['dist_coll_IS']
                                                             # mm
                                                             # Beam diameter
    # Input = smallest of beam or input port diameter

    assert stt['ap_input_port'] >= IS_aperture,'**WARNING** Input port smaller than beam **WARNING**' # WHAT TO DO??

    Area_input_port = pi * (IS_aperture * 1.E-3 / 2.)**2. # m2 (circular)
    Area_output_port = (stt['ap_output_port'] * 1.E-3)**2.      # m2 (squared)


    # APPLICATION OF THE "TRANSMISSION" CURVE OF THE IS

    iTrans_IS = LoadTrans(stt['Tpath'] + stt['ISTransn'], wave)

    # Two useful ratios:

    fj = (Area_input_port + Area_output_port) / Area_sphere
                                                              # Area efficiency
                                                              # of the IS
    fe = Area_output_port / Area_sphere


    # PREVIOUSLY (2.2<v<2.4)

    #ISFactor = iTrans_IS * fe / (pi * Area_output_port * \
    #           (1. - iTrans_IS * (1. - fj))) # Equivalent to the next expression

    # ISFactor = iTrans_IS / (pi * Area_sphere * # DEDUCED BY YANNIS ARGYRIOU
    #            (1. - iTrans_IS * (1. - fj)))   #

    ISFactor = 1. / (pi * Area_sphere *
               (1. - iTrans_IS * (1. - fj))) # DEDUCED BY RUYMAN AZZOLLINI
                                               # 5% boost wrt previous. (Comment by Yannis: what does "boost" signify (physically) in this context?)

    ISFactor *= stt['ISFctr']                  # scaling


    spcL = spcP * ISFactor                     # RADIANCE [W.m-2.um-1.sr-1]
    spcL *= stt['ro_optica']                   # due to the folding mirror.
                                               # Minor loss.

    # WARNING: Fudge factor introduced to account for discrepancy between the MRS PHOTOM CDP and the MRS PCE CDP values
    if stt['FUDGEFctr'] != 0.:
        print 'Fudge factor of {} introduced to account for discrepancy between the MRS PHOTOM CDP and the MRS PCE CDP values'.format(round(stt['FUDGEFctr'],2))
    spcL *= stt['FUDGEFctr']

    Vdata['SPC'] = spcL
    Vdata['spc_id'] = 'L'

    # BACKGROUND

    # BACKGROUND contribution of the "IS" deemed NEGLIGIBLE. (Comment by Yannis: It would appear that the background contribution of the "IS" is in fact not negligible in the longest MRS wavelengths..)

    assert Vdata['BKG_id'] == 'P','Expecting units of "P", but got {} instead.'.format(Vdata['BKG_id'])

    Vdata['BKG_BBFRAME'] *= ISFactor * stt['ro_optica'] 	# RADIANCE, L
    Vdata['BKG_BBSTRAP'] *= ISFactor * stt['ro_optica'] 	# RADIANCE, L

    Vdata['BKG_VAS'] *= ISFactor * stt['ro_optica'] 	  	# RADIANCE, L

    Vdata['BKG_id'] = 'L'

    return Vdata
