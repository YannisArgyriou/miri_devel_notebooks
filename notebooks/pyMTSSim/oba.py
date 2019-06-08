"""
NAME:

    OBA

PURPOSE:

    Include the effects of the Optical Bench Assembly on the beam:
            - Background contribution.

EXPLANATION:

    The OBA has been hypothesized as an additional source of background.

CALLING SEQUENCE:

    VDAT = OBA(VDATA, STT)

INPUTS:

    VDATA - Structure with the spectral information (See initialize.pro
            for more information).
    STT   - "internal" structure with a bunch of parameters which
            characterize the instruments and the code.

OUTPUTS:

    VDATA - (after modifications).

OPTIONAL INPUT KEYWORDS:

    None.
"""
def OBA(Vdata, stt):
    from numpy import cos,arctan
    from planck import planck

    wave = Vdata['wave']

    # BACKGROUND CONTRIBUTION

    # OBA Mechanical

    #AN_OBA = atan( stt['Side_OBA'] / (2. * stt['d_MIRI_OBA']) )

    # Irradiance from the OBA
    #E_OBA_IN = _get_Irr(wave, AN_OBA, stt['T_OBA'], stt['epsilon_OBA'], \
    #stt['theta_OBA'], stt['nF_MTS'])

    # ALTERNATIVE SHOTS!

    #E_OBA_IN = _get_Irr(wave, AN_OBA, stt['T_OBA'], stt['epsilon_OBA'], \
    #0., stt['nF_MTS'])


    M = planck( wave*1.E4, stt['T_OBA']) * stt['epsilon_OBA']
    M *= 1.E1 # units conversion
    THETA_MTS = arctan( 1./(2.*stt['nF_MTS']) )
    E_OBA_IN = M * 2. * (1. - cos(arctan(THETA_MTS)))

    Vdata['BKG_OBA'] = E_OBA_IN * stt['OBAFctr']

    return Vdata
