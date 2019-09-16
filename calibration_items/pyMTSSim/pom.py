"""
NAME:

    POM

PURPOSE:

    Include the effects of the Pick-Off-Mirror and its structure on
    the beam:
            - Background contribution.

EXPLANATION:

    The POM structure and mirror have been hypothesized as possible
    sources of background.

CALLING SEQUENCE:

    VDATA = POM(VDATA, STT)

INPUTS:

    VDATA - Structure with the spectral information (See initialize.pro
            for more information).
    CONFIG - structure with the configuration of the System (MTS [+ MIRI]).
    STT   - "internal" structure with a bunch of parameters which
            characterize the instruments and the code.

OUTPUTS:

    VDATA - (after modifications).

OPTIONAL INPUT KEYWORDS:

    None.
"""
def POM(Vdata, config, stt):
    from numpy import arctan
    from _get_irr import _get_irr

    wave = Vdata['wave']

    # ATTENUATION

    Vdata['SPC'] *= stt['ro_optica']
    if config['ispoint'] == stt['True']: Vdata['SPC_APER'] *= stt['ro_optica']
    Vdata['BKG_BBFRAME'] *= stt['ro_optica']
    Vdata['BKG_BBSTRAP'] *= stt['ro_optica']
    Vdata['BKG_VAS'] *= stt['ro_optica']
    Vdata['BKG_TRGT'] *= stt['ro_optica']
    Vdata['BKG_FM1'] *= stt['ro_optica']
    Vdata['BKG_FM2'] *= stt['ro_optica']
    Vdata['BKG_FM3'] *= stt['ro_optica']
    Vdata['BKG_FM4'] *= stt['ro_optica']
    Vdata['BKG_OBA'] *= stt['ro_optica']

    # BACKGROUND CONTRIBUTION

    # POM Mechanical

    AN_POM = arctan( stt['Side_POM'] / (2. * stt['d_MIRI_POM']))

    # Irradiance from the POM case.
    E_POM_IN = _get_irr(wave, AN_POM, stt['T_POM'], stt['epsilon_POM'], \
    stt['theta_POM'], stt['nF_MTS'])

    # POM Mirror

    AN_POM_m = AN_POM

    # Irradiance from the POM mirror.
    E_POM_m_IN = _get_irr(wave, AN_POM_m, stt['T_POM_M'], \
    stt['epsilon_POM_M'], stt['theta_POM_M'], stt['nF_MTS'])

    Vdata['BKG_POM'] = E_POM_IN * stt['POMFctr']
    Vdata['BKG_POM_M'] = E_POM_m_IN * stt['POMFctr']

    return Vdata
