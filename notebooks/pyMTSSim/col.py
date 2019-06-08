"""
NAME:

   COL

PURPOSE:

   To emulate the MTS-collimator.

EXPLANATION:

    ...

CALLING SEQUENCE:

    VDATA = COL(VDATA, CONFIG, STT)

INPUTS:

    VDATA  - Structure with the spectral information (See initialize.pro
             for more information).
    CONFIG - structure with the configuration of the System (MTS [+ MIRI]).

    STT    - "internal" structure with a bunch of parameters which
             characterize the instruments and the code.

OUTPUTS:

    VDATA  - (after modifications).

OPTIONAL INPUT KEYWORDS:

    None
"""
def COL(Vdata, config, stt):
    from numpy import cos,arctan

    divergence = config['divergence']        # Does the beam diverge?

    assert Vdata['spc_id'] == 'P','Expecting units of "P", but got {} instead.'.format(Vdata['spc_id'])

    spcP = Vdata['SPC']

    if divergence: beam_divergence = 0.
    else: beam_divergence = stt['dia_pinhole'] / stt['focal_collimator']
                                           # Collimation degree

    config['beam_divergence'] = arctan(beam_divergence)

    Teta2 = arctan(0.5 * stt['ap_collimator'] / stt['focal_collimator'])
                                        # Solid angle
                                        # subtended by the
                                        # collimator mirror

    # Radiation power arriving at the collimator aperture:

    spcP *= ( (2. * (1. - cos(Teta2)) * stt['ro_optica']) ) # W um-1

    spcP *= stt['ro_optica'] # Collimator's folding mirror


    Vdata['SPC'] = spcP
    Vdata['spc_id'] = 'P'

    # BACKGROUND

    assert Vdata['BKG_id'] == 'P','Expecting units of "P", but got {} instead.'.format(Vdata['BKG_id'])

    # BB
    # Kept for record:
    # Vdata['BKG_BBframe'] *= 4 * pi * stt['ap_collimator']**2. #-> v 1.9
    # Vdata['BKG_BBstrap'] *= 4 * pi * stt['ap_collimator']**2. #-> NOT DONE IN THE
    #                                                      #   v 1.9!
    #                                                      #   but should have been
    #                                                      #   done, for consistence.

    Vdata['BKG_BBFRAME'] *= (2. * (1. - cos(Teta2)) ) # W um-1
    Vdata['BKG_BBSTRAP'] *= (2. * (1. - cos(Teta2)) ) # W um-1
    Vdata['BKG_id'] = 'P'


    # BACKGROUND CONTRIBUTION OF THE COLLIMATOR: ASSUMED TO BE NEGLIGIBLE

    # ATTENUATION BY THE REFLECTION ON THE COLLIMATOR MIRROR and its FOLDING
    # MIRROR:


    Vdata['BKG_BBFRAME'] *= stt['ro_optica'] * stt['ro_optica']
    Vdata['BKG_BBSTRAP'] *= stt['ro_optica'] * stt['ro_optica']

    return Vdata
