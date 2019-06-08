"""
NAME:

    FW

PURPOSE:

    To emulate the MTS Filter Wheel and its effect on the beam:
            attenuation.
    The Etalons are treated apart, in "etalons.pro".

EXPLANATION:

    Unnecessary.

CALLING SEQUENCE:

    VDATA = FW(VDATA, CONFIG, STT)

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
def FW(Vdata, config, stt):
    from loadfilter_fw import LoadFilter_FW

    spcX = Vdata['SPC']
    Usefilter = config['Usefilter']

    # NOTE: this step does not involve change of units, and so spc_id
    # is not asked for confirmation.

    if Usefilter == 'HOLE': return Vdata
                                    # No action required!

    if Usefilter == 'BLOCK':
        Vdata['SPC'] *= 0. # block!
        Vdata['BKG_BBFRAME'] *= 0. # block!
        Vdata['BKG_BBSTRAP'] *= 0. # block!
        Vdata['BKG_VAS'] *= 0. # block!
        return Vdata

    if Usefilter == 'LWP': filtern = stt['Mpath']+ stt['LWP_fltr']
    if Usefilter == 'SWP': filtern = stt['Mpath']+ stt['SWP_fltr']
    assert (Usefilter == 'LWP') or (Usefilter == 'SWP') or (Usefilter[0:3] == 'ET_'),'**WARNING** Unknown filter **WARNING**'

    wave = Vdata['wave']

    iFilter_T = LoadFilter_FW(filtern, wave)

    # MULTIPLICATION OF spc BY INTERPOLATED TRANSMISSION CURVE

    spcX *= iFilter_T

    Vdata['SPC'] = spcX

    # BACKGROUND

    # Attenuation of BKG sources

    Vdata['BKG_BBFRAME'] *= iFilter_T
    Vdata['BKG_BBSTRAP'] *= iFilter_T
    Vdata['BKG_VAS'] *= iFilter_T

    # FW BACKGROUND: ASSUMED TO BE NEGLIGIBLE

    return Vdata
