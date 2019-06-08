"""
NAME:

    MOS

PURPOSE:

    Include the effects of the MOS on the beam:
            - attenuation.

EXPLANATION:

    This function multiplies the input spectrum by the attenuation
    corresponding to the MOS. It also applies this attenuation to
    the BACKGROUND sources previous to this element.

CALLING SEQUENCE:

    VDATA = MOS(VDATA, STT)

INPUTS:

    VDATA - Structure with the spectral information (See initialize.pro
            for more information).
    CONFIG - structure with the configuration of the System (MTS [+ MIRI]).
    STT   - "internal" structure with a bunch of parameters which
            characterize the instruments and the code.

OUTPUTS:

    VDATA - (after modifications).

OPTIONAL INPUT KEYWORDS:

    None
"""
def MOS(Vdata, config, stt):
    from loadtrans import LoadTrans

    spcX = Vdata['SPC']
    wave = Vdata['wave']

    # NOTE: this step does not involve change of units, and so spc_id is not
    # asked for confirmation.

    iTrans_MOS = LoadTrans(stt['Tpath']+stt['MOSTransn'], wave)

    spcX *= iTrans_MOS**4.        # 4 reflections
    spcX *= (1. - stt['MOSObsc'])    # Obscuration of the secondary in the
                                  # pupil mask!!

    Vdata['SPC'] = spcX

    if config['ispoint'] == stt['True']: \
        Vdata['SPC_APER'] *= iTrans_MOS**4. * (1.-stt['MOSObsc'])

    # BACKGROUND

    Vdata['BKG_BBFRAME'] *= iTrans_MOS**4. * (1.-stt['MOSObsc'])
    Vdata['BKG_BBSTRAP'] *= iTrans_MOS**4. * (1.-stt['MOSObsc'])

    Vdata['BKG_VAS'] *= iTrans_MOS**4. * (1.-stt['MOSObsc'])

    Vdata['BKG_TRGT'] *= iTrans_MOS**4. * (1.-stt['MOSObsc'])


    return Vdata
