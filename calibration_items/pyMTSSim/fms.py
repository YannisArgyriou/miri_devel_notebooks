"""
NAME:

   FMs

PURPOSE:

   To include the effects of the Folding Mirrors on the beam:
            - Background contribution and attenuation.

EXPLANATION:

   The 4 Folding Mirrors contribute to the Background,
   and also absorb some of the light they reflect.

CALLING SEQUENCE:

   VDATA = FMs(VDATA, STT)

INPUTS:

   VDATA - Structure with the spectral information (See initialize.pro
           for more information).
   CONFIG - structure with the configuration of the System (MTS [+ MIRI]).
   STT   - "internal" structure with a bunch of parameters which
           characterize the instruments and the code.

OUTPUTS:

   VDATA - (after modifications).

OPTIONAL INPUT KEYWORDS:
"""
def FMs(Vdata, config, stt):
    from numpy import arctan
    from loadtrans import LoadTrans
    from _get_irr import _get_irr

    wave = Vdata['wave']

    # APPLICATION OF THE TRANSMISSION CURVE OF THE FMs
    # NOTE: this step does not involve change of units, and so spc_id
    # is not asked for confirmation.

    iTrans_FMs = LoadTrans( stt['Tpath'] + stt['FMtransn'], wave)

    Vdata['SPC'] *= iTrans_FMs**4.                # 4 reflections
    if config['ispoint'] == stt['True']: Vdata['SPC_APER'] *= iTrans_FMs**4.

    # BACKGROUND

    # Note: a dubious "2" factor in the denominator has been eliminated
    # after conversation with T. Belenguer on 16th/09/2010 (Ry_FM1 is
    # a radius, and that factor seems to come from the wrong
    # interpretation of it being a diameter). The same has been done
    # with AN_FM2, AN_FM3 and AN_FM4

    AN_FM1 = arctan( stt['Ry_FM1'] / (stt['d_Image'] - stt['d_prim_FM1']))

    # Irradiance of FM1
    E_FM1_IN = _get_irr(wave, AN_FM1, stt['T_FM1'], stt['epsilon_FMs'], \
               stt['theta_FMs'], stt['nF_MTS'])


    AN_FM2 = arctan( stt['Ry_FM2'] / (stt['d_Image'] - stt['d_prim_FM1']) )

    # Irradiance of FM2
    E_FM2_IN = _get_irr(wave, AN_FM2, stt['T_FM2'], stt['epsilon_FMs'], \
               stt['theta_FMs'], stt['nF_MTS'])

    AN_FM3 = arctan( stt['Ry_FM3'] / (stt['d_Image'] - stt['d_prim_FM1'] - \
             stt['d_FM1_FM2']) )

    # Irradiance of FM3
    E_FM3_IN = _get_irr(wave, AN_FM3, stt['T_FM3'], stt['epsilon_FMs'], \
               stt['theta_FMs'], stt['nF_MTS'])

    AN_FM4 = arctan( stt['Ry_FM4'] / (stt['d_Image'] - stt['d_prim_FM1'] - \
             stt['d_FM1_FM2'] - stt['d_FM2_FM3'] - stt['d_FM3_FM4']) )

    # Irradiance of FM4

    E_FM4_IN = _get_irr(wave, AN_FM4, stt['T_FM4'], stt['epsilon_FMs'], \
               stt['theta_FMs'], stt['nF_MTS'])

    Vdata['BKG_FM1'] = E_FM1_IN * stt['FMFctr'] * iTrans_FMs**3.
    Vdata['BKG_FM2'] = E_FM2_IN * stt['FMFctr'] * iTrans_FMs**2.
    Vdata['BKG_FM3'] = E_FM3_IN * stt['FMFctr'] * iTrans_FMs
    Vdata['BKG_FM4'] = E_FM4_IN * stt['FMFctr']

    # BKG ATTENUATION

    # BB-BKG
    Vdata['BKG_BBFRAME'] *= iTrans_FMs**4.
    Vdata['BKG_BBSTRAP'] *= iTrans_FMs**4.

    # VAS-BKG

    Vdata['BKG_VAS'] *= iTrans_FMs**4.

    # TARGET-BKG

    Vdata['BKG_TRGT'] *= iTrans_FMs**4.

    return Vdata
