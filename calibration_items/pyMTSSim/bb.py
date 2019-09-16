"""
NAME:

    BB

PURPOSE:

    Provides the Radiation Energy Input from the "Black Body".
    Also includes the Background contribution from the BB-frame and the
    BB-strap.

EXPLANATION:

   This function provides the output emittance of the Black Body in the
   MTS, and the contributions to the background from the BB-Frame and
   the BB-strap.

CALLING SEQUENCE:

   VDATA = BB(VDATA, CONFIG, STT)

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
def BB (Vdata, config, stt):
    from numpy import pi
    from planck import planck
    # BB EMISSION

    T_BB = config['T_BB']

    wave = Vdata['wave'] # um

    spcM = planck(wave * 1.E4, T_BB) * stt['epsilon_BB']
                                                # Emittance, erg s-1 cm-2 A-1
    spcM = spcM * 1.E-7 * 1.E4 * 1.E4           # =*1.E1, Emittance, W m-2 um-1


    Area_BB = pi * (stt['dia_pinhole'] * 1.E-3 / 2.)**2.
                                                # Area (m2) of pin-hole in
                                                # front of BB.

    spcP = spcM * Area_BB

    Vdata['SPC'] = spcP
    Vdata['spc_id'] = 'P' # Spectral Power, W um-1

    # BB BKG
    #
    #            BB - frame

    spcM_bbframe = planck(wave*1.E4, stt['T_bbframe']) * stt['epsilon_bbframe']
                                                # erg s-1 cm-2 A-1
    spcM_bbframe *= 1.E1                        # W m-2 um-1


    spcM_bbframe *= stt['BBFrameFctr']          # Scaling to allow for unaccurate
                                                # knowledge of the design /
                                                # performance.

    #            BB - strap

    spcM_bbstrap = planck(wave*1.E4, stt['T_bbstrap']) * stt['epsilon_bbstrap']
    spcM_bbstrap *= 1.E1                        # W m-2 um-1
    spcM_bbstrap *= stt['BBStrapFctr']             # scaling

    Vdata['BKG_BBFRAME'] = spcM_bbframe
    Vdata['BKG_BBSTRAP'] = spcM_bbstrap
    Vdata['BKG_id'] = 'M'

    # THE BACKGROUNDS SHOULD BE MULTIPLIED BY
    # SOME AREA... !!!! VALUE??

    Vdata['BKG_BBFRAME'] *= stt['Area_bbframe']       # 1.E-4 m-2 (1 cm2) CHEAP TRICK!!
    Vdata['BKG_BBSTRAP'] *= stt['Area_bbstrap']       # 1.E-4 m-2 (1 cm2) CHEAP TRICK!!
    Vdata['BKG_id'] = 'P'                          # "CONSISTENT" WITH spc_id !!

    Vdata['units'] = 'MKS'

    return Vdata
