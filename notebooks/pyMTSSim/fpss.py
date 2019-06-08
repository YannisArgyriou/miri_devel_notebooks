"""
NAME:

    PSS

PURPOSE:

    Models the Pupil Scanning System.

EXPLANATION:

    Provides the radiation energy input from the PSS.

CALLING SEQUENCE:

    VDATA = PSS(VDATA, T_BB, CONFIG, STT)

INPUTS:

    VDATA  - Structure with the spectral information (See initialize.pro
             for more information).
    CONFIG - Structure with the configuration of the System (MTS [+ MIRI]).

    STT    - "internal" structure with a bunch of parameters which
             characterize the instruments and the code.

OUTPUTS:

    VDATA

OPTIONAL INPUT KEYWORDS:

    None
"""
def fPSS(Vdata, config, stt):
    from numpy import cos,arctan
    from planck import planck

    wave = Vdata['wave']
    d_Imagenpss = (stt['d_Image']  - stt['d_IR_Pupil']) * 1.E-3
                                            # m

    spcM = planck(wave*1.E4, stt['T_pss']) * stt['epsilon_pss']
                                            # erg s-1 cm-2 A-1
    spcM *= 1.E1                            # W m-2 um-1

    # IN PREVIOUS EPISODES...
    #Area_LED = pi * (stt['pssapert'] * 1.E-3 / 2.)**2.
                                            # m2
    #diam_FM4 = (1.E-3 * stt['Rx_FM4'] * 2. + 1.E-3 * stt['Ry_FM4'] * 2.) / 2.
                                            # Average of diameters, m
    # NOTE: in v.1.9, where it's stt['EP_dia'] was 151... I've just inferred
    # that was the corresponding parameter !!??

    #spcE = spcM * (sin( stt['EP_dia'] * 1.E-3 / 2./ d_Imagenpss))**2. * \
    #stt.epsilon_BB * Area_LED * ( sin( diam_FM4 / 2. / d_Imagenpss))^2.
    # END OF THE FLASHBACK

    phi_pss = arctan(stt['pssapert']*1.E-3/2./d_Imagenpss)
    spcE = spcM * 2.* (1. - cos(phi_pss)) # assuming lambertian emission
                                          # of the LED here.

    # Note that Vdata['SPC'] is 'reset' when using the PSS!!

    Vdata['SPC'] = spcE * stt['PssFctr']
    Vdata['spc_id'] = 'E'                         # Irradiance


    # I've suppossed that the PSS is ONLY used with FW = BLOCK,
    # and SO the following resetting is really redundant, but
    # it doesn't hurt... and just in case I forget, better leave it.
    # By the way, I really don't know if there's any GOOD
    # reason to use PSS=ON with FW != BLOCK!!

    if config['ispoint'] == stt['True']:
        Vdata['SPC_APER'] *= 0.
        Vdata['BKG_BBFRAME'] *= 0.
        Vdata['BKG_BBstrap'] *= 0.
        Vdata['BKG_VAS'] *= 0.

    return Vdata
