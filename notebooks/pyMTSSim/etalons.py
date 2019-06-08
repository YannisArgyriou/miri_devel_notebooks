"""
NAME:

    ETALONS

PURPOSE:

    Include the effect of the Etalons on the beam:
            attenuation.

EXPLANATION:

    ...

CALLING SEQUENCE:

    VDATA = ETALONS(VDATA, CONFIG, STT)

INPUTS:

    VDATA  - Structure with the spectral information (See initialize.pro
             for more information).
    CONFIG - structure with the configuration of the System (MTS [+ MIRI]).
           - filter : which Etalon is in the beam?
               One of {ET_1A, ET_2A, ET_1B, ET_2B}
    STT    - "internal" structure with a bunch of parameters which
             characterize the instruments and the code.

OUTPUTS:

    VDATA  - (after modifications).

OPTIONAL INPUT KEYWORDS:

    None
"""
def ETALONS(Vdata, config, stt):
    from numpy import pi,sin,cos,arcsin,sqrt,zeros,where
    from loadtrans import LoadTrans

    wave = Vdata['wave']
    nw = Vdata['nw']
    spcX = Vdata['SPC']
    Usefilter = config['Usefilter']
    beam_divergence = config['beam_divergence']

    # NOTE: this step does not involve change of units, and so spc_id is not
    # asked for confirmation.

    thetain_rad = stt['thetain_et'] * pi / 180.         # Incident angle in rad
    theta_etalon = arcsin(sin(thetain_rad) / stt['n_et'])    # Incident angle inside the etalon
                                                        # (rad)

    F_R = pi * sqrt(stt['R_et']) / (1. - stt['R_et'])      # Reflective finesse

    d = zeros(nw)

    A = zeros(nw)                                 # Absorption

    # THIS IS A TRICKY POINT, IN WHICH THE DISCRETIZATION OF CHANNELS IS
    # OVERCOME IN A RATHER CHEESY WAY... REALLY NOT SURE ABOUT THIS!!??

    for iCh in range(4):
        ixw = where((wave >= stt['waveminEt'][iCh]) & (wave <= stt['wavemaxEt'][iCh]))
        if len(ixw[0]) >= 0:
            d[ixw] = stt['d_ch'][iCh]
            A[ixw] = stt['A_Ch'][iCh]

    K_par = wave / stt['nm_par']
    K_flat = wave / stt['nm_flat']
    K_rms = wave / stt['nm_rms']

    F_def = ( (1/(K_par/3.**0.5))**2 + (1/(K_flat/2.))**2 + \
    (1/(K_rms/4.7))**2 )**(-0.5)

    dum = d * 1.E3                                      # um
    wave0 = 2. * stt['n_et'] * dum * cos(theta_etalon)

    int_ord = wave0 / wave                              # ??

    delta = 4. * pi / wave * stt['n_et'] * dum * cos(theta_etalon)

    Abss = ( 1. - A / ( 1. - stt['R_et'] ) )**2.

    Fa = 2. * pi * stt['n_et'] / (int_ord * pi * beam_divergence**2.)  # ??

    F_effect = 1. / sqrt(F_R**(-2.) + F_def**(-2.) + Fa**(-2.))
                                                                         # OK

    R_effect = 1. + 0.5 * (pi / F_effect)**2. *\
    ( 1. - sqrt( 1. + (2 * F_effect / pi )**2. ) )

    It_defect = Abss * (( 1. - stt['R_et'] ) / ( 1. + stt['R_et'] ))*\
    (( 1. + R_effect )/( 1. - R_effect ))/\
    ( 1. + ( 2. * F_effect / pi )**2. * ( sin( delta/2. ) )**2.)

    spcX *= It_defect

    # TRANSMISSION CURVE OF THE ETALON

    Etalons_tags = stt['Etalons_filtern'].keys()

    filtern = stt['Etalons_filtern'][Usefilter]

    # Note: some of the transmission curves of the Etalons are currently
    # unknown: ET_1A and ET_2A.

    filtern = stt['Tpath']+ filtern
    iTrans_et = LoadTrans(filtern,wave)
    spcX *= iTrans_et

    Vdata['SPC'] = spcX

    # EFFECT OVER BACKGROUND EMISSION

    Vdata['BKG_BBFRAME'] *= It_defect * iTrans_et
    Vdata['BKG_BBSTRAP'] *= It_defect * iTrans_et

    Vdata['BKG_VAS'] *= It_defect * iTrans_et

    return Vdata
