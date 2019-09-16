"""
NAME:
    ADDUP_BKG

PURPOSE:

    To add the spectra of the different background sources into a TOTAL
    Background.

EXPLANATION:

    Not needed.

CALLING SEQUENCE:

    VDATA  = ADDUP_BKG(VDATA, CONFIG, STT)

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

def ADDUP_BKG(Vdata, config, stt):
    from numpy import all,nan,pi,zeros

    nw = Vdata['nw']

    Vdata['BKG'] = zeros(nw)

    BKGtags = list(Vdata['BKGtags'])
    excl_list = ['BKG','BKG_PNH','BKG_EXT']
    for item in excl_list:
        BKGtags.remove(item)

    for BKGtag in BKGtags:
        subBKG = Vdata[BKGtag]
        Vdata['BKG'] += subBKG


    # EXTENDED VS. POINT-LIKE BACKGROUND

    if config['ispoint'] == stt['True']:

        Vdata['BKG_EXT'] = Vdata['BKG_TRGT'] + Vdata['BKG_POM'] + Vdata['BKG_POM_M'] + \
	    Vdata['BKG_FM1'] + Vdata['BKG_FM2'] + Vdata['BKG_FM3'] + Vdata['BKG_FM4'] + \
	    Vdata['BKG_OBA']

	r_Airy = 1.22 * (Vdata['wave']*1.E-6) * stt['nF_MTS']  # m
        Ap_Airy = pi * r_Airy**2.		   # m^2

	Vdata['BKG_PNH'] = (Vdata['BKG_BBFRAME'] + Vdata['BKG_BBSTRAP'] + \
	    Vdata['BKG_VAS']) * Ap_Airy

    else:
        Vdata['BKG_EXT'] = Vdata['BKG']
        Vdata['BKG_PNH'] = Vdata['BKG_PNH'] # left unchanged, as NaN!

    return Vdata
