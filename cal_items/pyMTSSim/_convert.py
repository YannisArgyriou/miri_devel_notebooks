def _switchmksconv(units, spc, wave, stt):
    from _mks2cgs import _MKS2cgs
    from _mks2ph import _MKS2ph
    from _mks2jy import _MKS2jy

    if units == 'MKS': return spc
    if units == 'CGS': return _MKS2cgs(spc)
    if units == 'PH': return _MKS2ph(spc, wave, stt)
    if units == 'JY': return _MKS2jy(spc, wave, stt)

def _CONVERT(Vdata, config, stt):
    """
    NAME:

        _CONVERT

    PURPOSE:

        To convert the Irradiances from MKS (W m-2 um-1) to other units.

    EXPLANATION:

        The spectra of the target and the background sources may be changed
        to other units.

    CALLING S==UENCE:

        VDATA = _CONVERT(VDATA, CONFIG, STT)

    INPUTS:

        VDATA  - Structure with the spectral information (See initialize.pro
                 for more information).
        CONFIG - structure with the configuration of the System (MTS [+ MIRI]).
                 The relevant parameter is "units":

                       'MKS' : W m-2 um-1
                       'CGS' : erg s-1 cm-2 um-1
                       'PH'  : ph s-1 cm-2 um-1
                       'JY'  : 1.E-26 W m-2 Hz-1

        STT    - "internal" structure with a bunch of parameters which
                 characterize the instruments and the code.

    OUTPUTS:

        VDATA  - (after modifications).

    OPTIONAL INPUTS KEYWORDS:

        None.
    """

    wave = Vdata['wave']

    # Default units of Irradiance: MKS, W.m-2.um-1
                                  # Actually, to be "MKS" it
                                  # should be W.m-3...
                                  # but don't be so "fastidious".

    units2conv = config['units']

    assert Vdata['spc_id'] == 'E', 'Expecting units of "E", but got {} instead.'.format(Vdata['spc_id'])

    assert Vdata['BKG_id'] == 'E', 'Expecting units of "E", but got instead.'.format(Vdata['BKG_id'])

    assert Vdata['units'] == 'MKS', 'Spectral Irradiance must be in units of W.m-2.um-1 ("MKS")'

    # Conversion of PHYSICAL UNITS of IRRADIANCE

    # Conversion of units of the TARGET spectrum & its UNCERTAINTY

    spc = Vdata['SPC']
    spc = _switchmksconv(units2conv, spc, wave, stt)

    Vdata['SPC'] = spc

    if config['ispoint'] == stt['True']:

        spc_aper = Vdata['spc_aper']

    	if (units2conv == 'CGS') or (units2conv == 'PH'): corrF = 1.e4
    	else: corrF = 1.

    	spc_aper = _switchmksconv(units2conv, spc_aper, wave, stt) * corrF
    	Vdata['spc_aper'] = spc_aper

    # Conversion of units of the BACKGROUND spectra

    BKGtags = list(Vdata['BKGtags'])
    excl_list = ['BKG_PNH']
    for item in excl_list:
        BKGtags.remove(item)

    for BKGtag in BKGtags:
        subBKG = Vdata[BKGtag]
        subBKG = _switchmksconv(units2conv, subBKG, wave, stt)
        Vdata[BKGtag] = subBKG


    if config['ispoint'] == stt['True']:
    	# We have to apply an ad-hoc correction to have the units right,
    	# as BKG_PNH is not an irradiance, as assumed by _switchmksconv,
    	# but a "power".

    	if (units2conv == 'CGS') or (units2conv == 'PH'): corrF = 1.e4
    	else: corrF = 1.

    	Vdata['BKG_PNH'] = _switchmksconv(units2conv, \
    	    Vdata['BKG_PNH'], wave, stt) * corrF

    Vdata['units'] = units2conv

    return Vdata
