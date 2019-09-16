"""
NAME:

    _MKS2JY

PURPOSE:

    Converts an [Spectral] Irradiance from "MKS" to "JY". This is,
    from W m-2 [um-1] to Janskys, i.e. units of 1.E-26 W m-2 Hz-1.

EXPLANATION:

    Used by _convert.pro and mirisim.pro

CALLING SEQUENCE:

    FLXJY = _MKS2JY(FLXMKS, WAVE, STT, [/REVERSE])

INPUTS:

    FLXMKS  - Input spectrum in "MKS" units.
    WAVE    - Wavelength[s] that correspond to FLXMKS, given in microns (um).
    STT     - "internal" structure with a bunch of parameters which
              characterize the instruments and the code.

OUTPUTS:

    FLXJY   - Output spectrum in "JY" units.

OPTIONAL INPUT KEYWORDS:

    REVERSE - If this keyword is set, then it converts from "JY" to "MKS".
"""
def _MKS2jy(flxMKS, wave, stt, REVERSE=False):
    # [wave] = um
    # input: W.m-2.um-1
    # output: Jy

    fctr = 1.E26 * wave**2. / (stt['clight']*1.E6)
                            # -> Jy (= 1.E-26 W.m-2.Hz-1)
                            # direct

    if REVERSE is True:
        # perform reverse transform
        fctr = fctr**(-1.)

    return flxMKS * fctr
