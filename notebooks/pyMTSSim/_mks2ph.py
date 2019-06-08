"""
NAME:

    _MKS2PH

PURPOSE:

    Converts an [Spectral] Irradiance from "MKS" to "PH". This is,
    from W m-2 [um-1] to photons s-1 cm-2 [um-1].

EXPLANATION:

    Used by _convert.pro and mirisim.pro

CALLING SEQUENCE:

    FLXPH = _MKS2PH(FLXMKS, WAVE, STT, [/REVERSE])

INPUTS:

    FLXMKS  - Input spectrum in "MKS" units.
    WAVE    - Wavelength[s] that correspond to FLXMKS, given in microns (um).
    STT     - "internal" structure with a bunch of parameters which
              characterize the instruments and the code.

OUTPUTS:

    FLXPH   - Output spectrum in "PH" units.

OPTIONAL INPUT KEYWORDS:

    REVERSE - If this keyword is set, then it converts from "PH" to "MKS".
"""
def _MKS2ph(flxMKS, wave, stt, REVERSE=False):
    # [wave] = um
    # input: W.m-2.um-1
    # output: ph.s-1.cm-2.um-1

    fctr = ((stt['hplanck'] * stt['clight'] / (wave * 1.E-6))**(-1.)) *\
                     1.E-4
                        # -> ph s-1 cm-2 um-1
                        # direct

    if REVERSE is True:
        # perform reverse transform
        fctr = fctr**(-1.)

    return flxMKS * fctr
