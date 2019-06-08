"""
NAME:

    _MKS2CGS

PURPOSE:

    Converts an [Spectral] Irradiance from "MKS" to "CGS". This is,
    from W m-2 [um-1] to erg s-1 cm-2 [um-1].

EXPLANATION:

    Used by _convert.pro and mirisim.pro

CALLING SEQUENCE:

    FLXCGS = _MKS2CGS(FLXMKS,[/REVERSE])

INPUTS:

    FLXMKS - Input spectrum in "MKS" units.


OUTPUTS:

    FLXCGS - Output spectrum in "CGS" units.

OPTIONAL INPUT KEYWORDS:

    REVERSE - If this keyword is set, then it converts from "CGS" to "MKS".
"""
def _MKS2cgs(flxMKS, REVERSE=False):

    # input: W.m-2.um-1
    # output: erg.s-1.cm-2.um-1
    fctr = 1.E7 * 1.E-4  # direct

    if REVERSE:
        # perform reverse transform
        fctr = fctr**(-1.)

    return flxMKS * fctr
