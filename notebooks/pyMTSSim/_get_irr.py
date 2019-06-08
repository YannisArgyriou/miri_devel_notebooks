"""
NAME:

    _GET_IRR

PURPOSE:

    It retrieves the Irradiance emitted from a generic element in the optical
    path.

EXPLANATION:

    It is an auxiliary task used by fms.pro, oba.pro and pom.pro.

CALLING SEQUENCE:

    E = _GET_IRR(WAVE, AN, T, EPSILON, THETA, NF_MTS)

INPUTS:

    WAVE    - wavelengths (um).
    AN      - Angle subtended by the element from the Input Plane of MIRI.
    T       - Temperature of the source.
    EPSILON - Emissivity of the source.
    THETA   - Angle that the surface of the source forms with the line
              of sight.
    NF_MTS  - Focal Ratio of the MTS in the image space.

OUTPUTS:

    E       - Irradiance, given in W m-2 um-1.

OPTIONAL INPUT KEYWORDS:

    None.
"""
def _get_irr(wave, AN, T, epsilon, theta, nF_MTS):
    from numpy import pi,cos,arctan
    from planck import planck

    M = planck(wave*1.E4, T) * epsilon      # erg s-1 cm-2 A-1
    M *= 1.E1                                # W m-2 um-1

    THETA_MTS = arctan( 1./(2.*nF_MTS) )

    E = M * cos(theta * pi / 180.) * 2. * (1. - cos(AN))
    E_IN = M * cos(theta * pi / 180.) * 2. * \
        (1. - cos(arctan(THETA_MTS)))

    if AN <= THETA_MTS: E_IN = E

    #if AN <= (1./(2.*nF_MTS)): print 'AN LT LO OTRO'
    #else: PRINT 'AN >= LO OTRO'

    return E_IN
