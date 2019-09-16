"""
NAME:

LOADTRANS

PURPOSE:

To read a transmission curve of some optical element from an ascii
file (organized in columns) and interpolate it to a desired set of
wavelengths.

EXPLANATION:

Only intended to be used within the scope of MTSSim, not for general
use.

CALLING SEQUENCE:

iTrans = LOADTRANS(CURVEN, WAVE)

INPUTS:

CURVE   - Name of the file with the transmission curve.
        It is expected that the information comes in columns, in number
        of 2. First would be the wavelength, in microns, the next is
        being the curve itself. See any of of the curves in
        $MTSSim/TRANSMIT/ for examples.

WAVE    - Vector with the wavelengths at which we want to interpolate
        the transmission curve. It must be in microns. If the spectral
        range is COMPLETELY out of bounds of that given by the vector
        of wavelengths in "FILTERN", the program will crash. If
        SOME requested wavelengths are outside the spectral range of
        CURVEN, then a transmission of "0" will be returned for them.
        The values in WAVE must be strictly ascending or descending.

OUTPUTS:

iTrans  - Vector with the transmission curve evaluated at the points
        listed in WAVE.

OPTIONAL INPUTS KEYWORDS:

None.
"""
def LoadTrans(curven, wave):
    from numpy import genfromtxt,where,zeros
    from scipy.interpolate import interp1d

    um, T = genfromtxt(curven, usecols=(0,1), delimiter = '',unpack=True)
    # INTERPOLATION OF TRANSMISSION CURVE

    assert (wave[0]>= min(um)) & (wave[-1] <= max(um)),'Wavelength extremes outside of interpolation range.'
    sub_iT = interp1d(um,T)(wave)

    iT = zeros(len(wave))
    iT = sub_iT

    iT /= 100.

    return iT
