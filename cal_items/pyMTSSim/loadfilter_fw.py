"""
NAME:

LOADFILTER_FW

PURPOSE:

To read a filter transmission curve from an ascii file (organized
in columns) and interpolate it to a desired set of wavelengths.

EXPLANATION:

Only intended to be used within the scope of MTSSim, not for general
use.

CALLING SEQUENCE:

iTrans = LOADFILTER_FW(FILTERN, WAVE)

INPUTS:

FILTERN - Name of the file with the transmission curve of the filter.
        It is expected that the information comes in columns, in number
        of 4. First would be the wavelength, in microns, the two
        next are ignored and the fourth is the curve itself. See any of
        of the curves in $MTSSim/MTS_FW/ for examples.

WAVE    - Vector with the wavelengths at which we want to interpolate
        the transmission curve. It must be in microns. If the spectral
        range is COMPLETELY out of bounds of that given by the vector
        of wavelengths in "FILTERN", the program will crash. If
        SOME requested wavelengths are outside the spectral range of
        FILTERN, then a transmission of "0" will be returned for them.
        The values in WAVE must be strictly ascending or descending.

OUTPUTS:

iTrans  - Vector with the transmission curve evaluated at the points
        listed in WAVE.

OPTIONAL INPUTS KEYWORDS:

None.
"""
def LoadFilter_FW(filtern, wave):
    from numpy import genfromtxt,where,zeros
    from scipy.interpolate import interp1d

    Filter_um,T1,T2,T3 = genfromtxt(filtern, skip_header = 14, skip_footer=1, usecols=(0,1,2,3), delimiter = '',unpack='True')
    Filter_T = T3         # it's the transmission at 35 K

    # INTERPOLATION OF TRANSMISSION CURVE

    selw = where((wave >= min(Filter_um)) & (wave <= max(Filter_um)))
    sub_iFilter_T = interp1d(Filter_um,Filter_T)(wave[selw])

    iFilter_T = zeros(len(wave))
    iFilter_T[selw] = sub_iFilter_T

    iFilter_T /= 100.

    return iFilter_T
