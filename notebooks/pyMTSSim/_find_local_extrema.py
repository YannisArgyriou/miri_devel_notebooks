"""
NAME:

       _FIND_LOCAL_EXTREMA

PURPOSE:

       This function retrieves the wavelengths of emission/absorption lines
       in a spectrum.

EXPLANATION:

       This is not a general purpose program, but it works on
       the spectra of etalons, as aimed. It is used by _find_lines.pro to
       find etalon lines where the flux is integrated (by sum_fluxes.pro).

CALLING SEQUENCE:

       LINES = _find_local_extrema(spc, wave)

INPUTS:

       SPC   - spectrum.

       WAVE  - wavelengths that correspond to "SPC".

OUTPUTS:

       LINES - List with the wavelengths of local maxima in SPC.

OPTIONAL INPUT KEYWORDS:

       NEG   - If this keyword is set the function will retrieve a list with
               the wavelengths of local minima, instead of maxima.
"""

def _FIND_LOCAL_EXTREMA(spc, wave, NEG=False):
    from numpy import where

    nw = len(wave)
    assert nw == len(spc),'_FIND_LINES: "spc" and "wave" do not have the same dimensions.'

    sub1plus = spc[1:nw-1]
    sub1min = spc[:nw-2]

    sub2 = spc[1:nw-2]

    difmin = sub2 - sub1min[:nw-3]
    difplus = sub2 - sub1plus[1:nw-2]

    if NEG is False:
        ixlines = where((difmin * difplus > 0.) & (difmin > 0.))[0]
    else:
        ixlines = where((difmin * difplus > 0.) & (difmin < 0.))[0]

    wavelines = wave[1:nw-2]
    wavelines = wavelines[ixlines]


    return wavelines
