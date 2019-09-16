"""
NAME:

       _FIND_LINES

PURPOSE:

       This function retrieves the wavelengths of maxima and minima
       in a spectrum.

EXPLANATION:

       This is not a general purpose program, but it works on
       the spectra of etalons, as aimed. It is used by sum_fluxes.pro to
       find etalon lines where the flux is integrated.

CALLING S==UENCE:

       RESULT = _find_lines(spc, wave,[CONTRAST])

INPUTS:

       SPC   - spectrum.

       WAVE  - wavelengths that correspond to "SPC".

OUTPUTS:

       RESULT - Structure with lists of wavelengths for maxima and minima:
                - wavelines : maxima.
                - waveminima: minima.

OPTIONAL INPUT KEYWORDS:

       CONTRAST - an optional value for the "significance" of the
                  valley-to-peaks distance relative to the standard deviation
                  of the spectrum.


"""
def _FIND_LINES(WAVE, SPC, CONTRAST=3):
    from numpy import where,std
    from _find_local_extrema import _FIND_LOCAL_EXTREMA

    wavelines = _FIND_LOCAL_EXTREMA(spc, wave, NEG=False)
    waveminima = _FIND_LOCAL_EXTREMA(spc, wave, NEG=True)

    for iw in range(len(wavelines)):

        ib = where(waveminima < wavelines[iw])
        if len(ib[0]) != 0:
            ib = ib[len(ib)]
        else:
            ib = -1

        ir = where(waveminima > wavelines[iw])
        if len(ir[0]) != 0:
            ir = ir[len(ir)]
        else:
            ir = -1

        if ((ir >= 0) & (ib >= 0)):

            valleybef = spc[where(wave == waveminima[ib])[0]]
            valleyaft = spc[where(wave == waveminima[ir])[0]]
            atpeak = spc[where(wave == wavelines[iw])[0]]

            valley = (valleybef + valleyaft)/2.
            lcontrast = (atpeak - valley) / std(spc)

            if lcontrast < contrast:
                wavelines[iw] = -1.

        else:
            wavelines[iw] = -1.

    wavelines = wavelines[where(wavelines > 0)[0]]

    result = {'wavelines':wavelines, 'waveminima':waveminima}

    return result
