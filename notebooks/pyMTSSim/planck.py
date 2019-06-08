"""
NAME:
   PLANCK()
PURPOSE:
   To calculate the Planck function in units of ergs/cm2/s/A

CALLING SEQUENCE:
   bbflux = PLANCK( wave, temp)

INPUT PARAMETERS:
   WAVE   Scalar or vector giving the wavelength(s) in **Angstroms**
           at which the Planck function is to be evaluated.
   TEMP   Scalar giving the temperature of the planck function in degree K

OUTPUT PARAMETERS:
   BBFLUX - Scalar or vector giving the blackbody flux (i.e. !pi*Intensity)
           in erg/cm^2/s/A in at the specified wavelength points.
"""
def planck(wave=None,temp=None):
    from numpy import exp
    assert wave is not None, 'Syntax - bbflux = planck( wave, temp)'

    assert temp is not None, 'Enter a blackbody temperature'

    # Gives the blackbody flux (i.e. PI*Intensity) ergs/cm2/s/a

    w = wave / 1.E8                              # Angstroms to cm
    # constants appropriate to cgs units.
    c1 =  3.7417749*1E-5                # =2*pi*h*c*c
    c2 =  1.4387687                     # =h*c/k
    val =  c2/(w*temp)

    bbflux =  c1 / ( w**5 * ( exp(val)-1. ) )
    return bbflux*1.E-8              # Convert to ergs/cm2/s/A
