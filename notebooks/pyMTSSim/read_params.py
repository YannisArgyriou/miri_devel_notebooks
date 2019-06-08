"""
NAME:

    READ_PARAMS

PURPOSE:

    To load the program with values for the "secondary" parameters.

EXPLANATION:

    See the acompanying manual for more information.

CALLING SEQUENCE:

    UNEXP = READ_PARAMS(PARFILE, STT)

INPUTS:

    PARFILE - File with secondary parameters, organized in two columns
              the first with the name of the parameter (lower case), and
              the second withe value (no quotation marks). Anything after
              a hash symbol (#) is ignored. See "inputs_test.txt" for
              an example.

            List of parameters:

                - outpath. Path to the "outputs" folder.
                - saveeps. Save graphs to .eps files? (yes/no).
                - epsroot. If "saveeps=yes", then the value of this
                  parameter is taken as the first part of the names
                  of the different plots. The second part is given by
                  the program and specifies the different plots.
                - savespc. Save the spectra to a plain text file?
                  (yes/no).
                - spcfile. Name of the file where the spectra will be
                  saved if "savespc" is "yes".
                - savesum. Save a hard-copy of the integrated flux
                  results? (yes/no).
                - spcsumfile. Name of the plain text where the
                  integrated fluxes will be stored if "savesum=yes".
                - Units. Units of the [spectral] irradiances, coded as
                  follows:
                       - MKS. W m-2 [um-1].
                       - CGS. erg s-1 cm-2 [um-1].
                       - PH. ph s-1 cm-2 [um-1].
                       - JY . Units of 10-26 W m-2 Hz-1 (only for
                         spectral irradiance).
                       - wmin. Minimum wavelength considered (>= 5 um).
                       - wmax. Maximum wavelength considered (<= 30 um).
                       - Divergence. Does the beam coming out of the BB
                         diverge?
                       - VASemit. Does the VAS emit? (yes/no).                        -
    STT         - "internal" structure with a bunch of parameters which
                  characterize the instruments and the code.

OUTPUTS:

    UNEXP       - Structure with the secondary parameters. Those parameters
                  not listed in PARFILE will be given default values (see the
                  code, first lines). If PARFILE is not found, again, the
                  default values will be used. So, be careful.

OPTIONAL INPUTS KEYWORDS:

    None
"""
def READ_PARAMS(parfile, stt):
    from numpy import genfromtxt
    unexp = {'OUTPATH':'OUTPUTS', 'SAVEEPS':stt['False'], 'EPSROOT':'None',\
    'SAVESPC':stt['False'], 'SPCFILE':'None', 'SAVESUM':stt['False'], 'SPCSUMFILE':'None',\
    'UNITS':'MKS', 'WMIN':stt['wmin'], 'WMAX':stt['wmax'], 'DIVERGENCE':stt['False'],\
    'VASEMIT':stt['True'],'CORRVAS':stt['False']}

    param_list = genfromtxt(parfile,skip_header=1, comments='#',dtype='|S10',usecols=(1), delimiter = '',unpack='True')
    unexp['OUTPATH'] = param_list[0]
    unexp['SAVEEPS'] = param_list[1]
    unexp['EPSROOT'] = param_list[2]
    unexp['SAVESPC'] = param_list[3]
    unexp['SPCFILE'] = param_list[4]
    unexp['SAVESUM'] = param_list[5]
    unexp['SPCSUMFILE'] = param_list[6]
    unexp['UNITS'] = param_list[7]
    unexp['WMIN'] = float(param_list[8])
    unexp['WMAX'] = float(param_list[9])
    unexp['DIVERGENCE'] = param_list[10]
    unexp['VASEMIT'] = param_list[11]
    unexp['CORRVAS'] = param_list[12]

    assert unexp['UNITS'] in stt['_units'],'units = {} ??'.format(unexp['UNITS'])

    return unexp
