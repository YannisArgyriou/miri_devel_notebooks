"""
NAME:

       _MKHDR

PURPOSE:

       To produce a header, in the form of a list of strings, with
       information about the configuration of MTS+MIRI.

EXPLANATION:

       It is used in "report" and "report_sum" to contextualize the
       information given in the screen and in the output text files.

CALLING SEQUENCE:

       PHDR = _MKHDR(VDATA, CONFIG, STT)

INPUTS:

       VDATA  - Structure with the spectral information (See initialize.pro
                for more information).
       CONFIG - structure with the configuration of the System (MTS [+ MIRI]).
                The relevant parameter is "units":
                      'MKS' : W m-2 um-1
                      'CGS' : erg s-1 cm-2 um-1
                      'PH'  : ph s-1 cm-2 um-1
                      'JY'  : 1.E-26 W m-2 Hz-1

       STT    - "internal" structure with a bunch of parameters which
                characterize the instruments and the code.

OUTPUTS:

       PHDR   - A list of "character strings" with the desired information.

OPTIONAL INPUTS KEYWORDS:

       None.
"""

def _MKHDR(Vdata, config, stt):
    from astropy.io import fits

    if config['units'] == 'MKS': beaunits = 'W m-2 um-1'
    if config['units'] == 'CGS': beaunits = 'erg s-1 cm-2 um-1'
    if config['units'] == 'PH': beaunits = 'ph s-1 cm-2 um-1'
    if config['units'] == 'JY': beaunits = 'JY'

    hdr = fits.PrimaryHDU()
    hdr.header['MTSSim version'] = (stt['version'],'Version')
    hdr.header['T_BB'] = (config['T_BB'],'Black Body Temperature')
    hdr.header['FW_MTS'] = (config['Usefilter'],'MTS-FW position')
    hdr.header['VASap'] = (int(config['VASap']),'MTS-VAS aperture')
    hdr.header['VASCOMM'] = (config['VASCOMM'],'commanded! MTS-VAS aperture')
    hdr.header['Target'] = (config['trgt'],'Position of the SSS')
    if config['pssON'] == stt['True']: spssON = 'True'
    else: spssON = 'False'
    hdr.header['PSSON'] = (spssON,'Is PSS on?')
    hdr.header['MIRI'] = (config['MIRI'],'MIRI Mode')
    hdr.header['nw'] = (Vdata['nw'],'Nr. of spectral bins')
    hdr.header['dw'] = (Vdata['dw']*1.E4,'Spectral sampling (AA)')
    hdr.header['units'] = (beaunits,'Units of irradiance')
    hdr.header['wmin'] = (config['wmin'],'Min. wavelength, um')
    hdr.header['wmax'] = (config['wmax'],'Max. wavelength, um')
    if config['divergence'] == stt['True']: sdivergence = 'True'
    else: sdivergence = 'False'
    hdr.header['divergence'] = (sdivergence,'Does the BB beam diverge?')

    return hdr
