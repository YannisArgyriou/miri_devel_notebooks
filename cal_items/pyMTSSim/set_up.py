"""
NAME:

      SET_UP

PURPOSE:

      Based on the inputs given to MTSSim, it retrieves a single structure
      that holds the configuration of the system (MTS [+MIRI]), that is used by
      different subroutines to yield their results.

EXPLANATION:

      Not required. This function has only interest within this program.

CALLING SEQUENCE:

      CONFIG = SET_UP(T_BB, FILTER, VASAP, TARGET, STT, pssON=pssON,
      MIRI=MIRI, PARFILE=PARFILE, INSTRUCT=INSTRUCT)

INPUTS:

      T_BB     - Black Body Temperature (0 < T_BB <= 800 K)
      FILTER   - Position of the Filter Wheel of the MTS. Possible values are:
                      "HOLE", "SWP", "LWP", "BLOCK", "ET_1A", "ET_2A",
                      "ET_1B", "ET_2B"
      VASAP    - Aperture of the VAS, within 0. <= VASap <= 100.
      TARGET   - Position of the "Source Scanning Subsystem" (SSS). One of:
                      "pnh", "pnh2" or "ext".
      STT      - "internal" structure with a bunch of parameters which
                 characterize the instruments (internal processing variables,
                 measures, names of files, etc.).

OUTPUTS:

      CONFIG   - structure with the configuration of the System.

OPTIONAL INPUT KEYWORDS:

      pssON    - is the Pupil Scanning System "ON"?
      MIRI     - status of MIRI. One of:
                   OFF, MRS, [one of the Imager Filters], [One of the
                   coronograph filters]
      PARFILE  - file with additional inputs (see "read_params.pro").

      INSTRUCT - Structure with additional inputs (see "read_params.pro").
"""
def SET_UP(T_BB=None, Usefilter=None, VASap=None, target=None, stt=None, pssON=None,\
        MIRI=None, parfile=None, instruct = None):
    from numpy import nan
    from read_params import READ_PARAMS
    # CHECK BB TEMPERATURE
    assert (T_BB >= 0) & (T_BB <= 1000.),'T_BB must be within (0,1000.)K & T_BB = {}'.format(T_BB)

    T_BB = float(T_BB)

    # MTS-FILTER WHEEL
    assert Usefilter in stt['FltrWheel'],'Unknown filter (MTS): {}'.format(Usefilter)

    # Is any of the Etalons in the optical path?

    if Usefilter[0:3] == 'ET_': EtalonON = stt['True']
    else: EtalonON = stt['False']

    # CHECK the VAS aperture is within range (0. - 100)
    assert (VASap >= -0.5) & (VASap <= 100.),'VASap must be within (-0.5,100.) & VASap = {}'.format(VASap)

    VASap = float( VASap )

    # Target : [pnh, pnh2, ext]
    assert target in ['PNH','PNH2','EXT'],'Unknown target: {}'.format(target)

    if target in ['PNH' 'PNH2']: ispoint = stt['True']
    else: ispoint = stt['False']

    # Is the Pupil Scanning System ON?
    if pssON is None : pssON = stt['False']
    else: pssON = stt['True']

    if pssON == stt['True']:
        assert Usefilter == 'BLOCK', 'If PSS is "ON", THEN filter MUST BE "BLOCK". QUITTING...'

    if pssON == stt['True']: ispoint = stt['False']

    # MIRI CONFIGURATION:
    #     OFF, [one of MIRIM filters], [one of CORON filters], 'MRS'
    if MIRI is None: MIRI = 'OFF'

    assert MIRI in ['OFF', stt['MIRIMfltrs'], stt['CORONfltrs'], 'MRS'],'Unknown MIRI configuration: {}. QUITTING...'.format(MIRI)

    # "instruct" is a structure with aditional parameters about the configuration
    # of the system. It is optional, and alternatively, a file with those
    # parameters may be passed, i.e. "parfile". If none of these, "instruct" or
    # "parfile" is passed to the program, default values, set up in read_params.py
    # are assumed.

    unexp = {'OUTPATH':'OUTPUTS', 'SAVEEPS':stt['False'], 'EPSROOT':'None',\
    'SAVESPC':stt['False'], 'SPCFILE':'None', 'SAVESUM':stt['False'], 'SPCSUMFILE':'None',\
    'UNITS':'MKS', 'WMIN':stt['wmin'], 'WMAX':stt['wmax'], 'DIVERGENCE':stt['False'],\
    'VASEMIT':stt['True'],'CORRVAS':stt['False']}

    if instruct is not None :
    	unexp['OUTPATH'] = instruct['OUTPATH']
    	unexp['SAVEEPS'] = instruct['SAVEEPS']
    	unexp['EPSROOT'] = instruct['EPSROOT']
    	unexp['SAVESPC'] = instruct['SAVESPC']
    	unexp['SPCFILE'] = instruct['SPCFILE']
        unexp['SAVESUM'] = instruct['SAVESUM']
    	unexp['SPCSUMFILE'] = instruct['SPCSUMFILE']
    	unexp['UNITS'] = instruct['UNITS']
    	unexp['WMIN'] = instruct['WMIN']
    	unexp['WMAX'] = instruct['WMAX']
    	unexp['DIVERGENCE'] = instruct['DIVERGENCE']
    	unexp['VASEMIT'] = instruct['VASEMIT']
    	unexp['CORRVAS'] = instruct['CORRVAS']

    # READING and PARSING the optional "parfile"

    if (instruct is None) & (parfile is not None):
        unexp = READ_PARAMS(parfile, stt)

    if (instruct is None) & (parfile is None):
        print 'Using defaults values in set_up.py for the "secondary" parameters...'
        unexp = read_params(parfile, stt)

    config = {'T_BB':T_BB, 'VASap':VASap, 'VASCOMM':VASap, 'Usefilter':Usefilter,\
    'beam_divergence':unexp['DIVERGENCE'], 'EtalonON':EtalonON, 'pssON':pssON,\
    'trgt':target, 'MIRI':MIRI,\
    'saveeps':unexp['SAVEEPS'], 'epsroot':unexp['EPSROOT'],\
    'savespc':unexp['SAVESPC'], 'spcfile':unexp['SPCFILE'],\
    'savesum':unexp['SAVESUM'], 'spcsumfile':unexp['SPCSUMFILE'],\
    'outpath':unexp['OUTPATH'], 'units':unexp['UNITS'], 'wmin':unexp['WMIN'],\
    'wmax':unexp['WMAX'], 'divergence':unexp['DIVERGENCE'],\
    'VASemit':unexp['VASEMIT'], 'ispoint':ispoint,\
    'corrVAS':unexp['CORRVAS']}

    return config
