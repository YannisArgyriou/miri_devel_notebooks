"""
NAME:

    MTSSim

PURPOSE:

    MIRI-Telescope-Simulator (MTS) Radiometric Simulator.
         (yep, it makes simulations of a simulator...)

EXPLANATION:

    This program provides estimates of the Irradiances at the input plane
    of MIRI given a configuration of the Miri Telescope Simulator (MTS).
    It is intended to serve as a reference for the planning and
    interpretation of Functional and Performance Tests of MIRI.

CALLING SEQUENCE:

    Vdata = MTSSim(T_BB, FILTER, VASAP, TARGET, [/pssON], [MIRI], [parfile],
            [instruct],[fluxes])

INPUTS:

    T_BB   - Temperature of the Black Body.
    FILTER - Filter set in the MTS.
    VASAP  - Aperture of the VAS (0 - 100%,
                  0% meaning fully closed).
    TARGET - Position of the Source Scanning System.

OUTPUTS:

    This function returns "VDATA", a structure with the spectral
    information (See initialize.pro for more information).

    It also produces graphs with the spectra (they may be saved as
    ".eps" files), and a summary of the configuration and the fluxes
    (which may be also stored as an ascii file). Optionally, it is
    also possible to store an ascii file with the spectra in columns.

OPTIONAL INPUT KEYWORDS:

    pssON    - If set, the Pupil Scanning System is "on". Otherwise
               it is "off".
    MIRI     - Status of MIRI. One of:
                 - OFF,
                 - {One of MIRIM filters},
                 - {One of the Coronograph filters},
                 - MRS.
    PARFILE  - Plain text file with additional inputs.
                      See "inputs_example.txt" for inspiration.
    INSTRUCT - Structure with additional inputs.
                      If "parfile" exists, "instruct" is ignored.
    FLUXES   - If named, it will retrieve a structure with fluxes.

REVISION HISTORY:

    Ruyman Azzollini - version 2.3.7 (April 2011).

    This a major upgrade, devised by Ruyman Azzollini (CSIC), over
    MTSSim vv. 0-1.9, which was written and maintained by Alvaro
    Labiano (ESAC) between 2006 and 2010.

    MTSSim vv. 0-1.9 were based on the Mathcad programs by INTA and the
    MIRI-MTS Optical and radiomentric analysis document, with other
    inputs from the INTA/MTS team.

    Alvaro thanked:

      - the INTA team, especially Tomas Belenguer and Raquel Lopez
        for their help and support with the radiometric analysis and
        design of the MTS.

      - Luis Colina, Almudena Alonso-Herrero for supervision and
        directions and Macarena Garcia Marin & Tanio Diaz Santos for their
        help with IDL.

      - Alejandro Bedregal, Bruno Merin and Eva Bauwens for
        testing and suffering the first versions of the program!

      - the Photometry Group for their continuous patience and
        suggestions.

    Ruyman wants to thank Alvaro Labiano for eating the "worst and largest
    part of the cake" (this is mostly a "remake"), the Miri Software Group
    for their support and valuable feedback, and also wants to extend his
    gratitude to those listed above.

    And big thanks to the Spanish National Football Team for bringing home
    our first shiny golden little World Cup on July 11th 2010... :-).

    VERSION HISTORY:

- See "version_history.txt" in the source directory.
"""
def MTSSim(T_BB, Usefilter, VASap, target, pssON=None,\
    MIRI=None, parfile=None, instruct=None, fluxes=None):
    # INITIALIZATION
    'Syntax -- Vdata = MTSSim(T_BB, Usefilter, VASap, target,  pssON=pssON, miri=miri, parfile=parfile, instruct=instruct,fluxes=fluxes)'
    """
    Some useful "definitions" before we start:

        - Spectral Radiation Power, P, W um-1, radiation power emitted/received, per
                wavelength interval.
        - Spectral Radiance, L, W m-2 sr-1 um-1
                radiation power EMITTED, per solid angle unit and unit area,
                and wavelength interval.
        - Spectral Emittance, M, W m-2 um-1
                radiation power EMITTED in a hemisphere, per unit area,
                and wavelength interval. For a lambertian source,
                M = L * pi
        - Spectral Irradiance, E, W m-2 um-1 (received)
                radiation power RECEIVED, per unit area and wavelength
                interval.


    The main data structure is "Vdata", which stores the following items:

        - vector of wavelengths ("wave"), its length ("nw") and sampling ("dw").
        - "spc": vector with the spectrum (its definition varies along the code)
          of the MTS calibration source.
        - "e_spc": vector of uncertainty of "spc".
        - "spc_id": physical amount that the spectrum of the target represents:
                M (emittance), L (radiance), P (power) or E (irradiance).
        - vectors with the spectra of the background sources ("BKG" and "BKG_*").
        - "spc_id": physical amount that the spectra of the background sources
                represent: M (emittance), L (radiance), P (power) or
                E (irradiance).
        - "BKGtags": list of the different Background Sources.
        - OPTIONAL (only if MIRI is not "OFF"):
                - vector with spectral irradiance in ADUs s-1 pix-1 ("spc_DNs").
                - vector with uncertainty in spc_DNs ("e_spc_DNs").
                - vector with background spectrum in ADUs s-1 pix-1 ("BKG_DNs").


    Give values to the "config" structure, which stores information relative
    to the configuration of the MTS and MIRI. It will halt the program if the
    inputs are not consistent with the capabilities and functionality of the
    system.
    """
    from internalvalues import stt
    from set_up import SET_UP
    from initialize import INITIALIZE
    from bb import BB
    from col import COL
    from vas import VAS
    from fw import FW
    from etalons import ETALONS
    from integrsphere import IS
    from sss import SSS
    from mos import MOS
    from fpss import fPSS
    from fms import FMs
    from oba import OBA
    from pom import POM
    from addup_bkg import ADDUP_BKG
    from _convert import _CONVERT
    from errors import ERRORS
    from report import REPORT
    from mirisim import MIRISIM
    from sum_fluxes import SUM_FLUXES


    # INITIALIZATION ----------------------------------------------------------
    config = SET_UP(T_BB, Usefilter, VASap, target, stt, pssON=pssON,\
    MIRI=MIRI, parfile=parfile,instruct=instruct)


    # INITIALIZE the main data structure
    # which holds the spectral vectors: Vdata

    Vdata = INITIALIZE(config, stt)

    # CALCULATIONS  HERE -------------------------------------------------

    # BB

    Vdata = BB(Vdata, config, stt) # P, [W um-1]

    # COLLIMATOR

    Vdata = COL(Vdata, config, stt)

    # VAS

    Vdata = VAS(Vdata, config, stt)

    # FILTER WHEEL (Etalons are worked out apart,
    # given that their calculations are far more involved).

    if config['EtalonON'] == stt['False']:  Vdata = FW(Vdata, config, stt)

    # ETALONs

    if config['EtalonON'] == stt['True']:  Vdata = ETALONS(Vdata, config, stt)

    # IS

    Vdata = IS(Vdata, config, stt)

    # TARGETS

    Vdata = SSS(Vdata, config, stt)

    # MOS

    Vdata = MOS(Vdata, config, stt)

    # PSS

    if config['pssON'] == stt['True']:  Vdata = fPSS(Vdata, config, stt)

    # FMs: add background and apply FM transmission curve
    Vdata = FMs(Vdata, config, stt)

    # OBA : add background
    Vdata = OBA(Vdata, stt)

    # POM and POM mirror: add background
    Vdata = POM(Vdata, config, stt)

    # ADDING-UP BACKGROUND SOURCES

    Vdata = ADDUP_BKG(Vdata, config, stt)

    # CONVERSION OF UNITS

    Vdata = _CONVERT(Vdata, config, stt)         # Converts irradiances to
                                                 # desired units.

    # UNCERTAINTIES

    Vdata = ERRORS(Vdata, config, stt)

    # REPORT (at MIRI input plane)

    if config['savespc']:
        REPORT(Vdata, config, stt)              # Dumps fits file with
                                                # results (spectra).

    # MIRI : IMAGER, MRS (Etalons)

    Vdata = MIRISIM(Vdata, config, stt)

    # INTEGRATED FLUXES (bands & lines)

    fluxes = SUM_FLUXES(Vdata, config, stt)

    return Vdata
