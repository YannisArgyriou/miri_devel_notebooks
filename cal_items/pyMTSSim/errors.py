"""
NAME:

    ERRORS

PURPOSE:

    Provide with at least some rough estimate of the uncertainty in the
    Irradiance from the MTS-calibration source at the MIRI input plane.

CALLING SEQUENCE:

    VDATA = ERRORS(VDATA, CONFIG, STT)

INPUTS:

    VDATA  - Structure with the spectral information (See initialize.pro
             for more information).
    CONFIG - structure with the configuration of the System (MTS [+ MIRI]).

    STT    - "internal" structure with a bunch of parameters which
             characterize the instruments and the code.

OUTPUT:

    VDATA - Modified version, with uncertainties added.
"""
def ERRORS(Vdata, config, stt):
    from numpy import floor

    wave = Vdata['wave']        # um

    # Calculates error in % (0-1)

    # [c2_PLANK] = cm K

    # delta_TBB : fractional uncertainty in the BB Temp
    # from 100 to 800K in 100K steps

    T_BB = config['T_BB']

    ixTBB = int(floor(T_BB/100.)-1.)
    if T_BB < 100.: ixTBB = 0
    assert T_BB <= 800.,'I did not expect a BB-TEMP = {}'.format(T_BB)

    delta_TBB = stt['delta_TBB'][ixTBB]

    D_BB_Irr = (delta_TBB * T_BB * stt['c2_PLANCK']) / (wave*1.E-4 * T_BB**2.)
                                            # Unc in BB-Irradiance

    D_r_pnhbb = 0.00525

    if config['ispoint'] == stt['True']: D_trgt = 0.005
    else: D_trgt = 0.
    if ((config['VASap'] > 0.) & (config['VASap'] < 100.)): D_VAS=0.005
    else: D_VAS = 0.
    D_IS = 0.03
    if ((config['Usefilter'] != 'HOLE') & (config['Usefilter'] != 'BLOCK')):D_fltr = 0.02
    else: D_fltr = 0.
    if (config['EtalonON'] == stt['True']): D_etal = 0.002
    else: D_etal = 0.

    D_MOS = 0.004   # (x4)
    D_gold = 0.0064 # (x4)

    if config['Usefilter'] != 'BLOCK':
        D_MTS_Irr_rel = D_BB_Irr + D_r_pnhbb + D_trgt + \
        D_VAS + D_IS + D_fltr + D_etal + 4. * D_gold + \
        4. * D_MOS
    else: D_MTS_Irr_rel = 0.

    Vdata['E_SPC'] = Vdata['SPC'] * D_MTS_Irr_rel
    if config['ispoint'] == stt['True']:
       Vdata['E_SPC_APER'] = Vdata['SPC_APER'] * D_MTS_Irr_rel

    return Vdata
