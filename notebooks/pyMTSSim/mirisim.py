def _switch2ph(units, spc, wave, stt):
    from _mks2cgs import _MKS2cgs
    from _mks2jy import _MKS2jy
    from _mks2ph import _MKS2ph

    # changes from different units to ph s-1 cm-2 um-1

    if units == 'PH': return spc
    if units == 'CGS': mks = _MKS2cgs(spc, REVERSE=True)
    if units == 'MKS': mks = spc
    if units == 'JY': mks = _MKS2jy(spc, wave, stt, REVERSE=True)

    return _MKS2ph(mks, wave, stt)

def _switch2ph_aper(units,spc,wave,stt):

    # changes from different units to ph s-1 um-1

    if units == 'PH': return spc

    if units == 'CGS': mks = _MKS2cgs(spc, REVERSE=True) * 1.E-4
    if units == 'MKS': mks = spc
    if units == 'JY': mks = _MKS2jy(spc, wave, stt, REVERSE=True)

    return _MKS2ph(mks, wave, stt) * 1.E4

"""
NAME:

       MIRISIM

PURPOSE:

       To provide with a first order aproximation of the count rates, in
       ADU s-1 pix-1, at the focal planes of MIRI.

EXPLANATION:

       WARNING: Only the quantum efficiency of the detectors, the
       electron-to-ADU conversion rate and the transmission curves of the
       filters are considered, but not the losses due to reflections within
       MIRI, slit-losses, or any other effects due to the design and/or
       imperfections of the MIRI optics. Thus, these predictions should be
       considered as what they are, mere proxies at best. For detailed
       and accurate predictions of the MIRI response to the MTS inputs, you
       must use MTSSim in combination with the simulators of the MIRI
       subsystems: MIRIMSim (for the MIRIM and the LRS), and SpecSim (for
       the MRS).

CALLING SEQUENCE:

       VDATA = MIRISIM(VDATA, CONFIG, STT)

INPUTS:

       VDATA  - Structure with the spectral information (See initialize.pro
                for more information).
       CONFIG - structure with the configuration of the System (MTS [+ MIRI]).
       STT    - "internal" structure with a bunch of parameters which
                characterize the instruments and the code.

OUTPUTS:

       VDATA  - (after modifications).

OPTIONAL INPUT KEYWORDS:

       None.

"""
def MIRISIM(Vdata, config, stt):
    from numpy import genfromtxt,where
    from scipy.interpolate import interp1d
    from loadtrans import LoadTrans

    MIRI = config['MIRI']

    if MIRI == 'OFF':
        return Vdata
    else:
        assert Vdata['spc_id'] == 'E','Expecting units of "E", but got {} instead.'.format(Vdata['units'])

        # Applying TRANSMISSION CURVES
        #            &
        # CONVERSION to DN s-1

        wave = Vdata['wave']

        # Reading quantum efficiency data:
        quantumf = stt['Opath'] + stt['quantumf']

        w_qe, qe_106, qe_105, qe_104 = genfromtxt(quantumf,skip_header=3,usecols=(0,2,4,6),unpack=True)

        if MIRI == 'MRS':
            qeSW = interp1d(w_qe,qe_105)(wave)
            qeLW = interp1d(w_qe,qe_104)(wave)

            ich12 = where(wave < stt['mrs_pbs']['CHAN2'][1])
            ich34 = where(wave >= stt['mrs_pbs']['CHAN2'][1])

            qe = wave * 0.
            if len(ich12[0]) != 0: qe[ich12] = qeSW[ich12]
            if len(ich34[0]) != 0: qe[ich34] = qeLW[ich34]

        isImager = stt['False']
        ixcheck = where(MIRI == stt['MIRIMfltrs'])
        if ixcheck[0] >= 0: isImager = stt['True']

        isCORON = stt['False']
        ixcheck = where(MIRI == stt['CORONfltrs'])
        if ixcheck[0] >= 0: isCORON = stt['True']

        if ((isImager == stt['True']) or (isCORON == stt['True'])):
            qe = interp1d(w_qe,qe_106)(wave)


        if isImager == stt['True']:
            ixMIRIMfltrs = where(MIRI == stt['MIRIMfltrs'])

            ixtag = str(ixMIRIMfltrs[0]+1)

            filtern = stt['Fpath'] + 'im' + ixtag + '_edit.dat'

            T = LoadTrans(filtern, wave)

            # Applying transmission curve to the target

            Vdata['SPC'] *= T
            Vdata['E_SPC'] *= T

            if config['ispoint'] == stt['True']:
                Vdata['SPC_APER'] *= T
                Vdata['E_SPC_APER'] *= T

            # Applying transmission curve to background

            for ixBKG in range(len(Vdata['BKGtags']) - 1):
                BKGtag = Vdata['BKGtags'][ixBKG]
                ixtag = where(Vdata.keys() == BKGtag)
                Vdata[ixtag] *= T  # applying transmission curve to backgrounds

        if isCORON == stt['True']:

            ixCORONfltrs = where(MIRI == stt['CORONfltrs'])

            filtern = stt['Cpath'] + 'cor0' + \
            split(str(ixCORONfltrs[0]+1),2) + '_edit.dat'

            T = LoadTrans(filtern, wave)

            # Applying transmission curve to the target

            Vdata['SPC'] *= T
            Vdata['E_SPC'] *= T

            if config['ispoint'] == stt['True']:
                    Vdata['SPC_APER'] *= T
                    Vdata['E_SPC_APER'] *= T

            # Applying transmission curve to background

            for ixBKG in range(len(Vdata['BKGtags']) - 1):
                BKGtag = Vdata['BKGtags'][ixBKG]
                ixtag = where(Vdata.keys() == BKGtag)
                Vdata[ixtag] *= T  # applying transmission curve to backgrounds

        if MIRI == 'MRS':
            wave = Vdata['wave']
            nchans = len(stt['mrs_chans'])
            nschans = len(stt['mrs_subch'])

            T = wave * 0.

            bolf = stt['Tpath'] + stt['mrs_bolf']
            ibT = LoadTrans(bolf,wave) * 100.

            for ischan in range(nschans):

                schan = stt['mrs_subch'][ischan][-2:]
                chan = schan[0]
                pb = stt['mrs_pbs_subch']['s'+schan]

                dicf = stt['Tpath']+stt['mrs_dicf']+schan+'.dat'
                mirf = stt['Tpath']+stt['mrs_mirf']+chan+'.dat'

                ixw = where((wave >= pb[0]) & (wave <= pb[1]))

                if len(ixw[0]) != 0:
                    idT = LoadTrans(dicf,wave[ixw]) * 100.
                    imT = LoadTrans(mirf,wave[ixw]) * 100.
                    T[ixw] = idT * imT * ibT[ixw]

            # Applying transmission curve to the target
            Vdata['SPC'] *= T
            Vdata['E_SPC'] *= T

    	if config['ispoint'] == stt['True']:
            Vdata['SPC_APER'] *= T
            Vdata['E_SPC_APER'] *= T

        # Applying transmission curve to background

        for ixBKG in range(len(Vdata['BKGtags']) - 1):
             BKGtag = Vdata['BKGtags'][ixBKG]
             Vdata[BKGtag] *= T  # applying transmission curve to backgrounds

        dnsfctr = qe / stt['dnphotrat']

        apix_cm2 = (Vdata['PX'][0,:] * Vdata['PX'][1,:]) * 1.E4

        spc = Vdata['SPC']
        e_spc = Vdata['E_SPC']
        BKG = Vdata['BKG']

        spc_ph = _switch2ph(Vdata['units'], spc, wave, stt)
        e_spc_ph = _switch2ph(Vdata['units'], e_spc, wave, stt)

        BKG_ph = _switch2ph(Vdata['units'], BKG, wave, stt)


        Vdata['BKG_ADUs'] = BKG_ph * dnsfctr * apix_cm2
                                                        # ADU s-1 pix-1 um-1
                                                        # BEWARE, VALUES PER um!
        Vdata['SPC_ADUs'] = spc_ph * dnsfctr * apix_cm2
                                                        # ADU s-1 pix-1 um-1
        Vdata['E_SPC_ADUs'] = e_spc_ph * dnsfctr * apix_cm2
                                                        # ADU s-1 pix-1 um-1


        if config['ispoint'] == stt['True']:

            spc_aper = Vdata['SPC_APER']
            e_spc_aper = Vdata['E_SPC_APER']

            spc_aper_ph = _switch2ph_aper(Vdata['units'], spc_aper, wave, stt)
            e_spc_aper_ph = _switch2ph_aper(Vdata['units'], e_spc_aper, wave, stt)

            Vdata['SPC_APER_ADUs'] = spc_aper_ph * dnsfctr
            Vdata['E_SPC_APER_ADUs'] = e_spc_aper_ph * dnsfctr


        return Vdata

"""
These estimates of the transmissions for MRS are no longer used. These
notes are only kept for the record.

Garbage Collector:

MRS peak tranmission.
(OBA, MIRI-DD-00001-AEU_Iss2_OBA-Design-Description.PDF, p.148-9)
CH1: 0.85 * 0.98^4 = 0.78
   CH2: 0.80 * 0.98^4 = 0.74
   CH3: 0.65 * 0.98^7 = 0.57
   CH4: 0.50 * 0.98^3 = 0.47
"""
