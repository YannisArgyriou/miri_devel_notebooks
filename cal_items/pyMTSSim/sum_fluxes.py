def _sum_asist(tag, stttags, sttcws, sttpbs, spc, BKG, e_spc, wave, spc_u, w_u, stt):
    from numpy import where
    from _flx_w_sum import _flx_w_sum

    _cw = sttcws[where(stttags == tag)[0]]
    _pb = sttpbs[where(stttags == tag)[0]]

    wmin = _pb[0]
    wmax = _pb[1]

    flux = _flx_w_sum(spc, e_spc, wave, wmin, wmax, spc_u, w_u, stt, e_flux=e_flux)

    fluxB = _flx_w_sum(BKG, e_spc, wave, wmin, wmax, spc_u, w_u,stt)

    data = {'flux':flux, 'fluxB':fluxB, 'e_flux':e_flux, 'wmin':wmin, 'wmax':wmax, 'cw':_cw}


    return data

"""
NAME:

      SUM_FLUXES

PURPOSE:

      This function retrieves the integrals of the spectra of the target and
      the background sources, over different spectral ranges, depending on
      the MTS and MIRI configurations.

EXPLANATION:

      Nothing else to declare.

CALLING S==UENCE:

      VDATA = SUM_FLUXES(VDATA, CONFIG, STT)

INPUTS:

      VDATA   - Structure with the spectral information (See initialize.pro
                for more information).
      CONFIG  - structure with the configuration of the System (MTS [+ MIRI]).
      STT     - "internal" structure with a bunch of parameters which
                characterize the instruments and the code.

OUTPUTS:

       FLUXES - Structure with information about selected fluxes, such as
                band name, flux, error in flux, wavelength, etc.


OPTIONAL INPUT KEYWorDS:

       None.
"""

def SUM_FLUXES(Vdata, config, stt):
    from numpy import nan,arange,where,isfinite
    from _find_lines import _FIND_LINES
    from _flx_w_sum import _flx_w_sum

    tags = ['None']
    flx_t = [nan]
    eflx_t = [nan]
    flx_bk = [nan]
    units = ['None']
    cw = [nan]
    pb = [nan]
    wunits = ['None']

    spc = Vdata['SPC']
    BKG = Vdata['BKG']
    e_spc = Vdata['E_SPC']
    wave = Vdata['wave']
    spc_u = Vdata['units']
    w_u = 'um'

    nspec = len(spc)


    if config['MIRI'] != 'OFF':
        spc_adus = Vdata['SPC_ADUs']
        e_spc_adus = Vdata['E_SPC_ADUs']
        BKG_adus = Vdata['BKG_ADUs']

    else:
        spc_adus = arange(nspec) * nan
        e_spc_adus = arange(nspec) * nan
        BKG_adus = arange(nspec) * nan

    if config['ispoint'] == stt['True']:
        spc_aper = Vdata['SPC_APER']
        e_spc_aper = Vdata['E_SPC_APER']
        bkg_aper = arange(nspec) * nan
    else:
        spc_aper = arange(nspec) * nan
        e_spc_aper = arange(nspec) * nan
        bkg_aper = arange(nspec) * nan

    if ((config['ispoint'] == stt['True']) & (config['MIRI'] != 'OFF')):
        spc_aper_adus = Vdata['SPC_APER_ADUs']
        e_spc_aper_adus = Vdata['E_SPC_APER_ADUs']
        bkg_aper_adus = arange(nspec) * nan
    else:
        spc_aper_adus = arange(nspec) * nan
        e_spc_aper_adus = arange(nspec) * nan
        bkg_aper_adus = arange(nspec) * nan

    # MIRIM

    ixMIRIM = where(config['MIRI'] == stt['MIRIMfltrs'])

    if len(ixMIRIM[0]) != 0:
        _pb = stt['MIRIMpbs'][ixmirim]
        _cw = stt['MIRIMcws'][ixmirim]

        wmin = config['wmin'] # BEWARE!
        wmax = config['wmax'] # BEWARE!

        # wmin = _pb[0]
        # wmax = _pb[1]

        flux,e_flux = _flx_w_sum(spc, e_spc, wave, wmin, wmax, spc_u, w_u, stt,e_flux=True)

        fluxB = _flx_w_sum(BKG, e_spc, wave, wmin, wmax, spc_u, w_u, stt)

        tags = tags+[config['MIRI']]
        flx_t = flx_t+[flux]
        eflx_t = eflx_t+[e_flux]
        flx_bk = flx_bk+[fluxB]
        units = units+[spc_u]
        cw = cw+[_cw]
        pb = pb+[_pb[1] - _pb[0]]
        wunits = wunits+[w_u]

    	if sum(isfinite(spc_adus)) != 0:

            flux_ADU,e_flux_ADUs = _flx_w_sum(spc_adus, e_spc_adus, wave, wmin, wmax, 'ADU', w_u, stt, e_flux=True)


            flux_BADU = _flx_w_sum(BKG_adus, e_spc_adus, wave, wmin, wmax, 'ADU', w_u,stt)

            tags = [tags,[config['MIRI']+'_ADU']]
            flx_t = [flx_t, [flux_ADU]]
            eflx_t = [eflx_t, [e_flux_ADUs]]
            flx_bk = [flx_bk, [flux_BADU]]
            units = [units, ['ADU']]
            cw = [cw, [_cw]]
            pb = [pb, [_pb[1] - _pb[0]]]
            wunits = [wunits, [w_u]]

    	if sum(isfinite(spc_aper)) != 0:

            flux_APER,e_flux_APER = _flx_w_sum(spc_aper, e_spc_aper, wave, wmin, wmax, spc_u, w_u, stt, e_flux = True)

            flux_BAPER = -99

            tags = [tags,[config['MIRI']+'_APER']]
            flx_t = [flx_t, [flux_APER]]
            eflx_t = [eflx_t, [e_flux_APER]]
            flx_bk = [flx_bk, [flux_BAPER]]
            units = [units, ['AP_'+spc_u]]
            cw = [cw, [_cw]]
            pb = [pb, [_pb[1] - _pb[0]]]
            wunits = [wunits, [w_u]]

    	if sum(isfinite(spc_aper_adus)) != 0:

            flux_APER_ADUs,e_flux_APER_ADUs = _flx_w_sum(spc_aper_adus, e_spc_aper_adus, wave, wmin, wmax, 'ADU', w_u, stt, e_flux = True)

            flux_BAPER_ADUs = -99

            tags = [tags,[config['MIRI']+'_AP-ADUs']]
            flx_t = [flx_t, [flux_APER_ADUs]]
            eflx_t = [eflx_t, [e_flux_APER_ADUs]]
            flx_bk = [flx_bk, [flux_BAPER_ADUs]]
            units = [units, ['AP_ADU']]
            cw = [cw, [_cw]]
            pb = [pb, [_pb[1] - _pb[0]]]
            wunits = [wunits, [w_u]]

    # CORONOGRAPH


    ixCorON = where(config['MIRI'] == stt['CORONfltrs'])

    if len(ixCorON[0]) != 0:

        _pb = stt['CorONpbs'][ixCorON]
        _cw = stt['CorONcws'][ixCorON]

        wmin = config['wmin'] # BEWARE!
        wmax = config['wmax'] # BEWARE!

        #wmin = _pb[0]
        #wmax = _pb[1]

        flux,e_flux = _flx_w_sum(spc, e_spc, wave, wmin, wmax, spc_u, w_u, stt, e_flux=True)

        fluxB = _flx_w_sum(BKG, e_spc, wave, wmin, wmax, spc_u, w_u, stt)

        tags = [tags,[config['MIRI']]]
        flx_t = [flx_t, [flux]]
        eflx_t = [eflx_t, [e_flux]]
        flx_bk = [flx_bk, [fluxB]]
        units = [units, [spc_u]]
        cw = [cw, [_cw]]
        pb = [pb, [_pb[1] - _pb[0]]]
        wunits = [wunits, [w_u]]


    	if sum(isfinite(spc_adus)) != 0:

            flux_ADU,e_flux_ADUs = _flx_w_sum(spc_adus, e_spc_adus, wave, wmin, wmax, 'ADU', w_u, stt, e_flux=True)

            flux_BADU = _flx_w_sum(BKG_adus, e_spc_adus, wave, wmin, wmax, 'ADU', w_u,stt)

            tags = [tags,[config['MIRI']+'_ADU']]
            flx_t = [flx_t, [flux_ADU]]
            eflx_t = [eflx_t, [e_flux_ADUs]]
            flx_bk = [flx_bk, [flux_BADU]]
            units = [units, ['ADU']]
            cw = [cw, [_cw]]
            pb = [pb, [_pb[1] - _pb[0]]]
            wunits = [wunits, [w_u]]

    	if sum(isfinite(spc_aper)) != 0:

            flux_APER,e_flux_APER = _flx_w_sum(spc_aper, e_spc_aper, wave, wmin, wmax, spc_u, w_u, stt, e_flux = True)

            flux_BAPER = -99

            tags = [tags,[config['MIRI']+'_APER']]
            flx_t = [flx_t, [flux_APER]]
            eflx_t = [eflx_t, [e_flux_APER]]
            flx_bk = [flx_bk, [flux_BAPER]]
            units = [units, ['AP_'+spc_u]]
            cw = [cw, [_cw]]
            pb = [pb, [_pb[1] - _pb[0]]]
            wunits = [wunits, [w_u]]

    	if sum(isfinite(spc_aper_adus)) != 0:

            flux_APER_ADUs,e_flux_APER_ADUs = _flx_w_sum(spc_aper_adus, e_spc_aper_adus, wave, wmin, wmax, 'ADU', w_u, stt, e_flux = True)

            flux_BAPER_ADUs = -99

            tags = [tags,[config['MIRI']+'_AP-ADUs']]
            flx_t = [flx_t, [flux_APER_ADUs]]
            eflx_t = [eflx_t, [e_flux_APER_ADUs]]
            flx_bk = [flx_bk, [flux_BAPER_ADUs]]
            units = [units, ['AP_ADU']]
            cw = [cw, [_cw]]
            pb = [pb, [_pb[1] - _pb[0]]]
            wunits = [wunits, [w_u]]

    # MIRI OFF/MRS, ETALON = OFF


    if (config['MIRI'] == 'OFF') or ((config['MIRI'] == 'MRS') & (config['EtalonON'] == stt['False'])):

        _pb = [config['wmin'], config['wmax']]
        _cw = (_pb[0]+_pb[1])/2.

        wmin = _pb[0]
        wmax = _pb[1]

        flux,e_flux = _flx_w_sum(spc, e_spc, wave, wmin, wmax, spc_u, w_u, stt,e_flux=True)

        fluxB = _flx_w_sum(BKG, e_spc, wave, wmin, wmax, spc_u, w_u, stt,e_flux=False)

        tags = tags+['WHOLE']
        flx_t = flx_t+ [flux]
        eflx_t = eflx_t+ [e_flux]
        flx_bk = flx_bk+ [fluxB]
        units = units+ [spc_u]
        cw = cw+ [_cw]
        pb = pb+ [wmax - wmin]
        wunits = wunits+ [w_u]

    	if sum(isfinite(spc_adus)) != 0:

            flux_ADU,e_flux_ADUs = _flx_w_sum(spc_adus, e_spc_adus, wave, wmin, wmax, 'ADU', w_u, stt, e_flux=True)

            flux_BADU = _flx_w_sum(BKG_adus, e_spc_adus, wave, wmin, wmax, 'ADU', w_u,stt)

            tags = tags+['WHOLE'+'_ADU']
            flx_t = flx_t+ [flux_ADU]
            eflx_t = eflx_t+ [e_flux_ADUs]
            flx_bk = flx_bk+ [flux_BADU]
            units = units+ ['ADU']
            cw = cw+ [_cw]
            pb = pb+ [wmax - wmin]
            wunits = wunits+ [w_u]

    	if sum(isfinite(spc_aper)) != 0:

            flux_APER,e_flux_APER = _flx_w_sum(spc_aper, e_spc_aper, wave, wmin, wmax, spc_u, w_u, stt, e_flux = True)

            flux_BAPER = -99

            tags = tags+['WHOLE'+'_APER']
            flx_t = flx_t+ [flux_APER]
            eflx_t = eflx_t+ [e_flux_APER]
            flx_bk = flx_bk+ [flux_BAPER]
            units = units+ ['AP_'+spc_u]
            cw = cw+ [_cw]
            pb = pb+ [wmax - wmin]
            wunits = wunits+ [w_u]

    	if sum(isfinite(spc_aper_adus)) != 0:

            flux_APER_ADUs,e_flux_APER_ADUs = _flx_w_sum(spc_aper_adus, e_spc_aper_adus, wave, wmin, wmax, 'ADU', w_u, stt, e_flux = True)

            flux_BAPER_ADUs = -99

            tags = tags+['WHOLE'+'_AP-ADUs']
            flx_t = flx_t+ [flux_APER_ADUs]
            eflx_t = eflx_t+ [e_flux_APER_ADUs]
            flx_bk = flx_bk+ [flux_BAPER_ADUs]
            units = units+ ['AP_ADU']
            cw = cw+ [_cw]
            pb = pb+ [wmax - wmin]
            wunits = wunits+ [w_u]

    # ETALONS (& MRS)


    if (config['EtalonON'] == stt['True']) & (config['MIRI'] == 'MRS'):

        wminEta = (stt['waveminEt'][where(config['Usefilter'] == stt['Etalons_names'])[0]])[0]
        wmaxEta = (stt['wavemaxEt'][where(config['Usefilter'] == stt['Etalons_names'])[0]])[0]

        chans = stt['mrs_chans']
        nchans = len(chans)
        subchs = stt['mrs_subch']
        nsubchs = len(subchs)

        # INTEGRALS OVER CHANNELS

        stttags = stt['mrs_chans']
        sttcws = stt['mrs_cws']
        sttpbs = stt['mrs_pbs']

        for ix in range(nchans):

            chan = chans[ix]

            tag = chan
            _pb = sttpbs[where(stttags == tag)[0]]

    	    if ((_pb[0] > wminEta) or (_pb[1] < wmaxEta)):

                _data = _sum_asist(tag, stttags, sttcws, sttpbs, spc, BKG, e_spc, wave, spc_u, w_u, stt)

                tags = tags+[tag]
                flx_t = flx_t+ [_data['flux']]
                eflx_t = eflx_t+ [_data['e_flux']]
                flx_bk = flx_bk+ [_data['fluxB']]
                units = units+ [spc_u]
                cw = cw+ [_data['cw']]
                pb = pb+ [_data['wmax'] - _data['wmin']]
                wunits = wunits+ [w_u]


    	    if sum(isfinite(spc_adus)) != 0:

    	        _data = _sum_asist(tag, stttags, sttcws, sttpbs, spc_adus, BKG_adus, e_spc_adus, wave, 'ADU', w_u, stt)

    	        tags = tags+[tag+'_ADU']
    	        flx_t = flx_t+ [_data['flux']]
                eflx_t = eflx_t+ [_data['e_flux']]
                flx_bk = flx_bk+ [_data['fluxB']]
                units = units+ ['ADU']
                cw = cw+ [_data['cw']]
                pb = pb+ [_data['wmax'] - _data['wmin']]
                wunits = wunits+ [w_u]

    	    if sum(isfinite(spc_aper)) != 0:

                _data = _sum_asist(tag, stttags, sttcws, sttpbs, spc_aper, BKG_APER, e_spc_aper, wave, 'AP_'+spc_u, w_u, stt)

                tags = tags+[tag + '_APER']
                flx_t = flx_t+ [_data['flux']]
                eflx_t = eflx_t+ [_data['e_flux']]
                flx_bk = flx_bk+ [-99.0]
                units = units+ ['AP_' + spc_u]
                cw = cw+ [_data['cw']]
                pb = pb+ [_data['wmax'] - _data['wmin']]
                wunits = wunits+ [w_u]

    	    if sum(isfinite(spc_aper_adus)) != 0:

                _data = _sum_asist(tag, stttags, sttcws, sttpbs, spc_aper_adus, BKG_APER_ADUs, e_spc_aper_adus, wave, 'AP_ADU', w_u, stt)

                tags = tags+[tag+'_AP-ADUs']
                flx_t = flx_t+ [_data['flux']]
                eflx_t = eflx_t+ [_data['e_flux']]
                flx_bk = flx_bk+ [-99]
                units = units+ ['AP_ADU']
                cw = cw+ [_data['cw']]
                pb = pb+ [_data['wmax'] - _data['wmin']]
                wunits = wunits+ [w_u]

        # INTEGRALS OVER SUBCHANNELS

        stttags = stt['mrs_subch']
        sttcws = stt['mrs_cws_subch']
        sttpbs = stt['mrs_pbs_subch']

        for jx in range(nsubchs):

            subchan = subchs[jx]
            tag = subchan

            _pb = sttpbs[where(stttags == tag)[0]]
    	    if (_pb[0] < wminEta) or (_pb[1] > wmaxEta):

                _data = _sum_asist(tag, stttags, sttcws, sttpbs, spc, BKG, e_spc, wave, spc_u, w_u, stt)

                tags = tags+[tag]
                flx_t = flx_t+ [_data['flux']]
                eflx_t = eflx_t+ [_data['e_flux']]
                flx_bk = flx_bk+ [_data['fluxB']]
                units = units+ [spc_u]
                cw = cw+ [_data['cw']]
                pb = pb+ [_data['wmax'] - _data['wmin']]
                wunits = wunits+ [w_u]


    	    if sum(isfinite(spc_adus)) != 0:

    	        _data = _sum_asist(tag, stttags, sttcws, sttpbs, spc_adus, BKG_adus, e_spc_adus, wave, 'ADU', w_u, stt)

    	        tags = tags+[tag+'_ADU']
    	        flx_t = flx_t+ [_data['flux']]
                eflx_t = eflx_t+ [_data['e_flux']]
                flx_bk = flx_bk+ [_data['fluxB']]
                units = units+ ['ADU']
                cw = cw+ [_data['cw']]
                pb = pb+ [_data['wmax'] - _data['wmin']]
                wunits = wunits+ [w_u]

    	    if sum(isfinite(spc_aper)) != 0:

                _data = _sum_asist(tag, stttags, sttcws, sttpbs, spc_aper, BKG_APER, e_spc_aper, wave, 'AP_'+spc_u, w_u, stt)

                tags = tags+[tag + '_APER']
                flx_t = flx_t+ [_data['flux']]
                eflx_t = eflx_t+ [_data['e_flux']]
                flx_bk = flx_bk+ [-99.]
                units = units+ ['AP_' + spc_u]
                cw = cw+ [_data['cw']]
                pb = pb+ [_data['wmax'] - _data['wmin']]
                wunits = wunits+ [w_u]

    	    if sum(isfinite(spc_aper_adus)) != 0:

                _data = _sum_asist(tag, stttags, sttcws, sttpbs, spc_aper_adus, BKG_APER_ADUs, e_spc_aper_adus, wave, 'AP_ADU', w_u, stt)

                tags = tags+[tag+'_AP-ADUs']
                flx_t = flx_t+ [_data['flux']]
                eflx_t = eflx_t+ [_data['e_flux']]
                flx_bk = flx_bk+ [-99.]
                units = units+ ['AP_ADU']
                cw = cw+ [_data['cw']]
                pb = pb+ [_data['wmax'] - _data['wmin']]
                wunits = wunits+ [w_u]

        # INTEGRALS OVER LINES

        stttags = stt['mrs_subch']
        sttcws = stt['mrs_cws_subch']
        sttpbs = stt['mrs_pbs_subch']

        for jx in range(nsubchs):

            subchan = subchs[jx]
            tag = subchan
            _pb = sttpbs[where(stttags == tag)[0]]

    	    if (_pb[0] < wminEta) or (_pb[1] > wmaxEta):

                ixwave = where((wave >= _pb[0]) & (wave <= _pb[1]))[0]

                disc = sum(spc[ixwave])

                if (sum(isfinite(disc)) == stt['True'] & disc > 0.):

                    swave = wave[ixwave]
                    sspc = spc[ixwave]
                    se_spc = e_spc[ixwave]
                    sBKG = BKG[ixwave]

                    extrema = _FIND_LINES(swave,sspc,CONTRAST = 3.)
                    wavelines = extrema['wavelines']
                    waveminima = extrema['waveminima']

                    # "BLUE" LINE

                    bline = wavelines[0]
                    ix_b = where(waveminima < bline)[0]
                    if ix_b[0] < 0:
                        bline = wavelines[1]
                        ix_b = where(waveminima < bline)[0]

		        ix_b = ix_b[len(ix_b)-1]
                ix_r = where(waveminima > bline)[0]

                brange = [waveminima[ix_b], waveminima[ix_r]]

                ltag = tag + '_B'
                lstttags = [ltag]

                lsttcws = {ltag:bline}
                lsttpbs = {ltag:[brange[0], brange[1]]}

                _linedata = _sum_asist(ltag, lstttags, lsttcws, lsttpbs, sspc, sBKG, se_spc, swave, spc_u, w_u, stt)

                tags = tags+ [ltag]
                flx_t = flx_t+ [_linedata['flux']]
                eflx_t = eflx_t+ [_linedata['e_flux']]
                flx_bk = flx_bk+ [_linedata['fluxB']]
                units = units+ [spc_u]
                cw = cw+ [_linedata['cw']]
                pb = pb+ [_linedata['wmax'] - _linedata['wmin']]
                wunits = wunits+ [w_u]

                if sum(isfinite(spc_adus)) != 0:

                    _linedata = _sum_asist(ltag, lstttags, lsttcws, lsttpbs, spc_adus[ixwave], BKG_adus[ixwave], e_spc_adus[ixwave], wave[ixwave], 'ADU', w_u, stt)

                    tags = tags+[ltag+'_ADU']
                    flx_t = flx_t+ [_linedata['flux']]
                    eflx_t = eflx_t+ [_linedata['e_flux']]
                    flx_bk = flx_bk+ [_linedata['fluxB']]
                    units = units+ ['ADU']
                    cw = cw+ [_linedata['cw']]
                    pb = pb+ [_linedata['wmax'] - _linedata['wmin']]
                    wunits = wunits+ [w_u]

    	        if sum(isfinite(spc_aper)) != 0:

                    _linedata = _sum_asist(ltag, lstttags, lsttcws, lsttpbs, spc_aper[ixwave], BKG_APER[ixwave], e_spc_aper[ixwave], wave[ixwave], 'AP_'+spc_u, w_u, stt)

                    tags = tags+[ltag + '_APER']
                    flx_t = flx_t+ [_linedata['flux']]
                    eflx_t = eflx_t+ [_linedata['e_flux']]
                    flx_bk = flx_bk+ [-99.]
                    units = units+ ['AP_' + spc_u]
                    cw = cw+ [_linedata['cw']]
                    pb = pb+ [_linedata['wmax'] - _linedata['wmin']]
                    wunits = wunits+ [w_u]

    	        if sum(isfinite(spc_aper_adus)) != 0:

                    _linedata = _sum_asist(ltag, lstttags, lsttcws, lsttpbs, spc_aper_adus[ixwave], BKG_APER_ADUs[ixwave], e_spc_aper_adus[ixwave], wave[ixwave], 'AP_ADU', w_u, stt)

                    tags = tags+[ltag + '_AP-ADUs']
                    flx_t = flx_t+ [_linedata['flux']]
                    eflx_t = eflx_t+ [_linedata['e_flux']]
                    flx_bk = flx_bk+ [-99.]
                    units = units+ ['AP_ADU']
                    cw = cw+ [_linedata['cw']]
                    pb = pb+ [_linedata['wmax'] - _linedata['wmin']]
                    wunits = wunits+ [w_u]


                # "CENTRAL" LINE


                cline = wavelines[len(wavelines)/2]

                ix_b = where(waveminima < cline)[0]
                ix_b = ix_b[len(ix_b)]
                ix_r = where(waveminima > cline)[0]

                crange = [waveminima[ix_b], waveminima[ix_r]]

                ctag = tag + '_C'
                lstttags = [ctag]

                lsttcws = {ctag:cline}
                lsttpbs = {ctag:[crange[0], crange[1]]}

                _linedata = _sum_asist(ctag, lstttags, lsttcws, lsttpbs, sspc, sBKG, se_spc, swave, spc_u, w_u, stt)

                tags = tags+ [ctag]
                flx_t = flx_t+ [_linedata['flux']]
                eflx_t = eflx_t+ [_linedata['e_flux']]
                flx_bk = flx_bk+ [_linedata['fluxB']]
                units = units+ [spc_u]
                cw = cw+ [_linedata['cw']]
                pb = pb+ [_linedata['wmax'] - _linedata['wmin']]
                wunits = wunits+ [w_u]


                if sum(isfinite(spc_adus)) != 0:

                    _linedata = _sum_asist(ctag, lstttags, lsttcws, lsttpbs, spc_adus[ixwave], BKG_adus[ixwave], e_spc_adus[ixwave], wave[ixwave], 'ADU', w_u, stt)

                    tags = tags+[ctag+'_ADU']
                    flx_t = flx_t+ [_linedata['flux']]
                    eflx_t = eflx_t+ [_linedata['e_flux']]
                    flx_bk = flx_bk+ [_linedata['fluxB']]
                    units = units+ ['ADU']
                    cw = cw+ [_linedata['cw']]
                    pb = pb+ [_linedata['wmax'] - _linedata['wmin']]
                    wunits = wunits+ [w_u]

    	        if sum(isfinite(spc_aper)) != 0:

                    _linedata = _sum_asist(ctag, lstttags, lsttcws, lsttpbs, spc_aper[ixwave], BKG_APER[ixwave], e_spc_aper[ixwave], wave[ixwave], 'AP_'+spc_u, w_u, stt)

                    tags = tags+[ctag + '_APER']
                    flx_t = flx_t+ [_linedata['flux']]
                    eflx_t = eflx_t+ [_linedata['e_flux']]
                    flx_bk = flx_bk+ [-99.]
                    units = units+ ['AP_' + spc_u]
                    cw = cw+ [_linedata['cw']]
                    pb = pb+ [_linedata['wmax'] - _linedata['wmin']]
                    wunits = wunits+ [w_u]

    	        if sum(isfinite(spc_aper_adus)) != 0:

                    _linedata = _sum_asist(ctag, lstttags, lsttcws, lsttpbs, spc_aper_adus[ixwave], BKG_APER_ADUs[ixwave], e_spc_aper_adus[ixwave], wave[ixwave], 'AP_ADU', w_u, stt)

                    tags = tags+[ctag + '_AP-ADUs']
                    flx_t = flx_t+ [_linedata['flux']]
                    eflx_t = eflx_t+ [_linedata['e_flux']]
                    flx_bk = flx_bk+ [-99]
                    units = units+ ['AP_ADU']
                    cw = cw+ [_linedata['cw']]
                    pb = pb+ [_linedata['wmax'] - _linedata['wmin']]
                    wunits = wunits+ [w_u]


                # "RED" LINE

                rline = wavelines[len(wavelines)-1]
                ix_r = (where(waveminima > rline))[0]
                if ix_r < 0:
                    rline = wavelines[len(wavelines)-2]
                    ix_r = (where(waveminima > rline))[0]

                ix_b = (where(waveminima < rline))
                ix_b = ix_b[len(ix_b)-1]

                rrange = [waveminima[ix_b], waveminima[ix_r]]

                rtag = tag + '_R'
                lstttags = [rtag]
                lsttcws = {rtag:rline}
                lsttpbs = {rtag:[rrange[0],rrange[1]]}

                _linedata = _sum_asist(rtag, lstttags, lsttcws, lsttpbs, sspc, sBKG, se_spc, swave, spc_u, w_u, stt)

                tags = tags+ [rtag]
                flx_t = flx_t+ [_linedata['flux']]
                eflx_t = eflx_t+ [_linedata['e_flux']]
                flx_bk = flx_bk+ [_linedata['fluxB']]
                units = units+ [spc_u]
                cw = cw+ [_linedata['cw']]
                pb = pb+ [_linedata['wmax'] - _linedata['wmin']]
                wunits = wunits+ [w_u]

                if sum(isfinite(spc_adus)) != 0:

                    _linedata = _sum_asist(rtag, lstttags, lsttcws, lsttpbs, spc_adus[ixwave], BKG_adus[ixwave], e_spc_adus[ixwave], wave[ixwave], 'ADU', w_u, stt)

                    tags = tags+[rtag+'_ADU']
                    flx_t = flx_t+ [_linedata['flux']]
                    eflx_t = eflx_t+ [_linedata['e_flux']]
                    flx_bk = flx_bk+ [_linedata['fluxB']]
                    units = units+ ['ADU']
                    cw = cw+ [_linedata['cw']]
                    pb = pb+ [_linedata['wmax'] - _linedata['wmin']]
                    wunits = wunits+ [w_u]

    	        if sum(isfinite(spc_aper)) != 0:

                    _linedata = _sum_asist(rtag, lstttags, lsttcws, lsttpbs, spc_aper[ixwave], BKG_APER[ixwave], e_spc_aper[ixwave], wave[ixwave], 'AP_'+spc_u, w_u, stt)

                    tags = tags+[rtag + '_APER']
                    flx_t = flx_t+ [_linedata['flux']]
                    eflx_t = eflx_t+ [_linedata['e_flux']]
                    flx_bk = flx_bk+ [-99.]
                    units = units+ ['AP_' + spc_u]
                    cw = cw+ [_linedata['cw']]
                    pb = pb+ [_linedata['wmax'] - _linedata['wmin']]
                    wunits = wunits+ [w_u]

    	        if sum(isfinite(spc_aper_adus)) != 0:

                    _linedata = _sum_asist(rtag, lstttags, lsttcws, lsttpbs, spc_aper_adus[ixwave], BKG_APER_ADUs[ixwave], e_spc_aper_adus[ixwave], wave[ixwave], 'AP_ADU', w_u, stt)

                    tags = tags+[rtag + '_AP-ADUs']
                    flx_t = flx_t+ [_linedata['flux']]
                    eflx_t = eflx_t+ [_linedata['e_flux']]
                    flx_bk = flx_bk+ [-99.]
                    units = units+ ['AP_ADU']
                    cw = cw+ [_linedata['cw']]
                    pb = pb+ [_linedata['wmax'] - _linedata['wmin']]
                    wunits = wunits+ [w_u]


    nflx = len(tags)

    fluxes = {'tags':tags[1:nflx], 'flx_t':flx_t[1:nflx], 'eflx_t':eflx_t[1:nflx],
    'flx_bk':flx_bk[1:nflx], 'units':units[1:nflx], 'cw':cw[1:nflx], 'pb':pb[1:nflx],
    'wunits':wunits[1:nflx]}


    return fluxes
