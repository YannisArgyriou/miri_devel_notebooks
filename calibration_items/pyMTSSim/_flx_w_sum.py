"""
NAME:

       _FLX_W_SUM

PURPOSE:

       To retrieve the wavelength integral of a given spectrum, within
       certain bounds.

EXPLANATION:

       It is used by sum_fluxes.pro to retrieve the fluxes in different
       wavelength intervals or spectral lines (etalons).
       In case of using a point source, it will retrieved spatially
       integrated fluxes, i.e. "radiation power".

CALLING SEQUENCE:

       FLUX = _FLX_W_SUM(SPC, E_SPC, WAVE, WMIN, WMAX, SPC_U, W_U, STT,$
               E_FLUX=E_FLUX)

INPUTS:

       SPC   - Spectrum
       E_SPC - Uncertainties in "SPC"
       WAVE  - wavelengths for SPC
       WMIN  - Lower integration limit, in wavelength.
       WMAX  - Upper integration limit, in wavelength.
       SPC_U - Units of "SPC". One of:
                       - MKS : W m-2 um-1
                       - CGS : ergs s-1 cm-2 um-1
                       - PH : photons s-1 cm-2 um-1
                       - JY : 1.E-26 W m-2 Hz-1
       W_U   - Units of "WAVE". By now, only microns ("um") are considered.
       STT   - "internal" structure with a bunch of parameters which
               characterize the instruments and the code.

OUTPUTS:

       FLUX - Integral of spc over the specified interval in wavelength.

OPTIONAL INPUT KEYWorDS:

       E_FLUX - If this keyword is used it will retrieve the uncertainty in
       in SPC, given E_SPC.
"""

def _flx_w_sum(spc, e_spc, wave, wmin, wmax, spc_u, w_u, stt, e_flux=False):
    from numpy import mean,ceil,arange,trapz,argsort
    from scipy.interpolate import interp1d

    fresample = 3. # The spectral sampling is multiplied by this factor,
                   # for the integration.

    if w_u == 'um': wf = 1.E-6

    assert w_u == 'um','Expecting wavelengths in [um]... but got {} ??'.format(w_u)

    # TRIMMING AND RESAMPLING

    if ((wmin < min(wave)) or (wmax > max(wave))):
         raise ValueError,'Integration range exceeds wavelength coverage.'


    nw = len(wave)
    dw = mean(wave[1:nw-1]-wave[:nw-2])
    ndw = dw / fresample
    nn_w = ceil((wmax-wmin)/ndw)
    ndw = ( wmax - wmin ) / (nn_w-1)
    nwave = arange(nn_w) * ndw + wmin

    nspc = interp1d(wave,spc)(nwave)
    ne_spc = interp1d(wave,e_spc)(nwave)

    nwave0 = nwave[0]

    if spc_u == 'JY':
        nx = nwave0 / nwave
    else:
        nx = nwave

    order = argsort(nx)
    nx = nx[order]
    nspc = nspc[order]
    ne_spc = ne_spc[order]

    flux = trapz(nspc,x=nx)

    midpoint = int(nn_w/2)
    err_flux = ne_spc[midpoint] / nspc[midpoint] * flux

    if spc_u == 'JY':
        fconv = 1.*1.E-26 * stt['clight'] / (nwave0 * wf)
        flux *= fconv         # W [m-2]
        err_flux *= fconv       # W [m-2]

    if e_flux:
        return flux,err_flux
    else:
        return flux
