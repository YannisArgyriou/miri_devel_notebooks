"""
NAME:

    SSS

PURPOSE:

    Include the effects of the target plate (Source Scanning Subsystem, SSS)
    on the beam:
            - Converts from Radiance (L) to Irradiances (E) at the Input
              Plane of MIRI.
            - Adds background.

EXPLANATION:


CALLING SEQUENCE:

    VDATA = SSS(VDATA, CONFIG, STT)

INPUTS:

    VDATA  - Structure with the spectral information (See initialize.pro
              for more information).
    CONFIG - Structure with the configuration of the System (MTS [+ MIRI]).

    STT    - "Internal" structure with a bunch of parameters which
              characterize the instruments and the code.

OUTPUTS:

    VDATA  - (after modifications).

OPTIONAL INPUT KEYWORDS:

    None
"""
def SSS(Vdata, config, stt):
    from numpy import pi,sin,cos,arctan,sqrt,exp,zeros
    from planck import planck

    spcL = Vdata['SPC']   # L
    wave = Vdata['wave']
    trgt = config['trgt']

    assert Vdata['spc_id'] == 'L','Expecting units of "L", but got {} instead.'.format(Vdata['spc_id'])

    # EP_dia : MTS exit pupil diameter (mm).
    # d_Image: Distance between exit pupil and image plane (mm).

    r_MTS = (stt['EP_dia'] * 1.E-3) / 2.                                 # m
    phi_MTS = arctan(r_MTS / (stt['d_Image'] * 1.E-3))                   # rad


    # Pinhole radius (negative for extended source, for "procedural" reasons.)
    if trgt == 'EXT': r_p = -1.
    if trgt == 'PNH': r_p = (stt['d_p_100'] * 1.E-3) / 2.        # m
    if trgt == 'PNH2': r_p = (stt['d_p_25'] * 1.E-3) / 2.        # m


    # Extended source

    # Irradiance at the MIRI input plane [E]

    # EXTENDED SOURCE

    if r_p < 0:
    	Efactor = 2. * pi * (1. - cos(phi_MTS)) # For YANNIS EDIT replace by: pi * (sin(phi_MTS))**2.
    	spcE = spcL * Efactor

        Vdata['SPC'] = spcE
        Vdata['spc_id'] = 'E'

        # BACKGROUND

        # BB

    	# Changed after conversation with T. Belenger on 16th/09/2010
            #Vdata['BKG_BBFRAME'] *= pi * (1. / (2. * stt['nF_MTS']))**2. #
            #Vdata['BKG_BBSTRAP'] *= pi * (1. / (2. * stt['nF_MTS']))**2. #

    	#Vdata['BKG_BBFRAME'] *= pi * (sin(phi_MTS))**2.
            #Vdata['BKG_BBSTRAP'] *= pi * (sin(phi_MTS))**2.

        # VAS

    	# In v. < 2.1 we had:
            # Vdata['BKG_VAS'] *= pi * (1. / (2. * stt['nF_MTS']))**2. #
            # In v>= 2.1:
    	#Vdata['BKG_VAS'] *= pi * (sin(phi_MTS))**2.

        # WHY NOT TO DO THE SAME AS THE SPECTRUM OF THE BB?
    	# In fact, both expressions are valid under a small-angle
    	# approximation... so, it has been changed by the "sin"
    	# version (21st/09/2010)!!

    	# ON TESTS

    	Vdata['BKG_BBFRAME'] *= Efactor
    	Vdata['BKG_BBSTRAP'] *= Efactor
    	Vdata['BKG_VAS'] *= Efactor

        Vdata['BKG_id'] = 'E'


    # PINHOLES

    if r_p >= 0:

        Ap = pi * r_p**2.   # Area of the pinhole, m2

    	PowFctr = Ap * 2. * pi * (1. - cos(phi_MTS))  * stt['Magnification']**2.
    	spcP = spcL * PowFctr

        # Value at the Peak

        r_Airy = 1.22 * (wave*1.E-6) * stt['nF_MTS']  # m

        Ap_Airy = pi * r_Airy**2.		   # m^2


        # Gaussian Approximation to an Airy Diffraction Pattern:
        # I = I0 * exp(-r^2 / (2*sigma^2))

        sigma = 0.42 * (wave*1.E-6) * stt['nF_MTS']   # m
        fwhm = sigma * 2.355                       # m

        # Irradiance at the peak of the PSF

        gaussareafactor = (2.* pi * sigma**2.)

        px = sqrt(Vdata['PX'][0,:] * Vdata['PX'][1,:]) # linear size of pixels,
	                                         # projected in the input plane!,
	                                         # in m.

	    # Irradiance at r = 0

        spcE0_aprox0 = spcP / gaussareafactor

    	# To compute the slope at the brightest pixel, the average irradiance
    	# within a pixel must be computed. This value shall be smaller than the
    	# irradiance at r=0. As a proxy, I get the average of the irradiance
    	# within a circle inscribed in a pixel, and a circle that circumscribes
    	# the pixel. The pixel size is an "effective" dimension, to account for
    	# the non squared pixels of MRS.


    	Rin = px / 2.         # Radius of a circle inscribed in
    	                      # a projected-pixel.
    	Rout = px / sqrt(2.)  # Radius of a circle that circumscribes
    	                      # a projected-pixel(pjpixel).

    	# Radiation power contained in Rin
    	P_in = 2. * pi * sigma**2. * spcE0_aprox0 * \
    	    (1. - exp(-Rin**2./(2.*sigma**2.)))

    	# Radiation power contained in Rout
    	P_out = 2. * pi * sigma**2. * spcE0_aprox0 * \
    	    (1. - exp(-Rout**2./(2.*sigma**2.)))

    	E0in = P_in / (pi*Rin**2.) # Average irradiance within Rin
    	E0out = P_out / (pi*Rout**2.) # Average irradiance within Rout

    	# Average Irradiance within 1 proj-pixel, given as the average of E0in
    	# and E0out... just a little bit "cheese".

    	# If you multiply this quantity by the area of 1 proj-pixel you would
    	# obtain the spectral power within 1 proj-pixel.

    	spcE0_aprox1 = (E0in+E0out) / 2.

        Vdata['SPC'] = spcE0_aprox1

        Vdata['SPC_APER'] = spcP

    	Vdata['spc_id'] = 'E'

    	Vdata['fwhm'] = fwhm                                        # m
    	Vdata['fwhm_pix'][0,:] = fwhm / Vdata['PX'][0,:]               # pixels
    	Vdata['fwhm_pix'][1,:] = fwhm / Vdata['PX'][1,:]               # pixels

        # BACKGROUND

        # BB

        # Input is RADIANCE, L, Output is E

        BKG_BBFrame_L = Vdata['BKG_BBFRAME']
        #BKG_BBFrame_P = BKG_BBFrame_L * Ap * Factor_dist
        #BKG_BBFrame_P *= pi * (1. /(2. * stt['nF_MTS']))**2.
                                                           # Power reaching pinhole

    	BKG_BBFrame_P = BKG_BBFrame_L * PowFctr

    	BKG_BBFrame_E = BKG_BBFrame_P / Ap_Airy    # !!
        Vdata['BKG_BBFRAME'] = BKG_BBFrame_E

        BKG_BBStrap_L = Vdata['BKG_BBSTRAP']
        BKG_BBStrap_P = BKG_BBStrap_L * PowFctr

    	#BKG_BBstrap_P = BKG_BBstrap_L * Ap * Factor_dist
        #BKG_BBstrap_P *= pi * (1. /(2. * stt['nF_MTS']))**2.
                                                           # Power reaching pinhole
        BKG_BBstrap_E = BKG_BBstrap_P /  Ap_Airy   # !!
        Vdata['BKG_BBstrap'] = BKG_BBstrap_E

        BKG_VAS_L = Vdata['BKG_VAS']
    	BKG_VAS_P = BKG_VAS_L * PowFctr

        #BKG_VAS_P = BKG_VAS_L * Ap * Factor_dist
        #BKG_VAS_P *= pi * (1. /(2. * stt['nF_MTS']))**2.
                                                           # Power reaching pinhole

    	BKG_VAS_E = BKG_VAS_P / Ap_Airy           # !!
        Vdata['BKG_VAS'] = BKG_VAS_E

        Vdata['BKG_id'] = 'E'

    # BACKGROUND CONTRIBUTION OF THE TARGET (only if a point source is used)

    if config['ispoint'] == stt['True']:

        BKG_TRGT_M = planck(wave*1.E4, stt['T_target']) * stt['epsilon_dif']
        BKG_TRGT_M *= 1.E1                                     # W m-2 um-1
        #BKG_TRGT_E = BKG_TRGT_M * (1. / (2. * stt['nF_MTS']))**2. # ??
        BKG_TRGT_E = BKG_TRGT_M * 2. * (1. - cos(phi_MTS)) # ??

        Vdata['BKG_TRGT'] = BKG_TRGT_E * stt['TargetFctr']

        Vdata['BKG_id'] = 'E'
    else:
        Vdata['BKG_TRGT'] = zeros(len(Vdata['wave']))

    return Vdata
