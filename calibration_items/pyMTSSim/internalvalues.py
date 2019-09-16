"""
 NAME:

        INTERNALVALUES

        This is not a procedure, nor a function, just a collection
        of variables.

 PURPOSE:

        Here are kept different parameters that characterize the MTS, and
        are used by MTSSim and its subroutines. They are all "bunched" together
        in a data structure, "stt". It also contains other accesory variables
        needed by the program, such as the value of some booleans.

 EXPLANATION:

        To make the simulations, a considerable number of parameters must be
        given as input, and it would be cumbersome to pass all of them as
        individual parameters or keywords to functions. Moreover, it is better
        to have them together in a single file, for easier inspection.

 CALLING SEQUENCE:

        @internalvalues.pro

 INPUTS:

        None.

 OUTPUTS:

        A data structure named "stt" is defined with all the parameters.

 OPTIONAL INPUT KEYWORDS:

        None.


-
 COMMON
"""
import numpy as np

version = '2.6.1'                        # version.

dw_lores = 0.05                          # um.
dw_hires = 5.E-5                         # um.
dw_hires = 5.E-3                         # only in tests... to run faster, like
                                         #  Fernando Alonso.
dw_print = 5.E-4                     # um.
wmin = 4.7                               # um..
wmax = 30.                               # um.
FltrWheel = ['HOLE','SWP','LWP','BLOCK','ET_1A','ET_1B','ET_2A','ET_2B']
root = ''                                # root path.
Tpath = 'TRANSMIT/'                       # Path to Transmission Curves.
Fpath = 'MIRIM_FILTERS_FM/'
Cpath = 'CORON_FILTERS_FM/'
Mpath = 'MTS_FW/'
Opath = 'OTHER_FILES/'
_units = ['MKS','CGS','PH','JY']
_flxunits = {'MKS':'W m-2','AP_MKS':'W','CGS':'erg s-1 cm-2',\
'AP_CGS':'erg s-1','PH':'ph s-1 cm-2','AP_PH':'ph s-1',\
'JY':'W m-2','AP_JY':'W', 'ADU':'ADU pix-1 s-1', 'AP_ADU':'ADU s-1'}

# Logical
BinaryTrue = 1
BinaryFalse = 0

# Physical constants
clight = 2.998*1.E8                    # Speed of light in Madrid
                                        #    (with fluid traffic), [m s-1].
hplanck = 6.6261*1.E-34                # Planck's constant (Mo. to Fr.) [J s].

d_Image = 3017.56                       # Distance exit pupil to image plane (mm).
ro_optica = 0.98                        # Reflectance of the optics, worst case scenario
d_IR_Pupil =11.826+5.5                  # Distance (mm) IR source to exit pupil plane.

nF_MTS = 19.98                          # F of the MTS in the image space.
Magnification = 1.887                   # MTS INVERSE magnification.
					  # (not used in v.>=2.1, but kept for record).
tau_opt = 0.85                           # Transmittance of the whole
                                         # MTS optical chain.
epsilon_mirror = 0.04                    # Emissivity of Mirrors.

# BB
dia_pinhole = 2.*0.4965                  # Diameter (mm) of pin-hole in-front of BB.
epsilon_BB = 0.95                        # Black body emissivity.

T_bbstrap = 39.2                         # T of the BB-strap.
epsilon_bbstrap = 0.1                    # Emissivity of the BB-strap
                                         #  PRIVATE COMUNICATION,
					  # TOMAS B. (16/09/2010)!!??
BBStrapFctr = 1.                         # To scale the emission, if needed.
Area_bbstrap = 1.E-4                     # 1 cm2 UNKNOWN!

T_bbframe = 53.3                         # T of the BB-frame.
epsilon_bbframe = 0.1                    # Emissivity of the BB-frame.
					  # PRIVATE COMUNICATION
					  # TOMAS B. (16/09/2010)!!??
BBFrameFctr = 1.                         # To scale the emission, if needed.
Area_bbframe = 1.E-4                     # 1cm2, UNKNOWN!!

delta_TBB = 0.00279                      # fractional uncertainty in the BB Temp

delta_TBB = np.array([0.3,0.5,0.8,1.1,1.4,1.6,1.9,2.2])**1.E-2

                                         # fractional uncertainty in the
					 # BB Temp for 100 to 800K in 100K
					 # steps.

c1_PLANCK = 3.7417749*1.E-5              #   = 2*pi*h*c*c, [erg cm2 s-1].
c2_PLANCK = 1.4387687                    # = h*c/k, [cm K].


# COLLIMATOR
focal_collimator = 80.                  # Focal length collimator -- BB (mm)
ap_collimator = 20.                     # Aperture (mm).

# VAS

epsilon_VAS = 0.6                       # Emissivity of the VAS.
                                        #  Assumed to be like that of the
                                        #  POM, in lack of a
					  # better "guess"...!!
VAScorrf = 'vas_corr.txt'		 # File with the relation between
					  # the commanded position of the VAS
					  # and its apparent real aperture.
T_VAS = 37.5                            # Temp of the VAS, for BKG emission.
# T_VAS = 100.                            # Temp of the VAS, for BKG emission!!
VASFctr = 1.                            # To scale the emission, if needed.

# FILTER WHEEL

LWP_fltr = 'LWP_20.dat'			 # LWP filter
SWP_fltr = 'SWP_20.dat'			 # SWP filter

# ETALONS
R_et = 0.93                             # Reflectance of coatings.
thetain_et = 0.                         #  Incidence angle (deg) of ray
                                        #   entering etalon.
n_et = 1.                               #  Refractive index for AIR-GAP.
d_ch = np.array([1.59,1.838,1.6,1.702]) #  Physical length of the cavity (mm),
                                        #   per channel.
Etalons_names = ['ET_1A','ET_1B','ET_2A','ET_2B']  # Etalons list.
waveminEt = [5.,7.71,11.89,18.35]         # Etalons' lambdas (low).
wavemaxEt = [7.71,11.89,18.35,28.3]       # Etalons' lambdas (high).
Aa = 0.                                   # Absorption coefficient (cm-1).
A_Ch = Aa * d_ch * 1.E-4                  # Absorption (0 in case of AIR-GAP
                                          #  etalon).
nm_par = 4.7E-2                           # Parallelism (micron).
nm_flat = 6.33E-2                         # Flatness (micron).
nm_rms = 1.62E-4                          # Roughness (micron).
Etalons_filtern = {'ET_1A':'etalon1a.dat','ET_1B':'etalon1b.dat',\
    'ET_2A':'etalon2a.dat','ET_2B':'etalon2b.dat'}
                                          # Etalon Transmission Curves.

# IS
R_sphere = 60.                            # Radius of the IS (mm).
dist_coll_IS = 80.                        # Distance (mm) from collimator to IS.
ap_input_port = 22.                       # Diameter of input port (mm).
ap_output_port = 77.                      # Side of the aperture port (mm).
                                          # MTSDESIGN!!??
ro_is = ro_optica                         # Throughput of the IS (now the curve
                                          # of transmission is used instead)
ISTransn = 'ISgold.dat'                   # Transmission Curve of the IS.
ISFctr = 1.                               # To scale the emission, if needed.


# TARGETS
EP_dia = 151.                             # MTS exit pupil diameter (mm).
dist_PH_IS = 46.182                       # Distance (mm) from IS to pinhole target.
T_target = 24.1                           # Temp of the target plate.
epsilon_dif = 0.6                         # Emissivity of the diffuse plate.
TargetFctr = 1.                           # To scale the emission if needed.
d_p_100 = 0.1018                          # 100 um pinhole diameter, mm.
d_p_25 = 0.025                            # 25 um pinhole diameter, mm.

# MOS

MOSTransn = 'MOScurve_capture_shift_blue0.06.txt' # CHANGED from  MOStransmit.dat (as used in MTSSim v2.6.1 of Ruyman Azzollini) to MOScurve_capture_shift_blue0.06.txt
                                          # Transmission Curve of each mirror
                                          #  of the MOS (there are 4)
MOSObsc = 0.053                           # Obscuration of the Pupil mask,
                                          #  equals 1/19.

# PSS
PssFctr = 1.                              # To calibrate PSS when real
                                          # measurements are available.
                                          # (CorrFact in MTSSim v1.9).
T_pss = 480.                              # Source temperature (K).
epsilon_pss = 0.8                         # PSS emissivity.
pssapert = 0.35                           # Aperture (mm).

# FMs
Ry_FM1 = 120.                             # y-Radius (mm).
Ry_FM2 = 119.                             # y-Radius (mm).
Ry_FM3 = 111.                             # y-Radius (mm).
Ry_FM4 = 112.                             # y-Radius (mm).
Rx_FM4=74.                                # x-Radius (mm).
d_prim_FM1 = 341.768                      # (mm).
d_FM1_FM2 = 583.599                       # (mm).
d_FM2_FM3 = 870.810                       # (mm).
d_FM3_FM4 = 237.627                       # (mm).
T_FM1 = 24.5                              # Temp of the FM1.
T_FM2 = 25.0                              # Temp of the FM2.
T_FM3 = 26.2                              # Temp of the FM3.
T_FM4 = 26.5                              # Temp of the FM4.
epsilon_FMs = 0.04                        # Emissivity of the FMs.
theta_FMs = 45.                           # Incidence angle on the FMs.
FMtransn = 'FM.dat'                       # Transmission Curve of each Folding Mirror.
FMFctr = 1.                               # To scale the emission, if needed.

# POM & POM mirror

d_MIRI_POM=450.                           # Distance from POM to image plane (mm).
Side_POM = 224.91                         # POM Side (mm).
T_POM = 25.1                              # Temp of the POM structure.
T_POM_M = T_POM                           # Temp of the POM mirror.
epsilon_POM = 0.6                         # Emissivity of the POM.
epsilon_POM_M = epsilon_mirror            # Emissivity of the POM mirror.
theta_POM=45.                             # Incidence angle of the surface of the
                                          # POM.
theta_POM_M= theta_POM                    # Incidence angle of the surface of
                                          # the POM mirror.
POMFctr = 1.                              # To scale the emission if needed.

# OBA
Side_OBA=1200.                            # Size at 45 deg (mm)
d_MIRI_OBA = 266.                         # Distance from OBA to image plane (mm).
epsilon_OBA = 0.6                         # Emissivity of the OBA.
T_OBA = 27.4                              # Temp of the MTS baseplate.
theta_OBA=89.427253                       # Incidence angle of the surface of the
                                          # OBA.
                                          #  THIS PARAMETER HAS BEEN CHANGED FROM
                                          #  0, TO ACCOUNT FOR A 0.01 FACTOR GIVEN
                                          #  IN THE ORIGINAL PROGRAM... TO MAKE IT
                                          #  FULLY CONSISTENT!!
OBAFctr = 1.                              # To scale the emission if needed.

# Fudge factor (introduced by Yannis Argyriou (KU Leuven) on Nov.12th 2018)
FUDGEFctr = 1./0.55

# MIRI

MIRIMfltrs = ['F560W','F770W','F1000W','F1130W','F1280W','F1500W','F1800W',\
'F2100W','F2550W','F2550WR']
MIRIMcws = {'F560W':5.6,'F770W':7.7,'F1000W':10.,'F1130W':11.3,'F1280W':12.8,'F1500W':15.0,\
'F1800W':18.0,'F2100W':21.0,'F2550W':25.5,'F2550WR':25.5}
MIRIMpbs = {'F560W':[5.,6.2],'F770W':[6.6,8.8], 'F1000W':[9.0,11.0],\
'F1130W':[10.95,11.65], 'F1280W':[11.6,14.0], 'F1500W':[13.5,16.5],\
'F1800W':[16.5,19.5], 'F2100W':[18.5,23.5], 'F2550W':[23.5,27.5],\
'F2550WR':[23.5,27.5]}

CORONfltrs = ['F1065C','F1140C','F1550C','F2300C']
CORONcws = {'F1065C':10.65,'F1140C':11.4,'F1550C':15.5,'F2300C':23.0}
CORONpbs = {'F1065C':[10.385,10.915],'F1140C':[11.115,11.685],\
    'F1550C':[15.11,15.89],'F2300C':[20.7,25.3]}

mrs_chans = ['CHAN1', 'CHAN2', 'CHAN3', 'CHAN4']
mrs_subch = ['S1A', 'S2A', 'S3A', 'S4A', 'S1B', 'S2B', 'S3B', 'S4B',\
'S1C', 'S2C', 'S3C', 'S4C']

mrs_cws_subch = {'s1A':5.37, 's2A':8.22, 's3A':12.575, 's4A':19.38,\
's1B':6.205, 's2B':9.49, 's3B':14.53, 's4B':22.59,'s1C':7.16,\
's2C':10.895, 's3C':16.78, 's4C':26.34}
mrs_pbs_subch = {'s1A':[4.92, 5.82], 's2A':[7.53,8.91], 's3A':[11.52,13.63],\
's4A':[17.65,21.11],'s1B':[5.68,6.73], 's2B':[8.69,10.29], 's3B':[13.31,15.75],\
's4B':[20.57,24.61],'s1C':[6.56,7.76], 's2C':[9.91,11.88], 's3C':[15.37,18.19],\
's4C':[23.99,28.69]}

mrs_cws = {'CHAN1':6.34, 'CHAN2':9.705, 'CHAN3':14.855, 'CHAN4':23.17}
mrs_pbs = {'CHAN1':[5.0,7.76], 'CHAN2':[7.53,11.88], 'CHAN3':[11.52,18.19],'CHAN4':[17.65,28.69]}
# Size of MRS pixels relative to imager pixels
mrs_pixrel = {'CHAN1':[1.614,1.620],'CHAN2':[1.564,2.586],'CHAN3':[2.038,3.656],'CHAN4':[2.263,6.062]}
mrs_pixrel = {'CHAN1':[1.61,1.62],'CHAN2':[1.56,2.59],'CHAN3':[2.04,3.66],'CHAN4':[2.26,6.06]}
mrs_peakTr = {'CHAN1': 0.78, 'CHAN2':0.74, 'CHAN3':0.57, 'CHAN4':0.47}  # NOT USED ANY MORE, in favor of more accurate values, and with finer dependance on wavelength.

mrs_dicf = 'Dichroic_Ch'
mrs_mirf = 'Mirrors_Ch'
mrs_bolf = 'BOL_Contam.dat'


# Taken from OBA-design report:
#   MIRI-DD-00001-AEU_Iss2_OBA-Design-Description.PDF, p.148-9)

MIRIstates =  [MIRIMfltrs,CORONfltrs,['MRS']]

dnphotrat = 5.95                         # Electron/DN ratio.
dnphotrat = 5.5				             # FM-VALUE!
FOV_Lineal = 72.                         # Lineal FoV (mm) At the aperture of MIRI.
num_px = 1024.                           # Number of pixels on a side of the detectors.
quantumf = 'quantum_eff.dat'             # quantum efficiency file.
quantumf = 'quantum_eff_FM.dat'          # (Model) FM quantum efficiencies.

#================================================================================


stt = {'version':version,\
'dw_lores':dw_lores, 'dw_hires':dw_hires, 'dw_print':dw_print,\
'wmin':wmin, 'wmax':wmax, 'FltrWheel':FltrWheel, 'root':root, 'Tpath':Tpath,\
'Opath':Opath, 'Fpath':Fpath, 'Cpath':Cpath, 'Mpath':Mpath, '_units':_units,\
'_flxunits':_flxunits, 'True':BinaryTrue, 'False':BinaryFalse,\
'clight':clight, 'hplanck':hplanck, 'd_Image':d_Image, 'ro_optica':ro_optica,\
'd_IR_Pupil':d_IR_Pupil, 'Magnification':Magnification,\
'nF_MTS':nF_MTS, 'tau_opt':tau_opt, 'epsilon_BB':epsilon_BB,\
'T_bbstrap':T_bbstrap, 'epsilon_bbstrap':epsilon_bbstrap,\
'T_bbframe':T_bbframe, 'epsilon_bbframe':epsilon_bbframe,\
'Area_bbstrap':Area_bbstrap, 'Area_bbframe':Area_bbframe,\
'BBStrapFctr':BBStrapFctr, 'BBFrameFctr':BBFrameFctr,\
'c1_PLANCK':c1_PLANCK, 'c2_PLANCK':c2_PLANCK, 'delta_TBB':delta_TBB,\
'dia_pinhole':dia_pinhole, 'focal_collimator':focal_collimator,\
'ap_collimator':ap_collimator, 'epsilon_VAS':epsilon_VAS,'VAScorrf':VAScorrf,\
'T_VAS':T_VAS, 'VASFctr':VASFctr,'LWP_fltr':LWP_fltr,'SWP_fltr':SWP_fltr,\
'R_et':R_et, 'thetain_et':thetain_et,\
'n_et':n_et, 'd_ch':d_ch, 'waveminEt':waveminEt, 'wavemaxEt':wavemaxEt,\
'A_Ch':A_Ch, 'nm_par':nm_par, 'nm_flat':nm_flat, 'nm_rms':nm_rms,\
'Etalons_filtern':Etalons_filtern,'Etalons_names':Etalons_names,\
'R_sphere':R_sphere, 'dist_coll_IS':dist_coll_IS, 'ap_input_port':ap_input_port,\
'ap_output_port':ap_output_port,'ro_is':ro_is, 'ISTransn':ISTransn, \
'ISFctr':ISFctr,'EP_dia':EP_dia, 'dist_PH_IS':dist_PH_IS,\
'T_target':T_target, 'epsilon_dif':epsilon_dif, 'TargetFctr':TargetFctr,\
'd_p_100':d_p_100, 'd_p_25':d_p_25,\
'MOSTransn':MOSTransn, 'MOSObsc':MOSObsc,'PssFctr':PssFctr,\
'T_pss':T_pss, 'epsilon_pss':epsilon_pss, 'pssapert':pssapert,\
'Ry_FM1':Ry_FM1, 'Ry_FM2':Ry_FM2, 'Ry_FM3':Ry_FM3, 'Ry_FM4':Ry_FM4,\
'Rx_FM4':Rx_FM4, 'd_prim_FM1':d_prim_FM1, 'd_FM1_FM2':d_FM1_FM2,\
'd_FM2_FM3':d_FM2_FM3, 'd_FM3_FM4':d_FM3_FM4,\
'T_FM1':T_FM1, 'T_FM2':T_FM2, 'T_FM3':T_FM3, 'T_FM4':T_FM4,\
'epsilon_FMs':epsilon_FMs, 'theta_FMs':theta_FMs, 'FMtransn':FMtransn, \
'FMFctr':FMFctr, 'd_MIRI_POM':d_MIRI_POM, 'Side_POM':Side_POM, 'T_POM':T_POM,\
'T_POM_M':T_POM_M, 'epsilon_POM':epsilon_POM, 'epsilon_POM_M':epsilon_POM_M,\
'theta_POM':theta_POM, 'theta_POM_M':theta_POM_M, 'POMFctr':POMFctr,\
'Side_OBA':Side_OBA, 'd_MIRI_OBA':d_MIRI_OBA, 'epsilon_OBA':epsilon_OBA,\
'T_OBA':T_OBA, 'theta_OBA':theta_OBA, 'OBAFctr': OBAFctr, 'FUDGEFctr': FUDGEFctr,\
'MIRIMfltrs':MIRIMfltrs, 'MIRIMcws':MIRIMcws, 'MIRIMpbs':MIRIMpbs,\
'CORONfltrs':CORONfltrs, 'CORONcws':CORONcws, 'CORONpbs':CORONpbs,\
'mrs_chans':mrs_chans, 'mrs_subch':mrs_subch,\
'mrs_cws_subch':mrs_cws_subch, 'mrs_pbs_subch':mrs_pbs_subch,\
'mrs_cws':mrs_cws, 'mrs_pbs':mrs_pbs, 'mrs_peakTr':mrs_peakTr,\
'mrs_pixrel':mrs_pixrel,'mrs_bolf':mrs_bolf,'mrs_dicf':mrs_dicf,\
'mrs_mirf':mrs_mirf,\
'dnphotrat':dnphotrat, 'FOV_Lineal':FOV_Lineal,\
'num_px':num_px, 'quantumf':quantumf}
