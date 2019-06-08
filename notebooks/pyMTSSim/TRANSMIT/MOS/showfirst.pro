
PRO showfirst
    ;
    fileA = 'MIRI_mirror_MIR_A.txt'
    fileB = 'MIRI_mirror_MIR_B.txt'
    
    oldMOSf = 'MOStransmit_oldcorr.txt'
    newMOSf = 'MOScurve_capture_shift_blue0.06.txt'
    
    READCOL,fileA,wnA,TA,format='(d,d)'
    READCOL,fileB,wnB,TB,format='(d,d)'
    
    waveA = 1.E4 / wnA
    waveB = 1.E4 / wnB
    
    READCOL,oldMOSf,waveO,TO,format='(d,d)'
    TO /= 100.
    
    READCOL,newMOSF,waveN,TN,format='(d,d)'
    TN /= 100.
    
    basic_colors, black ,white,red,green,blue,$
             yellow,cyan,magenta, orange,mint, $
             purple,pink,olive,lightb,gray,lgray,/LOAD

    window,0,retain=2,/COLOR
    
    plot,waveA,(TA-0.01)^4.,linestyle=0,color=0,xrange=[7.,9.5],$
    yrange=[0.8,1.1],xtitle='microns',ytitle='T',title='MOS REFLECTIVITY',$
    charsize=1.5,background=1,/NODATA
    
    oplot,waveA,(TA-0.01)^4.,linestyle=0,color=0,thick=2

    oplot,(waveO)/1.37,TO^4.,linestyle=2,color=3,thick=2
    oplot,waveO-3.0,TO^4.,linestyle=2,color=2,thick=2
    oplot,waveN,TN^4.,linestyle=0,color=4,thick=2
    
    ;psf = psf_gaussian(Npixel=500.,$
    ;	    FWHM=150.,/double,/normalize,NDIMEN=1)
	    
    ;convTN = CONVOLVE(TN,psf)   
    
    ;oplot,waveN-0.06,convTN^4.,linestyle=2,color=4,thick=2
    
    legend = ['RAL (T-0.01)^4',$
    'Old MTSSim, T^4, lambda / 1.37',$
    'Old MTSSim, T^4, lambda - 3.0 um',$
    'Manufacturer, T^4']
    
    AL_LEGEND_MTSSim,legend,linestyle=[0,2,2,0],thick=[2,2,2,2],$
    color=[0,3,2,4],/top,/left,textcolors=[0,0,0,0],outline=0,charsize=1.5
    
END
