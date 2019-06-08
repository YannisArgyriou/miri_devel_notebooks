
PRO MOScurve
    ; "Aim: to provide a correct transmission for the MOS"
    ; 
    fileA = 'MIRI_mirror_MIR_A.txt'
    fileB = 'MIRI_mirror_MIR_B.txt'
    
    oldMOSf = 'MOStransmit_oldcorr.txt'
    
    READCOL,fileA,wnA,TA,format='(d,d)'
    READCOL,fileB,wnB,TB,format='(d,d)'
    
    waveA = 1.E4 / wnA
    waveB = 1.E4 / wnB
    
    READCOL,oldMOSf,waveO,TO,format='(d,d)'
    TO /= 100.
    
    basic_colors, black ,white,red,green,blue,$
             yellow,cyan,magenta, orange,mint, $
             purple,pink,olive,lightb,gray,lgray,/LOAD

    window,0,retain=2,/COLOR
    
    stop
    
END
