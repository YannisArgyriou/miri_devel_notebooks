
PRO MOSresample
    
    curvef = 'MOScurve_capture_raw.txt'
    shift =  0.0 ; -0.06 ; um
    
        
    shifttag = '_shift_'
    IF shift LT 0 THEN shifttag += 'blue' 
    IF shift GT 0 THEN shifttag += 'red'
    shifttag += STRTRIM(STRING(ABS(shift),format='(f6.2)'),2)
    
    outf = 'MOScurve_capture'+shifttag+'.txt'
    
    READCOL,curvef,wave,T,format='(d,d)'
    
    T *= 100.
    minw = 2.0
    maxw = 30.
    dw = 0.005
    nw = (maxw-minw)/dw + 1
    newwave = findgen(nw)*dw + minw
    
    iT = INTERPOL(T,wave+shift,newwave)
    
    curve = {wave:newwave,T:iT}
    dumpfile,curve,outf

END
