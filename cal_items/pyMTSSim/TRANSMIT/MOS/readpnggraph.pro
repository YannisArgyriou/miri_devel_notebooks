
PRO readpnggraph
    
    png = 'MOSfeature_8um.png'
    
    img = READ_PNG(png)
    img = img[0,*,*]
    
    ;window,0,retain=2,xsize=1100,ysize=400
    tvframe,img
    
    ; Read mesh
    
    corners = {lowl:[0.,0.,0.,0.],lowr:[0.,0.,0.,0.],$
    upl:[0.,0.,0.,0.],upr:[0.,0.,0.,0.]}
    cornersnames = TAG_NAMES(corners)
    ncorners = n_elements(cornersnames)
    
    FOR i=0,ncorners-1 DO BEGIN
	
       name = cornersnames[i]
       PRINT,'mark corner '+name
       CURSOR,x,y,/DOWN,/DATA
       corners.(i)[0] = x
       corners.(i)[1] = y
       READ,name+'-x = ',xg
       READ,name+'-y = ',yg
       corners.(i)[2] = xg
       corners.(i)[3] = yg

   ENDFOR
    
    ; Linear Transform (no rotation!!)
    ; x = xlowl + (xlowr-xlowl)/(xplowr-xplowl) * xpx
    ; y = ylowl + (yupl-ylowl)/(ypxupl-ypxlowl) * ypx
    
    bx = (corners.lowr[2]-corners.lowl[2]) / (corners.lowr[0]-corners.lowl[0])
    ax = corners.lowl[2] - corners.lowl[0] * bx
    
    by = (corners.upl[3]-corners.lowl[3]) / (corners.upl[1]-corners.lowl[1])
    ay = corners.lowl[3] - corners.lowl[1] * by
    
    
    xpx = [-1.]
    ypx = [-1.]
    goon = 1
    
    PRINT,'Follow the curve you want to retrieve'
    
    WHILE goon EQ 1 DO BEGIN
        
	CURSOR,ixpx,iypx,/DOWN,/DATA
	xpx = [xpx,[ixpx]]
	ypx = [ypx,[iypx]]
	print,ixpx,iypx,ixpx * bx + ax, iypx * by + ay
	IF iypx LT 0. THEN goon = 0
    ENDWHILE
    
    xpx = xpx[1:n_elements(xpx)-2]
    ypx = ypx[1:n_elements(ypx)-2]    
    
    x = ax + bx * xpx
    y = ay + by * ypx
    
    wave = 1.E4 / x
    
    PLOT,wave,y
    curve = {wave:wave,T:y}
    curvef = 'MOSfeature_8um.txt'
    dumpfile,curve,curvef
    
    
    stop

END
