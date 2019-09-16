# REQUIRES rlim[0]=rlim[1]
# rlim is 3-element vector

# detx and dety are the actual detector x,y locations.  These are not
# needed by this code and are only for debugging, as are stopx and
# stopy which will stop at a particular cube x,y

import numpy as np
import matplotlib.pyplot as plt

def mmrs_cube(x, y, z, f, expnum, dim_out, rlim, scale=None, xsquash=None, ysquash=None, zsquash=None, maskcube=None,slice=None, wtype=None, expsig=None, detx=None, dety=None, detlam=None, stopx=None, stopy=None):

    # wtype:
    # 1: 1/d weighting
    # 2: 1/d^2 weighting
    # 3: gaussian weighting (define expsig)
    if wtype is None: wtype=2

    if xsquash is None: xsquash=1.
    if ysquash is None: ysquash=1.
    if zsquash is None: zsquash=1.
    if scale is None: scale=1.

    # Add some defaults for non-necessary arguments so code doesn't crash without
    if detx is None: detx=np.ones(1,len(x))
    if dety is None: dety=np.ones(1,len(y))
    if detlam is None: detlam=np.ones(len(z))
    if stopx is None: stopx=-1
    if stopy is None: stopy=-1

    # If the slice keyword is set, make only a single slice
    # and override the nominal output dimensions (ensure we
    # don't pass this backwards)
    thisdim_out=dim_out
    if slice is not None:
        if ((slice >= 0) & (slice < thisdim_out[2])):
            thisdim_out[2]=1
            arr_zcoord=[slice]
    else:
        arr_zcoord=None

    # Output cubes
    fcube=np.zeros(thisdim_out)
    maskcube=np.ones(thisdim_out)

    # XYZ output pixel coordinate arrays
    arr_xcoord = np.arange(thisdim_out[0])
    arr_ycoord = np.arange(thisdim_out[1])
    # Don't overwrite arr_zcoord it if was already set to one slice
    if arr_zcoord is None:
        arr_zcoord = np.arange(thisdim_out[2])

    # Loop over output cube building the cube slice by slice
    for k in range(thisdim_out[2]):
        # Print a message every 5% of the loop
        if (k % thisdim_out[2]/20 == 0):
            print('Constructing cube: {}/{} complete'.format(k,thisdim_out[2]))

        # First pass cut: trim to only stuff within rlim of this z location
        indexk=np.where(abs(z-arr_zcoord[k]-0.5) <= rlim[2])
        nindexk = len(indexk[0])

        # If nothing makes the cut,: do nothing.  Otherwise build the slice
        if (nindexk > 0):
            tempx=x[indexk]
            tempy=y[indexk]
            tempz=z[indexk]
            tempf=f[indexk]
            tempenum=expnum[indexk]
            temp_detx=detx[indexk]
            temp_dety=dety[indexk]
            temp_detlam=detlam[indexk]

            # QA plot
            # plt.figure()
            # plt.plot(tempx,tempy,'bo')
            # circxcen=np.median(tempx)
            # circycen=np.median(tempy)
            # circphi=np.arange(360.)*np.pi/180.
            # circx=circxcen+rlim[0]*np.cos(circphi)
            # circy=circycen+rlim[1]*np.sin(circphi)
            # plt.plot(circx,circy,'r')
            # plt.show()

            # Loop over output image, building the image row by row
            for j in range(thisdim_out[1]):
                # Second pass cut: trim to only stuff within rlim of this y location
                indexj=np.where(abs(tempy-arr_ycoord[j]) <= rlim[1])
                nindexj = len(indexj[0])

                # If nothing makes the cut, do nothing.  Otherwise
                # build the row
                if (nindexj > 0):
                    tempx2=tempx[indexj]
                    tempy2=tempy[indexj]
                    tempz2=tempz[indexj]
                    tempf2=tempf[indexj]
                    tempenum2=tempenum[indexj]
                    temp2_detx=temp_detx[indexj]
                    temp2_dety=temp_dety[indexj]
                    temp2_detlam=temp_detlam[indexj]

                    # plt.figure()
                    # plt.plot(tempx2,tempy2,'bo')
                    # circxcen=np.median(tempx2)
                    # circycen=np.median(tempy2)
                    # circphi=np.arange(360.)*np.pi/180.
                    # circx=circxcen+rlim[0]*np.cos(circphi)
                    # circy=circycen+rlim[1]*np.sin(circphi)
                    # plt.plot(circx,circy,'r')
                    # plt.show()

                    # Now do a 1d build within this slice, looping over input points
                    arr_weights=np.zeros((nindexj,thisdim_out[0]))

                    for q in range(nindexj):
                        arr_radius=(rlim[0]+1)*np.ones(thisdim_out[0])
                        arr_sradius=(rlim[0]+1)*np.ones(thisdim_out[0])

                        # Which output pixels are affected by input points, i.e.
                        # within rlim of this x location?
                        # Don't go outside output array boundaries
                        # sel = np.where((np.floor(tempx2-rlim[0]) > 0) & (np.ceil(tempx2+rlim[0]) < thisdim_out[0]-1))

                        # Number of points within box
                        nbox=len(arr_xcoord)

                        # Calculate physical spatial radius for ROI determination
                        rx=arr_xcoord-tempx2[q]
                        ry=arr_ycoord[j]-tempy2[q]
                        arr_radius=np.sqrt(rx**2 + np.ones(nbox)*ry**2)

                        # Determine points within the final circular ROI
                        tocalc=np.where(arr_radius <= rlim[0])
                        ncalc = len(tocalc[0])

                        # Squashed radii for weights
                        srx=rx/xsquash
                        sry=ry/ysquash
                        srz=(arr_zcoord[k]-tempz2[q]+0.5)/zsquash


                        # Combine normalized radii inside ROI
                        arr_sradius = np.sqrt( (srx**2) + np.ones(nbox)* sry**2 + np.ones(nbox)* srz**2  )

                        # Ensure no divide by zero
                        if (ncalc > 0):
                            if (wtype == 0):
                                arr_weights[tocalc+q*thisdim_out[0]]=1.
                            elif (wtype == 1):
                                arr_weights[tocalc+q*thisdim_out[0]]=1./arr_sradius[tocalc]
                            elif (wtype == 2):
                                arr_weights[q,tocalc]=1./arr_sradius[tocalc]**2
                            elif (wtype == 3):
                                arr_weights[tocalc+q*thisdim_out[0]]=np.exp(-0.5/expsig**2*arr_sradius[tocalc]**2)
                            elif (wtype == 4):
                                arr_weights[tocalc+q*thisdim_out[0]]=1./arr_sradius[tocalc]**4

                    # Normalization matrix
                    if (nindexj == 1):
                        matr_norm=arr_weights
                    else:
                        matr_norm = np.sum(arr_weights,0)

                    # Flag where the normalization matrix is zero; there is no good datum here
                    nodata=np.where(matr_norm == 0)
                    nnodata = len(nodata[0])
                    gooddata = np.where(matr_norm != 0)
                    ngood = len(gooddata[0])
                    # Mark good locations in output mask
                    if (ngood != 0):
                        maskcube[gooddata,j,k]=0

                    # We don't want to divide by zero where there is no data# set the normalization
                    # matrix to 1 in these cases
                    if (nnodata > 0):
                        matr_norm[nodata]=1.

                    # Apply the weights to calculate the output flux in this row
                    frow=np.zeros(thisdim_out[0])
                    for q in range (nindexj):
                        alpha=arr_weights[q,:]/matr_norm
                        frow+=tempf2[q]*alpha

                    # Put the row into the final cube
                    fcube[:,j,k]=frow*scale

            if ((j == stopy) & (slice is not None)):
                temp=arr_weights[stopx,:] # Cull the array weights for this x pixel
                thispix=np.where(temp != 0.)[0] # Identify where weights nonzero
                nthis = len(thispix)
                if (nthis > 0):
                    thispix_detx=temp2_detx[thispix]
                    thispix_dety=temp2_dety[thispix]
                    thispix_detlam=temp2_detlam[thispix]
                    thispix_dx=tempx2[thispix]-stopx
                    thispix_dy=tempy2[thispix]-stopy
                    thispix_dz=tempz2[thispix]-slice-0.5
                    thispix_enum=tempenum2[thispix]
                    thispix_flux=tempf2[thispix]
                    print('HARDCODED 0.1 arcsec and 0.002 micron spaxels')
                    print('Distances are in arcsec and microns')
                    print('Debug location is: {},{},{}'.format(stopx,stopy,slice))
                    print('Final value is: {}'.format(fcube[stopx,j,k]))
                    print('exp xdet ydet wave xdist ydist zdist rxy flux rweight nweight')
                    for r in range(nthis):
                        # NB- HARDCODING pixel size conversions to 0.1/0.1 arcsec and 0.002 micron
                        print(thispix_enum[r],thispix_detx[r],thispix_dety[r],thispix_detlam[r],thispix_dx[r]*0.1,thispix_dy[r]*0.1,thispix_dz[r]*0.002,np.sqrt(thispix_dx[r]**2+thispix_dy[r]**2)*0.1,thispix_flux[r],temp[thispix[r]],temp[thispix[r]]/matr_norm[stopx])
                else:
                  print('No non-zero weight contributing pixels!')
    print('{}/{} complete'.format(thisdim_out[2],thisdim_out[2]))

    return fcube
