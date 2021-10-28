#calculate the final correction autoslit makes to the x,y positions for the milling maching output file.
#This doesn't depend on WCS position, so I should just calculate a single table that the code reads in

import numpy as np
from astropy.table import Table,Column
def astrometry_calc(ra_center,dec_center):

    #create a grid of values to create the points for which the correction will be calculated
    x_range = np.linspace(174,436, num = 500)
    y_range = np.linspace(-174,174, num = 500)

    #set up arrays to hold results later
    x_values = []
    y_values = []
    x_offsets = []
    y_offsets = []
    offsets = []
#=================================================================================================================================
    #create the routine for the astrometry correction
    def astrometry_correction_calc(XIN,YIN):
        #set some constants
        CCD_SIDE = 2048.0
        CCD_SCALE = 0.15767
        X_CENTER = 305.0  #center of the mask in x
        Y_CENTER = 0.0    #center of the mask in y
        #set the fits values from autoslit that are used in the correction
        A=[2.0,  0.99476227,0.00728125,0.45412825E-6, -0.40955557E-5,0.96690643E-6]
        B=[1.0, -0.00001856, 0.99906324 ]

        #Convert mm to CCD pixels. (Approximate, but as long as you
        #use exactly the same formula for inversion, that is OK.)
        XCCD = CCD_SIDE / 2.0 + (X_CENTER - XIN) / CCD_SCALE
        YCCD = CCD_SIDE / 2.0 - (Y_CENTER - YIN) / CCD_SCALE

        XCCD_OUT = A[0] + A[1] * XCCD + A[2] * YCCD + A[3] * XCCD * XCCD + A[4] * YCCD * YCCD + A[5] * XCCD * YCCD
        YCCD_OUT = B[0] + B[1] * XCCD + B[2] * YCCD

        #Now convert back to mm on the slitmask.

        XOUT = X_CENTER - (XCCD_OUT - CCD_SIDE / 2.0) * CCD_SCALE
        YOUT = Y_CENTER + (YCCD_OUT - CCD_SIDE / 2.0) * CCD_SCALE

        return(XOUT,YOUT)
#=================================================================================================================================
    #set up more arrays to hold results
    corrected_x = np.zeros(shape=(len(x_range),len(y_range)))
    corrected_y = np.zeros(shape=(len(x_range),len(y_range)))

    #run through the values of the grid and put each point through the slit_astrometry calculation
    #appending results to the arrays from earlier in prep for creating an astropy table to return the results in
    for i in range(len(x_range)):
        for j in range(len(y_range)):
            temp = astrometry_correction_calc(x_range[i],y_range[j])
            corrected_x = temp[0]
            x_values.append(corrected_x)
            corrected_y= temp[1]
            y_values.append(corrected_y)
            x_offset = (corrected_x-x_range[i])
            x_offsets.append(x_offset)
            y_offset = (corrected_y-y_range[j])
            y_offsets.append(y_offset)
            offset = np.sqrt(x_offset**2+y_offset**2)
            offsets.append(offset)
#=================================================================================================================================            
    #create the astropy table and fill it in
    astrometry_results= Table()
    x_aparent = Column(x_values, name='Aparent X')
    y_aparent = Column(y_values, name='Aparent Y')
    xoffsets = Column(x_offsets, name='X Offset')
    yoffsets = Column(y_offsets, name='Y Offset')

    astrometry_results.add_column(x_aparent)
    astrometry_results.add_column(y_aparent)
    astrometry_results.add_column(xoffsets)
    astrometry_results.add_column(yoffsets)

    #return the astropy results table
    return(astrometry_results)
