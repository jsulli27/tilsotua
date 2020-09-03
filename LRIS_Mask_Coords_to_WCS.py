#code to calculate WCS positions from LRIS mask positions
#reads in file containing 

#import packages and set up some constants
import fileinput
import sys
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord, ICRS, Galactic, FK4, FK5
from astropy.table import Table,Column
from astropy.io import ascii,fits
from astropy.time import Time
from regions import DS9Parser, write_ds9
import astropy.units as u
import find_shifts as fs

def xytowcs(data_input_name,output_file):
        
    rpd = np.pi/180. #radians per degree
    rpas = rpd/3600. #radians per arcsec
    gap_angle = 0.09 * rpd #angle between the CCD chips in radians
    mask_angle = 8.06 * rpd #angle of mask to the focal plane
    bend = 1.94 * rpd #angle of the bend in the mask
    x_center = 305. #point treated as the center on the ccd
    #catalog_keyword = 'panstarrs'
    catalog_keyword = 'gaia'

#================================================================================================================================

#read in the data

#data_input_name = input('Input Data Filename:')

    hdu=fits.open(data_input_name)

#get the name for the output files

#    data_output_name = input('Names for Output Files:')

#obtain the rotation angle of the instrument

#rot_angle = float(input('LRIS Rotation Angle (in deg):'))
    rot_angle = hdu[2].data['PA_PNT'][0]
    
    print('Rotation Angle is selected at:',rot_angle)

#get the reference ra and dec
#assume this to be the center of the image

    ra0,dec0,equ = str(hdu[2].data['RA_PNT'][0]), str(hdu[2].data['DEC_PNT'][0]),str(hdu[2].data['EQUINPNT'][0])
    creation_date = hdu[2].data['DesDate'][0]
    
    if Time(creation_date) > Time('2007-08-01'):
        scale = 0.7253 *0.99857#scale of mask in mm/arcsec, corrected for ADC in use
        print('--------------ADC Correction Used-----------------')
    else:
        scale = 0.7253 #scale of mask in mm/arcsec, not corrected for ADC in use
        print('--------------NO ADC Correction Used-----------------')
    x0 = x_center / scale

#ra0 = input('Reference RA (in hms):')
#dec0 = input('Reference Dec (in dms):')
#equ = input('Equinox of Reference Coordinates:')

    ref_system = str(np.chararray.lower(str(hdu[2].data['RADEPNT'][0])))

#convert reference ra,dec to decimal degrees
#correct for precession of coordiantes based on the equinox they are given in

#set up SkyCoord object for the reference RA,Dec
    temp = SkyCoord(ra0+' '+dec0,frame=ref_system,unit=(u.deg,u.deg),equinox='J'+equ)
#set up J2000 coordinate frame
#fk5_2000 = FK5(equinox='J2000')
#take RA,Dec and transform to the J2000 frame
    #temp1=temp.transform_to(fk5_2000)
    ra0 = temp.ra.deg
    dec0 = temp.dec.deg
    print('Reference Coordinates set at: (',ra0,',',dec0,')')

#read in the slit data
    data = Table()
    x1col = Column(hdu[6].data['slitX1'][0:1], name='slitX1', dtype=float)
    y1col = Column(hdu[6].data['slitY1'][0:1], name='slitY1', dtype=float)
    data.add_column(x1col)
    data.add_column(y1col)
    data.add_row([hdu[6].data['slitX2'][0],hdu[6].data['slitY2'][0]])
    data.add_row([hdu[6].data['slitX3'][0],hdu[6].data['slitY3'][0]])
    data.add_row([hdu[6].data['slitX4'][0],hdu[6].data['slitY4'][0]])
    
    for i in range(1,len(hdu[6].data['slitX1'])):
        data.add_row([hdu[6].data['slitX1'][i],hdu[6].data['slitY1'][i]])
        data.add_row([hdu[6].data['slitX2'][i],hdu[6].data['slitY2'][i]])
        data.add_row([hdu[6].data['slitX3'][i],hdu[6].data['slitY3'][i]])
        data.add_row([hdu[6].data['slitX4'][i],hdu[6].data['slitY4'][i]])
    
    data.rename_column('slitX1','X')
    data.rename_column('slitY1','Y')
    
    for i in range(len(data['X'])):
        data['X'][i] = np.float(data['X'][i])
        data['Y'][i] = np.float(data['Y'][i])
        
#because the data from the html mask files is in the milling machine coordinate system, we have to convert to the mask coordinate
#system before working with the data
#x_mask = y_mill + 172.7
#y+mask = -x_mill +177.8

    x_mill = data['X']
    y_mill = data['Y']
    
    
    for i in range(len(data['X'])):
        a = y_mill[i]+172.7
        b = -x_mill[i]+177.8
        data['X'][i] = a
        data['Y'][i] = b
        
    x_input = data['X']
    y_input = data['Y']
    
#open the datafile and read in the data

    '''
    #if the file is an .html file we need to grab the data using pandas
    if data_input_name[-4:] == 'html':
        
        # This creates a list of dataframes:
        frame_list = pd.read_html(data_input_name)
        
        # Extract the dataframes for the mask design and desiSlits tables:
        maskDesign = frame_list[1]
        desiSlits = frame_list[5]
        bluSlits = frame_list[8]
    
        # By default, the table and column names are in the first two rows.
        # This assigns the appropriate column names and removes the junk.
        maskDesign = maskDesign.drop(index=[0,1])
    
        desiSlits.columns = desiSlits.iloc[1,:]
        desiSlits = desiSlits.drop(index=[0,1])
    
        bluSlits.columns = bluSlits.iloc[1,:]
        bluSlits = bluSlits.drop(index=[0,1])
        
        desiSlits = Table.from_pandas(desiSlits)
        bluSlits = Table.from_pandas(bluSlits)
        
        #take data and put it into the right format for my code
    
        data = Table()
        x1col = Column(bluSlits['slitX1'][0:1], name='slitX1', dtype=float)
        y1col = Column(bluSlits['slitY1'][0:1], name='slitY1', dtype=float)
        data.add_column(x1col)
        data.add_column(y1col)
        data.add_row([bluSlits['slitX2'][0],bluSlits['slitY2'][0]])
        data.add_row([bluSlits['slitX3'][0],bluSlits['slitY3'][0]])
        data.add_row([bluSlits['slitX4'][0],bluSlits['slitY4'][0]])
    
        for i in range(1,len(bluSlits)):
            data.add_row([bluSlits['slitX1'][i],bluSlits['slitY1'][i]])
            data.add_row([bluSlits['slitX2'][i],bluSlits['slitY2'][i]])
            data.add_row([bluSlits['slitX3'][i],bluSlits['slitY3'][i]])
            data.add_row([bluSlits['slitX4'][i],bluSlits['slitY4'][i]])
    
        data.rename_column('slitX1','X')
        data.rename_column('slitY1','Y')
    
    
        for i in range(len(data['X'])):
            data['X'][i] = np.float(data['X'][i])
            data['Y'][i] = np.float(data['Y'][i])
            
    
        #because the data from the html mask files is in the milling machine coordinate system, we have to convert to the mask coordinate
        #system before working with the data
        #x_mask = y_mill + 172.7
        #y+mask = -x_mill +177.8
    
        x_mill = data['X']
        y_mill = data['Y']
    
    
        for i in range(len(data['X'])):
            a = y_mill[i]+172.7
            b = -x_mill[i]+177.8
            data['X'][i] = a
            data['Y'][i] = b
        
    else:
        data = Table.read(data_input_name)
       
    x_input = data['X']
    y_input = data['Y']
    '''

#add columns to eventually hold the calculated RA and Dec for each object
    data.add_column(0.00000000000,name='Calc_RA')
    data.add_column(0.00000000000,name='Calc_Dec')
    print(data)

#=====================================================================================================================================

#calculate the RA,Dec from mask coordinates
    
#take theta to negative theta
    theta = (-rot_angle) * rpd
    
#have to correct the center coordiantes
    ra0 = ra0*rpd - (x0*np.cos(theta)*rpas/np.cos(dec0*rpd))
    dec0 = dec0*rpd - x0*np.sin(theta)*rpas
        
    for i in range(len(data['X'])):
    
#correct the x and y coords for the bend and tilt in the mask
        x_prime = np.cos(mask_angle)*(x_input[i])/scale
        y_prime = np.cos(bend)*(y_input[i])/scale
        
#take x,y to the eta and nu gnomic projection coordinates (rotation matrix)
        eta = (np.cos(theta)*x_prime-np.sin(theta)*y_prime)*rpas
        nu =  (np.sin(theta)*x_prime+np.cos(theta)*y_prime)*rpas
#take eta and nu to ra and dec (standard already calculated inversion)
        rho = np.sqrt(eta**2+nu**2)
        c = np.arctan2(rho,1.)
    
#calculate the final ra and dec using the inverse gnomic projection
        ra_t = eta*np.sin(c)
        ra_b = (rho*np.cos(dec0)*np.cos(c))-(nu*np.sin(dec0)*np.sin(c))
        dec_f = np.cos(c)*np.sin(dec0)
        dec_s = (nu*np.sin(c)*np.cos(dec0))/rho
        RA = (ra0+np.arctan2(ra_t/ra_b,1.))/rpd
        Dec = np.arcsin(dec_f+dec_s)/rpd
        
        temp = SkyCoord(str(RA)+' '+str(Dec),frame=ref_system,unit=(u.deg,u.deg),equinox='J'+equ)
        fk5_2000 = FK5(equinox='J2000')
    #take RA,Dec and transform to the J2000 frame
        temp1=temp.transform_to(fk5_2000)
        RA = temp1.ra.deg
        Dec = temp1.dec.deg
    
        data['Calc_RA'][i]=RA
        data['Calc_Dec'][i]=Dec
        
#====================================================================================================================================
#calculate the average offset for the dataset
    offset = fs.get_shift(data,theta,catalog_keyword,output_file,ref_system)
    print('offset is:',offset)
    
    data['Calc_RA']= data['Calc_RA']+offset[0]
    data['Calc_Dec']= data['Calc_Dec']+offset[1]
#=====================================================================================================================================


#write out the data and results to the output file
    ascii.write(data,output_file+'.csv',format='csv',overwrite=True)