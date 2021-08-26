#main routine to calculate WCS positions from LRIS mask positions
#reads in file containing

#import packages
import sys
import numpy as np
from astropy.coordinates import SkyCoord, ICRS, Galactic, FK4, FK5
from astropy.table import Table,Column
from astropy.io import ascii,fits
from astropy.time import Time
import astropy.units as u
from testreverseautoslitcode import find_shifts as fs
from testreverseautoslitcode import refractioncorrection as ref
from testreverseautoslitcode import astrometrycorrection as ac
from testreverseautoslitcode import create_err_plot
from testreverseautoslitcode import write_ds9_file as w9f
from testreverseautoslitcode import precessionroutine as pr
from shutil import copyfile

def xytowcs(data_input_name,output_file):

    #set up some constants
    rpd = np.pi/180. #radians per degree
    rph = rpd*15     #radians per hour
    rpas = rpd/3600. #radians per arcsec
    mask_angle = 8.06 * rpd #angle of mask to the focal plane
    bend = 1.94 * rpd #angle of the bend in the mask
    x_center = 305. #point treated as the center on the ccd
  #  catalog_keyword = 'panstarrs'   #uncomment the desired catalog to be used
    catalog_keyword = 'gaia'
 #   catalog_keyword = 'custom'
    if catalog_keyword == 'custom':
        catalog_file = input('Catalog filename:')
    else:
        catalog_file = ''

#================================================================================================================================

#read in the data\
#copy the original file to a new file that will be the one to which the results are added
    copyfile(data_input_name, output_file+'.fits')
#Open the copied file
    hdu=fits.open(output_file+'.fits',mode='update')
    maskbluID = str(hdu[5].data['BluID'][0])

#obtain the rotation angle of the instrument

    rot_angle = hdu[2].data['PA_PNT'][0]

    print('Rotation Angle is selected at:',rot_angle)

#get the reference ra and dec
#assume this to be the center of the image

    ra0,dec0,equ = hdu[2].data['RA_PNT'][0], hdu[2].data['DEC_PNT'][0],str(hdu[2].data['EQUINPNT'][0])
    epoch = hdu[2].data['EQUINPNT'][0]

#================================================================================================================================
    #Correct the center of the mask for refraction
    #Taken directly from autoslit
    HA = 0   #hour angle
    wave = 6000./10000  #wavelength in microns
    lat = 19.828 #keck latitude
    temperature = 0
    pres = 486.  #atmos pressure
    ra0 = ra0*rpd
    dec0 = dec0*rpd
    H = HA*rph + ra0 - ra0
    cosz = np.sin(lat*rpd)*np.sin(dec0)+np.cos(lat*rpd)*np.cos(dec0)*np.cos(H)
    sinz = np.sqrt(1-(cosz**2))
    tanz = sinz/cosz
    if sinz != 0:
        sinQ = np.cos(lat*rpd)*np.sin(H)/sinz
        cosQ = (np.cos(np.pi/2-lat*rpd)-cosz*np.cos(np.pi/2-dec0))/(sinz*np.sin(np.pi/2-dec0))
    else:
        sinQ = 0
        cosQ = 0
    WAVERS = 1. / (wave * wave)
    N = 64.328 + (29498.1/(146.-WAVERS)) + (255.4/(41.-WAVERS))
    N = 1.*10**-6 * N
    TCORR = 1. + 0.003661 * temperature
    NUM = 720.88 * TCORR
    N = N * ((pres * (1.+(1.049-0.0157*temperature) * 1.E-6 * pres)) / NUM)

    R = N * 206265. * tanz

    DA1 = R * sinQ * rpas / np.cos(dec0)
    DD1 = R * cosQ * rpas

    #set the center of the mask to the corrected value
    ra0 =  ra0 +  DA1
    dec0 = dec0 + DD1
    ra0 = str(ra0/rpd)
    dec0 = str(dec0/rpd)

    #Grab the date the mask was made
    creation_date = hdu[2].data['DesDate'][0]

    #Set the mask scale. This depends on whether the mask was made to be used with or without the ADC (ADC installed for B semester of 2007)
    if Time(creation_date) > Time('2007-08-01'):
        scale = 0.7253 *0.99857#scale of mask in mm/arcsec, corrected for ADC in use
        adcuse = 'post'
        print('--------------ADC Correction Used-----------------')
    else:
        scale = 0.7253 #scale of mask in mm/arcsec, not corrected for ADC in use
        print('--------------NO ADC Correction Used-----------------')
        adcuse = 'pre'
    x0 = x_center / scale

    ref_system = str(np.chararray.lower(str(hdu[2].data['RADEPNT'][0])))

#convert reference ra,dec to decimal degrees
#correct for precession of coordiantes based on the equinox they are given in

#set up SkyCoord object for the reference RA,Dec
    temp = SkyCoord(ra0+' '+dec0,frame=ref_system,unit=(u.deg,u.deg),equinox='J'+equ)
    ra0 = temp.ra.deg
    dec0 = temp.dec.deg
    racenter = temp.ra.deg
    deccenter = temp.dec.deg
    print('Reference Coordinates set at: (',ra0,',',dec0,')')

#read in the slit data
    print('------------Reading in Slit Data-----------------')
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

#The data from the html mask files is in the milling machine coordinate system. We have to convert to the mask coordinate
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

    #apply the astrometry correction to the x,y values in the mask frame
    print('-------------------Applying Astrometry Correction-------------------')
    astro_table = ac.astrometry_calc(ra0,dec0)
    for i in range(len(data['X'])):
        distance = np.zeros(len(astro_table))
        for j in range(len(astro_table)):
            distance[j] = np.sqrt(((data['X'][i]-astro_table['Aparent X'][j]))**2+(data['Y'][i]-astro_table['Aparent Y'][j])**2)
        val = np.argmin(distance)
        data['X'][i]= data['X'][i]-astro_table['X Offset'][val]
        data['Y'][i]= data['Y'][i]-astro_table['Y Offset'][val]


#add columns to data table to eventually hold the calculated RA and Dec for each object
    data.add_column(0.00000000000,name='Calc_RA')
    data.add_column(0.00000000000,name='Calc_Dec')

#=====================================================================================================================================

#calculate the RA,Dec from mask coordinates

#take theta to negative theta (internal angle opposite to recorded position angle)
    theta = (-rot_angle) * rpd

    print('--------------------Completing Inverse Gnomic Projection--------------------')
#have to correct the center coordiantes before inverse gnomic projectoin
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

        data['Calc_RA'][i]=RA
        data['Calc_Dec'][i]=Dec

#===================================================================================================================================
    #Now apply refraction correction to each of the points on the mask
    print('-----------------Calculating Refraction Lookup Table------------------------')
    refraction_table = ref.refraction_calc(racenter*rpd,deccenter*rpd)
    #apply correction for each point
    for i in range(len(data['Calc_RA'])):
        distance = np.zeros(len(refraction_table))
        for j in range(len(refraction_table)):
            distance[j] = np.sqrt(((data['Calc_RA'][i]-refraction_table['Aparent RA'][j])*np.cos(data['Calc_Dec'][i]))**2+(data['Calc_Dec'][i]-refraction_table['Aparent Dec'][j])**2)
        val = np.argmin(distance)
        data['Calc_RA'][i]= data['Calc_RA'][i]-refraction_table['RA Offset'][val]
        data['Calc_Dec'][i]= data['Calc_Dec'][i]-refraction_table['Dec Offset'][val]

    for i in range(len(data['Calc_RA'])):
         temp1 = pr.precession(data['Calc_RA'][i],data['Calc_Dec'][i],epoch,2000.)
         data['Calc_RA'][i] = temp1[0]
         data['Calc_Dec'][i] = temp1[1]
#        temp = SkyCoord(str(data['Calc_RA'][i])+' '+str(data['Calc_Dec'][i]),frame=ref_system,unit=(u.deg,u.deg),equinox='J'+equ)
#        fk5_2000 = FK5(equinox='J2000')
    #take RA,Dec and transform to the J2000 frame
#        temp1=temp.transform_to(fk5_2000)
#        data['Calc_RA'][i] = temp1.ra.deg
#        data['Calc_Dec'][i] = temp1.dec.deg

#====================================================================================================================================
#calculate the average offset for the dataset and refraction correction
    print('----------------Calculating Mask Shift-----------------')
    data,x_centers,y_centers,ra_shifted_centers,dec_shifted_centers,catalog_obj_ra,catalog_obj_dec,objects_ra,objects_dec = fs.get_shift(data,theta,catalog_keyword,output_file,ref_system,racenter,deccenter,catalog_file,maskbluID,adcuse)

    #create the error plot
    print('-----------------Creating Quick Look Plot---------------')
    create_err_plot(data,catalog_obj_ra,catalog_obj_dec,objects_ra,objects_dec,output_file)
    w9f.create_ds9_file(data,ra_shifted_centers,dec_shifted_centers,rot_angle,catalog_obj_ra,catalog_obj_dec,output_file)
#=====================================================================================================================================
#update the fits file extension to include the calculated center positions of the slits
    hdu[3].data['slitRA'] = data['RA_Center'][::4]
    hdu[3].data['slitDec'] = data['Dec_Center'][::4]
    hdu.flush()

#write out the data and results to the output file
    ascii.write(data,output_file+'.csv',format='csv',overwrite=True)
