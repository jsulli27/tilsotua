#code to query GAIA database for closest object to calculated coordinates

import numpy as np
from astropy.coordinates import SkyCoord, ICRS, Galactic, FK4, FK5
import astropy.units as u
from astroquery.gaia import Gaia
from astroquery.mast import Catalogs
from matplotlib import pyplot as plt
import create_err_plot as cep
import write_ds9_file as w9f

def get_shift(data,theta,catalog_keyword,output_file,ref_system):
    ra_shifts=[]
    dec_shifts = []
    ra_centers = []
    dec_centers = []
    x_centers = []
    y_centers = []
    ra_shift1_centers = []
    dec_shift1_centers = []
    ra_shift2_centers = []
    dec_shift2_centers = []
    catalog_obj_ra = []
    catalog_obj_dec = []

    i=0
    while i <= len(data)-4:
        #load in the corners for a given slit
        ra1,ra2,ra3,ra4 = data['Calc_RA'][i],data['Calc_RA'][i+1],data['Calc_RA'][i+2],data['Calc_RA'][i+3]
        dec1,dec2,dec3,dec4 = data['Calc_Dec'][i],data['Calc_Dec'][i+1],data['Calc_Dec'][i+2],data['Calc_Dec'][i+3]
        side1 = 3600*np.sqrt(((ra2-ra1)*np.cos(dec1*np.pi/180))**2+(dec1-dec2)**2)
        side2 = 3600*np.sqrt(((ra2-ra3)*np.cos(dec2*np.pi/180))**2+(dec3-dec2)**2)
        ratio = side1/side2
        print(side1,side2,ratio)
        if 0.9 < ratio < 1.1: #3.5 < side1 < 5.0 and 3.5 < side2 < 5.0:
            ra_arr=[ra1,ra2,ra3,ra4]
            dec_arr=[dec1,dec2,dec3,dec4]
            ra_avg = (np.max(ra_arr)+np.min(ra_arr))/2.
            dec_avg = (np.max(dec_arr)+np.min(dec_arr))/2.
            xpos = [data['X'][i],data['X'][i+1],data['X'][i+2],data['X'][i+3]]
            ypos = [data['Y'][i],data['Y'][i+1],data['Y'][i+2],data['Y'][i+3]]
            x_avg = (np.max(xpos)+np.min(xpos))/2.
            y_avg = (np.max(ypos)+np.min(ypos))/2.
            ra_centers.append(ra_avg)
            dec_centers.append(dec_avg)
            x_centers.append(x_avg)
            y_centers.append(y_avg)
            print('ra,dec of box =',ra_avg,dec_avg)
        
            #call astroquery
            
            coord = SkyCoord(ra=ra_avg, dec=dec_avg, unit=(u.deg,u.deg),frame=ref_system)
            if catalog_keyword == 'gaia':
                width = 200*u.arcsec
                height = 200*u.arcsec
                obj = Gaia.query_object_async(coordinate=coord, width=width,height=height)
                
                #find the object with the closest position
                diff = []
                for j in range(len(obj)):
                    diff.append(np.sqrt(((obj['ra'][j]-ra_avg)*np.cos(dec_avg))**2+(obj['dec'][j]-dec_avg)**2))
                min_diff = np.argmin(diff)
                
                #save catalog object to array
                catalog_obj_ra.append(obj[min_diff]['ra'])
                catalog_obj_dec.append(obj[min_diff]['dec'])
            
                print('GAIA object located at:',obj[min_diff]['ra'],obj[min_diff]['dec'])
                ra_shifts.append((ra_avg-obj['ra'][min_diff])*np.cos(dec_avg*np.pi/180.))
                dec_shifts.append(dec_avg-obj['dec'][min_diff])
            
            if catalog_keyword == 'panstarrs':
                obj = Catalogs.query_region(coord, radius=15*u.arcsec, catalog='Panstarrs', data_release='dr2')
        
                #find the object with the closest position
                diff = []
                for j in range(len(obj)):
                    diff.append(np.sqrt(((obj['raMean'][j]-ra_avg)*np.cos(dec_avg))**2+(obj['decMean'][j]-dec_avg)**2))
                min_diff = np.argmin(diff)

            
                #save catalog object to array
                catalog_obj_ra.append(obj[min_diff]['raMean'])
                catalog_obj_dec.append(obj[min_diff]['decMean'])
            
                print('PANSTARRS object located at:',obj[min_diff]['raMean'],obj[min_diff]['decMean'])
                ra_shifts.append((ra_avg-obj['raMean'][min_diff])*np.cos(dec_avg*np.pi/180.))
                dec_shifts.append(dec_avg-obj['decMean'][min_diff])
            i=i+4
        else:
            i=i+4
    j=0
    good_ra_shifts=[]
    good_dec_shifts=[]
    for i in range(len(ra_shifts)):
        if np.sqrt((ra_shifts[i]*np.cos(dec_shifts[i]*np.pi/180))**2+(dec_shifts[i]**2))*3600 < 4.0:
           good_ra_shifts.append(ra_shifts[i])
           good_dec_shifts.append(dec_shifts[i])
        else:
            j=j+1
    print(j,'elements removed from shift arrays')
        
    dec_shift_final = np.average(good_dec_shifts)
    ra_shift_final = np.average(good_ra_shifts)/np.cos(data['Calc_Dec'][0]*np.pi/180)
    ra_shift1_centers = ra_centers+ra_shift_final
    dec_shift1_centers = dec_centers+dec_shift_final
    ra_shift2_centers = ra_centers-ra_shift_final
    dec_shift2_centers = dec_centers-dec_shift_final
    print(ra_shifts)
    print('Shifted1 Positions:')
    print(ra_shift1_centers,dec_shift1_centers)
    print('Shifted2 Positions:')
    print(ra_shift2_centers,dec_shift2_centers)
        
    #create the error plot
    cep.create_err_plot(x_centers,y_centers,ra_shift2_centers,dec_shift2_centers,catalog_obj_ra,catalog_obj_dec,theta,output_file)
    w9f.create_ds9_file(data,ra_shift2_centers,dec_shift2_centers,catalog_obj_ra,catalog_obj_dec,output_file)
    
    return (-ra_shift_final ,-dec_shift_final)