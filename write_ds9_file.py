#function for creating the ds9 region files

#read in packages
import numpy as np
from astropy.coordinates import SkyCoord, ICRS, Galactic, FK4, FK5
from astropy.table import Table,Column
from regions import DS9Parser, Regions

def create_ds9_file(data,ra_centers,dec_centers,theta,catalog_object_ra,catalog_object_dec,output_file):
#create a region file associated with the coordinates given in the data file
    str_data = Table()
    ra_str = Column(data['RA_Center'], name='RA_Center', dtype=str)
    dec_str = Column(data['Dec_Center'], name='Dec_Center', dtype=str)
    str_data.add_column(ra_str)
    str_data.add_column(dec_str)
    pa = str(theta)
#=================================================================================================================================
#calculate the height and width (or grab them from the data?)
    width1 = str(np.abs(np.max(data['X'][0:4])-np.min(data['X'][0:4]))/3600)
    height1 = str(np.abs(np.max(data['Y'][0:4])-np.min(data['Y'][0:4]))/3600)
#=================================================================================================================================
#write out boxes to the reg_string
    reg_list ='icrs\nbox('+str_data['RA_Center'][0]+','+str_data['Dec_Center'][0]+','+width1+','+height1+','+pa+') # color=red'

    i=4
    while i in range(1,int(len(str_data))):#-4):
        width = str(np.abs(np.max(data['X'][i:i+4])-np.min(data['X'][i:i+4]))/(0.7253 *0.99857*3600))
        height = str(np.abs(np.max(data['Y'][i:i+4])-np.min(data['Y'][i:i+4]))/(0.7253 *0.99857*3600))
        reg_string='nbox('+str_data['RA_Center'][i]+','+str_data['Dec_Center'][i]+','+width+','+height+','+pa+') # color=red'
        reg_list = reg_list + reg_string
        i = i+4
#=================================================================================================================================

    #combine the two sections so the slits and vectors are both included
    reg_list = Regions.parse(reg_list,format='ds9')
    reg_filename = output_file+'.reg'
    reg_list.write(reg_filename,format='ds9',overwrite=True)
