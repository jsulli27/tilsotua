#function for creating the ds9 region files

#read in packages
import numpy as np
from astropy.coordinates import SkyCoord, ICRS, Galactic, FK4, FK5
from astropy.table import Table,Column

def create_ds9_file(data,ra_centers,dec_centers,theta,scale,ref_system,catalog_object_ra,catalog_object_dec,output_file):
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
#write out boxes to the reg_string\
    reg_filename = output_file+'.reg'
    with open(reg_filename,mode='wt') as f:
        f.write(ref_system+'\n')
        i=0
        while i in range(0,int(len(str_data))):#-4):
            width = str(np.abs(np.max(data['X'][i:i+4])-np.min(data['X'][i:i+4]))/(scale*3600)) #covert from pixels to arcsec with platescale
            height = str(np.abs(np.max(data['Y'][i:i+4])-np.min(data['Y'][i:i+4]))/(scale*3600))
            reg_string='box('+str_data['RA_Center'][i]+','+str_data['Dec_Center'][i]+','+width+','+height+','+pa+') # color=green\n'
            f.write(reg_string)
            i = i+4
    f.close()
