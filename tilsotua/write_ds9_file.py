"""
Module to write DS9 file of mask slit positions.
"""

#read in packages
import numpy as np
from astropy.coordinates import SkyCoord, ICRS, Galactic, FK4, FK5
from astropy.table import Table,Column

def create_ds9_file(data,theta,scale,ref_system,output_file):
    """
    Function to generated the DS9 file. Creates the string for lines to be written to the DS9 file for each slit.

    Args:
        data: table
            table of mask data

        theta: float
            mask position angle

        ref_system: str
            mask sky coordinate reference system

        output_file: str
            outfile file base name

    Returns:
        None
    """
    data_copy = data['RA_Center','Dec_Center','X','Y']
    #sort the data table by y position for labeling purpoes
    avg_y_values = np.zeros(len(data_copy))
    i=0
    while i in range(len(data_copy)):
        avg_y_values[i:i+4] = np.average(data_copy['Y'][i:i+4]),np.average(data_copy['Y'][i:i+4]),np.average(data_copy['Y'][i:i+4]),np.average(data_copy['Y'][i:i+4])
        i = i+4
    data_copy.add_column(Column(avg_y_values),name='Avg Y')
    vertex_order = np.tile([1,2,3,4],int(len(data)/4))
    data_copy.add_column(Column(vertex_order),name='Vertex')
    data_copy.sort(['Avg Y','Vertex'])

#create a region file associated with the coordinates given in the data file
    str_data = Table()
    ra_str = Column(data_copy['RA_Center'], name='RA_Center', dtype=str)
    dec_str = Column(data_copy['Dec_Center'], name='Dec_Center', dtype=str)
    str_data.add_column(ra_str)
    str_data.add_column(dec_str)
    pa = str(theta)
#=================================================================================================================================
#calculate the height and width (or grab them from the data?)
    width1 = str(np.abs(np.max(data_copy['X'][0:4])-np.min(data_copy['X'][0:4]))/3600)
    height1 = str(np.abs(np.max(data_copy['Y'][0:4])-np.min(data_copy['Y'][0:4]))/3600)
#=================================================================================================================================
#write out boxes to the reg_string\
    reg_filename = output_file+'.reg'
    with open(reg_filename,mode='wt') as f:
        f.write(ref_system+'\n')
        i=0
        while i in range(0,int(len(str_data))):#-4):
            width = str(np.abs(np.max(data_copy['X'][i:i+4])-np.min(data_copy['X'][i:i+4]))/(scale*3600)) #covert from pixels to arcsec with platescale
            height = str(np.abs(np.max(data_copy['Y'][i:i+4])-np.min(data_copy['Y'][i:i+4]))/(scale*3600))
            reg_string='box('+str_data['RA_Center'][i]+','+str_data['Dec_Center'][i]+','+width+','+height+','+pa+') # color=green text ={'+str(int(i/4+1))+'}\n'
            f.write(reg_string)
            i = i+4
    f.close()
