#read the file given in the command line and read 

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
import LRIS_Mask_Coords_to_WCS as lriscoords

#read in the input file from the command line
#reads in the file line by line
f = open(sys.argv[1],"r")
contents = f.readlines()
f.close()

#remove the "enter" from the end of lines if needed
for i in range(len(contents)):
    if contents[i][-1:] == '\n':
        contents[i] = contents[i][:-1]   

#split the command into the mask file name and the desired output file name
    names = contents[i].split(' ')
    mask_file= ''
    if len(names) > 2:
        for j in range(len(names)-1):
            mask_file = mask_file +' '+ names[j]
        mask_file = mask_file[1:]
        output_file = names[len(names)-1]
    else:
        mask_file = names[0]
        output_file = names[1]
    
    print('Opening file:',[i])
    print(output_file)
    lriscoords.xytowcs(mask_file,output_file)