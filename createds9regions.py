#create region files for ds9

from regions import DS9Parser, write_ds9
from astropy.table import Table,Column
from astropy.io import ascii
import numpy as np

#import the data
data = Table.read('outputmask1c_referencestarsonly.csv')
coords = Table()
ra = Column(data['Calc_RA'], name='Calc_RA', dtype=str)
dec = Column(data['Calc_Dec'], name='Calc_Dec', dtype=str)
coords.add_column(ra)
coords.add_column(dec)
#regs = list(np.zeros(len(coords)))
#convert to string? in the format the package needs
#reg_string = 'icrs\ncircle(35,-9,.003) # color=green'
#parser= DS9Parser(reg_string)
#regs = parser.shapes.to_regions()
'''
reg_string ='icrs\npolygon('+coords['Calc_RA'][0]+','+coords['Calc_Dec'][0]+','+coords['Calc_RA'][4]+','+coords['Calc_Dec'][4]+','+coords['Calc_RA'][8]+','+coords['Calc_Dec'][8]+','+coords['Calc_RA'][12]+','+coords['Calc_Dec'][12]+') # color=blue'
parser= DS9Parser(reg_string)
regs = parser.shapes.to_regions()
reg_list = regs
i=1
while i in range(1,int(len(coords)/4)):
    reg_string='icrs\npolygon('+coords['Calc_RA'][i]+','+coords['Calc_Dec'][i]+','+coords['Calc_RA'][i+4]+','+coords['Calc_Dec'][i+4]+','+coords['Calc_RA'][i+8]+','+coords['Calc_Dec'][i+8]+','+coords['Calc_RA'][i+12]+','+coords['Calc_Dec'][i+12]+') # color=blue'
    parser= DS9Parser(reg_string)
    regs = parser.shapes.to_regions()
    reg_list = reg_list + regs
    i = i+1
'''


reg_string='icrs\ncircle('+coords['Calc_RA'][0]+','+coords['Calc_Dec'][0]+',.001) # color=blue'
parser= DS9Parser(reg_string)
regs = parser.shapes.to_regions()
reg_list = regs
for i in range(1,len(coords)):
    reg_string='icrs\ncircle('+coords['Calc_RA'][i]+','+coords['Calc_Dec'][i]+',.001) # color=blue'
    parser= DS9Parser(reg_string)
    regs = parser.shapes.to_regions()
    reg_list = reg_list + regs
#convert the strings to regions?
#print(type(reg_list))


filename = 'mask1c_referencestars.reg'
write_ds9(reg_list, filename,radunit='deg')

