#Function to read in data from the autoslit files

from astropy.table import Table,vstack
from shutil import copyfile
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.time import Time
import pathlib

def gen_from_auto(autofile):
    """
    Generate a FITS file similar to the UCO/Lick archival masks from the autoslit
    output file.
    Args:
    autofile (str): Path to autoslit output file fed to AUTOSLIT

    Returns:
    None

    Note that the mask design date is a placeholder meant to trigger whether the
    ADC was used for caluclations in the main code.
    """

    dir_loc = str(pathlib.Path(__file__).parent.resolve())

    #read in the file
    auto_file = open(autofile+'.file3')

    content = auto_file.readlines()

    #need the center of the mask, epoch, position angle
    angle = content[26].split()[4]
    epoch = content[28].split()[-2]
    ra0 = content[58].split()[3]+content[58].split()[4]+content[58].split()[5]
    dec0 = content[59].split()[3]+content[59].split()[4]+content[59].split()[5]+content[59].split()[6]
    coords = SkyCoord(ra0,dec0,frame='fk5',equinox='J'+epoch)
    ra0 = coords.ra.deg
    dec0 = coords.dec.deg
    adc_use = content[41].split()[-2]

    #read in the X_mask,Y_mask slit positions
    def lines_that_equal(line_to_match, fp):
        #search for the lines in a file that
        matches = []
        for i in range(len(fp)):
            line = fp[i]
            if line == line_to_match:
                matches.append(i)
        return(matches)

    slit_start = lines_that_equal('newrow\n',content)
    slit_pos = Table(names=['slitX1','slitY1','slitX2','slitY2','slitX3','slitY3','slitX4','slitY4'], dtype=(float,float,float,float,float,float,float,float))
    for i in range(len(slit_start)-1):
        ind = slit_start[i]

        X,Y = content[ind+1].split()[0],content[ind+1].split()[1]
        X2,Y2 = content[ind+2].split()[0],content[ind+2].split()[1]
        X3,Y3 = content[ind+3].split()[0],content[ind+3].split()[1]
        X4,Y4 = content[ind+4].split()[0],content[ind+4].split()[1]

        slit_pos.add_row([X,Y,X2,Y2,X3,Y3,X4,Y4])

    #write the new fits file with the same structure as the UCO/Lick archive
    copyfile(dir_loc+'/blank_template.fits', autofile+'.fits')
    new_file = fits.open(autofile+'.fits',mode='update')

    #write the design information to the fits file
    mask_design = Table(new_file['MaskDesign'].data)
    mask_design.add_row()
    mask_design['RA_PNT'] = ra0
    mask_design['DEC_PNT'] = dec0
    mask_design['PA_PNT'] = angle
    mask_design['EQUINPNT'] = epoch
    mask_design['RADEPNT'] = 'fk5'
    if adc_use == 'TRUE':
        mask_design['DesDate'] = '2022-01-01'
    else:
        mask_design['DesDate'] = '2000-01-01'

    new_file['MaskDesign'] = fits.BinTableHDU(mask_design,header = new_file['MaskDesign'].header)

    #add the slit positions
    bluslits = Table(new_file['BluSlits'].data)
    bluslits = vstack([bluslits,slit_pos])
    new_file['BluSlits'] = fits.BinTableHDU(bluslits,header = new_file['BluSlits'].header)

    new_file.flush()
