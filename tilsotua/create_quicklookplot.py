#Create the quick look plot of the mask results and save as a pdf

import subprocess
import os
from astropy.io import fits
from astropy.wcs import WCS,utils
from astropy.table import Table,Column
from matplotlib import pyplot as plt
from astropy.coordinates import SkyCoord, ICRS, Galactic, FK4, FK5
import numpy as np
import astropy.units as u
from matplotlib.patches import Polygon as Pgon
from astropy.visualization import PercentileInterval, ImageNormalize
from glob import glob
from matplotlib.gridspec import GridSpec

#Function to create the quick look plots of the mask results
def create_qlp(data,catalog_obj_ra,catalog_obj_dec,objects_ra,objects_dec,output_file):
    data_copy = data['Calc_RA','Calc_Dec','X','Y']
    #Calculate the center of the mask to provide PanSTAMPS to download the image for the mask field  (should I just send over the actual center coords)
    x_center = np.average(data_copy['Calc_RA'])
    y_center = np.average(data_copy['Calc_Dec'])

    #Download the PanSTARRS field image for the mask location

    panstampscommand = 'panstamps --width=15 --filter=g stack '+str(x_center)+' '+str(y_center)
    process2 = subprocess.Popen(panstampscommand.split(), stdout=subprocess.PIPE)
    output, error = process2.communicate()
    process2.wait()

    #Open the field image. Just picks the first image on the list generated by glob
    field_file_path = str(x_center)+'*/''*.fits'
    all_field_files = glob(field_file_path)
    hdu = fits.open(all_field_files[0])[0]
#=================================================================================================================================
    #Set up the grid for the plot panels using GridSpec
    wcs = WCS(hdu.header)

    fig = plt.figure(constrained_layout=True,figsize=(7,4))
    gs = GridSpec(3, 3, figure=fig,width_ratios = [3,1,1],height_ratios=[1,1,1])
    ax1 = fig.add_subplot(gs[:, :-2],projection=wcs)
    ax2 = fig.add_subplot(gs[2, 1],projection=wcs)
    ax3 = fig.add_subplot(gs[2, 2],projection=wcs)
    ax4 = fig.add_subplot(gs[1, 1],projection=wcs)
    ax5 = fig.add_subplot(gs[1, 2],projection=wcs)
    ax6 = fig.add_subplot(gs[0, 1],projection=wcs)
    ax7 = fig.add_subplot(gs[0, 2],projection=wcs)


    #adjust tick label size on full field panel
    ax1.tick_params(axis='both', which='major', labelsize=8)

    #For clarity of the plot, set the ticks of the zoom in panels on the alignment boxes to be invisible
    lon2 = ax2.coords[0]
    lat2 = ax2.coords[1]
    lon2.set_ticks_visible(False)
    lat2.set_ticklabel_visible(False)
    lat2.set_ticks_visible(False)
    lon2.set_ticklabel_visible(False)
    lon3 = ax3.coords[0]
    lat3 = ax3.coords[1]
    lon3.set_ticks_visible(False)
    lat3.set_ticklabel_visible(False)
    lat3.set_ticks_visible(False)
    lon3.set_ticklabel_visible(False)
    lon4 = ax4.coords[0]
    lat4 = ax4.coords[1]
    lon4.set_ticks_visible(False)
    lat4.set_ticklabel_visible(False)
    lat4.set_ticks_visible(False)
    lon4.set_ticklabel_visible(False)
    lon5 = ax5.coords[0]
    lat5 = ax5.coords[1]
    lon5.set_ticks_visible(False)
    lat5.set_ticklabel_visible(False)
    lat5.set_ticks_visible(False)
    lon5.set_ticklabel_visible(False)
    lon6 = ax6.coords[0]
    lat6 = ax6.coords[1]
    lon6.set_ticks_visible(False)
    lat6.set_ticklabel_visible(False)
    lat6.set_ticks_visible(False)
    lon6.set_ticklabel_visible(False)
    lon7 = ax7.coords[0]
    lat7 = ax7.coords[1]
    lon7.set_ticks_visible(False)
    lat7.set_ticklabel_visible(False)
    lat7.set_ticks_visible(False)
    lon7.set_ticklabel_visible(False)

#=================================================================================================================================
    #Create the main plot. This is the full mask view.
    norm = ImageNormalize(hdu.data, interval=PercentileInterval(99))
    norm2 = ImageNormalize(hdu.data, interval=PercentileInterval(99))
    i=0
    avg_y_pos = []
    box_numbers = []

    #order the slits by y positions
    avg_y_values = np.zeros(len(data_copy))
    while i in range(len(data_copy)):
        avg_y_pos.append(np.average(data_copy['Y'][i:i+4]))
        avg_y_values[i:i+4] = np.average(data_copy['Y'][i:i+4]),np.average(data_copy['Y'][i:i+4]),np.average(data_copy['Y'][i:i+4]),np.average(data_copy['Y'][i:i+4])
        i = i+4

    avg_y_pos = np.sort(avg_y_pos)
    data_copy.add_column(Column(avg_y_values),name='Avg Y')
    vertex_order = np.tile([1,2,3,4],int(len(data_copy)/4))
    data_copy.add_column(Column(vertex_order),name='Vertex')
    data_copy.sort(['Avg Y','Vertex']) #trying to sort this to plot the boxes in order of the labels

    i=0
    while i in range(int(len(data_copy))):
        vertices = np.zeros(shape=(4,2))
        averages = np.zeros(shape=(2,2))
        vertices[0:4,0] = data_copy['Calc_RA'][i:i+4]
        vertices[0:4,1] = data_copy['Calc_Dec'][i:i+4]
        averages[0,0]=np.average((vertices[1,0],vertices[2,0]))
        averages[1,0] = np.average((vertices[0,0],vertices[3,0]))
        averages[0,1]=np.average((vertices[1,1],vertices[2,1]))
        averages[1,1]=np.average((vertices[0,1],vertices[3,1]))
        ra_center = np.average(vertices[:,0])
        dec_center = np.average(vertices[:,1])
        y_avg = np.average(data_copy['Y'][i:i+4])
        for k in range(len(avg_y_pos)):
            if y_avg == avg_y_pos[k]:
                ind = k+1
                ra1,ra2,ra3,ra4 = data_copy['Calc_RA'][i],data_copy['Calc_RA'][i+1],data_copy['Calc_RA'][i+2],data_copy['Calc_RA'][i+3]
                dec1,dec2,dec3,dec4 = data_copy['Calc_Dec'][i],data_copy['Calc_Dec'][i+1],data_copy['Calc_Dec'][i+2],data_copy['Calc_Dec'][i+3]
                side1 = 3600*np.sqrt(((ra2-ra1)*np.cos(dec1*np.pi/180))**2+(dec1-dec2)**2)
                side2 = 3600*np.sqrt(((ra2-ra3)*np.cos(dec2*np.pi/180))**2+(dec3-dec2)**2)
                ratio = side1/side2
                if 0.9 < ratio < 1.1:
                    box_numbers.append(ind)
            k=k+1
        pixelcoords = []
        for l in range(0,4):
            temp = SkyCoord(vertices[l,0],vertices[l,1], unit=(u.deg,u.deg),frame='fk5')
            temp2 = utils.skycoord_to_pixel(temp, wcs=wcs, origin=1, mode='all')
            pixelcoords.append(temp2)
        r = Pgon(pixelcoords, closed = True, edgecolor='green', facecolor='none',alpha=0.6,zorder=4)
        ax1.add_patch(r)
        ax1.text(ra_center,dec_center,str(ind),color='blue',transform=ax1.get_transform('fk5'))
        i = i+4


    im = ax1.imshow(hdu.data,origin='lower',zorder=0,cmap='gray_r',norm=norm)
    temp1a = SkyCoord(np.min(data_copy['Calc_RA'])-50/3600.,np.min(data_copy['Calc_Dec'])-50/3600.,unit=(u.deg,u.deg))
    temp1b = SkyCoord(np.max(data_copy['Calc_RA'])+50/3600.,np.max(data_copy['Calc_Dec'])+50/3600.,unit=(u.deg,u.deg))
    pixels1a=utils.skycoord_to_pixel(temp1a, wcs=wcs, origin=0, mode='all')
    pixels1b=utils.skycoord_to_pixel(temp1b, wcs=wcs, origin=0, mode='all')
    ax1.set_xlim([pixels1a[0],pixels1b[0]])
    ax1.set_ylim([pixels1a[1],pixels1b[1]])
    ax1.set_xlabel(r'$\alpha [J2000]$',fontsize='large')
    ax1.set_ylabel(r'$\delta [J2000]$',fontsize='large')
    ax1.invert_xaxis()
#=================================================================================================================================
    #Find the first six boxes and plot the cutout panels in the quick look plot.

    i=0
    boxes=0
    while i <= len(data_copy)-4 and boxes < 6:
        ra1,ra2,ra3,ra4 = data_copy['Calc_RA'][i],data_copy['Calc_RA'][i+1],data_copy['Calc_RA'][i+2],data_copy['Calc_RA'][i+3]
        dec1,dec2,dec3,dec4 = data_copy['Calc_Dec'][i],data_copy['Calc_Dec'][i+1],data_copy['Calc_Dec'][i+2],data_copy['Calc_Dec'][i+3]
        side1 = 3600*np.sqrt(((ra2-ra1)*np.cos(dec1*np.pi/180))**2+(dec1-dec2)**2)
        side2 = 3600*np.sqrt(((ra2-ra3)*np.cos(dec2*np.pi/180))**2+(dec3-dec2)**2)
        ratio = side1/side2
        if 0.9 < ratio < 1.1:
            boxes = boxes+1
            vertices = np.zeros(shape=(4,2))
            vertices[0:4,0] = data_copy['Calc_RA'][i:i+4]
            vertices[0:4,1] = data_copy['Calc_Dec'][i:i+4]
            ra_buffer = (5/(3600.*np.cos(vertices[0,1]*np.pi/180.)))
            dec_buffer = 5/3600.
            if boxes == 1:
                pixelcoords = []
                for l in range(len(vertices)):
                    temp = SkyCoord(vertices[l,0],vertices[l,1], unit=(u.deg,u.deg),frame='fk5')
                    temp2 = utils.skycoord_to_pixel(temp, wcs=wcs, origin=1, mode='all')
                    pixelcoords.append(temp2)
                r = Pgon(pixelcoords, closed = True, edgecolor='green', facecolor='none',alpha=0.6,zorder=4)
                ax2.add_patch(r)
                im2 = ax2.imshow(hdu.data,origin='lower',zorder=0,cmap='gray_r',norm=norm2)
                temp1a = SkyCoord(np.min(vertices[0:4,0])-ra_buffer,np.min(vertices[0:4,1])-dec_buffer,unit=(u.deg,u.deg))
                temp1b = SkyCoord(np.max(vertices[0:4,0])+ra_buffer,np.max(vertices[0:4,1])+dec_buffer,unit=(u.deg,u.deg))
                pixels1a=utils.skycoord_to_pixel(temp1a, wcs=wcs, origin=1, mode='all')
                pixels1b=utils.skycoord_to_pixel(temp1b, wcs=wcs, origin=1, mode='all')
                temp1c = SkyCoord(catalog_obj_ra,catalog_obj_dec,unit=(u.deg,u.deg))
                box_obj = np.argmin(temp1c.separation(temp1b).arcsec)
                pixels1c = utils.skycoord_to_pixel(temp1c[box_obj], wcs=wcs, origin=1, mode='all')
                ax2.scatter(pixels1c[0],pixels1c[1],zorder=5,color='yellow')#,transform=ax2.get_transform('fk5'))
                #ax2.scatter(catalog_obj_ra[boxes-1],catalog_obj_dec[boxes-1],zorder=5,color='yellow',transform=ax2.get_transform('fk5'))
                pixelcoordsra= []
                pixelcoordsdec= []
                for l in range(len(objects_ra[boxes-1])):
                    temp = SkyCoord(objects_ra[boxes-1],objects_dec[boxes-1], unit=(u.deg,u.deg),frame='fk5')
                    temp2 = utils.skycoord_to_pixel(temp, wcs=wcs, origin=1, mode='all')
                    pixelcoordsra.append(temp2[0])
                    pixelcoordsdec.append(temp2[1])
                ax2.scatter(pixelcoordsra,pixelcoordsdec,zorder=4,color='orange')
                ax2.text(np.min(vertices[0:4,0])-ra_buffer/3,np.min(vertices[0:4,1])-dec_buffer/3,str(box_numbers[boxes-1]),color='blue',transform=ax2.get_transform('fk5'))
                ax2.set_xlim([pixels1a[0],pixels1b[0]])
                ax2.set_ylim([pixels1a[1],pixels1b[1]])
                ax2.invert_xaxis()
            if boxes == 2:
                pixelcoords = []
                for l in range(len(vertices)):
                    temp = SkyCoord(vertices[l,0],vertices[l,1], unit=(u.deg,u.deg),frame='fk5')
                    temp2 = utils.skycoord_to_pixel(temp, wcs=wcs, origin=1, mode='all')
                    pixelcoords.append(temp2)
                r = Pgon(pixelcoords, closed = True, edgecolor='green', facecolor='none',alpha=0.6,zorder=4)
                ax3.add_patch(r)
                im3 = ax3.imshow(hdu.data,origin='lower',zorder=0,cmap='gray_r',norm=norm2)
                temp1a = SkyCoord(np.min(vertices[0:4,0])-ra_buffer,np.min(vertices[0:4,1])-dec_buffer,unit=(u.deg,u.deg))
                temp1b = SkyCoord(np.max(vertices[0:4,0])+ra_buffer,np.max(vertices[0:4,1])+dec_buffer,unit=(u.deg,u.deg))
                pixels1a=utils.skycoord_to_pixel(temp1a, wcs=wcs, origin=1, mode='all')
                pixels1b=utils.skycoord_to_pixel(temp1b, wcs=wcs, origin=1, mode='all')
                temp1c = SkyCoord(catalog_obj_ra,catalog_obj_dec,unit=(u.deg,u.deg))
                box_obj = np.argmin(temp1c.separation(temp1b).arcsec)
                pixels1c = utils.skycoord_to_pixel(temp1c[box_obj], wcs=wcs, origin=1, mode='all')
                ax3.scatter(pixels1c[0],pixels1c[1],zorder=5,color='yellow')
                #ax3.scatter(catalog_obj_ra[boxes-1],catalog_obj_dec[boxes-1],zorder=5,color='yellow',transform=ax3.get_transform('world'))
                pixelcoordsra= []
                pixelcoordsdec= []
                for l in range(len(objects_ra[boxes-1])):
                    temp = SkyCoord(objects_ra[boxes-1],objects_dec[boxes-1], unit=(u.deg,u.deg),frame='fk5')
                    temp2 = utils.skycoord_to_pixel(temp, wcs=wcs, origin=1, mode='all')
                    pixelcoordsra.append(temp2[0])
                    pixelcoordsdec.append(temp2[1])
                ax3.scatter(pixelcoordsra,pixelcoordsdec,zorder=4,color='orange')
                #ax3.scatter(objects_ra[boxes-1],objects_dec[boxes-1],zorder=4,color='purple',transform=ax3.get_transform('fk5'))
                ax3.text(np.min(vertices[0:4,0])-ra_buffer/3,np.min(vertices[0:4,1])-dec_buffer/3,str(box_numbers[boxes-1]),color='blue',transform=ax3.get_transform('fk5'))
                ax3.set_xlim([pixels1a[0],pixels1b[0]])
                ax3.set_ylim([pixels1a[1],pixels1b[1]])
                ax3.invert_xaxis()
            if boxes == 3:
                pixelcoords = []
                for l in range(len(vertices)):
                    temp = SkyCoord(vertices[l,0],vertices[l,1], unit=(u.deg,u.deg),frame='fk5')
                    temp2 = utils.skycoord_to_pixel(temp, wcs=wcs, origin=1, mode='all')
                    pixelcoords.append(temp2)
                r = Pgon(pixelcoords, closed = True, edgecolor='green', facecolor='none',alpha=0.6,zorder=4)
                ax4.add_patch(r)
                im4 = ax4.imshow(hdu.data,origin='lower',zorder=0,cmap='gray_r',norm=norm2)
                temp1a = SkyCoord(np.min(vertices[0:4,0])-ra_buffer,np.min(vertices[0:4,1])-dec_buffer,unit=(u.deg,u.deg))
                temp1b = SkyCoord(np.max(vertices[0:4,0])+ra_buffer,np.max(vertices[0:4,1])+dec_buffer,unit=(u.deg,u.deg))
                pixels1a=utils.skycoord_to_pixel(temp1a, wcs=wcs, origin=1, mode='all')
                pixels1b=utils.skycoord_to_pixel(temp1b, wcs=wcs, origin=1, mode='all')
                temp1c = SkyCoord(catalog_obj_ra,catalog_obj_dec,unit=(u.deg,u.deg))
                box_obj = np.argmin(temp1c.separation(temp1b).arcsec)
                pixels1c = utils.skycoord_to_pixel(temp1c[box_obj], wcs=wcs, origin=1, mode='all')
                ax4.scatter(pixels1c[0],pixels1c[1],zorder=5,color='yellow')
                #ax4.scatter(catalog_obj_ra[boxes-1],catalog_obj_dec[boxes-1],zorder=5,color='yellow',transform=ax4.get_transform('world'))
                pixelcoordsra= []
                pixelcoordsdec= []
                for l in range(len(objects_ra[boxes-1])):
                    temp = SkyCoord(objects_ra[boxes-1],objects_dec[boxes-1], unit=(u.deg,u.deg),frame='fk5')
                    temp2 = utils.skycoord_to_pixel(temp, wcs=wcs, origin=1, mode='all')
                    pixelcoordsra.append(temp2[0])
                    pixelcoordsdec.append(temp2[1])
                ax4.scatter(pixelcoordsra,pixelcoordsdec,zorder=4,color='orange')
              #  ax4.scatter(objects_ra[boxes-1],objects_dec[boxes-1],zorder=4,color='purple',transform=ax4.get_transform('fk5'))
                ax4.text(np.min(vertices[0:4,0])-ra_buffer/3,np.min(vertices[0:4,1])-dec_buffer/3,str(box_numbers[boxes-1]),color='blue',transform=ax4.get_transform('fk5'))
                ax4.set_xlim([pixels1a[0],pixels1b[0]])
                ax4.set_ylim([pixels1a[1],pixels1b[1]])
                ax4.invert_xaxis()
            if boxes == 4:
                pixelcoords = []
                for l in range(len(vertices)):
                    temp = SkyCoord(vertices[l,0],vertices[l,1], unit=(u.deg,u.deg),frame='fk5')
                    temp2 = utils.skycoord_to_pixel(temp, wcs=wcs, origin=1, mode='all')
                    pixelcoords.append(temp2)
                r = Pgon(pixelcoords, closed = True, edgecolor='green', facecolor='none',alpha=0.6,zorder=4)
                ax5.add_patch(r)
                im5 = ax5.imshow(hdu.data,origin='lower',zorder=0,cmap='gray_r',norm=norm2)
                temp1a = SkyCoord(np.min(vertices[0:4,0])-ra_buffer,np.min(vertices[0:4,1])-dec_buffer,unit=(u.deg,u.deg))
                temp1b = SkyCoord(np.max(vertices[0:4,0])+ra_buffer,np.max(vertices[0:4,1])+dec_buffer,unit=(u.deg,u.deg))
                pixels1a=utils.skycoord_to_pixel(temp1a, wcs=wcs, origin=1, mode='all')
                pixels1b=utils.skycoord_to_pixel(temp1b, wcs=wcs, origin=1, mode='all')
                temp1c = SkyCoord(catalog_obj_ra,catalog_obj_dec,unit=(u.deg,u.deg))
                box_obj = np.argmin(temp1c.separation(temp1b).arcsec)
                pixels1c = utils.skycoord_to_pixel(temp1c[box_obj], wcs=wcs, origin=1, mode='all')
                ax5.scatter(pixels1c[0],pixels1c[1],zorder=5,color='yellow')
                #ax5.scatter(catalog_obj_ra[boxes-1],catalog_obj_dec[boxes-1],zorder=5,color='yellow',transform=ax5.get_transform('world'))
                pixelcoordsra= []
                pixelcoordsdec= []
                for l in range(len(objects_ra[boxes-1])):
                    temp = SkyCoord(objects_ra[boxes-1],objects_dec[boxes-1], unit=(u.deg,u.deg),frame='fk5')
                    temp2 = utils.skycoord_to_pixel(temp, wcs=wcs, origin=1, mode='all')
                    pixelcoordsra.append(temp2[0])
                    pixelcoordsdec.append(temp2[1])
                ax5.scatter(pixelcoordsra,pixelcoordsdec,zorder=4,color='orange')
              #  ax5.scatter(objects_ra[boxes-1],objects_dec[boxes-1],zorder=4,color='purple',transform=ax5.get_transform('fk5'))
                ax5.text(np.min(vertices[0:4,0])-ra_buffer/3,np.min(vertices[0:4,1])-dec_buffer/3,str(box_numbers[boxes-1]),color='blue',transform=ax5.get_transform('fk5'))
                ax5.set_xlim([pixels1a[0],pixels1b[0]])
                ax5.set_ylim([pixels1a[1],pixels1b[1]])
                ax5.invert_xaxis()
            if boxes == 5:
                pixelcoords = []
                for l in range(len(vertices)):
                    temp = SkyCoord(vertices[l,0],vertices[l,1], unit=(u.deg,u.deg),frame='fk5')
                    temp2 = utils.skycoord_to_pixel(temp, wcs=wcs, origin=1, mode='all')
                    pixelcoords.append(temp2)
                r = Pgon(pixelcoords, closed = True, edgecolor='green', facecolor='none',alpha=0.6,zorder=4)
                ax6.add_patch(r)
                im6 = ax6.imshow(hdu.data,origin='lower',zorder=0,cmap='gray_r',norm=norm2)
                temp1a = SkyCoord(np.min(vertices[0:4,0])-ra_buffer,np.min(vertices[0:4,1])-dec_buffer,unit=(u.deg,u.deg))
                temp1b = SkyCoord(np.max(vertices[0:4,0])+ra_buffer,np.max(vertices[0:4,1])+dec_buffer,unit=(u.deg,u.deg))
                pixels1a=utils.skycoord_to_pixel(temp1a, wcs=wcs, origin=1, mode='all')
                pixels1b=utils.skycoord_to_pixel(temp1b, wcs=wcs, origin=1, mode='all')
                temp1c = SkyCoord(catalog_obj_ra,catalog_obj_dec,unit=(u.deg,u.deg))
                box_obj = np.argmin(temp1c.separation(temp1b).arcsec)
                pixels1c = utils.skycoord_to_pixel(temp1c[box_obj], wcs=wcs, origin=1, mode='all')
                ax6.scatter(pixels1c[0],pixels1c[1],zorder=5,color='yellow')
                #ax6.scatter(catalog_obj_ra[boxes-1],catalog_obj_dec[boxes-1],zorder=5,color='yellow',transform=ax6.get_transform('world'))
                pixelcoordsra= []
                pixelcoordsdec= []
                for l in range(len(objects_ra[boxes-1])):
                    temp = SkyCoord(objects_ra[boxes-1],objects_dec[boxes-1], unit=(u.deg,u.deg),frame='fk5')
                    temp2 = utils.skycoord_to_pixel(temp, wcs=wcs, origin=1, mode='all')
                    pixelcoordsra.append(temp2[0])
                    pixelcoordsdec.append(temp2[1])
                ax6.scatter(pixelcoordsra,pixelcoordsdec,zorder=4,color='orange')
            #    ax6.scatter(objects_ra[boxes-1],objects_dec[boxes-1],zorder=4,color='purple',transform=ax6.get_transform('fk5'))
                ax6.text(np.min(vertices[0:4,0])-ra_buffer/3,np.min(vertices[0:4,1])-dec_buffer/3,str(box_numbers[boxes-1]),color='blue',transform=ax6.get_transform('fk5'))
                ax6.set_xlim([pixels1a[0],pixels1b[0]])
                ax6.set_ylim([pixels1a[1],pixels1b[1]])
                ax6.invert_xaxis()
            if boxes == 6:
                pixelcoords = []
                for l in range(len(vertices)):
                    temp = SkyCoord(vertices[l,0],vertices[l,1], unit=(u.deg,u.deg),frame='fk5')
                    temp2 = utils.skycoord_to_pixel(temp, wcs=wcs, origin=1, mode='all')
                    pixelcoords.append(temp2)
                r = Pgon(pixelcoords, closed = True, edgecolor='green', facecolor='none',alpha=0.6,zorder=4)
                ax7.add_patch(r)
                im7 = ax7.imshow(hdu.data,origin='lower',zorder=0,cmap='gray_r',norm=norm2)
                temp1a = SkyCoord(np.min(vertices[0:4,0])-ra_buffer,np.min(vertices[0:4,1])-dec_buffer,unit=(u.deg,u.deg))
                temp1b = SkyCoord(np.max(vertices[0:4,0])+ra_buffer,np.max(vertices[0:4,1])+dec_buffer,unit=(u.deg,u.deg))
                pixels1a=utils.skycoord_to_pixel(temp1a, wcs=wcs, origin=1, mode='all')
                pixels1b=utils.skycoord_to_pixel(temp1b, wcs=wcs, origin=1, mode='all')
                temp1c = SkyCoord(catalog_obj_ra,catalog_obj_dec,unit=(u.deg,u.deg))
                box_obj = np.argmin(temp1c.separation(temp1b).arcsec)
                pixels1c = utils.skycoord_to_pixel(temp1c[box_obj], wcs=wcs, origin=1, mode='all')
                ax7.scatter(pixels1c[0],pixels1c[1],zorder=5,color='yellow')
                #ax7.scatter(catalog_obj_ra[boxes-1],catalog_obj_dec[boxes-1],zorder=5,color='yellow',transform=ax7.get_transform('world'))
                pixelcoordsra= []
                pixelcoordsdec= []
                for l in range(len(objects_ra[boxes-1])):
                    temp = SkyCoord(objects_ra[boxes-1],objects_dec[boxes-1], unit=(u.deg,u.deg),frame='fk5')
                    temp2 = utils.skycoord_to_pixel(temp, wcs=wcs, origin=1, mode='all')
                    pixelcoordsra.append(temp2[0])
                    pixelcoordsdec.append(temp2[1])
                ax7.scatter(pixelcoordsra,pixelcoordsdec,zorder=4,color='orange')
           #     ax7.scatter(objects_ra[boxes-1],objects_dec[boxes-1],zorder=4,color='purple',transform=ax7.get_transform('fk5'))
                ax7.text(np.min(vertices[0:4,0])-ra_buffer/3,np.min(vertices[0:4,1])-dec_buffer/3,str(box_numbers[boxes-1]),color='blue',transform=ax7.get_transform('fk5'))
                ax7.set_xlim([pixels1a[0],pixels1b[0]])
                ax7.set_ylim([pixels1a[1],pixels1b[1]])
                ax7.invert_xaxis()
        i = i+4
    plt.savefig(output_file+'quicklookplot.pdf',bbox_inches='tight',pad_inches=0.5)
#=================================================================================================================================
    #Delete the PanSTARRS image and directory
    #Comment this section out if you prefer to save the PanSTARRS image
    for i in range(len(all_field_files)):
        os.remove(all_field_files[i])
    if y_center < 0:
        os.rmdir(str(x_center)+'m'+str(y_center)[1:]+'/')
    if y_center > 0:
        os.rmdir(str(x_center)+'p'+str(y_center)+'/')
