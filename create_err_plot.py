#create the comparison plot between the calculated positions and the catalog positions
    
import numpy as np
from matplotlib import pyplot as plt
from astropy.table import Table,Column
from astropy.io import ascii,fits
from matplotlib import collections as mc
   
def create_err_plot(x_centers,y_centers,ra_centers,dec_centers,catalog_obj_ra,catalog_obj_dec,theta,output_file):
   ra_offsets = np.zeros(len(ra_centers))
   dec_offsets = np.zeros(len(ra_centers))

   for i in range(0,len(ra_centers)):
       ra_offsets[i] = (catalog_obj_ra[i]-ra_centers[i])*3600*np.cos(dec_centers[i]*np.pi/180)
       dec_offsets[i] = (catalog_obj_dec[i]-dec_centers[i])*3600
   x_offset = (np.cos(-theta)*ra_offsets-np.sin(-theta)*dec_offsets)*40/(0.7253 *0.99857)  # these offsets are scaled so they
   y_offset = (np.sin(-theta)*ra_offsets+np.cos(-theta)*dec_offsets)*40/(0.7253 *0.99857)  # are easily seen on the plot
   total_offset = np.sqrt(ra_offsets**2+dec_offsets**2)
   
   #create line segments from the actual to calculated positions
   lines = []
   for i in range(len(total_offset)):
       temp = [(x_centers[i],y_centers[i]),(x_centers[i]+x_offset[i],y_centers[i]+y_offset[i])]
       lines.append(temp)
   
   figtitle = input('Plot Title:')
   fig = plt.subplots(nrows=1,ncols=3,figsize=(20,8))
   plt.suptitle(figtitle)
   a1=plt.subplot(131)
   a1.title.set_text('Total Offsets in arcsec')
   plt.scatter(x_centers,y_centers,s=100,linestyle = '-',edgecolors='k',c=total_offset,cmap='Reds',vmin=0,vmax=1)
   lc = mc.LineCollection(lines,linewidths=2)
   a1.add_collection(lc)
   a1.set_xlabel('X Pos (mm)')
   a1.set_ylabel('Y Pos (mm)')
   a1.set_xlim([174,436])
   a1.set_ylim([-174,174])
   plt.colorbar()
   a2=plt.subplot(132)
   a2.title.set_text('RA Offsets in arcsec')
   plt.scatter(x_centers,y_centers,s=100,linestyle = '-',edgecolors='k',c=ra_offsets,cmap='bwr',vmin=-1,vmax=1)   
   plt.colorbar()
   a2.set_xlabel('X Pos (mm)')
   a2.set_ylabel('Y Pos (mm)')
   a2.set_xlim([174,436])
   a2.set_ylim([-174,174])

   a3=plt.subplot(133)
   a3.title.set_text('Dec Offsets in arcsec')
   plt.scatter(x_centers,y_centers,s=100,linestyle = '-',edgecolors='k',c=dec_offsets,cmap='bwr',vmin=-1,vmax=1)   
   plt.colorbar()
   a3.set_xlabel('X Pos (mm)')
   a3.set_ylabel('Y Pos (mm)')
   a3.set_xlim([174,436])
   a3.set_ylim([-174,174])

   plt.subplots_adjust(wspace=.2)
   plt.show()
   
   #save the figure
   figurefilename = output_file#input('Name for figure file:')
   plt.savefig('Figures/'+figurefilename+'.png')
   
   #print out the csv file
   xpositions = Column(x_centers, name='X')
   ypositions = Column(y_centers, name='Y')
   givenra = Column(catalog_obj_ra, name='Given_RA')
   givendec = Column(catalog_obj_dec, name='Given_Dec')
   calcra = Column(ra_centers, name='Calc_RA')
   calcdec = Column(dec_centers, name='Calc_Dec')
   xoffsets = Column(x_offset, name='X_Offset')
   yoffsets = Column(y_offset, name='Y_Offset')
   
   data_for_file = Table()
   data_for_file.add_column(xpositions)
   data_for_file.add_column(ypositions)
   data_for_file.add_column(givenra)
   data_for_file.add_column(givendec)
   data_for_file.add_column(calcra)
   data_for_file.add_column(calcdec)
   data_for_file.add_column(xoffsets)
   data_for_file.add_column(yoffsets)
  
   figuredatafilename = output_file+'comparison'#input('Name for figure creation file:')
   ascii.write(data_for_file,figuredatafilename+'.csv',format='csv',overwrite=True)
   
