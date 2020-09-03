# LRIS-X-Y-to-WCS

Code to calculate the RA and DEC positions of slits for LRIS masks based on the mask positions of those slits. The code then applies a shift to the calculated
positions based on querying either the PanSTARRS or GAIA catalog and comparing the postions of the alignment boxes to the objects that they should be centered on.
The following files are generated: a ds9 region file is generated which will allow the positions of the slits to be overlayed on top of an image of the 
corresponding field, a plot showing the errors between the centers of the alignment boxes for the mask and the objects picked out in either the PanSTARRS or GAIA 
catalogs that the box should contain, and two .csv files where one contains the results and the other contains all the informtation needed to reproduce the plot.

The code takes an input file from the command line with the names of the mask fits files as well as the desired output file names. From there, the code will run on
its own and will only require the manual input of a plot name.

All the code files should be in the same directory.

'callcode.py' is the file to call from the command line and passes the names of the files to the rest of the code.

'LRIS_Mask_Coords_to_WCS.py' contains the main routine for the code. The main calculation from LRIS mask position to WCS position is based on reversing the
calculations done by autoslit, the program used to create LRIS masks.

'find_shifts.py' queries the catalog of choice and calculates the shift to be applied to the slit positions based on the difference between the centers of the
alisgnment boxes and the catalog object determined to be the target of the box.

'create_err_plot.py' creates the plot and calculates the x and y offsets between the centers of the boxes and the catalog objects.

'createds9regions.py' creates the ds9 region file.
