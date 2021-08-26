## Reverse Autoslit Code

This code calculates the WCS positions of slits from slitmasks designed for the Low-Resolution Imaging Spectrograph (LRIS) on the Keck I telescope.

# Overview

This code takes data on the masks created for LRIS and stored in the LRIS archive and turn that data into useful information that allows . The code relies heavily on the routines developed for the autoslit3 code that is used to generate masks for LRIS.

The code takes input values for the locations of slits on a mask designed for LRIS in X,Y coordinates and reversing the major steps taken by the autoslit3 code. These steps include precessing the equinox and completing a reverse gnomic projection from X,Y coordinates to WCS coordinates.

The code then calculates a shift correction for the calculated slit positions based on the location of catalog objects.

These final WCS coordinates for the slits are recorded and output to a file, where they can then be used with the archive of LRIS observations.

# Using the Code

The code needs an input FITS file containing the mask information. That file must be in the format used by the LRIS archive. These FITS files have multiple extensions, each containing a specific table. The code uses this structure to grab the mask information. If the FITS files do not have the correct format, the code will not know where to grab the necessary information for its calculations.

The user gives the fits file name to the code and a name for the first CSV output file and the ds9 region file. The code then opens the file and reads in the center RA and Dec of the mask, the equinox of those coordinates, and the rotation angle.

The code then reads in the X,Y positions for the corners of each slit and completes a simple linear transformation between the mask milling machine coordinate system and the mask coordinate system (discussed in detail later) before completing the reverse gnomic projection.

With the WCS coordinates for the slits calculated in the equinox given by the mask at this point, the results are precessed to J2000 coordinates before being used for comparison against catalog objects. An average shift is calculated using the different between the raw slit positions of the guide star boxes and the catalog object positions selected at the targets meant to go in those boxes. This shift is applied to the raw slit positions and used as the final slit position.

The code then writes out the ds9 region file with the details of the slits and creates the error plot. The error plot shows the difference between the center of the box positions and the catalog object positions and contains a line for each point representing the direction and (not to scale) magnitude of the offset.

The user provides a file name for the figure file and the second CSV file that goes along with it.

# Parameters

The user is able to choose which catalog is queried, either PanSTARRS or GAIA, to calculate a shift. This is done with the catalog_keyword variable, which can be set to GAIA or Pan-STARRS.

Right now this variable is set at the beginning of the main routine, but in the future, it would be a good idea to make input parameters that can be selected by the user.

# Code Outputs

The code outputs four files: two CSV files, one PNG file, and one ds9 region file. These are named by the user.

The first CSV file contains the slit corner positions on the mask in X,Y coordinates and the calculated RA and Dec position. This means that for each box/slit, there are four rows in this file.

The region file creates polygon shapes using those four corners of the slits. These can be used in conjunction with an image of the field covered by the mask to see where the slits are calculated to lay and what objects the slits contain.

These files are created in the same directory the code is in.

The image of the error plot is saved and the second CSV file is written containing all the information necessary to reproduce the figure. This information includes the box center positions in X,Y, the WCS position of the box centers, the catalog object WCS position, and the offsets in X and Y directions.

The figure file is saved to a directory named ‘Figures’, while the CSV file is saved in the same directory as the code. (note: the code does not automatically create the ‘Figures’ directory at the moment.)

# Other Useful information

The code is made of three files. ‘LRIS_Mask_Coordiantes.py’ contains the main routine for the code. ‘Find_shifts.py’ contains the function to calculate the offsets between the initial box centers and the catalog objects. ‘Create_err_plot.py’ generates the diagnostic plots and the files that go with those plots.

For now these files need to be in the same directory.

# Limitations

There may also be fundamental limitations due to the accuracy of the input astrometry used to create a mask.
