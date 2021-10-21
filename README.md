## Reverse Autoslit Code

This code calculates the WCS positions of slits from slitmasks designed for the Low-Resolution Imaging Spectrograph (LRIS) on the Keck I telescope.

# Overview

This code takes archival FITS slitmask files for LRIS containing the positions of slits in those masks in the coordinate frame of the slitmask milling machine and reverses the steps taken by the [autoslit3](https://www2.keck.hawaii.edu/inst/lris/autoslit_WMKO.html) code (originally used to create the slitmasks) to obtain the WCS coordinates of the slits. A final astrometric correction (shift) is then applied to bring the mask results into agreement with the GAIA frame.

These final WCS coordinates for the slits are recorded and output to a file, where they can then be used with the archive of LRIS observations.

# Using the Code

The standard way to run the code on a given mask FITS file is:

from testreverseautoslitcode import xytowcs

xytowcs(maskfilename,outputfilename)

The code needs an input FITS file containing the mask information. That file must be in the format used by the UCO/Lick archive. These FITS files have multiple extensions, each containing a specific table. The code uses this structure to grab the mask information. If the FITS files do not have the correct format, the code will not know where to grab the necessary information for its calculations.

The given output filename will be used as the base for all four output files produced by the code.

# Note on Astrometric (Shift) Correction

The code defaults to using the GAIA catalog to calculate a correction to bring the slit positions into agreement with modern astrometry. There code does allow for use of the PanSTARRS catalog or a custom catalog of objects if the GAIA catalog


# Code Outputs

The code outputs four files: an updated FITS file, a CSV file, a PDF file, and a DS9 region file.

The updated FITS file is a copy of the original mask FITS file with the previously blank columns meant to record the slit center RA and slit center Dec positions filled in.

The PDF is the quick look plot generated by the code that shows the mask slit positions on the full mask field as well as zoomed in images around up to 6 alignment boxes with the selected corresponding GAIA catalog objects marked.

The CSV file contains the original coordinates of the slit corners in the mask coordinate frame, the calculated WCS positions of the slit corners, the slit center positions in both coordinate frames, and the coordinates of the GAIA catalog objects selected to go in the alignment boxes (zeros are recorded for non alignment box slits).

The DS9 region file contains the final slit positions for use in DS9.

# Other Useful information

'LRIS_Mask_Coords_to_WCS.py' contains the main function, which begins by reading the FITS file information and ends with writing the output files.
'astrometrycorrection.py' contains the routine for the distortion correction at the LRIS focal plane.
'precessionoutine.py' and 'refractioncorrection.py' contain the precession and refraction routines adapted from autoslit3.
'find_shifts.py' contains the function for finding the astrometry correction (shift) based on GAIA catalog object positions.
'create_err_plot.py' contains the function for writing the quick look plot.

The version of autoslit3 originally used to create the mask will affect the plate scale assumed when designing the mask due to the Atmospheric Dispersion Corrector (ADC) being installed. Our code makes a general cut on whether the pre or post-ADC plate scale was used based on the mask design date (pre August 2007 masks are assumed to be pre-ADC). This can be manually changed if needed.
