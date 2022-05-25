#Example script using tilsotua to process an archive mask FITS file
#This test mask is from the UCO/Lick archive, but tilsotua can also
#the ".file3" output file from autoslit.

#import the function that handles all of the processing from mask file to
#RA/Dec positions of slits.
from tilsotua import xytowcs

#Call function to process mask FITS file.
#"xytowcs" requires two arguments:the mask file name first, and a name for the
#output files that tilsotua will write.
xytowcs('test_mask.fits','tilsotua_example_mask')

#Execute "python test.py" from the command line to run script.
