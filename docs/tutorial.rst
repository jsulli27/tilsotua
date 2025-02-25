Tutorial
========

Example Application tilsotua
-----------------------------

The `Example <https://github.com/jsulli27/tilsotua/tree/master/Examples>`_ directory contains
a Python script and Jupyter Notebook with an example application of tilsotua.

Test mask files are provided along with the output files generated for the mask by tilsotua.

To run tilsotua on the archival mask file, call the ``xytowcs()`` function with:

>>> from tilsotua import xytowcs
>>> xytowcs(data_input_name:'test_Mask.fits',output_file:'tilsotua_example_mask')

To run tilsotua on the autoslit output file, call the ``xytowcs()`` function with:

>>> from tilsotua import xytowcs
>>> xytowcs(data_input_name:'test_Mask.file3',output_file:'tilsotua_example_mask')

This will generate the output files as provided in the Examples directory. You should
check the final alignment of the mask with the 'quick-look' plot.
