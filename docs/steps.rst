Steps
=====

Functions Handling Individual Steps in ``xytowcs()``
----------------------------------------------------

These functions are called while running ``xytowcs()``. They handle individual
steps required to read in LRIS slitmask information and transform the slitmask coordinates
in the milling machine frame to sky coordinates.

.. autofunction:: tilsotua.generate_object_cat

.. autofunction:: tilsotua.use_autoslit_input.gen_from_auto

.. autofunction:: tilsotua.refraction_correction

.. autofunction:: tilsotua.astrometrycorrection.astrometry_calc

.. autofunction:: tilsotua.find_shifts.get_shift

.. autofunction:: tilsotua.create_quicklookplot.create_qlp

.. autofunction:: tilsotua.write_ds9_file.create_ds9_file
