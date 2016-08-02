=======================
``validation_data_hsc``
=======================

``validation_data_hsc`` provides test data for exercising the LSST stack
through single frame and coadd processing based on engineering test data from
Hyper Suprime-Cam.

``raw`` contains 8 raw exposures (112 CCDs per exposure) from HSC:
* visit=903982 field=STRIPE82L filter=HSC-I expTime=30.0 ra=21:19:59.9840 dec=+00:00:00.060
* visit=904006 field=STRIPE82L filter=HSC-I expTime=30.0 ra=21:19:59.9840 dec=+00:15:00.130
* visit=904828 field=STRIPE82L filter=HSC-I expTime=60.0 ra=21:19:59.9800 dec=+00:00:00.040
* visit=904846 field=STRIPE82L filter=HSC-I expTime=300.0 ra=21:20:03.8470 dec=+00:00:15.480
* visit=903332 field=STRIPE82L filter=HSC-R expTime=30.0 ra=21:19:59.9840 dec=+00:00:00.010
* visit=903340 field=STRIPE82L filter=HSC-R expTime=30.0 ra=21:19:59.9820 dec=+00:15:00.090
* visit=904350 field=STRIPE82L filter=HSC-Y expTime=30.0 ra=21:19:59.9800 dec=+00:00:00.080
* visit=904378 field=STRIPE82L filter=HSC-Y expTime=30.0 ra=21:22:59.9830 dec=+00:15:00.020

* Visits 903982, 904006, 904828, 904846 are i-band from commissioning run 3 (November 2013),
  with exposure times ranging from 30 sec to 300 sec, and up to a 15 arcmin dither;
* Visits 903332, 903340 are r-band from commissioning run 2 (June 2013); and
* Visits 904350, 904378 are y-band from commissioning run 3.

``CALIB`` contains the calibration information sufficient to ISR the raw images.
No darks are available from commissioning run 2, so we use the ones from
commissioning run 3.  No fringe frames have been constructed from commissioning
run 3, so we use some from March 2014.

``DATA`` is a Butler repo containing both the ``raw`` data and the results of running
single-frame measurement, coadds, and multi-band forced photometry from the coadd detections.


Git LFS
=======

The data is stored using `Git LFS`_; refer to the `relevant
LSST documentation`_ for details on how to check out this repository.

.. _Git LFS: https://git-lfs.github.com
.. _relevant LSST documentation: http://developer.lsst.io/en/latest/tools/git_lfs.html

