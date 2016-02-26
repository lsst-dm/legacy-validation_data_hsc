=======================
``validation_data_hsc``
=======================

``validation_data_hsc`` provides test data for exercising the LSST stack
through single frame and coadd processing based on engineering test data from
Hyper Suprime-Cam.

``raw`` contains 33 raw images from HSC.

``CALIB`` contains the calibration information sufficient to ISR the raw images.

``DATA`` is a Butler repo containing both the ``raw`` data and the results running ``ci_hsc`` processing of single-frame measurement, coadds, and multi-band forced photometry from the coadd detections.

``scons.log`` is the output log from running ``ci_hsc``.


Git LFS
=======

The data used by ``ci_hsc`` is stored using `Git LFS`_; refer to the `relevant
LSST documentation`_ for details on how to check out this repository.

.. _Git LFS: https://git-lfs.github.com
.. _relevant LSST documentation: http://developer.lsst.io/en/latest/tools/git_lfs.html

