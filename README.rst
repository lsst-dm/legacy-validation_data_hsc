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

Note:

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

Processing
==========

Results of the processing are 'DATA'.

The stack version is based on v16.0.rc1

Here the a list of commands that were run on lsst-dev01.ncsa.illinois.edu to generate the processed data.

::
    #!/bin/bash

    source /software/lsstsw/stack/loadLSST.bash
    # Pick up pipe_drivers in lsst_distrib
    setup lsst_distrib -t w_2018_22  # which was v16.0.rc1
    # Set up our specific camera
    setup obs_subaru -t w_2018_22  # which was v16.0.rc1

    # Setup
    export OMP_NUM_THREADS=1

    # Fake out the pipeline about the origin of the reference catalog
    export SETUP_ASTROMETRY_NET_DATA="astrometry_net_data sdss-dr9-fink-v5b"
    export ASTROMETRY_NET_DATA_DIR=sdss-dr9-fink-v5b

    # Ingest raw data into repo
    DATA="DATA"
    mkdir "${DATA}"
    echo lsst.obs.hsc.HscMapper > "${DATA}"/_mapper
    ingestImages.py "${DATA}" --mode=link 'raw/*.fits'

    # Heavy lifting
    singleFrameDriver.py "${DATA}" --calib CALIB --output "${DATA}" -C config/singleFrameDriverConfig.py --job singleFrame --cores 16 --id ccd=0..103 visit=903982^904006^904828^904846^903332^903340^904350^904378
    makeDiscreteSkyMap.py "${DATA}" --output "${DATA}" --id ccd=0..103 visit=903982^904006^904828^904846^903332^903340^904350^904378
    # makeDiscreteSkyMap: tract 0 has corners (321.714, -1.294), (318.915, -1.294), (318.915, 1.504), (321.714, 1.504) (RA, Dec deg) and 15 x 15 patches
    coaddDriver.py "${DATA}" --calib CALIB --output "${DATA}" -C config/coaddConfig.py --job coadd --cores 16 --id tract=0 filter=HSC-I --selectId ccd=0..103 visit=903982^904006^904828^904846
    coaddDriver.py "${DATA}" --calib CALIB --output "${DATA}" -C config/coaddConfig.py --job coadd --cores 16 --id tract=0 filter=HSC-R --selectId ccd=0..103 visit=903332^903340
    coaddDriver.py "${DATA}" --calib CALIB --output "${DATA}" -C config/coaddConfig.py --job coadd --cores 16 --id tract=0 filter=HSC-Y --selectId ccd=0..103 visit=904350^904378
    multiBandDriver.py "${DATA}" --calib CALIB --output "${DATA}" -C config/multiBandConfig.py --job multiband --cores 16 --id tract=0 filter=HSC-R^HSC-I^HSC-Y


These commands are saved in `reprocess.sh` and were submitted using `sbatch reprocess.sl` which is a trivial SLURM wrapper for `reprocess.sh`.

Issues
======

Should update astrometric fitting to PS1 catalogs.

LOTS of warnings of the form (the particular measurement and coordinates change):

::

  multiBandDriver.measureCoaddSources.measurement WARNING: Error in base_PsfFlux.measure on record 704374639441: 
    File "src/CoaddPsf.cc", line 235, in virtual std::shared_ptr<lsst::afw::image::Image<double> > lsst::meas::algorithms::CoaddPsf::doComputeKernelImage(const Point2D&, const lsst::afw::image::Color&) const
      Cannot compute CoaddPsf at point (23775, 18940); no input images at that point. {0}
  lsst::pex::exceptions::InvalidParameterError: 'Cannot compute CoaddPsf at point (23775, 18940); no input images at that point.'

At the present, I don't believe these are important, but are just excess chatter
from the measurement framework.


Git LFS
=======

The data is stored using `Git LFS`_; refer to the `relevant
LSST documentation`_ for details on how to check out this repository.

.. _Git LFS: https://git-lfs.github.com
.. _relevant LSST documentation: http://developer.lsst.io/en/latest/tools/git_lfs.html

