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

Results of the processing are in the rerun "20160805".

The stack version is based on the current master, with the following changes:

* a version of GalSim patched for DM-7114 for the multiband stage, to avoid assertion failures.
* tickets/DM-7117 to prevent multiband failures in matching.
* tickets/DM-6612 to fix coadding.
* tickets/DM-7134 to fix singleFrameDriver parallelism.
* Disabled (unsetup -j) meas_modelfit because it makes everything run MUCH longer.

The versions actually used in the processing are recorded as the "packages" dataset.
Here's the result of `eups list -s`:

::

    activemqcpp           3.9.0.lsst2+3     b2153 b2152 b2151 b2150 setup
    afw                   12.0-16-ga7bfefa+2        b2153 setup
    afwdata               12.0+2            b2153 b2152 setup
    apr                   1.5.2             b2153 b2152 b2151 b2150 setup
    apr_util              1.5.4             b2153 b2152 b2151 b2150 setup
    astrometry_net        0.50.lsst3+1      b2153 b2152 setup
    astropy               0.0.1.lsst2       b2153 b2152 b2151 b2150 setup
    base                  12.0-4-g0d328d3+1         b2153 setup
    boost                 1.60+1            b2153 b2152 b2151 b2150 setup
    cat                   12.0-1-gbb14ecc+7         b2153 setup
    cfitsio               3360.lsst4        b2153 b2152 b2151 b2150 setup
    ci_hsc                LOCAL:/tigress/pprice/ci_hsc      setup
    coadd_chisquared      12.0-1-g10e615d+11        b2153 setup
    coadd_utils           12.0-1-g2dfa1ca+11        b2153 setup
    ctrl_events           12.0-2-gb6b4f2a   b2153 setup
    ctrl_execute          12.0+12           b2153 setup
    ctrl_orca             12.0+11           b2153 setup
    ctrl_platform_gordon  11.0+117          b2153 setup
    ctrl_platform_lsst    11.0+117          b2153 setup
    ctrl_pool             LOCAL:/home/pprice/LSST/ctrl/pool         setup
    ctrl_provenance       12.0+10           b2153 setup
    daf_base              12.0-1-g22b9da8+3         b2153 setup
    daf_butlerUtils       12.0-6-g459529e+2         b2153 setup
    daf_persistence       12.0-3-g711c365+1         b2153 setup
    datarel               12.0+28           b2153 setup
    db                    12.0+9            b2153 setup
    display_ds9           12.0-1-g2f8d434+7         b2153 setup
    doxygen               1.8.5.lsst1       b2153 b2152 b2151 b2150 setup
    eigen                 3.2.5.lsst2       b2153 b2152 b2151 b2150 setup
    esutil                0.5.3+1           b2153 setup
    fftw                  3.3.4.lsst2       b2153 b2152 b2151 b2150 setup
    galsim                LOCAL:/home/pprice/LSST/external/galsim/upstream/GalSim  setup
    geom                  10.0+63           b2153 setup
    gsl                   1.16.lsst3        b2153 b2152 b2151 b2150 setup
    ip_diffim             12.0-6-g840e2db+3         b2153 setup
    ip_isr                12.0-4-g9aeca5b+8         b2153 setup
    lmfit                 0.9.3+2           b2153 setup
    log                   12.0-3-g64b8957+1         b2153 setup
    log4cxx               0.10.0.lsst6+1    b2153 setup
    lsst                  12.0.rc1-2-g4137827+4     b2153 setup
    lsst_apps             12.0+27           b2153 setup
    lsst_build            LOCAL:/tigress/pprice/lsstsw/lsst_build   setup
    lsst_distrib          12.0-1-g51ba4ca+30        b2153 setup
    mariadb               10.1.11.lsst2     b2153 b2152 b2151 b2150 setup
    mariadbclient         10.1.11.lsst3     b2153 b2152 b2151 b2150 setup
    matplotlib            0.0.2-1-g49c793a  b2153 setup
    meas_algorithms       12.0-19-g23734c2+1        b2153 setup
    meas_astrom           LOCAL:/home/pprice/LSST/meas/astrom       setup
    meas_base             12.0-5-g9359e01   b2153 setup
    meas_deblender        12.0-2-gf5405dc+20        b2153 setup
    meas_extensions_photometryKron 7.3.1.0-14-gc1df0b4      b2153 setup
    meas_extensions_psfex 12.0-3-gc91cf5c+2         b2153 setup
    meas_extensions_shapeHSM 12.0+22        b2153 setup
    meas_extensions_simpleShape 12.0-1-g6a3db40+15  b2153 setup
    miniconda2            3.19.0.lsst4      current setup
    minuit2               5.34.14           b2153 b2152 b2151 b2150 setup
    mpi                   0.0.1+1           b2153 b2152 b2151 b2150 setup
    mpi4py                1.3.1.lsst2+1     b2153 b2152 setup
    mpich                 3.2               b2153 b2152 b2151 b2150 setup
    mysqlpython           1.2.3.lsst2+2     b2153 b2152 b2151 b2150 setup
    ndarray               12.0+6            b2153 setup
    numpy                 0.0.2+1           b2153 b2152 setup
    obs_lsstSim           12.0-5-gc3154c7+5         b2153 setup
    obs_sdss              12.0-4-g9494fec+13        b2153 setup
    obs_subaru            LOCAL:/home/pprice/LSST/obs/subaru        setup
    obs_test              12.0-5-gcdb69c4+2         b2153 setup
    pex_config            12.0-2-g55cb508+2         b2153 setup
    pex_exceptions        12.0+3            b2153 setup
    pex_logging           12.0+4            b2153 setup
    pex_policy            12.0-1-g298db87+3         b2153 setup
    pipe_base             12.0-5-g6728465+11        b2153 setup
    pipe_drivers          LOCAL:/home/pprice/LSST/pipe/drivers      setup
    pipe_tasks            LOCAL:/home/pprice/LSST/pipe/tasks        setup
    psfex                 12.0+2            b2153 setup
    pyfits                3.4.0+4           b2153 setup
    python                0.0.4             b2153 b2152 b2151 b2150 setup
    python_d2to1          0.2.12.lsst1      b2153 b2152 b2151 b2150 setup
    python_psutil         4.1.0+1           b2153 b2152 b2151 b2150 setup
    pyyaml                3.11.lsst1+1      b2153 b2152 b2151 b2150 setup
    scipy                 0.0.1.lsst1+1     b2153 setup
    scisql                0.3.5+14          b2153 b2152 setup
    scons                 2.5.0.lsst1       b2153 b2152 setup
    sconsUtils            12.0-2-ga8b26fe   b2153 b2152 setup
    shapelet              12.0-2-g4712fae+7         b2153 setup
    skymap                12.0-1-gff15314+3         b2153 setup
    skypix                10.0+422          b2153 setup
    sqlalchemy            1.0.8.lsst3+1     b2153 setup
    stsci_distutils       0.3.7.lsst1       b2153 b2152 b2151 b2150 setup
    swig                  3.0.2.lsst1       b2153 b2152 b2151 b2150 setup
    testdata_subaru       12.0              b2153 b2152 b2151 b2150 setup
    tmv                   0.73+2            b2153 setup
    utils                 12.0-3-gb6e8129   b2153 setup
    wcslib                5.13.lsst1        b2153 b2152 b2151 b2150 setup
    xpa                   2.1.15.lsst3      b2153 b2152 b2151 b2150 setup


Here is a list of commands that were run to generate the processed data.

::

    # Setup
    export OMP_NUM_THREADS=1

    # Fake out the pipeline about the origin of the reference catalog
    export SETUP_ASTROMETRY_NET_DATA="astrometry_net_data sdss-dr9-fink-v5b"
    export ASTROMETRY_NET_DATA_DIR=`pwd`/sdss-dr9-fink-v5b

    # Ingest raw data into repo
    mkdir DATA
    echo lsst.obs.hsc.HscMapper > DATA/_mapper
    ingestImages.py DATA --mode=link 'raw/*.fits'

    # Heavy lifting
    singleFrameDriver.py DATA --calib CALIB --rerun 20160805 --job singleFrame --cores 16 --id ccd=0..103 visit=903982^904006^904828^904846^903332^903340^904350^904378
    makeDiscreteSkyMap.py DATA --rerun 20160805 --id ccd=0..103 visit=903982^904006^904828^904846^903332^903340^904350^904378
    # makeDiscreteSkyMap: tract 0 has corners (321.714, -1.294), (318.915, -1.294), (318.915, 1.504), (321.714, 1.504) (RA, Dec deg) and 15 x 15 patches
    coaddDriver.py DATA --rerun 20160805 --job coadd --cores 16 --id tract=0 filter=HSC-I --selectId ccd=0..103 visit=903982^904006^904828^904846
    coaddDriver.py DATA --rerun 20160805 --job coadd --cores 16 --id tract=0 filter=HSC-R --selectId ccd=0..103 visit=903332^903340
    coaddDriver.py DATA --rerun 20160805 --job coadd --cores 16 --id tract=0 filter=HSC-Y --selectId ccd=0..103 visit=904350^904378
    multiBandDriver.py DATA --rerun 20160805 --job multiband --cores 16 --id tract=0 filter=HSC-R^HSC-I^HSC-Y -C multiband-config.py


Issues
======

Nothing serious known yet, besides the need to include local versions of products in the stack.

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

