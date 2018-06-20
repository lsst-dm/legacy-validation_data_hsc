#!/bin/bash

source /software/lsstsw/stack/loadLSST.bash
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
singleFrameDriver.py "${DATA}" --calib CALIB --output "${DATA}" -C config/singleFrameDriverconfig.py --job singleFrame --cores 16 --id ccd=0..103 visit=903982^904006^904828^904846^903332^903340^904350^904378
makeDiscreteSkyMap.py "${DATA}" --output "${DATA}" --id ccd=0..103 visit=903982^904006^904828^904846^903332^903340^904350^904378
# makeDiscreteSkyMap: tract 0 has corners (321.714, -1.294), (318.915, -1.294), (318.915, 1.504), (321.714, 1.504) (RA, Dec deg) and 15 x 15 patches
coaddDriver.py "${DATA}" --calib CALIB --output "${DATA}" -C config/coaddConfig.py --job coadd --cores 16 --id tract=0 filter=HSC-I --selectId ccd=0..103 visit=903982^904006^904828^904846
coaddDriver.py "${DATA}" --calib CALIB --output "${DATA}" -C config/coaddConfig.py --job coadd --cores 16 --id tract=0 filter=HSC-R --selectId ccd=0..103 visit=903332^903340
coaddDriver.py "${DATA}" --calib CALIB --output "${DATA}" -C config/coaddConfig.py --job coadd --cores 16 --id tract=0 filter=HSC-Y --selectId ccd=0..103 visit=904350^904378
multiBandDriver.py "${DATA}" --calib CALIB --output "${DATA}" -C config/multiBandConfig.py --job multiband --cores 16 --id tract=0 filter=HSC-R^HSC-I^HSC-Y
