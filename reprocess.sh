#!/bin/bash

source /software/lsstsw/stack/loadLSST.bash
# Set up our specific camera
# w_2018_35 was the weekly before the new HTM-based PS1 catalogs were merged
setup obs_subaru -t w_2018_35
# Add meas_extensions_photometryKron module to do Kron photometry.
setup meas_extensions_photometryKron -t w_2018_35
# Set up this specific repo just so we can generically refer to VALIDATION_DATA_HSC_DIR below.
setup -k -r validation_data_hsc
#
# Setup
export OMP_NUM_THREADS=1

# Ingest raw data into repo
DATA="data"
mkdir "${DATA}"
echo lsst.obs.hsc.HscMapper > "${DATA}"/_mapper
ingestImages.py "${DATA}" --mode=link "${VALIDATION_DATA_HSC_DIR}"/raw/*.fits

# Copy in the HSC transmission files (atmosphere+optics+sensors).
cp -pr transmission "${DATA}"

# Link in the reference catalogs
ln -s "${VALIDATION_DATA_HSC_DIR}"/ref_cats "${DATA}"/ref_cats

# Heavy lifting
singleFrameDriver.py "${DATA}" --calib CALIB --output "${DATA}" -C config/singleFrameDriverConfig.py --job singleFrame --cores 16 --id ccd=0..103 visit=903982^904006^904828^904846^903332^903340^904350^904378
makeDiscreteSkyMap.py "${DATA}" --output "${DATA}" --id ccd=0..103 visit=903982^904006^904828^904846^903332^903340^904350^904378
# makeDiscreteSkyMap: tract 0 has corners (321.714, -1.294), (318.915, -1.294), (318.915, 1.504), (321.714, 1.504) (RA, Dec deg) and 15 x 15 patches
coaddDriver.py "${DATA}" --calib CALIB --output "${DATA}" -C config/coaddConfig.py --job coadd --cores 16 --id tract=0 filter=HSC-I --selectId ccd=0..103 visit=903982^904006^904828^904846
coaddDriver.py "${DATA}" --calib CALIB --output "${DATA}" -C config/coaddConfig.py --job coadd --cores 16 --id tract=0 filter=HSC-R --selectId ccd=0..103 visit=903332^903340
coaddDriver.py "${DATA}" --calib CALIB --output "${DATA}" -C config/coaddConfig.py --job coadd --cores 16 --id tract=0 filter=HSC-Y --selectId ccd=0..103 visit=904350^904378
multiBandDriver.py "${DATA}" --calib CALIB --output "${DATA}" -C config/multiBandConfig.py --job multiband --cores 16 --id tract=0 filter=HSC-R^HSC-I^HSC-Y
