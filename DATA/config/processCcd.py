import lsst.pipe.tasks.processCcd
assert type(config)==lsst.pipe.tasks.processCcd.ProcessCcdConfig, 'config is of type %s.%s instead of lsst.pipe.tasks.processCcd.ProcessCcdConfig' % (type(config).__module__, type(config).__name__)
import filecmp
import eups.distrib.server
import eups.distrib.tarball
import eups.cmd
import lsst.pipe.tasks.setConfigFromEups
import eups.distrib.DistribFactory
import eups.db.Database
import eups.tags
import eups.Eups
import eups.distrib.eupspkg
import eups.distrib.Distrib
import eups
import eups.stack.ProductStack
import eups.exceptions
import eups.Uses
import ConfigParser
import eups.VersionParser
import lsst.obs.subaru.crosstalkYagi
import lsst.meas.extensions
import lsst.meas.extensions.psfex
import eups.distrib.Repositories
import eups.hooks
import eups.db
import lsst.obs.subaru.crosstalk
import lsst.meas.extensions.psfex.psfexLib
import eups.stack
import lsst.obs.hsc.vignette
import eups.app
import pipes
import lsst.obs.subaru.isr
import tarfile
import eups.VersionCompare
import eups.distrib
import lsst.afw.display.displayLib
import lsst.meas.extensions.psfex.psfexPsfDeterminer
import eups.distrib.Repository
import lsst.afw.display.rgb
import optparse
import eups.db.VersionFile
import eups.table
import eups.distrib.pacman
import eups.distrib.builder
import eups.lock
import eups.utils
import eups.Product
import eups.db.ChainFile
import eups.stack.ProductFamily
# Include HeavyFootprint data in source table?
config.doWriteHeavyFootprintsInSources=False

# Separate sources which are blended into distinct entities
config.doDeblend=True

# Write sources?
config.doWriteSources=True

import lsst.obs.subaru.isr
config.isr.retarget(target=lsst.obs.subaru.isr.SubaruIsrTask, ConfigClass=lsst.obs.subaru.isr.SubaruIsrConfig)
# Number of stdev below the background to set thumbnail minimum
config.isr.thumbnailStdev=3.0

# How to estimate the average value for BAD regions.
# Allowed values:
# 	MEDIAN	Correct using the median of the good data
# 	MEANCLIP	Correct using the (clipped) mean of good data
# 	None	Field is optional
# 
config.isr.badStatistic='MEANCLIP'

# Order of polynomial or to fit if overscan fit type is a polynomial, or number of spline knots if overscan fit type is a spline.
config.isr.overscanOrder=30

# Widen bleed trails based on their width?
config.isr.doWidenSaturationTrails=True

# Do fringe subtraction after flat-fielding?
config.isr.fringeAfterFlat=True

import lsst.obs.subaru.crosstalk
config.isr.crosstalk.retarget(target=lsst.obs.subaru.crosstalk.CrosstalkTask, ConfigClass=lsst.obs.subaru.crosstalk.CrosstalkConfig)
# Shape of coeffs array
config.isr.crosstalk.coeffs.shape=[4, 4]

# Crosstalk coefficients
config.isr.crosstalk.coeffs.values=[0.0, -0.000125, -0.000149, -0.000156, -0.000124, 0.0, -0.000132, -0.000157, -0.000171, -0.000134, 0.0, -0.000153, -0.000157, -0.000151, -0.000137, 0.0]

# Name for crosstalk mask plane
config.isr.crosstalk.crosstalkMaskPlane='CROSSTALK'

# Set crosstalk mask plane for pixels over this value
config.isr.crosstalk.minPixelToMask=45000.0

# Apply dark frame correction?
config.isr.doDark=True

# update exposure metadata in the assembled ccd to reflect the effective gain of the assembled chip
config.isr.setGainAssembledCcd=True

# Do overscan subtraction?
config.isr.doOverscan=True

# Correct for crosstalk
config.isr.doCrosstalk=True

# FWHM of PSF used when interpolating over bad columns (arcsec)
config.isr.fwhmForBadColumnInterpolation=1.0

# The read noise to use if no Detector is present in the Exposure
config.isr.readNoise=0.0

# FWHM of PSF (arcsec)
config.isr.fwhm=1.0

# Maximum number of iterations for the brighter fatter correction
config.isr.brighterFatterMaxIter=10

# If flatScalingType is 'USER' then scale flat by this amount; ignored otherwise
config.isr.flatUserScale=1.0

# Number of points to define the Vignette polygon
config.isr.numPolygonPoints=100

# Name of mask plane to use in saturation detection and interpolation
config.isr.saturatedMaskName='SAT'

# Should the gain be applied when applying the brighter fatter correction?
config.isr.brighterFatterApplyGain=True

# Should we set the level of all BAD patches of the chip to the chip's average value?
config.isr.doSetBadRegions=True

# Threshold used to stop iterating the brighter fatter correction.  It is the  absolute value of the difference between the current corrected image and the one from the previous iteration summed over all the pixels.
config.isr.brighterFatterThreshold=1000.0

# Correct for nonlinearity of the detector's response (ignored if coefficients are 0.0)
config.isr.doLinearize=True

# Kernel file used for the brighter fatter correction
config.isr.brighterFatterKernelFile=''

# Apply the brighter fatter correction
config.isr.doBrighterFatter=True

# Apply bias frame correction?
config.isr.doBias=True

# Normalize all the amplifiers in each CCD to have the same gain
# 
# This does not measure the gains, it simply forces the median of each amplifier to be equal
# after applying the nominal gain
# 
config.isr.normalizeGains=False

# Remove any PC cards in the header
config.isr.removePcCards=True

# Apply fringe correction?
config.isr.doFringe=True

# trim out non-data regions?
config.isr.assembleCcd.doTrim=True

# FITS headers to remove (in addition to DATASEC, BIASSEC, TRIMSEC and perhaps GAIN)
config.isr.assembleCcd.keysToRemove=[]

# renormalize to a gain of 1? (ignored if setGain false)
config.isr.assembleCcd.doRenorm=False

# set gain?
config.isr.assembleCcd.setGain=True

# The saturation level to use if no Detector is present in the Exposure (ignored if NaN)
config.isr.saturation=float('nan')

# Calculate variance?
config.isr.doVariance=True

# Assemble amp-level calibration exposures into ccd-level exposure?
config.isr.doAssembleIsrExposures=False

# Softening parameter for thumbnail mapping
config.isr.thumbnailQ=20.0

# fields to remove from the metadata of the assembled ccd.
config.isr.keysToRemoveFromAssembledCcd=[]

# Center of vignetting pattern, in x (focal plane coords)
config.isr.vignette.xCenter=-100.0

# Radius of vignetting pattern, in focal plane coords
config.isr.vignette.radius=17500.0

# Center of vignetting pattern, in y (focal plane coords)
config.isr.vignette.yCenter=100.0

# Apply flat field correction?
config.isr.doFlat=True

# The gain to use if no Detector is present in the Exposure (ignored if NaN)
config.isr.gain=float('nan')

# Rejection threshold (sigma) for collapsing overscan before fit
config.isr.overscanRej=3.0

# Border around saturated pixels for thumbnail
config.isr.thumbnailSatBorder=2

# Mask saturated pixels?
config.isr.doSaturation=True

# Trim guider shadow
config.isr.doGuider=False

# The method for scaling the flat on the fly.
# Allowed values:
# 	None	Field is optional
# 	MEDIAN	Scale by the inverse of the median
# 	USER	Scale by flatUserScale
# 	MEAN	Scale by the inverse of the mean
# 
config.isr.flatScalingType='USER'

# Number of pixels by which to grow the saturation footprints
config.isr.growSaturationFootprintSize=1

# Default value for fluxMag0T1 (for an unrecognised filter)
config.isr.defaultFluxMag0T1=158489319246.11172

# Correct the amplifiers for their gains
# 
# N.b. this is intended to be used *instead* of doFlat; it's useful if you're measuring system throughput
# 
config.isr.doApplyGains=False

# Persist Polygon used to define vignetted region?
config.isr.doWriteVignettePolygon=True

# Offset to the random number generator seed (full seed includes exposure ID)
config.isr.fringe.stats.rngSeedOffset=0

# Ignore pixels with these masks
config.isr.fringe.stats.badMaskPlanes=['SAT', 'NO_DATA']

# Statistic to use
config.isr.fringe.stats.stat=32

# Number of fitting iterations
config.isr.fringe.stats.iterations=3

# Sigma clip threshold
config.isr.fringe.stats.clip=3.0

# Only fringe-subtract these filters
config.isr.fringe.filters=['y', 'N921']

# Sigma clip threshold
config.isr.fringe.clip=3.0

# Half-size of large (background) measurements (pixels)
config.isr.fringe.large=30

# Number of fringe measurements
config.isr.fringe.num=30000

# Number of fitting iterations
config.isr.fringe.iterations=20

# Half-size of small (fringe) measurements (pixels)
config.isr.fringe.small=3

# Remove fringe pedestal?
config.isr.fringe.pedestal=False

# The method for fitting the overscan bias level.
# Allowed values:
# 	None	Field is optional
# 	LEG	Fit Legendre polynomial to the longest axis of the overscan region
# 	CUBIC_SPLINE	Fit cubic spline to the longest axis of the overscan region
# 	MEDIAN	Correct using the median of the overscan region
# 	POLY	Fit ordinary polynomial to the longest axis of the overscan region
# 	CHEB	Fit Chebyshev polynomial to the longest axis of the overscan region
# 	AKIMA_SPLINE	Fit Akima spline to the longest axis of the overscan region
# 	NATURAL_SPLINE	Fit natural spline to the longest axis of the overscan region
# 	MEAN	Correct using the mean of the overscan region
# 
config.isr.overscanFitType='AKIMA_SPLINE'

# Write OverScan-Subtracted thumbnail?
config.isr.qa.doThumbnailOss=True

# Mesh size in X (pix) to calculate count statistics
config.isr.qa.flatness.meshX=256

# How many times do we iterate clipping outliers in calculate count statistics?
config.isr.qa.flatness.nIter=3

# Do we clip outliers in calculate count statistics?
config.isr.qa.flatness.doClip=True

# How many sigma is used to clip outliers in calculate count statistics?
config.isr.qa.flatness.clipSigma=3.0

# Mesh size in Y (pix) to calculate count statistics
config.isr.qa.flatness.meshY=256

# Write flattened thumbnail?
config.isr.qa.doThumbnailFlattened=True

# Write OverScan-Subtracted image?
config.isr.qa.doWriteOss=False

# Write flattened image?
config.isr.qa.doWriteFlattened=False

# Range for thumbnail mapping
config.isr.thumbnailRange=5.0

# Persist postISRCCD?
config.isr.doWrite=False

# Mask defect pixels?
config.isr.doDefect=True

# The approximate flux of a zero-magnitude object in a one-second exposure, per filter
config.isr.fluxMag0T1={'g': 398107170553.49854, 'N816': 15848931924.611174, 'i': 275422870333.81744, 'r': 398107170553.49854, 'N921': 19054607179.632523, 'N515': 20892961308.54041, 'y': 91201083935.59116, 'z': 120226443461.74132}

# Maximum deviation from the median for overscan
config.isr.overscanMaxDev=1000.0

# Assemble amp-level exposures into a ccd-level exposure?
config.isr.doAssembleCcd=True

# Tweak flats to match observed amplifier ratios?
config.isr.doTweakFlat=False

# Binning factor for thumbnail
config.isr.thumbnailBinning=4

# Detect sources?
config.doDetection=True

# Run fake sources injection task
config.doFakes=False

# Perform ISR?
config.doIsr=True

# Mask plane to set on pixels affected by fakes.  Will be added if not already present.
config.fakes.maskPlaneName='FAKE'

# Estimate the background again after final source detection?
config.detection.reEstimateBackground=True

# Grow detections by nSigmaToGrow * sigma; if 0 then do not grow
config.detection.nSigmaToGrow=2.4

# detected sources with fewer than the specified number of pixels will be ignored
# 	Valid Range = [0,inf)
config.detection.minPixels=1

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.detection.tempLocalBackground.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.detection.tempLocalBackground.approxOrderY=-1

# Names of mask planes to ignore while estimating the background
config.detection.tempLocalBackground.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# behaviour if there are too few points in grid for requested interpolation style
# Allowed values:
# 	None	Field is optional
# 	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
# 	THROW_EXCEPTION	throw an exception if there are too few points
# 	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
# 
config.detection.tempLocalBackground.undersampleStyle='REDUCE_INTERP_ORDER'

# how to interpolate the background values. This maps to an enum; see afw::math::Background
# Allowed values:
# 	NONE	No background estimation is to be attempted
# 	CONSTANT	Use a single constant value
# 	None	Field is optional
# 	LINEAR	Use linear interpolation
# 	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
# 	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
# 
config.detection.tempLocalBackground.algorithm='AKIMA_SPLINE'

# Ignore NaNs when estimating the background
config.detection.tempLocalBackground.isNanSafe=False

# how large a region of the sky should be used for each background point
# 	Valid Range = [10,inf)
config.detection.tempLocalBackground.binSize=64

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.detection.tempLocalBackground.weighting=True

# type of statistic to use for grid points
# Allowed values:
# 	None	Field is optional
# 	MEDIAN	median
# 	MEANCLIP	clipped mean
# 	MEAN	unclipped mean
# 
config.detection.tempLocalBackground.statisticsProperty='MEANCLIP'

# Use Approximate (Chebyshev) to model background.
config.detection.tempLocalBackground.useApprox=False

# Include threshold relative to thresholdValue
# 	Valid Range = [0.0,inf)
config.detection.includeThresholdMultiplier=1.0

# Pixels should be grown as isotropically as possible (slower)
config.detection.isotropicGrow=True

# Fiddle factor to add to the background; debugging only
config.detection.adjustBackground=0.0

# Do temporary interpolated background subtraction before footprint detection?
config.detection.doTempLocalBackground=False

# specifies whether to detect positive, or negative sources, or both
# Allowed values:
# 	positive	detect only positive sources
# 	negative	detect only negative sources
# 	both	detect both positive and negative sources
# 
config.detection.thresholdPolarity='positive'

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.detection.background.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.detection.background.approxOrderY=-1

# Names of mask planes to ignore while estimating the background
config.detection.background.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# behaviour if there are too few points in grid for requested interpolation style
# Allowed values:
# 	None	Field is optional
# 	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
# 	THROW_EXCEPTION	throw an exception if there are too few points
# 	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
# 
config.detection.background.undersampleStyle='REDUCE_INTERP_ORDER'

# how to interpolate the background values. This maps to an enum; see afw::math::Background
# Allowed values:
# 	NONE	No background estimation is to be attempted
# 	CONSTANT	Use a single constant value
# 	None	Field is optional
# 	LINEAR	Use linear interpolation
# 	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
# 	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
# 
config.detection.background.algorithm='NATURAL_SPLINE'

# Ignore NaNs when estimating the background
config.detection.background.isNanSafe=False

# how large a region of the sky should be used for each background point
# 	Valid Range = [10,inf)
config.detection.background.binSize=128

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.detection.background.weighting=True

# type of statistic to use for grid points
# Allowed values:
# 	None	Field is optional
# 	MEDIAN	median
# 	MEANCLIP	clipped mean
# 	MEAN	unclipped mean
# 
config.detection.background.statisticsProperty='MEANCLIP'

# Use Approximate (Chebyshev) to model background.
config.detection.background.useApprox=False

# Grow detections to set the image mask bits, but return the original (not-grown) footprints
config.detection.returnOriginalFootprints=False

# specifies the desired flavor of Threshold
# Allowed values:
# 	pixel_stdev	threshold applied to per-pixel std deviation
# 	variance	threshold applied to image variance
# 	value	threshold applied to image value
# 	stdev	threshold applied to image std deviation
# 
config.detection.thresholdType='stdev'

# Threshold for footprints
# 	Valid Range = [0.0,inf)
config.detection.thresholdValue=5.0

# Write icSrc to reference matches?
config.doWriteCalibrateMatches=True

# Compute and write src to reference matches?
config.doWriteSourceMatches=True

# Pixel size (arcsec).  Only needed if no Wcs is provided
config.calibrate.initialPsf.pixelScale=0.25

# PSF model type
# Allowed values:
# 	DoubleGaussian	Double Gaussian model
# 	None	Field is optional
# 	SingleGaussian	Single Gaussian model
# 
config.calibrate.initialPsf.model='SingleGaussian'

# FWHM of PSF model (arcsec)
config.calibrate.initialPsf.fwhm=1.0

# Size of PSF model (pixels)
config.calibrate.initialPsf.size=15

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.calibrate.background.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.calibrate.background.approxOrderY=-1

# Names of mask planes to ignore while estimating the background
config.calibrate.background.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# behaviour if there are too few points in grid for requested interpolation style
# Allowed values:
# 	None	Field is optional
# 	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
# 	THROW_EXCEPTION	throw an exception if there are too few points
# 	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
# 
config.calibrate.background.undersampleStyle='REDUCE_INTERP_ORDER'

# how to interpolate the background values. This maps to an enum; see afw::math::Background
# Allowed values:
# 	NONE	No background estimation is to be attempted
# 	CONSTANT	Use a single constant value
# 	None	Field is optional
# 	LINEAR	Use linear interpolation
# 	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
# 	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
# 
config.calibrate.background.algorithm='NATURAL_SPLINE'

# Ignore NaNs when estimating the background
config.calibrate.background.isNanSafe=False

# how large a region of the sky should be used for each background point
# 	Valid Range = [10,inf)
config.calibrate.background.binSize=128

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.calibrate.background.weighting=True

# type of statistic to use for grid points
# Allowed values:
# 	None	Field is optional
# 	MEDIAN	median
# 	MEANCLIP	clipped mean
# 	MEAN	unclipped mean
# 
config.calibrate.background.statisticsProperty='MEANCLIP'

# Use Approximate (Chebyshev) to model background.
config.calibrate.background.useApprox=False

# Interpolate over defects? (ignored unless you provide a list of defects)
config.calibrate.repair.doInterpolate=True

# Type of statistic to calculate edge fallbackValue for interpolation
# Allowed values:
# 	None	Field is optional
# 	MEAN	mean
# 	MEDIAN	median
# 	USER	user value set in fallbackUserValue config
# 	MEANCLIP	clipped mean
# 
config.calibrate.repair.interp.fallbackValueType='MEANCLIP'

# Add a Gaussian to represent wings?
config.calibrate.repair.interp.modelPsf.addWing=True

# wing amplitude, as a multiple of core amplitude (dimensionless); ignored if addWing false
config.calibrate.repair.interp.modelPsf.wingAmplitude=0.1

# Minimum kernel size if using sizeFactor (pixels); ignored if size is not None
config.calibrate.repair.interp.modelPsf.minSize=5

# Maximum kernel size if using sizeFactor (pixels); ignored if size is not None
config.calibrate.repair.interp.modelPsf.maxSize=None

# Kernel size as a factor of fwhm (dimensionless); size = sizeFactor * fwhm; ignored if size is not None
config.calibrate.repair.interp.modelPsf.sizeFactor=3.0

# Default FWHM of Gaussian model of core of star (pixels)
config.calibrate.repair.interp.modelPsf.defaultFwhm=3.0

# wing width, as a multiple of core width (dimensionless); ignored if addWing false
config.calibrate.repair.interp.modelPsf.wingFwhmFactor=2.5

# Kernel size (width and height) (pixels); if None then sizeFactor is used
config.calibrate.repair.interp.modelPsf.size=None

# If fallbackValueType is 'USER' then use this as the fallbackValue; ignored otherwise
config.calibrate.repair.interp.fallbackUserValue=0.0

# Smoothly taper to the fallback value at the edge of the image?
config.calibrate.repair.interp.useFallbackValueAtEdge=True

# Allow negative values for egde interpolation fallbackValue?  If False, set fallbackValue to max(fallbackValue, 0.0)
config.calibrate.repair.interp.negativeFallbackAllowed=True

# Find and mask out cosmic rays?
config.calibrate.repair.doCosmicRay=True

# Don't interpolate over CR pixels
config.calibrate.repair.cosmicray.keepCRs=False

# used in condition 3 for CR; see CR.cc code
config.calibrate.repair.cosmicray.cond3_fac=2.5

# used in condition 3 for CR; see CR.cc code
config.calibrate.repair.cosmicray.cond3_fac2=0.4

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.calibrate.repair.cosmicray.background.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.calibrate.repair.cosmicray.background.approxOrderY=-1

# Names of mask planes to ignore while estimating the background
config.calibrate.repair.cosmicray.background.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# behaviour if there are too few points in grid for requested interpolation style
# Allowed values:
# 	None	Field is optional
# 	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
# 	THROW_EXCEPTION	throw an exception if there are too few points
# 	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
# 
config.calibrate.repair.cosmicray.background.undersampleStyle='REDUCE_INTERP_ORDER'

# how to interpolate the background values. This maps to an enum; see afw::math::Background
# Allowed values:
# 	NONE	No background estimation is to be attempted
# 	CONSTANT	Use a single constant value
# 	None	Field is optional
# 	LINEAR	Use linear interpolation
# 	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
# 	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
# 
config.calibrate.repair.cosmicray.background.algorithm='AKIMA_SPLINE'

# Ignore NaNs when estimating the background
config.calibrate.repair.cosmicray.background.isNanSafe=False

# how large a region of the sky should be used for each background point
# 	Valid Range = [10,inf)
config.calibrate.repair.cosmicray.background.binSize=100000

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.calibrate.repair.cosmicray.background.weighting=True

# type of statistic to use for grid points
# Allowed values:
# 	None	Field is optional
# 	MEDIAN	median
# 	MEANCLIP	clipped mean
# 	MEAN	unclipped mean
# 
config.calibrate.repair.cosmicray.background.statisticsProperty='MEDIAN'

# Use Approximate (Chebyshev) to model background.
config.calibrate.repair.cosmicray.background.useApprox=False

# number of times to look for contaminated pixels near known CR pixels
config.calibrate.repair.cosmicray.niteration=3

# maximum number of contaminated pixels
config.calibrate.repair.cosmicray.nCrPixelMax=1000000

# CRs must be > this many sky-sig above sky
config.calibrate.repair.cosmicray.minSigma=6.0

# CRs must have > this many DN (== electrons/gain) in initial detection
config.calibrate.repair.cosmicray.min_DN=150.0

# Field name prefix for the flux other measurements should be aperture corrected to match
config.calibrate.measureApCorr.refFluxName='base_CircularApertureFlux_17_0'

# Minimum number of degrees of freedom (# of valid data points - # of parameters); if this is exceeded, the order of the fit is decreased (in both dimensions), and if we can't decrease it enough, we'll raise ValueError.
# 	Valid Range = [1,inf)
config.calibrate.measureApCorr.minDegreesOfFreedom=1

# Name of a flag field that indicates that a source should be used to constrain the aperture corrections
config.calibrate.measureApCorr.inputFilterFlag='calib_psfUsed'

# Number of standard devisations to clip at
config.calibrate.measureApCorr.numSigmaClip=3.0

# if true, only include terms where the sum of the x and y order is less than or equal to max(orderX, orderY)
config.calibrate.measureApCorr.fitConfig.triangular=True

# maximum Chebyshev function order in x
config.calibrate.measureApCorr.fitConfig.orderX=2

# maximum Chebyshev function order in y
config.calibrate.measureApCorr.fitConfig.orderY=2

# Number of iterations for sigma clipping
config.calibrate.measureApCorr.numIter=4

# Perform PSF fitting?
config.calibrate.doPsf=True

# Maximum allowed shift of WCS, due to matching (pixel)
# 	Valid Range = [-inf,4000)
config.calibrate.astrometry.matcher.maxOffsetPix=750

# Type of source flux; typically one of Ap or Psf
config.calibrate.astrometry.matcher.sourceFluxType='Ap'

# Number of bright stars to use
# 	Valid Range = [2,inf)
config.calibrate.astrometry.matcher.numBrightStars=50

# number of points to define a shape for matching
config.calibrate.astrometry.matcher.numPointsForShape=6

# Allowed non-perpendicularity of x and y (degree)
# 	Valid Range = [-inf,45.0)
config.calibrate.astrometry.matcher.allowedNonperpDeg=3.0

# Maximum separation between reference objects and sources beyond which they will not be considered a match (arcsec)
# 	Valid Range = [0,inf)
config.calibrate.astrometry.matcher.maxMatchDistArcSec=2.0

# Minimum number of matched pairs as a fraction of the smaller of the number of reference stars or the number of good sources; the actual minimum is the smaller of this value or minMatchedPairs
# 	Valid Range = [0,1)
config.calibrate.astrometry.matcher.minFracMatchedPairs=0.3

# maximum determinant of linear transformation matrix for a usable solution
config.calibrate.astrometry.matcher.maxDeterminant=0.02

# Rotation angle allowed between sources and position reference objects (degrees)
# 	Valid Range = [-inf,6.0)
config.calibrate.astrometry.matcher.maxRotationDeg=1.0

# Minimum number of matched pairs; see also minFracMatchedPairs
# 	Valid Range = [2,inf)
config.calibrate.astrometry.matcher.minMatchedPairs=30

# the match distance below which further iteration is pointless (arcsec); ignored if not fitting a WCS
# 	Valid Range = [0,inf)
config.calibrate.astrometry.minMatchDistanceArcSec=0.001

# If True then load reference objects and match sources but do not fit a WCS;  this simply controls whether 'run' calls 'solve' or 'loadAndMatch'
config.calibrate.astrometry.forceKnownWcs=False

# Default reference catalog filter to use if filter not specified in exposure; if blank then filter must be specified in exposure
config.calibrate.astrometry.refObjLoader.defaultFilter=''

# Padding to add to 4 all edges of the bounding box (pixels)
# 	Valid Range = [0,inf)
config.calibrate.astrometry.refObjLoader.pixelMargin=50

# Mapping of camera filter name: reference catalog filter name; each reference filter must exist
config.calibrate.astrometry.refObjLoader.filterMap={'y': 'z'}

# maximum number of iterations of match sources and fit WCSignored if not fitting a WCS
# 	Valid Range = [1,inf)
config.calibrate.astrometry.maxIter=3

# the maximum match distance is set to  mean_match_distance + matchDistanceSigma*std_dev_match_distance; ignored if not fitting a WCS
# 	Valid Range = [0,inf)
config.calibrate.astrometry.matchDistanceSigma=2.0

# Number of standard deviations for clipping level
# 	Valid Range = [0.0,inf)
config.calibrate.astrometry.wcsFitter.rejSigma=3.0

# maximum median scatter of a WCS fit beyond which the fit fails (arcsec); be generous, as this is only intended to catch catastrophic failures
# 	Valid Range = [0,inf)
config.calibrate.astrometry.wcsFitter.maxScatterArcsec=10.0

# number of rejection iterations
# 	Valid Range = [0,inf)
config.calibrate.astrometry.wcsFitter.numRejIter=1

# number of iterations of fitter (which fits X and Y separately, and so benefits from a few iterations
# 	Valid Range = [1,inf)
config.calibrate.astrometry.wcsFitter.numIter=3

# order of SIP polynomial
# 	Valid Range = [0,inf)
config.calibrate.astrometry.wcsFitter.order=3

# Subtract background (after computing it, if not supplied)?
config.calibrate.doBackground=True

# Estimate the background again after final source detection?
config.calibrate.detection.reEstimateBackground=True

# Grow detections by nSigmaToGrow * sigma; if 0 then do not grow
config.calibrate.detection.nSigmaToGrow=2.4

# detected sources with fewer than the specified number of pixels will be ignored
# 	Valid Range = [0,inf)
config.calibrate.detection.minPixels=1

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.calibrate.detection.tempLocalBackground.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.calibrate.detection.tempLocalBackground.approxOrderY=-1

# Names of mask planes to ignore while estimating the background
config.calibrate.detection.tempLocalBackground.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# behaviour if there are too few points in grid for requested interpolation style
# Allowed values:
# 	None	Field is optional
# 	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
# 	THROW_EXCEPTION	throw an exception if there are too few points
# 	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
# 
config.calibrate.detection.tempLocalBackground.undersampleStyle='REDUCE_INTERP_ORDER'

# how to interpolate the background values. This maps to an enum; see afw::math::Background
# Allowed values:
# 	NONE	No background estimation is to be attempted
# 	CONSTANT	Use a single constant value
# 	None	Field is optional
# 	LINEAR	Use linear interpolation
# 	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
# 	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
# 
config.calibrate.detection.tempLocalBackground.algorithm='AKIMA_SPLINE'

# Ignore NaNs when estimating the background
config.calibrate.detection.tempLocalBackground.isNanSafe=False

# how large a region of the sky should be used for each background point
# 	Valid Range = [10,inf)
config.calibrate.detection.tempLocalBackground.binSize=64

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.calibrate.detection.tempLocalBackground.weighting=True

# type of statistic to use for grid points
# Allowed values:
# 	None	Field is optional
# 	MEDIAN	median
# 	MEANCLIP	clipped mean
# 	MEAN	unclipped mean
# 
config.calibrate.detection.tempLocalBackground.statisticsProperty='MEANCLIP'

# Use Approximate (Chebyshev) to model background.
config.calibrate.detection.tempLocalBackground.useApprox=False

# Include threshold relative to thresholdValue
# 	Valid Range = [0.0,inf)
config.calibrate.detection.includeThresholdMultiplier=10.0

# Pixels should be grown as isotropically as possible (slower)
config.calibrate.detection.isotropicGrow=False

# Fiddle factor to add to the background; debugging only
config.calibrate.detection.adjustBackground=0.0

# Do temporary interpolated background subtraction before footprint detection?
config.calibrate.detection.doTempLocalBackground=False

# specifies whether to detect positive, or negative sources, or both
# Allowed values:
# 	positive	detect only positive sources
# 	negative	detect only negative sources
# 	both	detect both positive and negative sources
# 
config.calibrate.detection.thresholdPolarity='positive'

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.calibrate.detection.background.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.calibrate.detection.background.approxOrderY=-1

# Names of mask planes to ignore while estimating the background
config.calibrate.detection.background.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# behaviour if there are too few points in grid for requested interpolation style
# Allowed values:
# 	None	Field is optional
# 	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
# 	THROW_EXCEPTION	throw an exception if there are too few points
# 	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
# 
config.calibrate.detection.background.undersampleStyle='REDUCE_INTERP_ORDER'

# how to interpolate the background values. This maps to an enum; see afw::math::Background
# Allowed values:
# 	NONE	No background estimation is to be attempted
# 	CONSTANT	Use a single constant value
# 	None	Field is optional
# 	LINEAR	Use linear interpolation
# 	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
# 	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
# 
config.calibrate.detection.background.algorithm='NATURAL_SPLINE'

# Ignore NaNs when estimating the background
config.calibrate.detection.background.isNanSafe=False

# how large a region of the sky should be used for each background point
# 	Valid Range = [10,inf)
config.calibrate.detection.background.binSize=128

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.calibrate.detection.background.weighting=True

# type of statistic to use for grid points
# Allowed values:
# 	None	Field is optional
# 	MEDIAN	median
# 	MEANCLIP	clipped mean
# 	MEAN	unclipped mean
# 
config.calibrate.detection.background.statisticsProperty='MEANCLIP'

# Use Approximate (Chebyshev) to model background.
config.calibrate.detection.background.useApprox=False

# Grow detections to set the image mask bits, but return the original (not-grown) footprints
config.calibrate.detection.returnOriginalFootprints=True

# specifies the desired flavor of Threshold
# Allowed values:
# 	pixel_stdev	threshold applied to per-pixel std deviation
# 	variance	threshold applied to image variance
# 	value	threshold applied to image value
# 	stdev	threshold applied to image std deviation
# 
config.calibrate.detection.thresholdType='stdev'

# Threshold for footprints
# 	Valid Range = [0.0,inf)
config.calibrate.detection.thresholdValue=5.0

# Compute photometric zeropoint?
config.calibrate.doPhotoCal=True

# Apply aperture corrections? Silently ignored if endOrder <= lsst.meas.base.APCORR_ORDER when calling run
# Allowed values:
# 	None	Field is optional
# 	yes	apply aperture corrections; fail if data not available
# 	noButWarn	do not apply aperture corrections, but warn if data available (since aperture corrections could have been applied)
# 	yesOrWarn	apply aperture corrections if data available, else warn
# 	no	do not apply aperture corrections
# 
config.calibrate.initialMeasurement.doApplyApCorr='no'

# When measuring, replace other detected footprints with noise?
config.calibrate.initialMeasurement.doReplaceWithNoise=True

# flux measurement algorithms in getApCorrNameSet() to ignore; if a name is listed that does not appear in getApCorrNameSet() then a warning is logged
config.calibrate.initialMeasurement.applyApCorr.ignoreList=[]

# set the general failure flag for a flux when it cannot be aperture-corrected?
config.calibrate.initialMeasurement.applyApCorr.doFlagApCorrFailures=True

# The seed multiplier value to use for random number generation.  0 will not set seed.
config.calibrate.initialMeasurement.noiseReplacer.noiseSeedMultiplier=1

# Add ann offset to the generated noise.
config.calibrate.initialMeasurement.noiseReplacer.noiseOffset=0.0

# How to choose mean and variance of the Gaussian noise we generate?
# Allowed values:
# 	variance	Mean = 0, variance = the image's variance
# 	meta	Mean = 0, variance = the "BGMEAN" metadata entry
# 	measure	Measure clipped mean and variance from the whole image
# 
config.calibrate.initialMeasurement.noiseReplacer.noiseSource='measure'

# whether to run this plugin in single-object mode
config.calibrate.initialMeasurement.plugins['base_CircularApertureFlux'].doMeasure=True

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.calibrate.initialMeasurement.plugins['base_CircularApertureFlux'].shiftKernel='lanczos5'

# Radius (in pixels) of apertures.
config.calibrate.initialMeasurement.plugins['base_CircularApertureFlux'].radii=[3.0, 4.5, 6.0, 9.0, 12.0, 17.0, 25.0, 35.0, 50.0, 70.0]

# Maximum radius (in pixels) for which the sinc algorithm should be used instead of the faster naive algorithm.  For elliptical apertures, this is the minor axis radius.
config.calibrate.initialMeasurement.plugins['base_CircularApertureFlux'].maxSincRadius=10.0

# whether to run this plugin in single-object mode
config.calibrate.initialMeasurement.plugins['base_SkyCoord'].doMeasure=True

# whether to run this plugin in single-object mode
config.calibrate.initialMeasurement.plugins['base_FPPosition'].doMeasure=True

# whether to run this plugin in single-object mode
config.calibrate.initialMeasurement.plugins['base_NaiveCentroid'].doMeasure=True

# Value to subtract from the image pixel values
config.calibrate.initialMeasurement.plugins['base_NaiveCentroid'].background=0.0

# whether to run this plugin in single-object mode
config.calibrate.initialMeasurement.plugins['base_PixelFlags'].doMeasure=True

# List of mask planes to be searched for which occur in the center of a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.calibrate.initialMeasurement.plugins['base_PixelFlags'].masksFpCenter=[]

# List of mask planes to be searched for which occur anywhere within a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.calibrate.initialMeasurement.plugins['base_PixelFlags'].masksFpAnywhere=[]

# whether to run this plugin in single-object mode
config.calibrate.initialMeasurement.plugins['base_InputCount'].doMeasure=True

# whether to run this plugin in single-object mode
config.calibrate.initialMeasurement.plugins['base_ScaledApertureFlux'].doMeasure=True

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.calibrate.initialMeasurement.plugins['base_ScaledApertureFlux'].shiftKernel='lanczos5'

# Scaling factor of PSF FWHM for aperture radius.
config.calibrate.initialMeasurement.plugins['base_ScaledApertureFlux'].scale=3.14

# whether to run this plugin in single-object mode
config.calibrate.initialMeasurement.plugins['base_PsfFlux'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.calibrate.initialMeasurement.plugins['base_PsfFlux'].badMaskPlanes=[]

# whether to run this plugin in single-object mode
config.calibrate.initialMeasurement.plugins['base_Variance'].doMeasure=True

# Scale factor to apply to shape for aperture
config.calibrate.initialMeasurement.plugins['base_Variance'].scale=5.0

# Mask planes to ignore
config.calibrate.initialMeasurement.plugins['base_Variance'].mask=['DETECTED', 'DETECTED_NEGATIVE', 'BAD', 'SAT']

# whether to run this plugin in single-object mode
config.calibrate.initialMeasurement.plugins['base_Jacobian'].doMeasure=True

# Nominal pixel size (arcsec)
config.calibrate.initialMeasurement.plugins['base_Jacobian'].pixelScale=0.5

# whether to run this plugin in single-object mode
config.calibrate.initialMeasurement.plugins['base_GaussianCentroid'].doMeasure=True

# whether to run this plugin in single-object mode
config.calibrate.initialMeasurement.plugins['base_PeakLikelihoodFlux'].doMeasure=True

# Name of warping kernel (e.g. "lanczos4") used to compute the peak
config.calibrate.initialMeasurement.plugins['base_PeakLikelihoodFlux'].warpingKernelName='lanczos4'

# whether to run this plugin in single-object mode
config.calibrate.initialMeasurement.plugins['base_SdssCentroid'].doMeasure=True

# if the peak's less than this insist on binning at least once
config.calibrate.initialMeasurement.plugins['base_SdssCentroid'].peakMin=-1.0

# fiddle factor for adjusting the binning
config.calibrate.initialMeasurement.plugins['base_SdssCentroid'].wfac=1.5

# maximum allowed binning
config.calibrate.initialMeasurement.plugins['base_SdssCentroid'].binmax=16

# whether to run this plugin in single-object mode
config.calibrate.initialMeasurement.plugins['base_PeakCentroid'].doMeasure=True

# whether to run this plugin in single-object mode
config.calibrate.initialMeasurement.plugins['base_GaussianFlux'].doMeasure=True

# FIXME! NEVER DOCUMENTED!
config.calibrate.initialMeasurement.plugins['base_GaussianFlux'].background=0.0

# Maximum centroid shift, limited to 2-10
config.calibrate.initialMeasurement.plugins['base_SdssShape'].maxShift=0.0

# Convergence tolerance for FWHM
config.calibrate.initialMeasurement.plugins['base_SdssShape'].tol2=9.999999747378752e-05

# Convergence tolerance for e1,e2
config.calibrate.initialMeasurement.plugins['base_SdssShape'].tol1=9.999999747378752e-06

# whether to run this plugin in single-object mode
config.calibrate.initialMeasurement.plugins['base_SdssShape'].doMeasure=True

# Additional value to add to background
config.calibrate.initialMeasurement.plugins['base_SdssShape'].background=0.0

# Maximum number of iterations
config.calibrate.initialMeasurement.plugins['base_SdssShape'].maxIter=100

# whether to run this plugin in single-object mode
config.calibrate.initialMeasurement.plugins['base_ClassificationExtendedness'].doMeasure=True

# correction factor for psfFlux error
config.calibrate.initialMeasurement.plugins['base_ClassificationExtendedness'].psfErrFactor=0.0

# correction factor for modelFlux error
config.calibrate.initialMeasurement.plugins['base_ClassificationExtendedness'].modelErrFactor=0.0

# critical ratio of model to psf flux
config.calibrate.initialMeasurement.plugins['base_ClassificationExtendedness'].fluxRatio=0.925

config.calibrate.initialMeasurement.plugins.names=['base_CircularApertureFlux', 'base_NaiveCentroid', 'base_PixelFlags', 'base_SkyCoord', 'base_PsfFlux', 'base_Variance', 'base_GaussianCentroid', 'base_SdssCentroid', 'base_GaussianFlux', 'base_SdssShape']
# the name of the flux measurement algorithm used for calibration
config.calibrate.initialMeasurement.slots.calibFlux='base_CircularApertureFlux_12_0'

# the name of the algorithm used to set the source aperture flux slot
config.calibrate.initialMeasurement.slots.apFlux='base_CircularApertureFlux_3_0'

# the name of the algorithm used to set the source inst flux slot
config.calibrate.initialMeasurement.slots.instFlux='base_GaussianFlux'

# the name of the algorithm used to set source moments parameters
config.calibrate.initialMeasurement.slots.shape='base_SdssShape'

# the name of the centroiding algorithm used to set source x,y
config.calibrate.initialMeasurement.slots.centroid='base_SdssCentroid'

# the name of the algorithm used to set the source model flux slot
config.calibrate.initialMeasurement.slots.modelFlux='base_GaussianFlux'

# the name of the algorithm used to set the source psf flux slot
config.calibrate.initialMeasurement.slots.psfFlux='base_PsfFlux'

# floor for variance is lam*data
config.calibrate.measurePsf.psfDeterminer['pca'].lam=0.05

# for psf candidate evaluation
config.calibrate.measurePsf.psfDeterminer['pca'].reducedChi2ForPsfCandidates=2.0

# Use non-linear fitter for spatial variation of Kernel
config.calibrate.measurePsf.psfDeterminer['pca'].nonLinearSpatialFit=False

# Mask blends in image?
config.calibrate.measurePsf.psfDeterminer['pca'].doMaskBlends=True

# size of cell used to determine PSF (pixels, column direction)
config.calibrate.measurePsf.psfDeterminer['pca'].sizeCellX=256

# size of cell used to determine PSF (pixels, row direction)
config.calibrate.measurePsf.psfDeterminer['pca'].sizeCellY=256

# Rejection threshold (stdev) for candidates based on spatial fit
config.calibrate.measurePsf.psfDeterminer['pca'].spatialReject=3.0

# radius of the kernel to create, relative to the square root of the stellar quadrupole moments
config.calibrate.measurePsf.psfDeterminer['pca'].kernelSize=10.0

# Reject candidates that are blended?
config.calibrate.measurePsf.psfDeterminer['pca'].doRejectBlends=False

# number of iterations of PSF candidate star list
config.calibrate.measurePsf.psfDeterminer['pca'].nIterForPsf=3

# Should each PSF candidate be given the same weight, independent of magnitude?
config.calibrate.measurePsf.psfDeterminer['pca'].constantWeight=True

# Maximum radius of the kernel
config.calibrate.measurePsf.psfDeterminer['pca'].kernelSizeMax=45

# number of stars per psf Cell for spatial fitting
config.calibrate.measurePsf.psfDeterminer['pca'].nStarPerCellSpatialFit=5

# Number of pixels to ignore around the edge of PSF candidate postage stamps
config.calibrate.measurePsf.psfDeterminer['pca'].borderWidth=0

# number of eigen components for PSF kernel creation
config.calibrate.measurePsf.psfDeterminer['pca'].nEigenComponents=4

# Threshold (stdev) for rejecting extraneous pixels around candidate; applied if positive
config.calibrate.measurePsf.psfDeterminer['pca'].pixelThreshold=0.0

# specify spatial order for PSF kernel creation
config.calibrate.measurePsf.psfDeterminer['pca'].spatialOrder=2

# tolerance of spatial fitting
config.calibrate.measurePsf.psfDeterminer['pca'].tolerance=0.01

# Minimum radius of the kernel
config.calibrate.measurePsf.psfDeterminer['pca'].kernelSizeMin=25

# number of stars per psf cell for PSF kernel creation
config.calibrate.measurePsf.psfDeterminer['pca'].nStarPerCell=3

# floor for variance is lam*data
config.calibrate.measurePsf.psfDeterminer['psfex'].lam=0.05

# for psf candidate evaluation
config.calibrate.measurePsf.psfDeterminer['psfex'].reducedChi2ForPsfCandidates=2.0

# specify spatial order for PSF kernel creation
config.calibrate.measurePsf.psfDeterminer['psfex'].spatialOrder=2

# number of eigen components for PSF kernel creation
config.calibrate.measurePsf.psfDeterminer['psfex']._PsfexPsfDeterminerConfig__nEigenComponents=4

# Should PSFEX be permitted to recentroid PSF candidates?
config.calibrate.measurePsf.psfDeterminer['psfex'].recentroid=False

# Resolution of the internal PSF model relative to the pixel size; e.g. 0.5 is equal to 2x oversampling
config.calibrate.measurePsf.psfDeterminer['psfex'].samplingSize=0.5

# size of cell used to determine PSF (pixels, row direction)
config.calibrate.measurePsf.psfDeterminer['psfex'].sizeCellY=256

# List of mask bits which cause a source to be rejected as bad
# N.b. INTRP is used specially in PsfCandidateSet; it means "Contaminated by neighbour"
# 
config.calibrate.measurePsf.psfDeterminer['psfex'].badMaskBits=['INTRP', 'SAT']

# Rejection threshold (stdev) for candidates based on spatial fit
config.calibrate.measurePsf.psfDeterminer['psfex'].spatialReject=3.0

# Number of pixels to ignore around the edge of PSF candidate postage stamps
config.calibrate.measurePsf.psfDeterminer['psfex']._PsfexPsfDeterminerConfig__borderWidth=0

# radius of the kernel to create, relative to the square root of the stellar quadrupole moments
config.calibrate.measurePsf.psfDeterminer['psfex'].kernelSize=81.0

# size of cell used to determine PSF (pixels, column direction)
config.calibrate.measurePsf.psfDeterminer['psfex'].sizeCellX=256

# Should each PSF candidate be given the same weight, independent of magnitude?
config.calibrate.measurePsf.psfDeterminer['psfex']._PsfexPsfDeterminerConfig__constantWeight=True

# Maximum radius of the kernel
config.calibrate.measurePsf.psfDeterminer['psfex'].kernelSizeMax=45

# number of stars per psf Cell for spatial fitting
config.calibrate.measurePsf.psfDeterminer['psfex']._PsfexPsfDeterminerConfig__nStarPerCellSpatialFit=5

# number of stars per psf cell for PSF kernel creation
config.calibrate.measurePsf.psfDeterminer['psfex']._PsfexPsfDeterminerConfig__nStarPerCell=3

# number of iterations of PSF candidate star list
config.calibrate.measurePsf.psfDeterminer['psfex']._PsfexPsfDeterminerConfig__nIterForPsf=3

# tolerance of spatial fitting
config.calibrate.measurePsf.psfDeterminer['psfex'].tolerance=0.01

# Minimum radius of the kernel
config.calibrate.measurePsf.psfDeterminer['psfex'].kernelSizeMin=25

config.calibrate.measurePsf.psfDeterminer.name='psfex'
# maximum width to include in histogram
config.calibrate.measurePsf.starSelector['objectSize'].widthMax=10.0

# Standard deviation of width allowed to be interpreted as good stars
config.calibrate.measurePsf.starSelector['objectSize'].widthStdAllowed=0.15

# minimum width to include in histogram
config.calibrate.measurePsf.starSelector['objectSize'].widthMin=0.9

# size of the Psf kernel to create
config.calibrate.measurePsf.starSelector['objectSize'].kernelSize=21

# specify the minimum psfFlux for good Psf Candidates
config.calibrate.measurePsf.starSelector['objectSize'].fluxMin=4000.0

# number of pixels to ignore around the edge of PSF candidate postage stamps
config.calibrate.measurePsf.starSelector['objectSize'].borderWidth=0

# List of flags which cause a source to be rejected as bad
config.calibrate.measurePsf.starSelector['objectSize'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter', 'base_PixelFlags_flag_bad']

# specify the maximum psfFlux for good Psf Candidates (ignored if == 0)
config.calibrate.measurePsf.starSelector['objectSize'].fluxMax=0.0

# Name of field in Source to use for flux measurement
config.calibrate.measurePsf.starSelector['objectSize'].sourceFluxField='base_PsfFlux_flux'

# Keep objects within this many sigma of cluster 0's median
config.calibrate.measurePsf.starSelector['objectSize'].nSigmaClip=2.0

# specify the minimum psfFlux for good Psf Candidates
config.calibrate.measurePsf.starSelector['catalog'].fluxLim=0.0

# specify the maximum psfFlux for good Psf Candidates (ignored if == 0)
config.calibrate.measurePsf.starSelector['catalog'].fluxMax=0.0

# PSF candidate objects may not have any of these bits set
config.calibrate.measurePsf.starSelector['catalog'].badStarPixelFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter']

# size of the kernel to create
config.calibrate.measurePsf.starSelector['catalog'].kernelSize=21

# number of pixels to ignore around the edge of PSF candidate postage stamps
config.calibrate.measurePsf.starSelector['catalog'].borderWidth=0

# size of the kernel to create
config.calibrate.measurePsf.starSelector['secondMoment'].kernelSize=21

# Multiplier of mean for minimum moments histogram range
config.calibrate.measurePsf.starSelector['secondMoment'].histMomentMinMultiplier=2.0

# Number of bins in moment histogram
config.calibrate.measurePsf.starSelector['secondMoment'].histSize=64

# Clipping threshold for moments histogram range
config.calibrate.measurePsf.starSelector['secondMoment'].histMomentClip=5.0

# number of pixels to ignore around the edge of PSF candidate postage stamps
config.calibrate.measurePsf.starSelector['secondMoment'].borderWidth=0

# List of flags which cause a source to be rejected as bad
config.calibrate.measurePsf.starSelector['secondMoment'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter']

# Multiplier of mean for maximum moments histogram range
config.calibrate.measurePsf.starSelector['secondMoment'].histMomentMaxMultiplier=5.0

# specify the maximum psfFlux for good Psf Candidates (ignored if == 0)
config.calibrate.measurePsf.starSelector['secondMoment'].fluxMax=0.0

# candidate PSF's shapes must lie within this many sigma of the average shape
config.calibrate.measurePsf.starSelector['secondMoment'].clumpNSigma=2.0

# specify the minimum psfFlux for good Psf Candidates
config.calibrate.measurePsf.starSelector['secondMoment'].fluxLim=12500.0

# Maximum moment to consider
config.calibrate.measurePsf.starSelector['secondMoment'].histMomentMax=100.0

# Fraction of objects to use in first pass
config.calibrate.measurePsf.starSelector['sizeMagnitude'].startn1=0.1

# Minimum size to use
config.calibrate.measurePsf.starSelector['sizeMagnitude'].minsize=0.0

# Order of polynomial of fit of size(x,y)
config.calibrate.measurePsf.starSelector['sizeMagnitude'].fitorder=1

# Minimum magnitude to use
config.calibrate.measurePsf.starSelector['sizeMagnitude'].minmag=0.0

# What fraction of objects are likely stars?
config.calibrate.measurePsf.starSelector['sizeMagnitude'].starfrac=0.5

# Maximum magnitude to use
config.calibrate.measurePsf.starSelector['sizeMagnitude'].maxmag=1e+100

# Are sizes already log(size)?
config.calibrate.measurePsf.starSelector['sizeMagnitude'].logsize=False

# Perform size(x,y) fit with fitStars brightest stars
config.calibrate.measurePsf.starSelector['sizeMagnitude'].starsperbin=30

# nSigma to reject a star as an outlier
config.calibrate.measurePsf.starSelector['sizeMagnitude'].fitsigclip=4.0

# Maximum size to use
config.calibrate.measurePsf.starSelector['sizeMagnitude'].maxsize=1e+100

# nSigma to reject a star as an outlier
config.calibrate.measurePsf.starSelector['sizeMagnitude'].aperture=5.0

# Smaller = purer smaple of stars, larger = more stars
config.calibrate.measurePsf.starSelector['sizeMagnitude'].purityratio=0.05

config.calibrate.measurePsf.starSelector.name='objectSize'
# This number will be multplied by the exposure ID to set the random seed for reserving candidates
config.calibrate.measurePsf.reserveSeed=1

# Fraction PSF candidates to reserve from fitting
config.calibrate.measurePsf.reserveFraction=-1.0

# Require photometric calibration to succeed?
config.calibrate.requirePhotoCal=True

# Compute aperture corrections?
config.calibrate.doMeasureApCorr=True

# number of iterations
config.calibrate.photocal.nIter=20

# Write a field name astrom_usedByPhotoCal to the schema
config.calibrate.photocal.doWriteOutput=True

config.calibrate.photocal.colorterms.data={}
config.calibrate.photocal.colorterms.data['hsc*']=lsst.pipe.tasks.colorterms.ColortermDict()
config.calibrate.photocal.colorterms.data['hsc*'].data={}
config.calibrate.photocal.colorterms.data['hsc*'].data['i']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.calibrate.photocal.colorterms.data['hsc*'].data['i'].c2=0.0

# First-order parameter
config.calibrate.photocal.colorterms.data['hsc*'].data['i'].c1=0.0

# Constant parameter
config.calibrate.photocal.colorterms.data['hsc*'].data['i'].c0=0.0

# name of primary filter
config.calibrate.photocal.colorterms.data['hsc*'].data['i'].primary='i'

# name of secondary filter
config.calibrate.photocal.colorterms.data['hsc*'].data['i'].secondary='i'

config.calibrate.photocal.colorterms.data['hsc*'].data['y']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.calibrate.photocal.colorterms.data['hsc*'].data['y'].c2=0.0

# First-order parameter
config.calibrate.photocal.colorterms.data['hsc*'].data['y'].c1=0.0

# Constant parameter
config.calibrate.photocal.colorterms.data['hsc*'].data['y'].c0=0.0

# name of primary filter
config.calibrate.photocal.colorterms.data['hsc*'].data['y'].primary='y'

# name of secondary filter
config.calibrate.photocal.colorterms.data['hsc*'].data['y'].secondary='y'

config.calibrate.photocal.colorterms.data['hsc*'].data['r']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.calibrate.photocal.colorterms.data['hsc*'].data['r'].c2=0.0

# First-order parameter
config.calibrate.photocal.colorterms.data['hsc*'].data['r'].c1=0.0

# Constant parameter
config.calibrate.photocal.colorterms.data['hsc*'].data['r'].c0=0.0

# name of primary filter
config.calibrate.photocal.colorterms.data['hsc*'].data['r'].primary='r'

# name of secondary filter
config.calibrate.photocal.colorterms.data['hsc*'].data['r'].secondary='r'

config.calibrate.photocal.colorterms.data['hsc*'].data['z']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.calibrate.photocal.colorterms.data['hsc*'].data['z'].c2=0.0

# First-order parameter
config.calibrate.photocal.colorterms.data['hsc*'].data['z'].c1=0.0

# Constant parameter
config.calibrate.photocal.colorterms.data['hsc*'].data['z'].c0=0.0

# name of primary filter
config.calibrate.photocal.colorterms.data['hsc*'].data['z'].primary='z'

# name of secondary filter
config.calibrate.photocal.colorterms.data['hsc*'].data['z'].secondary='z'

config.calibrate.photocal.colorterms.data['hsc*'].data['g']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.calibrate.photocal.colorterms.data['hsc*'].data['g'].c2=0.0

# First-order parameter
config.calibrate.photocal.colorterms.data['hsc*'].data['g'].c1=0.0

# Constant parameter
config.calibrate.photocal.colorterms.data['hsc*'].data['g'].c0=0.0

# name of primary filter
config.calibrate.photocal.colorterms.data['hsc*'].data['g'].primary='g'

# name of secondary filter
config.calibrate.photocal.colorterms.data['hsc*'].data['g'].secondary='g'

config.calibrate.photocal.colorterms.data['sdss*']=lsst.pipe.tasks.colorterms.ColortermDict()
config.calibrate.photocal.colorterms.data['sdss*'].data={}
config.calibrate.photocal.colorterms.data['sdss*'].data['g']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.calibrate.photocal.colorterms.data['sdss*'].data['g'].c2=-0.00726883

# First-order parameter
config.calibrate.photocal.colorterms.data['sdss*'].data['g'].c1=-0.08366937

# Constant parameter
config.calibrate.photocal.colorterms.data['sdss*'].data['g'].c0=-0.00816446

# name of primary filter
config.calibrate.photocal.colorterms.data['sdss*'].data['g'].primary='g'

# name of secondary filter
config.calibrate.photocal.colorterms.data['sdss*'].data['g'].secondary='r'

config.calibrate.photocal.colorterms.data['sdss*'].data['N816']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.calibrate.photocal.colorterms.data['sdss*'].data['N816'].c2=-0.05474862

# First-order parameter
config.calibrate.photocal.colorterms.data['sdss*'].data['N816'].c1=-0.63558358

# Constant parameter
config.calibrate.photocal.colorterms.data['sdss*'].data['N816'].c0=0.00927133

# name of primary filter
config.calibrate.photocal.colorterms.data['sdss*'].data['N816'].primary='i'

# name of secondary filter
config.calibrate.photocal.colorterms.data['sdss*'].data['N816'].secondary='z'

config.calibrate.photocal.colorterms.data['sdss*'].data['i']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.calibrate.photocal.colorterms.data['sdss*'].data['i'].c2=-0.01374245

# First-order parameter
config.calibrate.photocal.colorterms.data['sdss*'].data['i'].c1=-0.16922042

# Constant parameter
config.calibrate.photocal.colorterms.data['sdss*'].data['i'].c0=0.00130204

# name of primary filter
config.calibrate.photocal.colorterms.data['sdss*'].data['i'].primary='i'

# name of secondary filter
config.calibrate.photocal.colorterms.data['sdss*'].data['i'].secondary='z'

config.calibrate.photocal.colorterms.data['sdss*'].data['i2']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.calibrate.photocal.colorterms.data['sdss*'].data['i2'].c2=-0.01067212

# First-order parameter
config.calibrate.photocal.colorterms.data['sdss*'].data['i2'].c1=-0.20739606

# Constant parameter
config.calibrate.photocal.colorterms.data['sdss*'].data['i2'].c0=0.00124676

# name of primary filter
config.calibrate.photocal.colorterms.data['sdss*'].data['i2'].primary='i'

# name of secondary filter
config.calibrate.photocal.colorterms.data['sdss*'].data['i2'].secondary='z'

config.calibrate.photocal.colorterms.data['sdss*'].data['r']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.calibrate.photocal.colorterms.data['sdss*'].data['r'].c2=-0.03068248

# First-order parameter
config.calibrate.photocal.colorterms.data['sdss*'].data['r'].c1=0.01284177

# Constant parameter
config.calibrate.photocal.colorterms.data['sdss*'].data['r'].c0=0.0023181

# name of primary filter
config.calibrate.photocal.colorterms.data['sdss*'].data['r'].primary='r'

# name of secondary filter
config.calibrate.photocal.colorterms.data['sdss*'].data['r'].secondary='i'

config.calibrate.photocal.colorterms.data['sdss*'].data['N921']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.calibrate.photocal.colorterms.data['sdss*'].data['N921'].c2=-0.05451118

# First-order parameter
config.calibrate.photocal.colorterms.data['sdss*'].data['N921'].c1=0.0986353

# Constant parameter
config.calibrate.photocal.colorterms.data['sdss*'].data['N921'].c0=0.00752972

# name of primary filter
config.calibrate.photocal.colorterms.data['sdss*'].data['N921'].primary='z'

# name of secondary filter
config.calibrate.photocal.colorterms.data['sdss*'].data['N921'].secondary='i'

config.calibrate.photocal.colorterms.data['sdss*'].data['y']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.calibrate.photocal.colorterms.data['sdss*'].data['y'].c2=0.00574408

# First-order parameter
config.calibrate.photocal.colorterms.data['sdss*'].data['y'].c1=0.35652971

# Constant parameter
config.calibrate.photocal.colorterms.data['sdss*'].data['y'].c0=0.01739708

# name of primary filter
config.calibrate.photocal.colorterms.data['sdss*'].data['y'].primary='z'

# name of secondary filter
config.calibrate.photocal.colorterms.data['sdss*'].data['y'].secondary='i'

config.calibrate.photocal.colorterms.data['sdss*'].data['z']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.calibrate.photocal.colorterms.data['sdss*'].data['z'].c2=0.01479369

# First-order parameter
config.calibrate.photocal.colorterms.data['sdss*'].data['z'].c1=0.01353969

# Constant parameter
config.calibrate.photocal.colorterms.data['sdss*'].data['z'].c0=-0.0068062

# name of primary filter
config.calibrate.photocal.colorterms.data['sdss*'].data['z'].primary='z'

# name of secondary filter
config.calibrate.photocal.colorterms.data['sdss*'].data['z'].secondary='i'

config.calibrate.photocal.colorterms.data['ps1*']=lsst.pipe.tasks.colorterms.ColortermDict()
config.calibrate.photocal.colorterms.data['ps1*'].data={}
config.calibrate.photocal.colorterms.data['ps1*'].data['g']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.calibrate.photocal.colorterms.data['ps1*'].data['g'].c2=-0.0151057

# First-order parameter
config.calibrate.photocal.colorterms.data['ps1*'].data['g'].c1=0.06508481

# Constant parameter
config.calibrate.photocal.colorterms.data['ps1*'].data['g'].c0=0.00730066

# name of primary filter
config.calibrate.photocal.colorterms.data['ps1*'].data['g'].primary='g'

# name of secondary filter
config.calibrate.photocal.colorterms.data['ps1*'].data['g'].secondary='r'

config.calibrate.photocal.colorterms.data['ps1*'].data['N816']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.calibrate.photocal.colorterms.data['ps1*'].data['N816'].c2=-0.10781564

# First-order parameter
config.calibrate.photocal.colorterms.data['ps1*'].data['N816'].c1=-0.68757034

# Constant parameter
config.calibrate.photocal.colorterms.data['ps1*'].data['N816'].c0=0.01191062

# name of primary filter
config.calibrate.photocal.colorterms.data['ps1*'].data['N816'].primary='i'

# name of secondary filter
config.calibrate.photocal.colorterms.data['ps1*'].data['N816'].secondary='z'

config.calibrate.photocal.colorterms.data['ps1*'].data['i']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.calibrate.photocal.colorterms.data['ps1*'].data['i'].c2=-0.03034094

# First-order parameter
config.calibrate.photocal.colorterms.data['ps1*'].data['i'].c1=-0.13944659

# Constant parameter
config.calibrate.photocal.colorterms.data['ps1*'].data['i'].c0=0.00166891

# name of primary filter
config.calibrate.photocal.colorterms.data['ps1*'].data['i'].primary='i'

# name of secondary filter
config.calibrate.photocal.colorterms.data['ps1*'].data['i'].secondary='z'

config.calibrate.photocal.colorterms.data['ps1*'].data['i2']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.calibrate.photocal.colorterms.data['ps1*'].data['i2'].c2=-0.02675511

# First-order parameter
config.calibrate.photocal.colorterms.data['ps1*'].data['i2'].c1=-0.18483562

# Constant parameter
config.calibrate.photocal.colorterms.data['ps1*'].data['i2'].c0=0.00180361

# name of primary filter
config.calibrate.photocal.colorterms.data['ps1*'].data['i2'].primary='i'

# name of secondary filter
config.calibrate.photocal.colorterms.data['ps1*'].data['i2'].secondary='z'

config.calibrate.photocal.colorterms.data['ps1*'].data['r']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.calibrate.photocal.colorterms.data['ps1*'].data['r'].c2=-0.01877566

# First-order parameter
config.calibrate.photocal.colorterms.data['ps1*'].data['r'].c1=0.02093734

# Constant parameter
config.calibrate.photocal.colorterms.data['ps1*'].data['r'].c0=0.00279757

# name of primary filter
config.calibrate.photocal.colorterms.data['ps1*'].data['r'].primary='r'

# name of secondary filter
config.calibrate.photocal.colorterms.data['ps1*'].data['r'].secondary='i'

config.calibrate.photocal.colorterms.data['ps1*'].data['N921']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.calibrate.photocal.colorterms.data['ps1*'].data['N921'].c2=-0.25059679

# First-order parameter
config.calibrate.photocal.colorterms.data['ps1*'].data['N921'].c1=-0.59278367

# Constant parameter
config.calibrate.photocal.colorterms.data['ps1*'].data['N921'].c0=0.00142051

# name of primary filter
config.calibrate.photocal.colorterms.data['ps1*'].data['N921'].primary='z'

# name of secondary filter
config.calibrate.photocal.colorterms.data['ps1*'].data['N921'].secondary='y'

config.calibrate.photocal.colorterms.data['ps1*'].data['y']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.calibrate.photocal.colorterms.data['ps1*'].data['y'].c2=0.02880125

# First-order parameter
config.calibrate.photocal.colorterms.data['ps1*'].data['y'].c1=0.14747401

# Constant parameter
config.calibrate.photocal.colorterms.data['ps1*'].data['y'].c0=-0.00156858

# name of primary filter
config.calibrate.photocal.colorterms.data['ps1*'].data['y'].primary='y'

# name of secondary filter
config.calibrate.photocal.colorterms.data['ps1*'].data['y'].secondary='z'

config.calibrate.photocal.colorterms.data['ps1*'].data['z']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.calibrate.photocal.colorterms.data['ps1*'].data['z'].c2=-0.00316369

# First-order parameter
config.calibrate.photocal.colorterms.data['ps1*'].data['z'].c1=-0.28840221

# Constant parameter
config.calibrate.photocal.colorterms.data['ps1*'].data['z'].c0=-0.00907517

# name of primary filter
config.calibrate.photocal.colorterms.data['ps1*'].data['z'].primary='z'

# name of secondary filter
config.calibrate.photocal.colorterms.data['ps1*'].data['z'].secondary='y'

# Name of photometric reference catalog; used to select a color term dict in colorterms. see also applyColorTerms
config.calibrate.photocal.photoCatName='sdss-dr9-fink-v5b'

# Don't use objects fainter than this magnitude
config.calibrate.photocal.magLimit=22.0

# maximum sigma to use when clipping
config.calibrate.photocal.sigmaMax=0.25

# Name of the source flux field to use.  The associated flag field
# ('<name>_flags') will be implicitly included in badFlags.
config.calibrate.photocal.fluxField='slot_CalibFlux_flux'

# List of source flag fields that must be set for a source to be used.
config.calibrate.photocal.goodFlags=[]

# Additional magnitude uncertainty to be added in quadrature with measurement errors.
# 	Valid Range = [0.0,inf)
config.calibrate.photocal.magErrFloor=0.0

# clip at nSigma
config.calibrate.photocal.nSigma=3.0

# List of source flag fields that will cause a source to be rejected when they are set.
config.calibrate.photocal.badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolated', 'base_PixelFlags_flag_saturated']

# use median instead of mean to compute zeropoint
config.calibrate.photocal.useMedian=True

# Apply photometric color terms to reference stars? One of:
# None: apply if colorterms and photoCatName are not None;
#       fail if color term data is not available for the specified ref catalog and filter.
# True: always apply colorterms; fail if color term data is not available for the
#       specified reference catalog and filter.
# False: do not apply.
config.calibrate.photocal.applyColorTerms=True

# Compute astrometric solution?
config.calibrate.doAstrometry=True

# Apply aperture corrections? Silently ignored if endOrder <= lsst.meas.base.APCORR_ORDER when calling run
# Allowed values:
# 	None	Field is optional
# 	yes	apply aperture corrections; fail if data not available
# 	noButWarn	do not apply aperture corrections, but warn if data available (since aperture corrections could have been applied)
# 	yesOrWarn	apply aperture corrections if data available, else warn
# 	no	do not apply aperture corrections
# 
config.calibrate.measurement.doApplyApCorr='yes'

# When measuring, replace other detected footprints with noise?
config.calibrate.measurement.doReplaceWithNoise=True

# flux measurement algorithms in getApCorrNameSet() to ignore; if a name is listed that does not appear in getApCorrNameSet() then a warning is logged
config.calibrate.measurement.applyApCorr.ignoreList=[]

# set the general failure flag for a flux when it cannot be aperture-corrected?
config.calibrate.measurement.applyApCorr.doFlagApCorrFailures=True

# The seed multiplier value to use for random number generation.  0 will not set seed.
config.calibrate.measurement.noiseReplacer.noiseSeedMultiplier=1

# Add ann offset to the generated noise.
config.calibrate.measurement.noiseReplacer.noiseOffset=0.0

# How to choose mean and variance of the Gaussian noise we generate?
# Allowed values:
# 	variance	Mean = 0, variance = the image's variance
# 	meta	Mean = 0, variance = the "BGMEAN" metadata entry
# 	measure	Measure clipped mean and variance from the whole image
# 
config.calibrate.measurement.noiseReplacer.noiseSource='measure'

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_CircularApertureFlux'].doMeasure=True

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.calibrate.measurement.plugins['base_CircularApertureFlux'].shiftKernel='lanczos5'

# Radius (in pixels) of apertures.
config.calibrate.measurement.plugins['base_CircularApertureFlux'].radii=[3.0, 4.5, 6.0, 9.0, 12.0, 17.0, 25.0, 35.0, 50.0, 70.0]

# Maximum radius (in pixels) for which the sinc algorithm should be used instead of the faster naive algorithm.  For elliptical apertures, this is the minor axis radius.
config.calibrate.measurement.plugins['base_CircularApertureFlux'].maxSincRadius=12.0

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_SkyCoord'].doMeasure=True

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_FPPosition'].doMeasure=True

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_NaiveCentroid'].doMeasure=True

# Value to subtract from the image pixel values
config.calibrate.measurement.plugins['base_NaiveCentroid'].background=0.0

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_PixelFlags'].doMeasure=True

# List of mask planes to be searched for which occur in the center of a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.calibrate.measurement.plugins['base_PixelFlags'].masksFpCenter=[]

# List of mask planes to be searched for which occur anywhere within a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.calibrate.measurement.plugins['base_PixelFlags'].masksFpAnywhere=[]

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_InputCount'].doMeasure=True

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_ScaledApertureFlux'].doMeasure=True

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.calibrate.measurement.plugins['base_ScaledApertureFlux'].shiftKernel='lanczos5'

# Scaling factor of PSF FWHM for aperture radius.
config.calibrate.measurement.plugins['base_ScaledApertureFlux'].scale=3.14

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_PsfFlux'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.calibrate.measurement.plugins['base_PsfFlux'].badMaskPlanes=[]

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_Variance'].doMeasure=True

# Scale factor to apply to shape for aperture
config.calibrate.measurement.plugins['base_Variance'].scale=5.0

# Mask planes to ignore
config.calibrate.measurement.plugins['base_Variance'].mask=['DETECTED', 'DETECTED_NEGATIVE', 'BAD', 'SAT']

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_Jacobian'].doMeasure=True

# Nominal pixel size (arcsec)
config.calibrate.measurement.plugins['base_Jacobian'].pixelScale=0.5

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_GaussianCentroid'].doMeasure=True

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_PeakLikelihoodFlux'].doMeasure=True

# Name of warping kernel (e.g. "lanczos4") used to compute the peak
config.calibrate.measurement.plugins['base_PeakLikelihoodFlux'].warpingKernelName='lanczos4'

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_SdssCentroid'].doMeasure=True

# if the peak's less than this insist on binning at least once
config.calibrate.measurement.plugins['base_SdssCentroid'].peakMin=-1.0

# fiddle factor for adjusting the binning
config.calibrate.measurement.plugins['base_SdssCentroid'].wfac=1.5

# maximum allowed binning
config.calibrate.measurement.plugins['base_SdssCentroid'].binmax=16

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_PeakCentroid'].doMeasure=True

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_GaussianFlux'].doMeasure=True

# FIXME! NEVER DOCUMENTED!
config.calibrate.measurement.plugins['base_GaussianFlux'].background=0.0

# Maximum centroid shift, limited to 2-10
config.calibrate.measurement.plugins['base_SdssShape'].maxShift=0.0

# Convergence tolerance for FWHM
config.calibrate.measurement.plugins['base_SdssShape'].tol2=9.999999747378752e-05

# Convergence tolerance for e1,e2
config.calibrate.measurement.plugins['base_SdssShape'].tol1=9.999999747378752e-06

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_SdssShape'].doMeasure=True

# Additional value to add to background
config.calibrate.measurement.plugins['base_SdssShape'].background=0.0

# Maximum number of iterations
config.calibrate.measurement.plugins['base_SdssShape'].maxIter=100

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_ClassificationExtendedness'].doMeasure=True

# correction factor for psfFlux error
config.calibrate.measurement.plugins['base_ClassificationExtendedness'].psfErrFactor=0.0

# correction factor for modelFlux error
config.calibrate.measurement.plugins['base_ClassificationExtendedness'].modelErrFactor=0.0

# critical ratio of model to psf flux
config.calibrate.measurement.plugins['base_ClassificationExtendedness'].fluxRatio=0.95

config.calibrate.measurement.plugins.names=['base_CircularApertureFlux', 'base_NaiveCentroid', 'base_PixelFlags', 'base_SkyCoord', 'base_PsfFlux', 'base_Variance', 'base_GaussianCentroid', 'base_SdssCentroid', 'base_GaussianFlux', 'base_SdssShape', 'base_ClassificationExtendedness']
# the name of the flux measurement algorithm used for calibration
config.calibrate.measurement.slots.calibFlux='base_CircularApertureFlux_12_0'

# the name of the algorithm used to set the source aperture flux slot
config.calibrate.measurement.slots.apFlux='base_CircularApertureFlux_3_0'

# the name of the algorithm used to set the source inst flux slot
config.calibrate.measurement.slots.instFlux='base_GaussianFlux'

# the name of the algorithm used to set source moments parameters
config.calibrate.measurement.slots.shape='base_SdssShape'

# the name of the centroiding algorithm used to set source x,y
config.calibrate.measurement.slots.centroid='base_SdssCentroid'

# the name of the algorithm used to set the source model flux slot
config.calibrate.measurement.slots.modelFlux='base_GaussianFlux'

# the name of the algorithm used to set the source psf flux slot
config.calibrate.measurement.slots.psfFlux='base_PsfFlux'

# Require astrometry to succeed, if activated?
config.calibrate.requireAstrometry=True

# Apply aperture corrections? Silently ignored if endOrder <= lsst.meas.base.APCORR_ORDER when calling run
# Allowed values:
# 	None	Field is optional
# 	yes	apply aperture corrections; fail if data not available
# 	noButWarn	do not apply aperture corrections, but warn if data available (since aperture corrections could have been applied)
# 	yesOrWarn	apply aperture corrections if data available, else warn
# 	no	do not apply aperture corrections
# 
config.measurement.doApplyApCorr='yes'

# When measuring, replace other detected footprints with noise?
config.measurement.doReplaceWithNoise=True

# flux measurement algorithms in getApCorrNameSet() to ignore; if a name is listed that does not appear in getApCorrNameSet() then a warning is logged
config.measurement.applyApCorr.ignoreList=[]

# set the general failure flag for a flux when it cannot be aperture-corrected?
config.measurement.applyApCorr.doFlagApCorrFailures=True

# The seed multiplier value to use for random number generation.  0 will not set seed.
config.measurement.noiseReplacer.noiseSeedMultiplier=1

# Add ann offset to the generated noise.
config.measurement.noiseReplacer.noiseOffset=0.0

# How to choose mean and variance of the Gaussian noise we generate?
# Allowed values:
# 	variance	Mean = 0, variance = the image's variance
# 	meta	Mean = 0, variance = the "BGMEAN" metadata entry
# 	measure	Measure clipped mean and variance from the whole image
# 
config.measurement.noiseReplacer.noiseSource='measure'

# whether to run this plugin in single-object mode
config.measurement.plugins['base_CircularApertureFlux'].doMeasure=True

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.measurement.plugins['base_CircularApertureFlux'].shiftKernel='lanczos5'

# Radius (in pixels) of apertures.
config.measurement.plugins['base_CircularApertureFlux'].radii=[3.0, 4.5, 6.0, 9.0, 12.0, 17.0, 25.0, 35.0, 50.0, 70.0]

# Maximum radius (in pixels) for which the sinc algorithm should be used instead of the faster naive algorithm.  For elliptical apertures, this is the minor axis radius.
config.measurement.plugins['base_CircularApertureFlux'].maxSincRadius=12.0

# whether to run this plugin in single-object mode
config.measurement.plugins['base_SkyCoord'].doMeasure=True

# whether to run this plugin in single-object mode
config.measurement.plugins['base_FPPosition'].doMeasure=True

# whether to run this plugin in single-object mode
config.measurement.plugins['base_NaiveCentroid'].doMeasure=True

# Value to subtract from the image pixel values
config.measurement.plugins['base_NaiveCentroid'].background=0.0

# whether to run this plugin in single-object mode
config.measurement.plugins['base_PixelFlags'].doMeasure=True

# List of mask planes to be searched for which occur in the center of a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.measurement.plugins['base_PixelFlags'].masksFpCenter=[]

# List of mask planes to be searched for which occur anywhere within a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.measurement.plugins['base_PixelFlags'].masksFpAnywhere=[]

# whether to run this plugin in single-object mode
config.measurement.plugins['base_InputCount'].doMeasure=True

# whether to run this plugin in single-object mode
config.measurement.plugins['base_ScaledApertureFlux'].doMeasure=True

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.measurement.plugins['base_ScaledApertureFlux'].shiftKernel='lanczos5'

# Scaling factor of PSF FWHM for aperture radius.
config.measurement.plugins['base_ScaledApertureFlux'].scale=3.14

# whether to run this plugin in single-object mode
config.measurement.plugins['base_PsfFlux'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.measurement.plugins['base_PsfFlux'].badMaskPlanes=[]

# whether to run this plugin in single-object mode
config.measurement.plugins['base_Variance'].doMeasure=True

# Scale factor to apply to shape for aperture
config.measurement.plugins['base_Variance'].scale=5.0

# Mask planes to ignore
config.measurement.plugins['base_Variance'].mask=['DETECTED', 'DETECTED_NEGATIVE', 'BAD', 'SAT']

# whether to run this plugin in single-object mode
config.measurement.plugins['base_Jacobian'].doMeasure=True

# Nominal pixel size (arcsec)
config.measurement.plugins['base_Jacobian'].pixelScale=0.5

# whether to run this plugin in single-object mode
config.measurement.plugins['base_GaussianCentroid'].doMeasure=True

# whether to run this plugin in single-object mode
config.measurement.plugins['base_PeakLikelihoodFlux'].doMeasure=True

# Name of warping kernel (e.g. "lanczos4") used to compute the peak
config.measurement.plugins['base_PeakLikelihoodFlux'].warpingKernelName='lanczos4'

# whether to run this plugin in single-object mode
config.measurement.plugins['base_SdssCentroid'].doMeasure=True

# if the peak's less than this insist on binning at least once
config.measurement.plugins['base_SdssCentroid'].peakMin=-1.0

# fiddle factor for adjusting the binning
config.measurement.plugins['base_SdssCentroid'].wfac=1.5

# maximum allowed binning
config.measurement.plugins['base_SdssCentroid'].binmax=16

# whether to run this plugin in single-object mode
config.measurement.plugins['base_PeakCentroid'].doMeasure=True

# whether to run this plugin in single-object mode
config.measurement.plugins['base_GaussianFlux'].doMeasure=True

# FIXME! NEVER DOCUMENTED!
config.measurement.plugins['base_GaussianFlux'].background=0.0

# Maximum centroid shift, limited to 2-10
config.measurement.plugins['base_SdssShape'].maxShift=0.0

# Convergence tolerance for FWHM
config.measurement.plugins['base_SdssShape'].tol2=9.999999747378752e-05

# Convergence tolerance for e1,e2
config.measurement.plugins['base_SdssShape'].tol1=9.999999747378752e-06

# whether to run this plugin in single-object mode
config.measurement.plugins['base_SdssShape'].doMeasure=True

# Additional value to add to background
config.measurement.plugins['base_SdssShape'].background=0.0

# Maximum number of iterations
config.measurement.plugins['base_SdssShape'].maxIter=100

# whether to run this plugin in single-object mode
config.measurement.plugins['base_ClassificationExtendedness'].doMeasure=True

# correction factor for psfFlux error
config.measurement.plugins['base_ClassificationExtendedness'].psfErrFactor=0.0

# correction factor for modelFlux error
config.measurement.plugins['base_ClassificationExtendedness'].modelErrFactor=0.0

# critical ratio of model to psf flux
config.measurement.plugins['base_ClassificationExtendedness'].fluxRatio=0.95

config.measurement.plugins.names=['base_CircularApertureFlux', 'base_NaiveCentroid', 'base_PixelFlags', 'base_SkyCoord', 'base_PsfFlux', 'base_Variance', 'base_GaussianCentroid', 'base_SdssCentroid', 'base_GaussianFlux', 'base_SdssShape', 'base_ClassificationExtendedness']
# the name of the flux measurement algorithm used for calibration
config.measurement.slots.calibFlux='base_CircularApertureFlux_12_0'

# the name of the algorithm used to set the source aperture flux slot
config.measurement.slots.apFlux='base_CircularApertureFlux_3_0'

# the name of the algorithm used to set the source inst flux slot
config.measurement.slots.instFlux='base_GaussianFlux'

# the name of the algorithm used to set source moments parameters
config.measurement.slots.shape='base_SdssShape'

# the name of the centroiding algorithm used to set source x,y
config.measurement.slots.centroid='base_SdssCentroid'

# the name of the algorithm used to set the source model flux slot
config.measurement.slots.modelFlux='base_GaussianFlux'

# the name of the algorithm used to set the source psf flux slot
config.measurement.slots.psfFlux='base_PsfFlux'

# Measure sources?
config.doMeasurement=True

# What to do when a peak to be deblended is close to the edge of the image
# Allowed values:
# 	ramp	Ramp down flux at the image edge by the PSF
# 	noclip	Ignore the edge when building the symmetric template.
# 	clip	Clip the template at the edge AND the mirror of the edge.
# 	None	Field is optional
# 
config.deblend.edgeHandling='ramp'

# Assign stray flux to deblend children.  Implies findStrayFlux.
config.deblend.assignStrayFlux=True

# Mask planes with the corresponding limit on the fraction of masked pixels. Sources violating this limit will not be deblended.
config.deblend.maskLimits={'NO_DATA': 0.25}

# How to split flux among peaks
# Allowed values:
# 	trim	Shrink the parent footprint to pixels that are not assigned to children
# 	r-to-peak	~ 1/(1+R^2) to the peak
# 	r-to-footprint	~ 1/(1+R^2) to the closest pixel in the footprint.  CAUTION: this can be computationally expensive on large footprints!
# 	nearest-footprint	Assign 100% to the nearest footprint (using L-1 norm aka Manhattan distance)
# 	None	Field is optional
# 
config.deblend.strayFluxRule='r-to-peak'

# If True, catch exceptions thrown by the deblender, log them, and set a flag on the parent, instead of letting them propagate up
config.deblend.catchFailures=False

# When the deblender should attribute stray flux to point sources
# Allowed values:
# 	always	Always
# 	never	Never; stray flux will not be attributed to any deblended child if the deblender thinks all peaks look like point sources
# 	necessary	When there is not an extended object in the footprint
# 	None	Field is optional
# 
config.deblend.strayFluxToPointSources='necessary'

# Chi-squared per DOF cut for deciding a source is a PSF during deblending (shifted PSF model #2)
config.deblend.psfChisq2b=1.5

# Mask planes to ignore when performing statistics
config.deblend.maskPlanes=['SAT', 'INTRP', 'NO_DATA']

# Maximum area for footprints before they are ignored as large; non-positive means no threshold applied
config.deblend.maxFootprintArea=10000

# Minimum axis ratio for footprints before they are ignored as large; non-positive means no threshold applied
config.deblend.minFootprintAxisRatio=0.0

# Maximum linear dimension for footprints before they are ignored as large; non-positive means no threshold applied
config.deblend.maxFootprintSize=0

# Chi-squared per DOF cut for deciding a source is a PSF during deblending (un-shifted PSF model)
config.deblend.psfChisq1=1.5

# Find stray flux---flux not claimed by any child in the deblender.
config.deblend.findStrayFlux=True

# Footprints smaller in width or height than this value will be ignored; 0 to never ignore.
config.deblend.tinyFootprintSize=2

# Chi-squared per DOF cut for deciding a source is PSF during deblending (shifted PSF model)
config.deblend.psfChisq2=1.5

# Guarantee that all peaks produce a child source.
config.deblend.propagateAllPeaks=False

# Mask name for footprints not deblended, or None
config.deblend.notDeblendedMask='NOT_DEBLENDED'

# Only deblend the brightest maxNumberOfPeaks peaks in the parent (<= 0: unlimited)
config.deblend.maxNumberOfPeaks=0

# When splitting stray flux, clip fractions below this value to zero.
config.deblend.clipStrayFluxFraction=0.01

# Perform calibration?
config.doCalibrate=True

# If True persist background model with background subtracted calexp.  If False persist calexp with the background included.
config.persistBackgroundModel=True

# Write calibration results?
config.doWriteCalibrate=True

