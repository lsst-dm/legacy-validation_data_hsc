import lsst.pipe.tasks.assembleCoadd
assert type(config)==lsst.pipe.tasks.assembleCoadd.SafeClipAssembleCoaddConfig, 'config is of type %s.%s instead of lsst.pipe.tasks.assembleCoadd.SafeClipAssembleCoaddConfig' % (type(config).__module__, type(config).__name__)
# Sigma for outlier rejection; ignored if doSigmaClip false.
config.sigmaClip=1.5

# desired photometric zero point
config.scaleZeroPoint.zeroPoint=27.0

# Automatically select the coadd temp exposure to use as a reference for background matching? Ignored if doMatchBackgrounds false. If False you must specify the reference temp exposure as the data Id
config.autoReference=True

# Match backgrounds of coadd temp exposures before coadding them? If False, the coadd temp expsosures must already have been background subtracted or matched
config.doMatchBackgrounds=False

# Minimum fractional overlap of clipped footprint with visit DETECTED to be clipped
config.minClipFootOverlap=0.65

# Minimum fractional overlap of clipped footprints with visit DETECTED to be clipped when two visits overlap
config.minClipFootOverlapDouble=0.45

# Estimate the background again after final source detection?
config.clipDetection.reEstimateBackground=False

# Grow detections by nSigmaToGrow * sigma; if 0 then do not grow
config.clipDetection.nSigmaToGrow=4.0

# detected sources with fewer than the specified number of pixels will be ignored
# 	Valid Range = [0,inf)
config.clipDetection.minPixels=4

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.clipDetection.tempLocalBackground.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.clipDetection.tempLocalBackground.approxOrderY=-1

# Names of mask planes to ignore while estimating the background
config.clipDetection.tempLocalBackground.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# behaviour if there are too few points in grid for requested interpolation style
# Allowed values:
# 	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
# 	None	Field is optional
# 	THROW_EXCEPTION	throw an exception if there are too few points
# 	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
# 
config.clipDetection.tempLocalBackground.undersampleStyle='REDUCE_INTERP_ORDER'

# how to interpolate the background values. This maps to an enum; see afw::math::Background
# Allowed values:
# 	NONE	No background estimation is to be attempted
# 	CONSTANT	Use a single constant value
# 	None	Field is optional
# 	LINEAR	Use linear interpolation
# 	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
# 	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
# 
config.clipDetection.tempLocalBackground.algorithm='AKIMA_SPLINE'

# Ignore NaNs when estimating the background
config.clipDetection.tempLocalBackground.isNanSafe=False

# how large a region of the sky should be used for each background point
# 	Valid Range = [10,inf)
config.clipDetection.tempLocalBackground.binSize=64

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.clipDetection.tempLocalBackground.weighting=True

# type of statistic to use for grid points
# Allowed values:
# 	None	Field is optional
# 	MEDIAN	median
# 	MEANCLIP	clipped mean
# 	MEAN	unclipped mean
# 
config.clipDetection.tempLocalBackground.statisticsProperty='MEANCLIP'

# Use Approximate (Chebyshev) to model background.
config.clipDetection.tempLocalBackground.useApprox=False

# Include threshold relative to thresholdValue
# 	Valid Range = [0.0,inf)
config.clipDetection.includeThresholdMultiplier=1.0

# Pixels should be grown as isotropically as possible (slower)
config.clipDetection.isotropicGrow=True

# Fiddle factor to add to the background; debugging only
config.clipDetection.adjustBackground=0.0

# Do temporary interpolated background subtraction before footprint detection?
config.clipDetection.doTempLocalBackground=False

# specifies whether to detect positive, or negative sources, or both
# Allowed values:
# 	positive	detect only positive sources
# 	negative	detect only negative sources
# 	both	detect both positive and negative sources
# 
config.clipDetection.thresholdPolarity='both'

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.clipDetection.background.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.clipDetection.background.approxOrderY=-1

# Names of mask planes to ignore while estimating the background
config.clipDetection.background.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# behaviour if there are too few points in grid for requested interpolation style
# Allowed values:
# 	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
# 	None	Field is optional
# 	THROW_EXCEPTION	throw an exception if there are too few points
# 	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
# 
config.clipDetection.background.undersampleStyle='REDUCE_INTERP_ORDER'

# how to interpolate the background values. This maps to an enum; see afw::math::Background
# Allowed values:
# 	NONE	No background estimation is to be attempted
# 	CONSTANT	Use a single constant value
# 	None	Field is optional
# 	LINEAR	Use linear interpolation
# 	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
# 	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
# 
config.clipDetection.background.algorithm='NATURAL_SPLINE'

# Ignore NaNs when estimating the background
config.clipDetection.background.isNanSafe=False

# how large a region of the sky should be used for each background point
# 	Valid Range = [10,inf)
config.clipDetection.background.binSize=256

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.clipDetection.background.weighting=True

# type of statistic to use for grid points
# Allowed values:
# 	None	Field is optional
# 	MEDIAN	median
# 	MEANCLIP	clipped mean
# 	MEAN	unclipped mean
# 
config.clipDetection.background.statisticsProperty='MEANCLIP'

# Use Approximate (Chebyshev) to model background.
config.clipDetection.background.useApprox=False

# Grow detections to set the image mask bits, but return the original (not-grown) footprints
config.clipDetection.returnOriginalFootprints=False

# specifies the desired flavor of Threshold
# Allowed values:
# 	pixel_stdev	threshold applied to per-pixel std deviation
# 	variance	threshold applied to image variance
# 	value	threshold applied to image value
# 	stdev	threshold applied to image std deviation
# 
config.clipDetection.thresholdType='pixel_stdev'

# Threshold for footprints
# 	Valid Range = [0.0,inf)
config.clipDetection.thresholdValue=2.0

# Type of statistic to calculate edge fallbackValue for interpolation
# Allowed values:
# 	None	Field is optional
# 	MEAN	mean
# 	MEDIAN	median
# 	USER	user value set in fallbackUserValue config
# 	MEANCLIP	clipped mean
# 
config.interpImage.fallbackValueType='MEDIAN'

# Add a Gaussian to represent wings?
config.interpImage.modelPsf.addWing=True

# wing amplitude, as a multiple of core amplitude (dimensionless); ignored if addWing false
config.interpImage.modelPsf.wingAmplitude=0.1

# Minimum kernel size if using sizeFactor (pixels); ignored if size is not None
config.interpImage.modelPsf.minSize=5

# Maximum kernel size if using sizeFactor (pixels); ignored if size is not None
config.interpImage.modelPsf.maxSize=None

# Kernel size as a factor of fwhm (dimensionless); size = sizeFactor * fwhm; ignored if size is not None
config.interpImage.modelPsf.sizeFactor=3.0

# Default FWHM of Gaussian model of core of star (pixels)
config.interpImage.modelPsf.defaultFwhm=3.0

# wing width, as a multiple of core width (dimensionless); ignored if addWing false
config.interpImage.modelPsf.wingFwhmFactor=2.5

# Kernel size (width and height) (pixels); if None then sizeFactor is used
config.interpImage.modelPsf.size=None

# If fallbackValueType is 'USER' then use this as the fallbackValue; ignored otherwise
config.interpImage.fallbackUserValue=0.0

# Smoothly taper to the fallback value at the edge of the image?
config.interpImage.useFallbackValueAtEdge=True

# Allow negative values for egde interpolation fallbackValue?  If False, set fallbackValue to max(fallbackValue, 0.0)
config.interpImage.negativeFallbackAllowed=False

# Mask planes that, if set, the associated pixel should not be included in the coaddTempExp.
config.badMaskPlanes=['BAD', 'EDGE', 'SAT', 'INTRP', 'NO_DATA']

# Add a Gaussian to represent wings?
config.modelPsf.addWing=True

# wing amplitude, as a multiple of core amplitude (dimensionless); ignored if addWing false
config.modelPsf.wingAmplitude=0.1

# Minimum kernel size if using sizeFactor (pixels); ignored if size is not None
config.modelPsf.minSize=5

# Maximum kernel size if using sizeFactor (pixels); ignored if size is not None
config.modelPsf.maxSize=None

# Kernel size as a factor of fwhm (dimensionless); size = sizeFactor * fwhm; ignored if size is not None
config.modelPsf.sizeFactor=3.0

# Default FWHM of Gaussian model of core of star (pixels)
config.modelPsf.defaultFwhm=3.0

# wing width, as a multiple of core width (dimensionless); ignored if addWing false
config.modelPsf.wingFwhmFactor=2.5

# Kernel size (width and height) (pixels); if None then sizeFactor is used
config.modelPsf.size=None

# Maximum RMS of residuals of the background offset fit in matchBackgrounds.
config.maxMatchResidualRMS=1.0

# Minimum number of pixels in footprint to use DETECTED mask from the single visits when labeling clipped footprints
config.minBigOverlap=100

# Minimum fractional overlap of clipped footprint with visit DETECTED to be clipped when only one visit overlaps
config.minClipFootOverlapSingle=0.5

# Type of statistic to estimate pixel value for the grid points
# Allowed values:
# 	None	Field is optional
# 	MEDIAN	median
# 	MEANCLIP	clipped mean
# 	MEAN	mean
# 
config.matchBackgrounds.gridStatistic='MEAN'

# Algorithm to interpolate the background values; ignored if usePolynomial is TrueMaps to an enum; see afw.math.Background
# Allowed values:
# 	NONE	No background estimation is to be attempted
# 	CONSTANT	Use a single constant value
# 	None	Field is optional
# 	LINEAR	Use linear interpolation
# 	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
# 	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
# 
config.matchBackgrounds.interpStyle='AKIMA_SPLINE'

# Weight given to mean background level when calculating best reference exposure. Higher weight prefers exposures with low mean background level. Ignored when reference visit is supplied.
# 	Valid Range = [0.0,1.0)
config.matchBackgrounds.bestRefWeightLevel=0.2

# Behaviour if there are too few points in grid for requested interpolation style. Note: INCREASE_NXNYSAMPLE only allowed for usePolynomial=True.
# Allowed values:
# 	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
# 	None	Field is optional
# 	THROW_EXCEPTION	throw an exception if there are too few points
# 	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
# 
config.matchBackgrounds.undersampleStyle='REDUCE_INTERP_ORDER'

# Bin size for gridding the difference image and fitting a spatial model
config.matchBackgrounds.binSize=256

# Order of Chebyshev polynomial background model. Ignored if usePolynomial False
config.matchBackgrounds.order=8

# Names of mask planes to ignore while estimating the background
config.matchBackgrounds.badMaskPlanes=['NO_DATA', 'DETECTED', 'DETECTED_NEGATIVE', 'SAT', 'BAD', 'INTRP', 'CR']

# Fit background difference with Chebychev polynomial interpolation (using afw.math.Approximate)? If False, fit with spline interpolation using afw.math.Background
config.matchBackgrounds.usePolynomial=False

# Sigma for outlier rejection; ignored if gridStatistic != 'MEANCLIP'.
config.matchBackgrounds.numSigmaClip=3

# Tolerance on almost zero standard deviation in a background-offset grid bin. If all bins have a standard deviation below this value, the background offset model is approximated without inverse-variance weighting. (usePolynomial=True)
# 	Valid Range = [0.0,inf)
config.matchBackgrounds.gridStdevEpsilon=1e-08

# Weight given to image variance when calculating best reference exposure. Higher weight prefers exposures with low image variance. Ignored when reference visit is supplied
# 	Valid Range = [0.0,1.0)
config.matchBackgrounds.bestRefWeightVariance=0.4

# Weight given to coverage (number of pixels that overlap with patch), when calculating best reference exposure. Higher weight prefers exposures with high coverage.Ignored when reference visit is supplied
# 	Valid Range = [0.0,1.0)
config.matchBackgrounds.bestRefWeightCoverage=0.4

# Use inverse-variance weighting when approximating background offset model? This will fail when the background offset is constant (this is usually only the case in testing with artificial images).(usePolynomial=True)
config.matchBackgrounds.approxWeighting=True

# Number of iterations of outlier rejection; ignored if gridStatistic != 'MEANCLIP'.
config.matchBackgrounds.numIter=2

# Coadd name: typically one of deep or goodSeeing.
config.coaddName='deep'

# Match to modelPsf?
config.doPsfMatch=False

# Threshold (in fractional weight) of rejection at which we propagate a mask plane to the coadd; that is, we set the mask bit on the coadd if the fraction the rejected frames would have contributed exceeds this value.
config.maskPropagationThresholds={'SAT': 0.1}

# Apply meas_mosaic ubercal results to input calexps?
config.doApplyUberCal=False

# Width, height of stack subregion size; make small enough that a full stack of images will fit into memory at once.
config.subregionSize=[2000, 2000]

# Maximum ratio of the mean squared error of the background matching model to the variance of the difference in backgrounds
config.maxMatchResidualRatio=1.1

# Number of iterations of outlier rejection; ignored if doSigmaClip false.
config.clipIter=3

# Save weights in the CCDs table as well as the visits table? (This is necessary for easy construction of CoaddPsf, but otherwise duplicate information.)
config.inputRecorder.saveCcdWeights=True

# Save the total number of good pixels in each coaddTempExp (redundant with a sum of good pixels in associated CCDs)
config.inputRecorder.saveVisitGoodPix=True

# Add records for CCDs we iterated over but did not add a coaddTempExp due to a lack of unmasked pixels in the coadd footprint.
config.inputRecorder.saveEmptyCcds=False

# Add records for CCDs we iterated over but did not add a coaddTempExp due to an exception (often due to the calexp not being found on disk).
config.inputRecorder.saveErrorCcds=False

# Interpolate over NaN pixels? Also extrapolate, if necessary, but the results are ugly.
config.doInterp=True

# Persist coadd?
config.doWrite=True

# Perform sigma clipped outlier rejection? If False then compute a simple mean.
config.doSigmaClip=True

# Maximum fractional overlap of clipped footprints with visit DETECTED when considering two visits
config.maxClipFootOverlapDouble=0.15

