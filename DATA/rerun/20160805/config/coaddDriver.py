import lsst.pipe.drivers.coaddDriver
assert type(config)==lsst.pipe.drivers.coaddDriver.CoaddDriverConfig, 'config is of type %s.%s instead of lsst.pipe.drivers.coaddDriver.CoaddDriverConfig' % (type(config).__module__, type(config).__name__)
# Estimate the background again after final source detection?
config.detectCoaddSources.detection.reEstimateBackground=True

# Grow detections by nSigmaToGrow * sigma; if 0 then do not grow
config.detectCoaddSources.detection.nSigmaToGrow=2.4

# detected sources with fewer than the specified number of pixels will be ignored
# 	Valid Range = [0,inf)
config.detectCoaddSources.detection.minPixels=1

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.detectCoaddSources.detection.tempLocalBackground.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.detectCoaddSources.detection.tempLocalBackground.approxOrderY=-1

# Names of mask planes to ignore while estimating the background
config.detectCoaddSources.detection.tempLocalBackground.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# behaviour if there are too few points in grid for requested interpolation style
# Allowed values:
# 	None	Field is optional
# 	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
# 	THROW_EXCEPTION	throw an exception if there are too few points
# 	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
# 
config.detectCoaddSources.detection.tempLocalBackground.undersampleStyle='REDUCE_INTERP_ORDER'

# how to interpolate the background values. This maps to an enum; see afw::math::Background
# Allowed values:
# 	NONE	No background estimation is to be attempted
# 	CONSTANT	Use a single constant value
# 	LINEAR	Use linear interpolation
# 	None	Field is optional
# 	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
# 	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
# 
config.detectCoaddSources.detection.tempLocalBackground.algorithm='AKIMA_SPLINE'

# Ignore NaNs when estimating the background
config.detectCoaddSources.detection.tempLocalBackground.isNanSafe=False

# how large a region of the sky should be used for each background point
# 	Valid Range = [1,inf)
config.detectCoaddSources.detection.tempLocalBackground.binSize=64

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.detectCoaddSources.detection.tempLocalBackground.weighting=True

# type of statistic to use for grid points
# Allowed values:
# 	None	Field is optional
# 	MEDIAN	median
# 	MEANCLIP	clipped mean
# 	MEAN	unclipped mean
# 
config.detectCoaddSources.detection.tempLocalBackground.statisticsProperty='MEANCLIP'

# Use Approximate (Chebyshev) to model background.
config.detectCoaddSources.detection.tempLocalBackground.useApprox=False

# Include threshold relative to thresholdValue
# 	Valid Range = [0.0,inf)
config.detectCoaddSources.detection.includeThresholdMultiplier=1.0

# Pixels should be grown as isotropically as possible (slower)
config.detectCoaddSources.detection.isotropicGrow=True

# Fiddle factor to add to the background; debugging only
config.detectCoaddSources.detection.adjustBackground=0.0

# Do temporary interpolated background subtraction before footprint detection?
config.detectCoaddSources.detection.doTempLocalBackground=False

# specifies whether to detect positive, or negative sources, or both
# Allowed values:
# 	positive	detect only positive sources
# 	negative	detect only negative sources
# 	both	detect both positive and negative sources
# 
config.detectCoaddSources.detection.thresholdPolarity='positive'

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.detectCoaddSources.detection.background.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.detectCoaddSources.detection.background.approxOrderY=-1

# Names of mask planes to ignore while estimating the background
config.detectCoaddSources.detection.background.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# behaviour if there are too few points in grid for requested interpolation style
# Allowed values:
# 	None	Field is optional
# 	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
# 	THROW_EXCEPTION	throw an exception if there are too few points
# 	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
# 
config.detectCoaddSources.detection.background.undersampleStyle='REDUCE_INTERP_ORDER'

# how to interpolate the background values. This maps to an enum; see afw::math::Background
# Allowed values:
# 	NONE	No background estimation is to be attempted
# 	CONSTANT	Use a single constant value
# 	LINEAR	Use linear interpolation
# 	None	Field is optional
# 	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
# 	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
# 
config.detectCoaddSources.detection.background.algorithm='NATURAL_SPLINE'

# Ignore NaNs when estimating the background
config.detectCoaddSources.detection.background.isNanSafe=False

# how large a region of the sky should be used for each background point
# 	Valid Range = [1,inf)
config.detectCoaddSources.detection.background.binSize=4096

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.detectCoaddSources.detection.background.weighting=True

# type of statistic to use for grid points
# Allowed values:
# 	None	Field is optional
# 	MEDIAN	median
# 	MEANCLIP	clipped mean
# 	MEAN	unclipped mean
# 
config.detectCoaddSources.detection.background.statisticsProperty='MEANCLIP'

# Use Approximate (Chebyshev) to model background.
config.detectCoaddSources.detection.background.useApprox=False

# Grow detections to set the image mask bits, but return the original (not-grown) footprints
config.detectCoaddSources.detection.returnOriginalFootprints=False

# specifies the desired flavor of Threshold
# Allowed values:
# 	pixel_stdev	threshold applied to per-pixel std deviation
# 	variance	threshold applied to image variance
# 	value	threshold applied to image value
# 	stdev	threshold applied to image std deviation
# 
config.detectCoaddSources.detection.thresholdType='pixel_stdev'

# Threshold for footprints
# 	Valid Range = [0.0,inf)
config.detectCoaddSources.detection.thresholdValue=5.0

# Mask planes for pixels to ignore when scaling variance
config.detectCoaddSources.mask=['DETECTED', 'BAD', 'SAT', 'NO_DATA', 'INTRP']

# Scale variance plane using empirical noise?
config.detectCoaddSources.doScaleVariance=True

# Name of coadd
config.detectCoaddSources.coaddName='deep'

# Sigma for outlier rejection; ignored if doSigmaClip false.
config.assembleCoadd.sigmaClip=1.5

# desired photometric zero point
config.assembleCoadd.scaleZeroPoint.zeroPoint=27.0

# Apply meas_mosaic ubercal results to input calexps?
config.assembleCoadd.doApplyUberCal=False

# Match backgrounds of coadd temp exposures before coadding them? If False, the coadd temp expsosures must already have been background subtracted or matched
config.assembleCoadd.doMatchBackgrounds=False

# Minimum fractional overlap of clipped footprint with visit DETECTED to be clipped
config.assembleCoadd.minClipFootOverlap=0.65

import lsst.pipe.drivers.utils
import lsst.pex.config.config
config.assembleCoadd.select.retarget(target=lsst.pipe.drivers.utils.NullSelectImagesTask, ConfigClass=lsst.pex.config.config.Config)
# Minimum fractional overlap of clipped footprints with visit DETECTED to be clipped when two visits overlap
config.assembleCoadd.minClipFootOverlapDouble=0.45

# Estimate the background again after final source detection?
config.assembleCoadd.clipDetection.reEstimateBackground=False

# Grow detections by nSigmaToGrow * sigma; if 0 then do not grow
config.assembleCoadd.clipDetection.nSigmaToGrow=4.0

# detected sources with fewer than the specified number of pixels will be ignored
# 	Valid Range = [0,inf)
config.assembleCoadd.clipDetection.minPixels=4

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.assembleCoadd.clipDetection.tempLocalBackground.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.assembleCoadd.clipDetection.tempLocalBackground.approxOrderY=-1

# Names of mask planes to ignore while estimating the background
config.assembleCoadd.clipDetection.tempLocalBackground.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# behaviour if there are too few points in grid for requested interpolation style
# Allowed values:
# 	None	Field is optional
# 	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
# 	THROW_EXCEPTION	throw an exception if there are too few points
# 	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
# 
config.assembleCoadd.clipDetection.tempLocalBackground.undersampleStyle='REDUCE_INTERP_ORDER'

# how to interpolate the background values. This maps to an enum; see afw::math::Background
# Allowed values:
# 	NONE	No background estimation is to be attempted
# 	CONSTANT	Use a single constant value
# 	LINEAR	Use linear interpolation
# 	None	Field is optional
# 	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
# 	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
# 
config.assembleCoadd.clipDetection.tempLocalBackground.algorithm='AKIMA_SPLINE'

# Ignore NaNs when estimating the background
config.assembleCoadd.clipDetection.tempLocalBackground.isNanSafe=False

# how large a region of the sky should be used for each background point
# 	Valid Range = [1,inf)
config.assembleCoadd.clipDetection.tempLocalBackground.binSize=64

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.assembleCoadd.clipDetection.tempLocalBackground.weighting=True

# type of statistic to use for grid points
# Allowed values:
# 	None	Field is optional
# 	MEDIAN	median
# 	MEANCLIP	clipped mean
# 	MEAN	unclipped mean
# 
config.assembleCoadd.clipDetection.tempLocalBackground.statisticsProperty='MEANCLIP'

# Use Approximate (Chebyshev) to model background.
config.assembleCoadd.clipDetection.tempLocalBackground.useApprox=False

# Include threshold relative to thresholdValue
# 	Valid Range = [0.0,inf)
config.assembleCoadd.clipDetection.includeThresholdMultiplier=1.0

# Pixels should be grown as isotropically as possible (slower)
config.assembleCoadd.clipDetection.isotropicGrow=True

# Fiddle factor to add to the background; debugging only
config.assembleCoadd.clipDetection.adjustBackground=0.0

# Do temporary interpolated background subtraction before footprint detection?
config.assembleCoadd.clipDetection.doTempLocalBackground=False

# specifies whether to detect positive, or negative sources, or both
# Allowed values:
# 	positive	detect only positive sources
# 	negative	detect only negative sources
# 	both	detect both positive and negative sources
# 
config.assembleCoadd.clipDetection.thresholdPolarity='both'

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.assembleCoadd.clipDetection.background.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.assembleCoadd.clipDetection.background.approxOrderY=-1

# Names of mask planes to ignore while estimating the background
config.assembleCoadd.clipDetection.background.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# behaviour if there are too few points in grid for requested interpolation style
# Allowed values:
# 	None	Field is optional
# 	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
# 	THROW_EXCEPTION	throw an exception if there are too few points
# 	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
# 
config.assembleCoadd.clipDetection.background.undersampleStyle='REDUCE_INTERP_ORDER'

# how to interpolate the background values. This maps to an enum; see afw::math::Background
# Allowed values:
# 	NONE	No background estimation is to be attempted
# 	CONSTANT	Use a single constant value
# 	LINEAR	Use linear interpolation
# 	None	Field is optional
# 	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
# 	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
# 
config.assembleCoadd.clipDetection.background.algorithm='NATURAL_SPLINE'

# Ignore NaNs when estimating the background
config.assembleCoadd.clipDetection.background.isNanSafe=False

# how large a region of the sky should be used for each background point
# 	Valid Range = [1,inf)
config.assembleCoadd.clipDetection.background.binSize=256

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.assembleCoadd.clipDetection.background.weighting=True

# type of statistic to use for grid points
# Allowed values:
# 	None	Field is optional
# 	MEDIAN	median
# 	MEANCLIP	clipped mean
# 	MEAN	unclipped mean
# 
config.assembleCoadd.clipDetection.background.statisticsProperty='MEANCLIP'

# Use Approximate (Chebyshev) to model background.
config.assembleCoadd.clipDetection.background.useApprox=False

# Grow detections to set the image mask bits, but return the original (not-grown) footprints
config.assembleCoadd.clipDetection.returnOriginalFootprints=False

# specifies the desired flavor of Threshold
# Allowed values:
# 	pixel_stdev	threshold applied to per-pixel std deviation
# 	variance	threshold applied to image variance
# 	value	threshold applied to image value
# 	stdev	threshold applied to image std deviation
# 
config.assembleCoadd.clipDetection.thresholdType='pixel_stdev'

# Threshold for footprints
# 	Valid Range = [0.0,inf)
config.assembleCoadd.clipDetection.thresholdValue=2.0

# Type of statistic to calculate edge fallbackValue for interpolation
# Allowed values:
# 	None	Field is optional
# 	MEAN	mean
# 	MEDIAN	median
# 	USER	user value set in fallbackUserValue config
# 	MEANCLIP	clipped mean
# 
config.assembleCoadd.interpImage.fallbackValueType='MEDIAN'

# Add a Gaussian to represent wings?
config.assembleCoadd.interpImage.modelPsf.addWing=True

# wing amplitude, as a multiple of core amplitude (dimensionless); ignored if addWing false
config.assembleCoadd.interpImage.modelPsf.wingAmplitude=0.1

# Minimum kernel size if using sizeFactor (pixels); ignored if size is not None
config.assembleCoadd.interpImage.modelPsf.minSize=5

# Maximum kernel size if using sizeFactor (pixels); ignored if size is not None
config.assembleCoadd.interpImage.modelPsf.maxSize=None

# Kernel size as a factor of fwhm (dimensionless); size = sizeFactor * fwhm; ignored if size is not None
config.assembleCoadd.interpImage.modelPsf.sizeFactor=3.0

# Default FWHM of Gaussian model of core of star (pixels)
config.assembleCoadd.interpImage.modelPsf.defaultFwhm=3.0

# wing width, as a multiple of core width (dimensionless); ignored if addWing false
config.assembleCoadd.interpImage.modelPsf.wingFwhmFactor=2.5

# Kernel size (width and height) (pixels); if None then sizeFactor is used
config.assembleCoadd.interpImage.modelPsf.size=None

# If fallbackValueType is 'USER' then use this as the fallbackValue; ignored otherwise
config.assembleCoadd.interpImage.fallbackUserValue=0.0

# Smoothly taper to the fallback value at the edge of the image?
config.assembleCoadd.interpImage.useFallbackValueAtEdge=True

# Allow negative values for egde interpolation fallbackValue?  If False, set fallbackValue to max(fallbackValue, 0.0)
config.assembleCoadd.interpImage.negativeFallbackAllowed=False

# Mask planes that, if set, the associated pixel should not be included in the coaddTempExp.
config.assembleCoadd.badMaskPlanes=['BAD', 'EDGE', 'SAT', 'INTRP', 'NO_DATA']

# Add a Gaussian to represent wings?
config.assembleCoadd.modelPsf.addWing=True

# wing amplitude, as a multiple of core amplitude (dimensionless); ignored if addWing false
config.assembleCoadd.modelPsf.wingAmplitude=0.1

# Minimum kernel size if using sizeFactor (pixels); ignored if size is not None
config.assembleCoadd.modelPsf.minSize=5

# Maximum kernel size if using sizeFactor (pixels); ignored if size is not None
config.assembleCoadd.modelPsf.maxSize=None

# Kernel size as a factor of fwhm (dimensionless); size = sizeFactor * fwhm; ignored if size is not None
config.assembleCoadd.modelPsf.sizeFactor=3.0

# Default FWHM of Gaussian model of core of star (pixels)
config.assembleCoadd.modelPsf.defaultFwhm=3.0

# wing width, as a multiple of core width (dimensionless); ignored if addWing false
config.assembleCoadd.modelPsf.wingFwhmFactor=2.5

# Kernel size (width and height) (pixels); if None then sizeFactor is used
config.assembleCoadd.modelPsf.size=None

# Name of mask bit used for bright objects
config.assembleCoadd.brightObjectMaskName='BRIGHT_OBJECT'

# Maximum RMS of residuals of the background offset fit in matchBackgrounds.
config.assembleCoadd.maxMatchResidualRMS=1.0

# Minimum number of pixels in footprint to use DETECTED mask from the single visits when labeling clipped footprints
config.assembleCoadd.minBigOverlap=100

# Minimum fractional overlap of clipped footprint with visit DETECTED to be clipped when only one visit overlaps
config.assembleCoadd.minClipFootOverlapSingle=0.5

# Type of statistic to estimate pixel value for the grid points
# Allowed values:
# 	None	Field is optional
# 	MEDIAN	median
# 	MEANCLIP	clipped mean
# 	MEAN	mean
# 
config.assembleCoadd.matchBackgrounds.gridStatistic='MEAN'

# Algorithm to interpolate the background values; ignored if usePolynomial is TrueMaps to an enum; see afw.math.Background
# Allowed values:
# 	NONE	No background estimation is to be attempted
# 	CONSTANT	Use a single constant value
# 	LINEAR	Use linear interpolation
# 	None	Field is optional
# 	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
# 	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
# 
config.assembleCoadd.matchBackgrounds.interpStyle='AKIMA_SPLINE'

# Weight given to mean background level when calculating best reference exposure. Higher weight prefers exposures with low mean background level. Ignored when reference visit is supplied.
# 	Valid Range = [0.0,1.0)
config.assembleCoadd.matchBackgrounds.bestRefWeightLevel=0.2

# Behaviour if there are too few points in grid for requested interpolation style. Note: INCREASE_NXNYSAMPLE only allowed for usePolynomial=True.
# Allowed values:
# 	None	Field is optional
# 	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
# 	THROW_EXCEPTION	throw an exception if there are too few points
# 	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
# 
config.assembleCoadd.matchBackgrounds.undersampleStyle='REDUCE_INTERP_ORDER'

# Bin size for gridding the difference image and fitting a spatial model
config.assembleCoadd.matchBackgrounds.binSize=256

# Order of Chebyshev polynomial background model. Ignored if usePolynomial False
config.assembleCoadd.matchBackgrounds.order=8

# Names of mask planes to ignore while estimating the background
config.assembleCoadd.matchBackgrounds.badMaskPlanes=['NO_DATA', 'DETECTED', 'DETECTED_NEGATIVE', 'SAT', 'BAD', 'INTRP', 'CR']

# Fit background difference with Chebychev polynomial interpolation (using afw.math.Approximate)? If False, fit with spline interpolation using afw.math.Background
config.assembleCoadd.matchBackgrounds.usePolynomial=False

# Sigma for outlier rejection; ignored if gridStatistic != 'MEANCLIP'.
config.assembleCoadd.matchBackgrounds.numSigmaClip=3

# Tolerance on almost zero standard deviation in a background-offset grid bin. If all bins have a standard deviation below this value, the background offset model is approximated without inverse-variance weighting. (usePolynomial=True)
# 	Valid Range = [0.0,inf)
config.assembleCoadd.matchBackgrounds.gridStdevEpsilon=1e-08

# Weight given to image variance when calculating best reference exposure. Higher weight prefers exposures with low image variance. Ignored when reference visit is supplied
# 	Valid Range = [0.0,1.0)
config.assembleCoadd.matchBackgrounds.bestRefWeightVariance=0.4

# Weight given to coverage (number of pixels that overlap with patch), when calculating best reference exposure. Higher weight prefers exposures with high coverage.Ignored when reference visit is supplied
# 	Valid Range = [0.0,1.0)
config.assembleCoadd.matchBackgrounds.bestRefWeightCoverage=0.4

# Use inverse-variance weighting when approximating background offset model? This will fail when the background offset is constant (this is usually only the case in testing with artificial images).(usePolynomial=True)
config.assembleCoadd.matchBackgrounds.approxWeighting=True

# Number of iterations of outlier rejection; ignored if gridStatistic != 'MEANCLIP'.
config.assembleCoadd.matchBackgrounds.numIter=2

# Coadd name: typically one of deep or goodSeeing.
config.assembleCoadd.coaddName='deep'

# Match to modelPsf?
config.assembleCoadd.doPsfMatch=False

# Set mask and flag bits for bright objects?
config.assembleCoadd.doMaskBrightObjects=True

# Threshold (in fractional weight) of rejection at which we propagate a mask plane to the coadd; that is, we set the mask bit on the coadd if the fraction the rejected frames would have contributed exceeds this value.
config.assembleCoadd.maskPropagationThresholds={'SAT': 0.1}

# Automatically select the coadd temp exposure to use as a reference for background matching? Ignored if doMatchBackgrounds false. If False you must specify the reference temp exposure as the data Id
config.assembleCoadd.autoReference=True

# Width, height of stack subregion size; make small enough that a full stack of images will fit into memory at once.
config.assembleCoadd.subregionSize=[10000, 200]

# Maximum ratio of the mean squared error of the background matching model to the variance of the difference in backgrounds
config.assembleCoadd.maxMatchResidualRatio=1.1

# Number of iterations of outlier rejection; ignored if doSigmaClip false.
config.assembleCoadd.clipIter=3

# Save weights in the CCDs table as well as the visits table? (This is necessary for easy construction of CoaddPsf, but otherwise duplicate information.)
config.assembleCoadd.inputRecorder.saveCcdWeights=True

# Save the total number of good pixels in each coaddTempExp (redundant with a sum of good pixels in associated CCDs)
config.assembleCoadd.inputRecorder.saveVisitGoodPix=True

# Add records for CCDs we iterated over but did not add a coaddTempExp due to a lack of unmasked pixels in the coadd footprint.
config.assembleCoadd.inputRecorder.saveEmptyCcds=False

# Add records for CCDs we iterated over but did not add a coaddTempExp due to an exception (often due to the calexp not being found on disk).
config.assembleCoadd.inputRecorder.saveErrorCcds=False

# Interpolate over NaN pixels? Also extrapolate, if necessary, but the results are ugly.
config.assembleCoadd.doInterp=True

# Persist coadd?
config.assembleCoadd.doWrite=True

# Perform sigma clipped outlier rejection? If False then compute a simple mean.
config.assembleCoadd.doSigmaClip=False

# Mask planes to remove before coadding
config.assembleCoadd.removeMaskPlanes=['CROSSTALK', 'NOT_DEBLENDED']

# Maximum fractional overlap of clipped footprints with visit DETECTED when considering two visits
config.assembleCoadd.maxClipFootOverlapDouble=0.15

# Name for coadd
config.coaddName='deep'

# Build background reference?
config.doBackgroundReference=False

# Overwrite coadd?
config.doOverwriteCoadd=False

# overwrite <coaddName>Coadd_tempExp; If False, continue if the file exists on disk
config.makeCoaddTempExp.doOverwrite=False

# Scale kernelSize, alardGaussians by input Fwhm
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].scaleByFwhm=True

# Type of spatial functions for kernel and background
# Allowed values:
# 	chebyshev1	Chebyshev polynomial of the first kind
# 	polynomial	Standard x,y polynomial
# 	None	Field is optional
# 
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].spatialModelType='chebyshev1'

# Calculate kernel and background uncertainties for each kernel candidate?
#                  This comes from the inverse of the covariance matrix.
#                  Warning: regularization can cause problems for this step.
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].calculateKernelUncertainty=False

# Include terms (including kernel cross terms) for background in ip_diffim
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].fitForBackground=False

# Number of Gaussians in AL basis during deconvolution
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].alardNGaussDeconv=3

# Maximum Kernel Size
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].kernelSizeMax=35

# Number of Gaussians in alard-lupton basis
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].alardNGauss=3

# Do sigma clipping on each raw kernel candidate
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].singleKernelClipping=False

# Maximum condition number for a well conditioned spatial matrix
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].maxSpatialConditionNumber=10000000000.0

# Minimum Kernel Size
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].kernelSizeMin=11

# Radius for calculation of stats in 'core' of KernelCandidate diffim.
#                  Total number of pixels used will be (2*radius)**2.
#                  This is used both for 'core' diffim quality as well as ranking of
#                  KernelCandidates by their total flux in this core
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].candidateCoreRadius=3

# Size (rows) in pixels of each SpatialCell for spatial modeling
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].sizeCellX=128

# Size (columns) in pixels of each SpatialCell for spatial modeling
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].sizeCellY=128

# Test for maximum condition number when inverting a kernel matrix.
#                  Anything above maxConditionNumber is not used and the candidate is set as BAD.
#                  Also used to truncate inverse matrix in estimateBiasedRisk.  However,
#                  if you are doing any deconvolution you will want to turn this off, or use
#                  a large maxConditionNumber
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].checkConditionNumber=False

# Minimum Sigma (pixels) for Gaussians during deconvolution;
#         make smaller than alardMinSig as this is only indirectly used
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].alardMinSigDeconv=0.4

# Mask planes to ignore when calculating diffim statistics
#                  Options: NO_DATA EDGE SAT BAD CR INTRP
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].badMaskPlanes=['NO_DATA', 'EDGE', 'SAT']

# Degree of spatial modification of ALL gaussians in AL basis during deconvolution
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].alardDegGaussDeconv=3

# Use Pca to reduce the dimensionality of the kernel basis sets.
#                  This is particularly useful for delta-function kernels.
#                  Functionally, after all Cells have their raw kernels determined, we run
#                  a Pca on these Kernels, re-fit the Cells using the eigenKernels and then
#                  fit those for spatial variation using the same technique as for Alard-Lupton kernels.
#                  If this option is used, the first term will have no spatial variation and the
#                  kernel sum will be conserved.
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].usePcaForSpatialKernel=False

# Subtract off the mean feature before doing the Pca
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].subtractMeanForPca=True

# Number of KernelCandidates in each SpatialCell to use in the spatial fitting
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].nStarPerCell=3

# Rejects KernelCandidates yielding bad difference image quality.
#                  Used by BuildSingleKernelVisitor, AssessSpatialKernelVisitor.
#                  Represents average over pixels of (image/sqrt(variance)).
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].candidateResidualMeanMax=0.25

# Use the core of the footprint for the quality statistics, instead of the entire footprint.
#                  WARNING: if there is deconvolution we probably will need to turn this off
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].useCoreStats=False

# Use constant variance weighting in single kernel fitting?
#                  In some cases this is better for bright star residuals.
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].constantVarianceWeighting=True

# Type of basis set for PSF matching kernel.
# Allowed values:
# 	alard-lupton	Alard-Lupton sum-of-gaussians basis set,
#                            * The first term has no spatial variation
#                            * The kernel sum is conserved
#                            * You may want to turn off 'usePcaForSpatialKernel'
# 	None	Field is optional
# 	delta-function	Delta-function kernel basis set,
#                            * You may enable the option useRegularization
#                            * You should seriously consider usePcaForSpatialKernel, which will also
#                              enable kernel sum conservation for the delta function kernels
# 
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].kernelBasisSet='alard-lupton'

# Sigma in pixels of Gaussians (FWHM = 2.35 sigma).  Must in number equal alardNGauss
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].alardSigGauss=[0.7, 1.5, 3.0]

# Default scale factor between Gaussian sigmas 
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].alardGaussBeta=2.0

# Use Bayesian Information Criterion to select the number of bases going into the kernel
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].useBicForKernelBasis=False

# Number of rows/columns in the convolution kernel; should be odd-valued.
#                  Modified by kernelSizeFwhmScaling if scaleByFwhm = true
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].kernelSize=11

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].afwBackgroundConfig.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].afwBackgroundConfig.approxOrderY=-1

# Names of mask planes to ignore while estimating the background
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].afwBackgroundConfig.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# behaviour if there are too few points in grid for requested interpolation style
# Allowed values:
# 	None	Field is optional
# 	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
# 	THROW_EXCEPTION	throw an exception if there are too few points
# 	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
# 
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].afwBackgroundConfig.undersampleStyle='REDUCE_INTERP_ORDER'

# how to interpolate the background values. This maps to an enum; see afw::math::Background
# Allowed values:
# 	NONE	No background estimation is to be attempted
# 	CONSTANT	Use a single constant value
# 	LINEAR	Use linear interpolation
# 	None	Field is optional
# 	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
# 	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
# 
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].afwBackgroundConfig.algorithm='NATURAL_SPLINE'

# Ignore NaNs when estimating the background
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].afwBackgroundConfig.isNanSafe=False

# how large a region of the sky should be used for each background point
# 	Valid Range = [1,inf)
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].afwBackgroundConfig.binSize=256

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].afwBackgroundConfig.weighting=True

# type of statistic to use for grid points
# Allowed values:
# 	None	Field is optional
# 	MEDIAN	median
# 	MEANCLIP	clipped mean
# 	MEAN	unclipped mean
# 
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].afwBackgroundConfig.statisticsProperty='MEANCLIP'

# Use Approximate (Chebyshev) to model background.
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].afwBackgroundConfig.useApprox=False

# Do sigma clipping after building the spatial model
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].spatialKernelClipping=False

# Maximum allowed sigma for outliers from kernel sum distribution.
#                  Used to reject variable objects from the kernel model
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].maxKsumSigma=3.0

# Spatial order of convolution kernel variation
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].spatialKernelOrder=2

# Minimum number of pixels in an acceptable Footprint
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].detectionConfig.fpNpixMin=5

# Type of detection threshold
# Allowed values:
# 	pixel_stdev	Use stdev derived from variance plane
# 	variance	Use variance of image plane
# 	None	Field is optional
# 	value	Use counts as the detection threshold type
# 	stdev	Use standard deviation of image plane
# 
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].detectionConfig.detThresholdType='pixel_stdev'

# If true run detection on the template (image to convolve);
#                  if false run detection on the science image
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].detectionConfig.detOnTemplate=True

# Value of footprint detection threshold
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].detectionConfig.detThreshold=10.0

# Mask planes that lead to an invalid detection.
#                  Options: NO_DATA EDGE SAT BAD CR INTRP
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].detectionConfig.badMaskPlanes=['NO_DATA', 'EDGE', 'SAT']

# Maximum number of pixels in an acceptable Footprint;
#                  too big and the subsequent convolutions become unwieldy
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].detectionConfig.fpNpixMax=500

# Growing radius (in pixels) for each raw detection
#                  footprint.  The smaller the faster; however the
#                  kernel sum does not converge if the stamp is too
#                  small; and the kernel is not constrained at all if
#                  the stamp is the size of the kernel.  The grown stamp
#                  is 2 * fpGrowPix pixels larger in each dimension.
#                  This is overridden by fpGrowKernelScaling if scaleByFwhm
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].detectionConfig.fpGrowPix=30

# If config.scaleByFwhm, grow the footprint based on
#                  the final kernelSize.  Each footprint will be
#                  2*fpGrowKernelScaling*kernelSize x
#                  2*fpGrowKernelScaling*kernelSize.  With the value
#                  of 1.0, the remaining pixels in each KernelCandiate
#                  after convolution by the basis functions will be
#                  eqaul to the kernel size iteslf.
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].detectionConfig.fpGrowKernelScaling=1.0

# Scale fpGrowPix by input Fwhm?
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].detectionConfig.scaleByFwhm=True

# Maximum number of iterations for rejecting bad KernelCandidates in spatial fitting
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].maxSpatialIterations=3

# How much to scale the kernel size based on the largest AL Sigma
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].kernelSizeFwhmScaling=6.0

# Polynomial order of spatial modification of Gaussians.  Must in number equal alardNGauss
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].alardDegGauss=[4, 2, 2]

# Spatial order of differential background variation
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].spatialBgOrder=1

# Rejects KernelCandidates yielding bad difference image quality.
#                  Used by BuildSingleKernelVisitor, AssessSpatialKernelVisitor.
#                  Represents stddev over pixels of (image/sqrt(variance)).
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].candidateResidualStdMax=1.5

# Do sigma clipping on the ensemble of kernel sums
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].kernelSumClipping=False

# Number of principal components to use for Pca basis, including the
#                  mean kernel if requested.
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].numPrincipalComponents=5

# Warping kernel
# Allowed values:
# 	bilinear	bilinear interpolation
# 	lanczos3	Lanczos kernel of order 3
# 	lanczos4	Lanczos kernel of order 4
# 	lanczos5	Lanczos kernel of order 5
# 	None	Field is optional
# 
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].warpingConfig.warpingKernelName='lanczos3'

# interpLength argument to lsst.afw.math.warpExposure
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].warpingConfig.interpLength=10

# cacheSize argument to lsst.afw.math.SeparableKernel.computeCache
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].warpingConfig.cacheSize=1000000

# use GPU acceleration?
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].warpingConfig.devicePreference=0

# Warping kernel for mask (use warpingKernelName if '')
# Allowed values:
# 		use the regular warping kernel for the mask plane, as well as the image and variance planes
# 	bilinear	bilinear interpolation
# 	lanczos3	Lanczos kernel of order 3
# 	lanczos4	Lanczos kernel of order 4
# 	lanczos5	Lanczos kernel of order 5
# 	None	Field is optional
# 
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].warpingConfig.maskWarpingKernelName='bilinear'

# mask bits to grow to full width of image/variance kernel,
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].warpingConfig.growFullMask=16

# Remake KernelCandidate using better variance estimate after first pass?
#                  Primarily useful when convolving a single-depth image, otherwise not necessary.
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].iterateSingleKernel=False

# Use singular values (SVD) or eigen values (EIGENVALUE) to determine condition number
# Allowed values:
# 	EIGENVALUE	Use eigen values (faster)
# 	SVD	Use singular values
# 	None	Field is optional
# 
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].conditionNumberType='EIGENVALUE'

# Maximum condition number for a well conditioned matrix
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].maxConditionNumber=50000000.0

# Minimum Sigma (pixels) for Gaussians
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].alardMinSig=0.7

# Use afw background subtraction instead of ip_diffim
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].useAfwBackground=False

config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel.name='AL'
# Warping kernel
# Allowed values:
# 	bilinear	bilinear interpolation
# 	lanczos3	Lanczos kernel of order 3
# 	lanczos4	Lanczos kernel of order 4
# 	lanczos5	Lanczos kernel of order 5
# 	None	Field is optional
# 
config.makeCoaddTempExp.warpAndPsfMatch.warp.warpingKernelName='lanczos3'

# interpLength argument to lsst.afw.math.warpExposure
config.makeCoaddTempExp.warpAndPsfMatch.warp.interpLength=10

# cacheSize argument to lsst.afw.math.SeparableKernel.computeCache
config.makeCoaddTempExp.warpAndPsfMatch.warp.cacheSize=1000000

# use GPU acceleration?
config.makeCoaddTempExp.warpAndPsfMatch.warp.devicePreference=0

# Warping kernel for mask (use warpingKernelName if '')
# Allowed values:
# 		use the regular warping kernel for the mask plane, as well as the image and variance planes
# 	bilinear	bilinear interpolation
# 	lanczos3	Lanczos kernel of order 3
# 	lanczos4	Lanczos kernel of order 4
# 	lanczos5	Lanczos kernel of order 5
# 	None	Field is optional
# 
config.makeCoaddTempExp.warpAndPsfMatch.warp.maskWarpingKernelName='bilinear'

# mask bits to grow to full width of image/variance kernel,
config.makeCoaddTempExp.warpAndPsfMatch.warp.growFullMask=16

# Save weights in the CCDs table as well as the visits table? (This is necessary for easy construction of CoaddPsf, but otherwise duplicate information.)
config.makeCoaddTempExp.inputRecorder.saveCcdWeights=True

# Save the total number of good pixels in each coaddTempExp (redundant with a sum of good pixels in associated CCDs)
config.makeCoaddTempExp.inputRecorder.saveVisitGoodPix=True

# Add records for CCDs we iterated over but did not add a coaddTempExp due to a lack of unmasked pixels in the coadd footprint.
config.makeCoaddTempExp.inputRecorder.saveEmptyCcds=False

# Add records for CCDs we iterated over but did not add a coaddTempExp due to an exception (often due to the calexp not being found on disk).
config.makeCoaddTempExp.inputRecorder.saveErrorCcds=False

# Coadd name: typically one of deep or goodSeeing.
config.makeCoaddTempExp.coaddName='deep'

# Apply meas_mosaic ubercal results to input calexps?
config.makeCoaddTempExp.doApplyUberCal=False

# Add a Gaussian to represent wings?
config.makeCoaddTempExp.modelPsf.addWing=True

# wing amplitude, as a multiple of core amplitude (dimensionless); ignored if addWing false
config.makeCoaddTempExp.modelPsf.wingAmplitude=0.1

# Minimum kernel size if using sizeFactor (pixels); ignored if size is not None
config.makeCoaddTempExp.modelPsf.minSize=5

# Maximum kernel size if using sizeFactor (pixels); ignored if size is not None
config.makeCoaddTempExp.modelPsf.maxSize=None

# Kernel size as a factor of fwhm (dimensionless); size = sizeFactor * fwhm; ignored if size is not None
config.makeCoaddTempExp.modelPsf.sizeFactor=3.0

# Default FWHM of Gaussian model of core of star (pixels)
config.makeCoaddTempExp.modelPsf.defaultFwhm=3.0

# wing width, as a multiple of core width (dimensionless); ignored if addWing false
config.makeCoaddTempExp.modelPsf.wingFwhmFactor=2.5

# Kernel size (width and height) (pixels); if None then sizeFactor is used
config.makeCoaddTempExp.modelPsf.size=None

# Match to modelPsf?
config.makeCoaddTempExp.doPsfMatch=False

import lsst.pipe.drivers.utils
import lsst.pex.config.config
config.makeCoaddTempExp.select.retarget(target=lsst.pipe.drivers.utils.NullSelectImagesTask, ConfigClass=lsst.pex.config.config.Config)
# persist <coaddName>Coadd_tempExp
config.makeCoaddTempExp.doWrite=True

# Mask planes that, if set, the associated pixel should not be included in the coaddTempExp.
config.makeCoaddTempExp.badMaskPlanes=['NO_DATA']

# Work with a background subtracted calexp?
config.makeCoaddTempExp.bgSubtracted=True

