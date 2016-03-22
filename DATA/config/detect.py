import lsst.pipe.tasks.multiBand
assert type(config)==lsst.pipe.tasks.multiBand.DetectCoaddSourcesConfig, 'config is of type %s.%s instead of lsst.pipe.tasks.multiBand.DetectCoaddSourcesConfig' % (type(config).__module__, type(config).__name__)
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
# 	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
# 	None	Field is optional
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
# 	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
# 	None	Field is optional
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
config.detection.background.binSize=256

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
config.detection.thresholdType='pixel_stdev'

# Threshold for footprints
# 	Valid Range = [0.0,inf)
config.detection.thresholdValue=5.0

# Mask planes for pixels to ignore when scaling variance
config.mask=['DETECTED', 'BAD', 'SAT', 'NO_DATA', 'INTRP']

# Scale variance plane using empirical noise?
config.doScaleVariance=True

# Name of coadd
config.coaddName='deep'

