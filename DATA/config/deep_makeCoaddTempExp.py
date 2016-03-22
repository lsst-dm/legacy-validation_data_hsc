import lsst.pipe.tasks.makeCoaddTempExp
assert type(config)==lsst.pipe.tasks.makeCoaddTempExp.MakeCoaddTempExpConfig, 'config is of type %s.%s instead of lsst.pipe.tasks.makeCoaddTempExp.MakeCoaddTempExpConfig' % (type(config).__module__, type(config).__name__)
# overwrite <coaddName>Coadd_tempExp; If False, continue if the file exists on disk
config.doOverwrite=True

# Scale kernelSize, alardGaussians by input Fwhm
config.warpAndPsfMatch.psfMatch.kernel['AL'].scaleByFwhm=True

# Type of spatial functions for kernel and background
# Allowed values:
# 	chebyshev1	Chebyshev polynomial of the first kind
# 	polynomial	Standard x,y polynomial
# 	None	Field is optional
# 
config.warpAndPsfMatch.psfMatch.kernel['AL'].spatialModelType='chebyshev1'

# Calculate kernel and background uncertainties for each kernel candidate?
#                  This comes from the inverse of the covariance matrix.
#                  Warning: regularization can cause problems for this step.
config.warpAndPsfMatch.psfMatch.kernel['AL'].calculateKernelUncertainty=False

# Include terms (including kernel cross terms) for background in ip_diffim
config.warpAndPsfMatch.psfMatch.kernel['AL'].fitForBackground=False

# Number of Gaussians in AL basis during deconvolution
config.warpAndPsfMatch.psfMatch.kernel['AL'].alardNGaussDeconv=3

# Maximum Kernel Size
config.warpAndPsfMatch.psfMatch.kernel['AL'].kernelSizeMax=35

# Number of Gaussians in alard-lupton basis
config.warpAndPsfMatch.psfMatch.kernel['AL'].alardNGauss=3

# Do sigma clipping on each raw kernel candidate
config.warpAndPsfMatch.psfMatch.kernel['AL'].singleKernelClipping=False

# Maximum condition number for a well conditioned spatial matrix
config.warpAndPsfMatch.psfMatch.kernel['AL'].maxSpatialConditionNumber=10000000000.0

# Minimum Kernel Size
config.warpAndPsfMatch.psfMatch.kernel['AL'].kernelSizeMin=11

# Radius for calculation of stats in 'core' of KernelCandidate diffim.
#                  Total number of pixels used will be (2*radius)**2.
#                  This is used both for 'core' diffim quality as well as ranking of
#                  KernelCandidates by their total flux in this core
config.warpAndPsfMatch.psfMatch.kernel['AL'].candidateCoreRadius=3

# Size (rows) in pixels of each SpatialCell for spatial modeling
config.warpAndPsfMatch.psfMatch.kernel['AL'].sizeCellX=128

# Size (columns) in pixels of each SpatialCell for spatial modeling
config.warpAndPsfMatch.psfMatch.kernel['AL'].sizeCellY=128

# Test for maximum condition number when inverting a kernel matrix.
#                  Anything above maxConditionNumber is not used and the candidate is set as BAD.
#                  Also used to truncate inverse matrix in estimateBiasedRisk.  However,
#                  if you are doing any deconvolution you will want to turn this off, or use
#                  a large maxConditionNumber
config.warpAndPsfMatch.psfMatch.kernel['AL'].checkConditionNumber=False

# Minimum Sigma (pixels) for Gaussians during deconvolution; 
#         make smaller than alardMinSig as this is only indirectly used
config.warpAndPsfMatch.psfMatch.kernel['AL'].alardMinSigDeconv=0.4

# Mask planes to ignore when calculating diffim statistics
#                  Options: NO_DATA EDGE SAT BAD CR INTRP
config.warpAndPsfMatch.psfMatch.kernel['AL'].badMaskPlanes=['NO_DATA', 'EDGE', 'SAT']

# Degree of spatial modification of ALL gaussians in AL basis during deconvolution
config.warpAndPsfMatch.psfMatch.kernel['AL'].alardDegGaussDeconv=3

# Use Pca to reduce the dimensionality of the kernel basis sets.
#                  This is particularly useful for delta-function kernels.
#                  Functionally, after all Cells have their raw kernels determined, we run
#                  a Pca on these Kernels, re-fit the Cells using the eigenKernels and then
#                  fit those for spatial variation using the same technique as for Alard-Lupton kernels.
#                  If this option is used, the first term will have no spatial variation and the
#                  kernel sum will be conserved.
config.warpAndPsfMatch.psfMatch.kernel['AL'].usePcaForSpatialKernel=False

# Subtract off the mean feature before doing the Pca
config.warpAndPsfMatch.psfMatch.kernel['AL'].subtractMeanForPca=True

# Number of KernelCandidates in each SpatialCell to use in the spatial fitting
config.warpAndPsfMatch.psfMatch.kernel['AL'].nStarPerCell=3

# Rejects KernelCandidates yielding bad difference image quality.
#                  Used by BuildSingleKernelVisitor, AssessSpatialKernelVisitor.
#                  Represents average over pixels of (image/sqrt(variance)).
config.warpAndPsfMatch.psfMatch.kernel['AL'].candidateResidualMeanMax=0.25

# Use the core of the footprint for the quality statistics, instead of the entire footprint.
#                  WARNING: if there is deconvolution we probably will need to turn this off
config.warpAndPsfMatch.psfMatch.kernel['AL'].useCoreStats=False

# Use constant variance weighting in single kernel fitting?
#                  In some cases this is better for bright star residuals.
config.warpAndPsfMatch.psfMatch.kernel['AL'].constantVarianceWeighting=True

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
config.warpAndPsfMatch.psfMatch.kernel['AL'].kernelBasisSet='alard-lupton'

# Sigma in pixels of Gaussians (FWHM = 2.35 sigma).  Must in number equal alardNGauss
config.warpAndPsfMatch.psfMatch.kernel['AL'].alardSigGauss=[0.7, 1.5, 3.0]

# Default scale factor between Gaussian sigmas 
config.warpAndPsfMatch.psfMatch.kernel['AL'].alardGaussBeta=2.0

# Use Bayesian Information Criterion to select the number of bases going into the kernel
config.warpAndPsfMatch.psfMatch.kernel['AL'].useBicForKernelBasis=False

# Number of rows/columns in the convolution kernel; should be odd-valued.
#                  Modified by kernelSizeFwhmScaling if scaleByFwhm = true
config.warpAndPsfMatch.psfMatch.kernel['AL'].kernelSize=11

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.warpAndPsfMatch.psfMatch.kernel['AL'].afwBackgroundConfig.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.warpAndPsfMatch.psfMatch.kernel['AL'].afwBackgroundConfig.approxOrderY=-1

# Names of mask planes to ignore while estimating the background
config.warpAndPsfMatch.psfMatch.kernel['AL'].afwBackgroundConfig.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# behaviour if there are too few points in grid for requested interpolation style
# Allowed values:
# 	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
# 	THROW_EXCEPTION	throw an exception if there are too few points
# 	None	Field is optional
# 	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
# 
config.warpAndPsfMatch.psfMatch.kernel['AL'].afwBackgroundConfig.undersampleStyle='REDUCE_INTERP_ORDER'

# how to interpolate the background values. This maps to an enum; see afw::math::Background
# Allowed values:
# 	NONE	No background estimation is to be attempted
# 	CONSTANT	Use a single constant value
# 	None	Field is optional
# 	LINEAR	Use linear interpolation
# 	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
# 	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
# 
config.warpAndPsfMatch.psfMatch.kernel['AL'].afwBackgroundConfig.algorithm='NATURAL_SPLINE'

# Ignore NaNs when estimating the background
config.warpAndPsfMatch.psfMatch.kernel['AL'].afwBackgroundConfig.isNanSafe=False

# how large a region of the sky should be used for each background point
# 	Valid Range = [10,inf)
config.warpAndPsfMatch.psfMatch.kernel['AL'].afwBackgroundConfig.binSize=256

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.warpAndPsfMatch.psfMatch.kernel['AL'].afwBackgroundConfig.weighting=True

# type of statistic to use for grid points
# Allowed values:
# 	None	Field is optional
# 	MEDIAN	median
# 	MEANCLIP	clipped mean
# 	MEAN	unclipped mean
# 
config.warpAndPsfMatch.psfMatch.kernel['AL'].afwBackgroundConfig.statisticsProperty='MEANCLIP'

# Use Approximate (Chebyshev) to model background.
config.warpAndPsfMatch.psfMatch.kernel['AL'].afwBackgroundConfig.useApprox=False

# Do sigma clipping after building the spatial model
config.warpAndPsfMatch.psfMatch.kernel['AL'].spatialKernelClipping=False

# Maximum allowed sigma for outliers from kernel sum distribution.
#                  Used to reject variable objects from the kernel model
config.warpAndPsfMatch.psfMatch.kernel['AL'].maxKsumSigma=3.0

# Spatial order of convolution kernel variation
config.warpAndPsfMatch.psfMatch.kernel['AL'].spatialKernelOrder=2

# Minimum number of pixels in an acceptable Footprint
config.warpAndPsfMatch.psfMatch.kernel['AL'].detectionConfig.fpNpixMin=5

# Type of detection threshold
# Allowed values:
# 	pixel_stdev	Use stdev derived from variance plane
# 	variance	Use variance of image plane
# 	None	Field is optional
# 	value	Use counts as the detection threshold type
# 	stdev	Use standard deviation of image plane
# 
config.warpAndPsfMatch.psfMatch.kernel['AL'].detectionConfig.detThresholdType='pixel_stdev'

# If true run detection on the template (image to convolve);
#                  if false run detection on the science image
config.warpAndPsfMatch.psfMatch.kernel['AL'].detectionConfig.detOnTemplate=True

# Value of footprint detection threshold
config.warpAndPsfMatch.psfMatch.kernel['AL'].detectionConfig.detThreshold=10.0

# Mask planes that lead to an invalid detection.
#                  Options: NO_DATA EDGE SAT BAD CR INTRP
config.warpAndPsfMatch.psfMatch.kernel['AL'].detectionConfig.badMaskPlanes=['NO_DATA', 'EDGE', 'SAT']

# Maximum number of pixels in an acceptable Footprint;
#                  too big and the subsequent convolutions become unwieldy
config.warpAndPsfMatch.psfMatch.kernel['AL'].detectionConfig.fpNpixMax=500

# Growing radius (in pixels) for each raw detection
#                  footprint.  The smaller the faster; however the
#                  kernel sum does not converge if the stamp is too
#                  small; and the kernel is not constrained at all if
#                  the stamp is the size of the kernel.  The grown stamp
#                  is 2 * fpGrowPix pixels larger in each dimension.
#                  This is overridden by fpGrowKernelScaling if scaleByFwhm
config.warpAndPsfMatch.psfMatch.kernel['AL'].detectionConfig.fpGrowPix=30

# If config.scaleByFwhm, grow the footprint based on
#                  the final kernelSize.  Each footprint will be
#                  2*fpGrowKernelScaling*kernelSize x
#                  2*fpGrowKernelScaling*kernelSize.  With the value
#                  of 1.0, the remaining pixels in each KernelCandiate
#                  after convolution by the basis functions will be
#                  eqaul to the kernel size iteslf.
config.warpAndPsfMatch.psfMatch.kernel['AL'].detectionConfig.fpGrowKernelScaling=1.0

# Scale fpGrowPix by input Fwhm?
config.warpAndPsfMatch.psfMatch.kernel['AL'].detectionConfig.scaleByFwhm=True

# Maximum number of iterations for rejecting bad KernelCandidates in spatial fitting
config.warpAndPsfMatch.psfMatch.kernel['AL'].maxSpatialIterations=3

# How much to scale the kernel size based on the largest AL Sigma
config.warpAndPsfMatch.psfMatch.kernel['AL'].kernelSizeFwhmScaling=6.0

# Polynomial order of spatial modification of Gaussians.  Must in number equal alardNGauss
config.warpAndPsfMatch.psfMatch.kernel['AL'].alardDegGauss=[4, 2, 2]

# Spatial order of differential background variation
config.warpAndPsfMatch.psfMatch.kernel['AL'].spatialBgOrder=1

# Rejects KernelCandidates yielding bad difference image quality.
#                  Used by BuildSingleKernelVisitor, AssessSpatialKernelVisitor.
#                  Represents stddev over pixels of (image/sqrt(variance)).
config.warpAndPsfMatch.psfMatch.kernel['AL'].candidateResidualStdMax=1.5

# Do sigma clipping on the ensemble of kernel sums
config.warpAndPsfMatch.psfMatch.kernel['AL'].kernelSumClipping=False

# Number of principal components to use for Pca basis, including the
#                  mean kernel if requested.
config.warpAndPsfMatch.psfMatch.kernel['AL'].numPrincipalComponents=5

# Warping kernel
# Allowed values:
# 	bilinear	bilinear interpolation
# 	lanczos3	Lanczos kernel of order 3
# 	lanczos4	Lanczos kernel of order 4
# 	lanczos5	Lanczos kernel of order 5
# 	None	Field is optional
# 
config.warpAndPsfMatch.psfMatch.kernel['AL'].warpingConfig.warpingKernelName='lanczos3'

# interpLength argument to lsst.afw.math.warpExposure
config.warpAndPsfMatch.psfMatch.kernel['AL'].warpingConfig.interpLength=10

# cacheSize argument to lsst.afw.math.SeparableKernel.computeCache
config.warpAndPsfMatch.psfMatch.kernel['AL'].warpingConfig.cacheSize=1000000

# use GPU acceleration?
config.warpAndPsfMatch.psfMatch.kernel['AL'].warpingConfig.devicePreference=0

# Warping kernel for mask (use warpingKernelName if '')
# Allowed values:
# 		use the regular warping kernel for the mask plane, as well as the image and variance planes
# 	bilinear	bilinear interpolation
# 	None	Field is optional
# 	lanczos3	Lanczos kernel of order 3
# 	lanczos4	Lanczos kernel of order 4
# 	lanczos5	Lanczos kernel of order 5
# 
config.warpAndPsfMatch.psfMatch.kernel['AL'].warpingConfig.maskWarpingKernelName='bilinear'

# mask bits to grow to full width of image/variance kernel,
config.warpAndPsfMatch.psfMatch.kernel['AL'].warpingConfig.growFullMask=16

# Remake KernelCandidate using better variance estimate after first pass?
#                  Primarily useful when convolving a single-depth image, otherwise not necessary.
config.warpAndPsfMatch.psfMatch.kernel['AL'].iterateSingleKernel=False

# Use singular values (SVD) or eigen values (EIGENVALUE) to determine condition number
# Allowed values:
# 	EIGENVALUE	Use eigen values (faster)
# 	SVD	Use singular values
# 	None	Field is optional
# 
config.warpAndPsfMatch.psfMatch.kernel['AL'].conditionNumberType='EIGENVALUE'

# Maximum condition number for a well conditioned matrix
config.warpAndPsfMatch.psfMatch.kernel['AL'].maxConditionNumber=50000000.0

# Minimum Sigma (pixels) for Gaussians
config.warpAndPsfMatch.psfMatch.kernel['AL'].alardMinSig=0.7

# Use afw background subtraction instead of ip_diffim
config.warpAndPsfMatch.psfMatch.kernel['AL'].useAfwBackground=False

config.warpAndPsfMatch.psfMatch.kernel.name='AL'
# Warping kernel
# Allowed values:
# 	bilinear	bilinear interpolation
# 	lanczos3	Lanczos kernel of order 3
# 	lanczos4	Lanczos kernel of order 4
# 	lanczos5	Lanczos kernel of order 5
# 	None	Field is optional
# 
config.warpAndPsfMatch.warp.warpingKernelName='lanczos3'

# interpLength argument to lsst.afw.math.warpExposure
config.warpAndPsfMatch.warp.interpLength=10

# cacheSize argument to lsst.afw.math.SeparableKernel.computeCache
config.warpAndPsfMatch.warp.cacheSize=1000000

# use GPU acceleration?
config.warpAndPsfMatch.warp.devicePreference=0

# Warping kernel for mask (use warpingKernelName if '')
# Allowed values:
# 		use the regular warping kernel for the mask plane, as well as the image and variance planes
# 	bilinear	bilinear interpolation
# 	None	Field is optional
# 	lanczos3	Lanczos kernel of order 3
# 	lanczos4	Lanczos kernel of order 4
# 	lanczos5	Lanczos kernel of order 5
# 
config.warpAndPsfMatch.warp.maskWarpingKernelName='bilinear'

# mask bits to grow to full width of image/variance kernel,
config.warpAndPsfMatch.warp.growFullMask=16

# Save weights in the CCDs table as well as the visits table? (This is necessary for easy construction of CoaddPsf, but otherwise duplicate information.)
config.inputRecorder.saveCcdWeights=True

# Save the total number of good pixels in each coaddTempExp (redundant with a sum of good pixels in associated CCDs)
config.inputRecorder.saveVisitGoodPix=True

# Add records for CCDs we iterated over but did not add a coaddTempExp due to a lack of unmasked pixels in the coadd footprint.
config.inputRecorder.saveEmptyCcds=False

# Add records for CCDs we iterated over but did not add a coaddTempExp due to an exception (often due to the calexp not being found on disk).
config.inputRecorder.saveErrorCcds=False

# Coadd name: typically one of deep or goodSeeing.
config.coaddName='deep'

# Apply meas_mosaic ubercal results to input calexps?
config.doApplyUberCal=False

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

# Match to modelPsf?
config.doPsfMatch=False

# persist <coaddName>Coadd_tempExp
config.doWrite=True

# Mask planes that, if set, the associated pixel should not be included in the coaddTempExp.
config.badMaskPlanes=['NO_DATA']

# Work with a background subtracted calexp?
config.bgSubtracted=True

