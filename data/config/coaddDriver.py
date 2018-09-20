import lsst.pipe.drivers.coaddDriver
assert type(config)==lsst.pipe.drivers.coaddDriver.CoaddDriverConfig, 'config is of type %s.%s instead of lsst.pipe.drivers.coaddDriver.CoaddDriverConfig' % (type(config).__module__, type(config).__name__)
# Name for coadd
config.coaddName='deep'

import lsst.pipe.tasks.selectImages
config.select.retarget(target=lsst.pipe.tasks.selectImages.PsfWcsSelectImagesTask, ConfigClass=lsst.pipe.tasks.selectImages.PsfWcsSelectImagesConfig)

# Maximum median ellipticity residual
config.select.maxEllipResidual=0.007

# Maximum scatter in the size residuals
config.select.maxSizeScatter=None

# Maximum scatter in the size residuals, scaled by the median size
config.select.maxScaledSizeScatter=0.009

# select star with this field
config.select.starSelection='calib_psf_used'

# name of star shape
config.select.starShape='base_SdssShape'

# name of psf shape
config.select.psfShape='base_SdssShape_psf'

# Coadd name: typically one of deep or goodSeeing.
config.makeCoaddTempExp.coaddName='deep'

import lsst.pipe.drivers.utils
import lsst.pex.config.config
config.makeCoaddTempExp.select.retarget(target=lsst.pipe.drivers.utils.NullSelectImagesTask, ConfigClass=lsst.pex.config.config.Config)

# Mask planes that, if set, the associated pixel should not be included in the coaddTempExp.
config.makeCoaddTempExp.badMaskPlanes=['NO_DATA']

# Add records for CCDs we iterated over but did not add a coaddTempExp due to a lack of unmasked pixels in the coadd footprint.
config.makeCoaddTempExp.inputRecorder.saveEmptyCcds=False

# Add records for CCDs we iterated over but did not add a coaddTempExp due to an exception (often due to the calexp not being found on disk).
config.makeCoaddTempExp.inputRecorder.saveErrorCcds=False

# Save the total number of good pixels in each coaddTempExp (redundant with a sum of good pixels in associated CCDs)
config.makeCoaddTempExp.inputRecorder.saveVisitGoodPix=True

# Save weights in the CCDs table as well as the visits table? (This is necessary for easy construction of CoaddPsf, but otherwise duplicate information.)
config.makeCoaddTempExp.inputRecorder.saveCcdWeights=True

# Match to modelPsf? Deprecated. Sets makePsfMatched=True, makeDirect=False
config.makeCoaddTempExp.doPsfMatch=False

# Kernel size (width and height) (pixels); if None then sizeFactor is used
config.makeCoaddTempExp.modelPsf.size=None

# Kernel size as a factor of fwhm (dimensionless); size = sizeFactor * fwhm; ignored if size is not None
config.makeCoaddTempExp.modelPsf.sizeFactor=3.0

# Minimum kernel size if using sizeFactor (pixels); ignored if size is not None
config.makeCoaddTempExp.modelPsf.minSize=5

# Maximum kernel size if using sizeFactor (pixels); ignored if size is not None
config.makeCoaddTempExp.modelPsf.maxSize=None

# Default FWHM of Gaussian model of core of star (pixels)
config.makeCoaddTempExp.modelPsf.defaultFwhm=7.7

# Add a Gaussian to represent wings?
config.makeCoaddTempExp.modelPsf.addWing=True

# wing width, as a multiple of core width (dimensionless); ignored if addWing false
config.makeCoaddTempExp.modelPsf.wingFwhmFactor=2.5

# wing amplitude, as a multiple of core amplitude (dimensionless); ignored if addWing false
config.makeCoaddTempExp.modelPsf.wingAmplitude=0.1

# Apply meas_mosaic ubercal results to input calexps?
config.makeCoaddTempExp.doApplyUberCal=False

# Size in pixels of matching kernel. Must be odd.
config.makeCoaddTempExp.matchingKernelSize=29

# Warping kernel
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].warpingConfig.warpingKernelName='lanczos3'

# Warping kernel for mask (use warpingKernelName if '')
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].warpingConfig.maskWarpingKernelName='bilinear'

# interpLength argument to lsst.afw.math.warpExposure
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].warpingConfig.interpLength=10

# cacheSize argument to lsst.afw.math.SeparableKernel.computeCache
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].warpingConfig.cacheSize=1000000

# mask bits to grow to full width of image/variance kernel,
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].warpingConfig.growFullMask=16

# Value of footprint detection threshold
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].detectionConfig.detThreshold=10.0

# Type of detection threshold
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].detectionConfig.detThresholdType='pixel_stdev'

# If true run detection on the template (image to convolve);
#                  if false run detection on the science image
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].detectionConfig.detOnTemplate=True

# Mask planes that lead to an invalid detection.
#                  Options: NO_DATA EDGE SAT BAD CR INTRP
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].detectionConfig.badMaskPlanes=['NO_DATA', 'EDGE', 'SAT']

# Minimum number of pixels in an acceptable Footprint
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].detectionConfig.fpNpixMin=5

# Maximum number of pixels in an acceptable Footprint;
#                  too big and the subsequent convolutions become unwieldy
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].detectionConfig.fpNpixMax=500

# If config.scaleByFwhm, grow the footprint based on
#                  the final kernelSize.  Each footprint will be
#                  2*fpGrowKernelScaling*kernelSize x
#                  2*fpGrowKernelScaling*kernelSize.  With the value
#                  of 1.0, the remaining pixels in each KernelCandiate
#                  after convolution by the basis functions will be
#                  equal to the kernel size itself.
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].detectionConfig.fpGrowKernelScaling=1.0

# Growing radius (in pixels) for each raw detection
#                  footprint.  The smaller the faster; however the
#                  kernel sum does not converge if the stamp is too
#                  small; and the kernel is not constrained at all if
#                  the stamp is the size of the kernel.  The grown stamp
#                  is 2 * fpGrowPix pixels larger in each dimension.
#                  This is overridden by fpGrowKernelScaling if scaleByFwhm
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].detectionConfig.fpGrowPix=30

# Scale fpGrowPix by input Fwhm?
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].detectionConfig.scaleByFwhm=True

# type of statistic to use for grid points
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].afwBackgroundConfig.statisticsProperty='MEANCLIP'

# behaviour if there are too few points in grid for requested interpolation style
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].afwBackgroundConfig.undersampleStyle='REDUCE_INTERP_ORDER'

# how large a region of the sky should be used for each background point
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].afwBackgroundConfig.binSize=128

# Sky region size to be used for each background point in X direction. If 0, the binSize config is used.
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].afwBackgroundConfig.binSizeX=0

# Sky region size to be used for each background point in Y direction. If 0, the binSize config is used.
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].afwBackgroundConfig.binSizeY=0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].afwBackgroundConfig.algorithm='AKIMA_SPLINE'

# Names of mask planes to ignore while estimating the background
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].afwBackgroundConfig.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# Ignore NaNs when estimating the background
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].afwBackgroundConfig.isNanSafe=False

# Use Approximate (Chebyshev) to model background.
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].afwBackgroundConfig.useApprox=True

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].afwBackgroundConfig.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].afwBackgroundConfig.approxOrderY=-1

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].afwBackgroundConfig.weighting=True

# Use afw background subtraction instead of ip_diffim
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].useAfwBackground=False

# Include terms (including kernel cross terms) for background in ip_diffim
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].fitForBackground=False

# Type of basis set for PSF matching kernel.
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].kernelBasisSet='alard-lupton'

# Number of rows/columns in the convolution kernel; should be odd-valued.
#                  Modified by kernelSizeFwhmScaling if scaleByFwhm = true
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].kernelSize=21

# Scale kernelSize, alardGaussians by input Fwhm
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].scaleByFwhm=False

# How much to scale the kernel size based on the largest AL Sigma
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].kernelSizeFwhmScaling=6.0

# Minimum Kernel Size
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].kernelSizeMin=21

# Maximum Kernel Size
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].kernelSizeMax=35

# Type of spatial functions for kernel and background
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].spatialModelType='chebyshev1'

# Spatial order of convolution kernel variation
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].spatialKernelOrder=2

# Spatial order of differential background variation
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].spatialBgOrder=1

# Size (rows) in pixels of each SpatialCell for spatial modeling
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].sizeCellX=128

# Size (columns) in pixels of each SpatialCell for spatial modeling
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].sizeCellY=128

# Number of KernelCandidates in each SpatialCell to use in the spatial fitting
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].nStarPerCell=3

# Maximum number of iterations for rejecting bad KernelCandidates in spatial fitting
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].maxSpatialIterations=3

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

# Number of principal components to use for Pca basis, including the
#                  mean kernel if requested.
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].numPrincipalComponents=5

# Do sigma clipping on each raw kernel candidate
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].singleKernelClipping=False

# Do sigma clipping on the ensemble of kernel sums
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].kernelSumClipping=False

# Do sigma clipping after building the spatial model
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].spatialKernelClipping=False

# Test for maximum condition number when inverting a kernel matrix.
#                  Anything above maxConditionNumber is not used and the candidate is set as BAD.
#                  Also used to truncate inverse matrix in estimateBiasedRisk.  However,
#                  if you are doing any deconvolution you will want to turn this off, or use
#                  a large maxConditionNumber
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].checkConditionNumber=False

# Mask planes to ignore when calculating diffim statistics
#                  Options: NO_DATA EDGE SAT BAD CR INTRP
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].badMaskPlanes=['NO_DATA', 'EDGE', 'SAT']

# Rejects KernelCandidates yielding bad difference image quality.
#                  Used by BuildSingleKernelVisitor, AssessSpatialKernelVisitor.
#                  Represents average over pixels of (image/sqrt(variance)).
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].candidateResidualMeanMax=0.25

# Rejects KernelCandidates yielding bad difference image quality.
#                  Used by BuildSingleKernelVisitor, AssessSpatialKernelVisitor.
#                  Represents stddev over pixels of (image/sqrt(variance)).
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].candidateResidualStdMax=1.5

# Use the core of the footprint for the quality statistics, instead of the entire footprint.
#                  WARNING: if there is deconvolution we probably will need to turn this off
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].useCoreStats=False

# Radius for calculation of stats in 'core' of KernelCandidate diffim.
#                  Total number of pixels used will be (2*radius)**2.
#                  This is used both for 'core' diffim quality as well as ranking of
#                  KernelCandidates by their total flux in this core
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].candidateCoreRadius=3

# Maximum allowed sigma for outliers from kernel sum distribution.
#                  Used to reject variable objects from the kernel model
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].maxKsumSigma=3.0

# Maximum condition number for a well conditioned matrix
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].maxConditionNumber=50000000.0

# Use singular values (SVD) or eigen values (EIGENVALUE) to determine condition number
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].conditionNumberType='EIGENVALUE'

# Maximum condition number for a well conditioned spatial matrix
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].maxSpatialConditionNumber=10000000000.0

# Remake KernelCandidate using better variance estimate after first pass?
#                  Primarily useful when convolving a single-depth image, otherwise not necessary.
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].iterateSingleKernel=False

# Use constant variance weighting in single kernel fitting?
#                  In some cases this is better for bright star residuals.
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].constantVarianceWeighting=True

# Calculate kernel and background uncertainties for each kernel candidate?
#                  This comes from the inverse of the covariance matrix.
#                  Warning: regularization can cause problems for this step.
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].calculateKernelUncertainty=False

# Use Bayesian Information Criterion to select the number of bases going into the kernel
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].useBicForKernelBasis=False

# Number of Gaussians in alard-lupton basis
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].alardNGauss=3

# Polynomial order of spatial modification of Gaussians.  Must in number equal alardNGauss
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].alardDegGauss=[4, 2, 2]

# Sigma in pixels of Gaussians (FWHM = 2.35 sigma).  Must in number equal alardNGauss
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].alardSigGauss=[1.0, 2.0, 4.5]

# Default scale factor between Gaussian sigmas 
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].alardGaussBeta=2.0

# Minimum Sigma (pixels) for Gaussians
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].alardMinSig=0.7

# Degree of spatial modification of ALL gaussians in AL basis during deconvolution
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].alardDegGaussDeconv=3

# Minimum Sigma (pixels) for Gaussians during deconvolution;
#         make smaller than alardMinSig as this is only indirectly used
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].alardMinSigDeconv=0.4

# Number of Gaussians in AL basis during deconvolution
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel['AL'].alardNGaussDeconv=3

config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.kernel.name='AL'
# If too small, automatically pad the science Psf? Pad to smallest dimensions appropriate for the matching kernel dimensions, as specified by autoPadPsfTo. If false, pad by the padPsfBy config.
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.doAutoPadPsf=True

# Minimum Science Psf dimensions as a fraction of matching kernel dimensions. If the dimensions of the Psf to be matched are less than the matching kernel dimensions * autoPadPsfTo, pad Science Psf to this size. Ignored if doAutoPadPsf=False.
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.autoPadPsfTo=1.4

# Pixels (even) to pad Science Psf by before matching. Ignored if doAutoPadPsf=True
config.makeCoaddTempExp.warpAndPsfMatch.psfMatch.padPsfBy=0

# Warping kernel
config.makeCoaddTempExp.warpAndPsfMatch.warp.warpingKernelName='lanczos3'

# Warping kernel for mask (use warpingKernelName if '')
config.makeCoaddTempExp.warpAndPsfMatch.warp.maskWarpingKernelName='bilinear'

# interpLength argument to lsst.afw.math.warpExposure
config.makeCoaddTempExp.warpAndPsfMatch.warp.interpLength=10

# cacheSize argument to lsst.afw.math.SeparableKernel.computeCache
config.makeCoaddTempExp.warpAndPsfMatch.warp.cacheSize=1000000

# mask bits to grow to full width of image/variance kernel,
config.makeCoaddTempExp.warpAndPsfMatch.warp.growFullMask=16

# persist <coaddName>Coadd_<warpType>Warp
config.makeCoaddTempExp.doWrite=True

# Work with a background subtracted calexp?
config.makeCoaddTempExp.bgSubtracted=True

# Warping kernel cache size
config.makeCoaddTempExp.coaddPsf.cacheSize=10000

# Name of warping kernel; choices: lanczos3,lanczos4,lanczos5,bilinear,nearest
config.makeCoaddTempExp.coaddPsf.warpingKernelName='lanczos3'

# Make direct Warp/Coadds
config.makeCoaddTempExp.makeDirect=True

# Make Psf-Matched Warp/Coadd?
config.makeCoaddTempExp.makePsfMatched=True

# Apply sky correction?
config.makeCoaddTempExp.doApplySkyCorr=False

# Build background reference?
config.doBackgroundReference=False

import lsst.pipe.tasks.assembleCoadd
config.assembleCoadd.retarget(target=lsst.pipe.tasks.assembleCoadd.CompareWarpAssembleCoaddTask, ConfigClass=lsst.pipe.tasks.assembleCoadd.CompareWarpAssembleCoaddConfig)

# Coadd name: typically one of deep or goodSeeing.
config.assembleCoadd.coaddName='deep'

import lsst.pipe.drivers.utils
import lsst.pex.config.config
config.assembleCoadd.select.retarget(target=lsst.pipe.drivers.utils.NullSelectImagesTask, ConfigClass=lsst.pex.config.config.Config)

# Mask planes that, if set, the associated pixel should not be included in the coaddTempExp.
config.assembleCoadd.badMaskPlanes=['NO_DATA', 'BAD', 'SAT', 'SUSPECT']

# Add records for CCDs we iterated over but did not add a coaddTempExp due to a lack of unmasked pixels in the coadd footprint.
config.assembleCoadd.inputRecorder.saveEmptyCcds=False

# Add records for CCDs we iterated over but did not add a coaddTempExp due to an exception (often due to the calexp not being found on disk).
config.assembleCoadd.inputRecorder.saveErrorCcds=False

# Save the total number of good pixels in each coaddTempExp (redundant with a sum of good pixels in associated CCDs)
config.assembleCoadd.inputRecorder.saveVisitGoodPix=True

# Save weights in the CCDs table as well as the visits table? (This is necessary for easy construction of CoaddPsf, but otherwise duplicate information.)
config.assembleCoadd.inputRecorder.saveCcdWeights=True

# Match to modelPsf? Deprecated. Sets makePsfMatched=True, makeDirect=False
config.assembleCoadd.doPsfMatch=False

# Kernel size (width and height) (pixels); if None then sizeFactor is used
config.assembleCoadd.modelPsf.size=None

# Kernel size as a factor of fwhm (dimensionless); size = sizeFactor * fwhm; ignored if size is not None
config.assembleCoadd.modelPsf.sizeFactor=3.0

# Minimum kernel size if using sizeFactor (pixels); ignored if size is not None
config.assembleCoadd.modelPsf.minSize=5

# Maximum kernel size if using sizeFactor (pixels); ignored if size is not None
config.assembleCoadd.modelPsf.maxSize=None

# Default FWHM of Gaussian model of core of star (pixels)
config.assembleCoadd.modelPsf.defaultFwhm=3.0

# Add a Gaussian to represent wings?
config.assembleCoadd.modelPsf.addWing=True

# wing width, as a multiple of core width (dimensionless); ignored if addWing false
config.assembleCoadd.modelPsf.wingFwhmFactor=2.5

# wing amplitude, as a multiple of core amplitude (dimensionless); ignored if addWing false
config.assembleCoadd.modelPsf.wingAmplitude=0.1

# Apply meas_mosaic ubercal results to input calexps?
config.assembleCoadd.doApplyUberCal=False

# Size in pixels of matching kernel. Must be odd.
config.assembleCoadd.matchingKernelSize=29

# Warp name: one of 'direct' or 'psfMatched'
config.assembleCoadd.warpType='direct'

# Width, height of stack subregion size; make small enough that a full stack of images will fit into memory at once.
config.assembleCoadd.subregionSize=[10000, 200]

# Main stacking statistic for aggregating over the epochs.
config.assembleCoadd.statistic='MEAN'

# Perform sigma clipped outlier rejection with MEANCLIP statistic? (DEPRECATED)
config.assembleCoadd.doSigmaClip=False

# Sigma for outlier rejection; ignored if non-clipping statistic selected.
config.assembleCoadd.sigmaClip=3.0

# Number of iterations of outlier rejection; ignored if non-clipping statistic selected.
config.assembleCoadd.clipIter=2

# Calculate coadd variance from input variance by stacking statistic.Passed to StatisticsControl.setCalcErrorFromInputVariance()
config.assembleCoadd.calcErrorFromInputVariance=True

# desired photometric zero point
config.assembleCoadd.scaleZeroPoint.zeroPoint=27.0

# Interpolate over NaN pixels? Also extrapolate, if necessary, but the results are ugly.
config.assembleCoadd.doInterp=True

# Kernel size (width and height) (pixels); if None then sizeFactor is used
config.assembleCoadd.interpImage.modelPsf.size=None

# Kernel size as a factor of fwhm (dimensionless); size = sizeFactor * fwhm; ignored if size is not None
config.assembleCoadd.interpImage.modelPsf.sizeFactor=3.0

# Minimum kernel size if using sizeFactor (pixels); ignored if size is not None
config.assembleCoadd.interpImage.modelPsf.minSize=5

# Maximum kernel size if using sizeFactor (pixels); ignored if size is not None
config.assembleCoadd.interpImage.modelPsf.maxSize=None

# Default FWHM of Gaussian model of core of star (pixels)
config.assembleCoadd.interpImage.modelPsf.defaultFwhm=3.0

# Add a Gaussian to represent wings?
config.assembleCoadd.interpImage.modelPsf.addWing=True

# wing width, as a multiple of core width (dimensionless); ignored if addWing false
config.assembleCoadd.interpImage.modelPsf.wingFwhmFactor=2.5

# wing amplitude, as a multiple of core amplitude (dimensionless); ignored if addWing false
config.assembleCoadd.interpImage.modelPsf.wingAmplitude=0.1

# Smoothly taper to the fallback value at the edge of the image?
config.assembleCoadd.interpImage.useFallbackValueAtEdge=True

# Type of statistic to calculate edge fallbackValue for interpolation
config.assembleCoadd.interpImage.fallbackValueType='MEDIAN'

# If fallbackValueType is 'USER' then use this as the fallbackValue; ignored otherwise
config.assembleCoadd.interpImage.fallbackUserValue=0.0

# Allow negative values for egde interpolation fallbackValue?  If False, set fallbackValue to max(fallbackValue, 0.0)
config.assembleCoadd.interpImage.negativeFallbackAllowed=False

# Persist coadd?
config.assembleCoadd.doWrite=True

# Create image of number of contributing exposures for each pixel
config.assembleCoadd.doNImage=True

# Use ValidPolygons from shrunk Psf-Matched Calexps? Should be set to True by CompareWarp only.
config.assembleCoadd.doUsePsfMatchedPolygons=True

# Threshold (in fractional weight) of rejection at which we propagate a mask plane to the coadd; that is, we set the mask bit on the coadd if the fraction the rejected frames would have contributed exceeds this value.
config.assembleCoadd.maskPropagationThresholds={'SAT': 0.1}

# Mask planes to remove before coadding
config.assembleCoadd.removeMaskPlanes=['NOT_DEBLENDED', 'EDGE', 'CROSSTALK']

# Set mask and flag bits for bright objects?
config.assembleCoadd.doMaskBrightObjects=True

# Name of mask bit used for bright objects
config.assembleCoadd.brightObjectMaskName='BRIGHT_OBJECT'

# Warping kernel cache size
config.assembleCoadd.coaddPsf.cacheSize=10000

# Name of warping kernel; choices: lanczos3,lanczos4,lanczos5,bilinear,nearest
config.assembleCoadd.coaddPsf.warpingKernelName='lanczos3'

# Attach a piecewise TransmissionCurve for the coadd? (requires all input Exposures to have TransmissionCurves).
config.assembleCoadd.doAttachTransmissionCurve=False

# Coadd name: typically one of deep or goodSeeing.
config.assembleCoadd.assembleStaticSkyModel.coaddName='deep'

# Mask planes that, if set, the associated pixel should not be included in the coaddTempExp.
config.assembleCoadd.assembleStaticSkyModel.badMaskPlanes=['NO_DATA']

# Add records for CCDs we iterated over but did not add a coaddTempExp due to a lack of unmasked pixels in the coadd footprint.
config.assembleCoadd.assembleStaticSkyModel.inputRecorder.saveEmptyCcds=False

# Add records for CCDs we iterated over but did not add a coaddTempExp due to an exception (often due to the calexp not being found on disk).
config.assembleCoadd.assembleStaticSkyModel.inputRecorder.saveErrorCcds=False

# Save the total number of good pixels in each coaddTempExp (redundant with a sum of good pixels in associated CCDs)
config.assembleCoadd.assembleStaticSkyModel.inputRecorder.saveVisitGoodPix=True

# Save weights in the CCDs table as well as the visits table? (This is necessary for easy construction of CoaddPsf, but otherwise duplicate information.)
config.assembleCoadd.assembleStaticSkyModel.inputRecorder.saveCcdWeights=True

# Match to modelPsf? Deprecated. Sets makePsfMatched=True, makeDirect=False
config.assembleCoadd.assembleStaticSkyModel.doPsfMatch=False

# Kernel size (width and height) (pixels); if None then sizeFactor is used
config.assembleCoadd.assembleStaticSkyModel.modelPsf.size=None

# Kernel size as a factor of fwhm (dimensionless); size = sizeFactor * fwhm; ignored if size is not None
config.assembleCoadd.assembleStaticSkyModel.modelPsf.sizeFactor=3.0

# Minimum kernel size if using sizeFactor (pixels); ignored if size is not None
config.assembleCoadd.assembleStaticSkyModel.modelPsf.minSize=5

# Maximum kernel size if using sizeFactor (pixels); ignored if size is not None
config.assembleCoadd.assembleStaticSkyModel.modelPsf.maxSize=None

# Default FWHM of Gaussian model of core of star (pixels)
config.assembleCoadd.assembleStaticSkyModel.modelPsf.defaultFwhm=3.0

# Add a Gaussian to represent wings?
config.assembleCoadd.assembleStaticSkyModel.modelPsf.addWing=True

# wing width, as a multiple of core width (dimensionless); ignored if addWing false
config.assembleCoadd.assembleStaticSkyModel.modelPsf.wingFwhmFactor=2.5

# wing amplitude, as a multiple of core amplitude (dimensionless); ignored if addWing false
config.assembleCoadd.assembleStaticSkyModel.modelPsf.wingAmplitude=0.1

# Apply meas_mosaic ubercal results to input calexps?
config.assembleCoadd.assembleStaticSkyModel.doApplyUberCal=False

# Size in pixels of matching kernel. Must be odd.
config.assembleCoadd.assembleStaticSkyModel.matchingKernelSize=21

# Warp name: one of 'direct' or 'psfMatched'
config.assembleCoadd.assembleStaticSkyModel.warpType='psfMatched'

# Width, height of stack subregion size; make small enough that a full stack of images will fit into memory at once.
config.assembleCoadd.assembleStaticSkyModel.subregionSize=[10000, 200]

# Main stacking statistic for aggregating over the epochs.
config.assembleCoadd.assembleStaticSkyModel.statistic='MEANCLIP'

# Perform sigma clipped outlier rejection with MEANCLIP statistic? (DEPRECATED)
config.assembleCoadd.assembleStaticSkyModel.doSigmaClip=False

# Sigma for outlier rejection; ignored if non-clipping statistic selected.
config.assembleCoadd.assembleStaticSkyModel.sigmaClip=2.5

# Number of iterations of outlier rejection; ignored if non-clipping statistic selected.
config.assembleCoadd.assembleStaticSkyModel.clipIter=3

# Calculate coadd variance from input variance by stacking statistic.Passed to StatisticsControl.setCalcErrorFromInputVariance()
config.assembleCoadd.assembleStaticSkyModel.calcErrorFromInputVariance=False

# desired photometric zero point
config.assembleCoadd.assembleStaticSkyModel.scaleZeroPoint.zeroPoint=27.0

# Interpolate over NaN pixels? Also extrapolate, if necessary, but the results are ugly.
config.assembleCoadd.assembleStaticSkyModel.doInterp=True

# Kernel size (width and height) (pixels); if None then sizeFactor is used
config.assembleCoadd.assembleStaticSkyModel.interpImage.modelPsf.size=None

# Kernel size as a factor of fwhm (dimensionless); size = sizeFactor * fwhm; ignored if size is not None
config.assembleCoadd.assembleStaticSkyModel.interpImage.modelPsf.sizeFactor=3.0

# Minimum kernel size if using sizeFactor (pixels); ignored if size is not None
config.assembleCoadd.assembleStaticSkyModel.interpImage.modelPsf.minSize=5

# Maximum kernel size if using sizeFactor (pixels); ignored if size is not None
config.assembleCoadd.assembleStaticSkyModel.interpImage.modelPsf.maxSize=None

# Default FWHM of Gaussian model of core of star (pixels)
config.assembleCoadd.assembleStaticSkyModel.interpImage.modelPsf.defaultFwhm=3.0

# Add a Gaussian to represent wings?
config.assembleCoadd.assembleStaticSkyModel.interpImage.modelPsf.addWing=True

# wing width, as a multiple of core width (dimensionless); ignored if addWing false
config.assembleCoadd.assembleStaticSkyModel.interpImage.modelPsf.wingFwhmFactor=2.5

# wing amplitude, as a multiple of core amplitude (dimensionless); ignored if addWing false
config.assembleCoadd.assembleStaticSkyModel.interpImage.modelPsf.wingAmplitude=0.1

# Smoothly taper to the fallback value at the edge of the image?
config.assembleCoadd.assembleStaticSkyModel.interpImage.useFallbackValueAtEdge=True

# Type of statistic to calculate edge fallbackValue for interpolation
config.assembleCoadd.assembleStaticSkyModel.interpImage.fallbackValueType='MEDIAN'

# If fallbackValueType is 'USER' then use this as the fallbackValue; ignored otherwise
config.assembleCoadd.assembleStaticSkyModel.interpImage.fallbackUserValue=0.0

# Allow negative values for egde interpolation fallbackValue?  If False, set fallbackValue to max(fallbackValue, 0.0)
config.assembleCoadd.assembleStaticSkyModel.interpImage.negativeFallbackAllowed=False

# Persist coadd?
config.assembleCoadd.assembleStaticSkyModel.doWrite=False

# Create image of number of contributing exposures for each pixel
config.assembleCoadd.assembleStaticSkyModel.doNImage=False

# Use ValidPolygons from shrunk Psf-Matched Calexps? Should be set to True by CompareWarp only.
config.assembleCoadd.assembleStaticSkyModel.doUsePsfMatchedPolygons=False

# Threshold (in fractional weight) of rejection at which we propagate a mask plane to the coadd; that is, we set the mask bit on the coadd if the fraction the rejected frames would have contributed exceeds this value.
config.assembleCoadd.assembleStaticSkyModel.maskPropagationThresholds={'SAT': 0.1}

# Mask planes to remove before coadding
config.assembleCoadd.assembleStaticSkyModel.removeMaskPlanes=['NOT_DEBLENDED']

# Set mask and flag bits for bright objects?
config.assembleCoadd.assembleStaticSkyModel.doMaskBrightObjects=False

# Name of mask bit used for bright objects
config.assembleCoadd.assembleStaticSkyModel.brightObjectMaskName='BRIGHT_OBJECT'

# Warping kernel cache size
config.assembleCoadd.assembleStaticSkyModel.coaddPsf.cacheSize=10000

# Name of warping kernel; choices: lanczos3,lanczos4,lanczos5,bilinear,nearest
config.assembleCoadd.assembleStaticSkyModel.coaddPsf.warpingKernelName='lanczos3'

# Attach a piecewise TransmissionCurve for the coadd? (requires all input Exposures to have TransmissionCurves).
config.assembleCoadd.assembleStaticSkyModel.doAttachTransmissionCurve=False

# detected sources with fewer than the specified number of pixels will be ignored
config.assembleCoadd.detect.minPixels=4

# Pixels should be grown as isotropically as possible (slower)
config.assembleCoadd.detect.isotropicGrow=True

# Grow all footprints at the same time? This allows disconnected footprints to merge.
config.assembleCoadd.detect.combinedGrow=True

# Grow detections by nSigmaToGrow * [PSF RMS width]; if 0 then do not grow
config.assembleCoadd.detect.nSigmaToGrow=2.0

# Grow detections to set the image mask bits, but return the original (not-grown) footprints
config.assembleCoadd.detect.returnOriginalFootprints=False

# Threshold for footprints
config.assembleCoadd.detect.thresholdValue=5.0

# Include threshold relative to thresholdValue
config.assembleCoadd.detect.includeThresholdMultiplier=1.0

# specifies the desired flavor of Threshold
config.assembleCoadd.detect.thresholdType='pixel_stdev'

# specifies whether to detect positive, or negative sources, or both
config.assembleCoadd.detect.thresholdPolarity='both'

# Fiddle factor to add to the background; debugging only
config.assembleCoadd.detect.adjustBackground=0.0

# Estimate the background again after final source detection?
config.assembleCoadd.detect.reEstimateBackground=False

# type of statistic to use for grid points
config.assembleCoadd.detect.background.statisticsProperty='MEANCLIP'

# behaviour if there are too few points in grid for requested interpolation style
config.assembleCoadd.detect.background.undersampleStyle='REDUCE_INTERP_ORDER'

# how large a region of the sky should be used for each background point
config.assembleCoadd.detect.background.binSize=128

# Sky region size to be used for each background point in X direction. If 0, the binSize config is used.
config.assembleCoadd.detect.background.binSizeX=0

# Sky region size to be used for each background point in Y direction. If 0, the binSize config is used.
config.assembleCoadd.detect.background.binSizeY=0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
config.assembleCoadd.detect.background.algorithm='AKIMA_SPLINE'

# Names of mask planes to ignore while estimating the background
config.assembleCoadd.detect.background.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# Ignore NaNs when estimating the background
config.assembleCoadd.detect.background.isNanSafe=False

# Use Approximate (Chebyshev) to model background.
config.assembleCoadd.detect.background.useApprox=True

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.assembleCoadd.detect.background.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.assembleCoadd.detect.background.approxOrderY=-1

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.assembleCoadd.detect.background.weighting=True

# type of statistic to use for grid points
config.assembleCoadd.detect.tempLocalBackground.statisticsProperty='MEANCLIP'

# behaviour if there are too few points in grid for requested interpolation style
config.assembleCoadd.detect.tempLocalBackground.undersampleStyle='REDUCE_INTERP_ORDER'

# how large a region of the sky should be used for each background point
config.assembleCoadd.detect.tempLocalBackground.binSize=64

# Sky region size to be used for each background point in X direction. If 0, the binSize config is used.
config.assembleCoadd.detect.tempLocalBackground.binSizeX=0

# Sky region size to be used for each background point in Y direction. If 0, the binSize config is used.
config.assembleCoadd.detect.tempLocalBackground.binSizeY=0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
config.assembleCoadd.detect.tempLocalBackground.algorithm='AKIMA_SPLINE'

# Names of mask planes to ignore while estimating the background
config.assembleCoadd.detect.tempLocalBackground.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# Ignore NaNs when estimating the background
config.assembleCoadd.detect.tempLocalBackground.isNanSafe=False

# Use Approximate (Chebyshev) to model background.
config.assembleCoadd.detect.tempLocalBackground.useApprox=False

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.assembleCoadd.detect.tempLocalBackground.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.assembleCoadd.detect.tempLocalBackground.approxOrderY=-1

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.assembleCoadd.detect.tempLocalBackground.weighting=True

# Enable temporary local background subtraction? (see tempLocalBackground)
config.assembleCoadd.detect.doTempLocalBackground=False

# type of statistic to use for grid points
config.assembleCoadd.detect.tempWideBackground.statisticsProperty='MEANCLIP'

# behaviour if there are too few points in grid for requested interpolation style
config.assembleCoadd.detect.tempWideBackground.undersampleStyle='REDUCE_INTERP_ORDER'

# how large a region of the sky should be used for each background point
config.assembleCoadd.detect.tempWideBackground.binSize=512

# Sky region size to be used for each background point in X direction. If 0, the binSize config is used.
config.assembleCoadd.detect.tempWideBackground.binSizeX=0

# Sky region size to be used for each background point in Y direction. If 0, the binSize config is used.
config.assembleCoadd.detect.tempWideBackground.binSizeY=0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
config.assembleCoadd.detect.tempWideBackground.algorithm='AKIMA_SPLINE'

# Names of mask planes to ignore while estimating the background
config.assembleCoadd.detect.tempWideBackground.ignoredPixelMask=['BAD', 'EDGE', 'NO_DATA']

# Ignore NaNs when estimating the background
config.assembleCoadd.detect.tempWideBackground.isNanSafe=False

# Use Approximate (Chebyshev) to model background.
config.assembleCoadd.detect.tempWideBackground.useApprox=False

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.assembleCoadd.detect.tempWideBackground.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.assembleCoadd.detect.tempWideBackground.approxOrderY=-1

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.assembleCoadd.detect.tempWideBackground.weighting=True

# Do temporary wide (large-scale) background subtraction before footprint detection?
config.assembleCoadd.detect.doTempWideBackground=False

# The maximum number of peaks in a Footprint before trying to replace its peaks using the temporary local background
config.assembleCoadd.detect.nPeaksMaxSimple=1

# Multiple of PSF RMS size to use for convolution kernel bounding box size; note that this is not a half-size. The size will be rounded up to the nearest odd integer
config.assembleCoadd.detect.nSigmaForKernel=7.0

# Mask planes to ignore when calculating statistics of image (for thresholdType=stdev)
config.assembleCoadd.detect.statsMask=['BAD', 'SAT', 'EDGE', 'NO_DATA']

# detected sources with fewer than the specified number of pixels will be ignored
config.assembleCoadd.detectTemplate.minPixels=1

# Pixels should be grown as isotropically as possible (slower)
config.assembleCoadd.detectTemplate.isotropicGrow=False

# Grow all footprints at the same time? This allows disconnected footprints to merge.
config.assembleCoadd.detectTemplate.combinedGrow=True

# Grow detections by nSigmaToGrow * [PSF RMS width]; if 0 then do not grow
config.assembleCoadd.detectTemplate.nSigmaToGrow=2.0

# Grow detections to set the image mask bits, but return the original (not-grown) footprints
config.assembleCoadd.detectTemplate.returnOriginalFootprints=False

# Threshold for footprints
config.assembleCoadd.detectTemplate.thresholdValue=5.0

# Include threshold relative to thresholdValue
config.assembleCoadd.detectTemplate.includeThresholdMultiplier=1.0

# specifies the desired flavor of Threshold
config.assembleCoadd.detectTemplate.thresholdType='stdev'

# specifies whether to detect positive, or negative sources, or both
config.assembleCoadd.detectTemplate.thresholdPolarity='positive'

# Fiddle factor to add to the background; debugging only
config.assembleCoadd.detectTemplate.adjustBackground=0.0

# Estimate the background again after final source detection?
config.assembleCoadd.detectTemplate.reEstimateBackground=False

# type of statistic to use for grid points
config.assembleCoadd.detectTemplate.background.statisticsProperty='MEANCLIP'

# behaviour if there are too few points in grid for requested interpolation style
config.assembleCoadd.detectTemplate.background.undersampleStyle='REDUCE_INTERP_ORDER'

# how large a region of the sky should be used for each background point
config.assembleCoadd.detectTemplate.background.binSize=128

# Sky region size to be used for each background point in X direction. If 0, the binSize config is used.
config.assembleCoadd.detectTemplate.background.binSizeX=0

# Sky region size to be used for each background point in Y direction. If 0, the binSize config is used.
config.assembleCoadd.detectTemplate.background.binSizeY=0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
config.assembleCoadd.detectTemplate.background.algorithm='AKIMA_SPLINE'

# Names of mask planes to ignore while estimating the background
config.assembleCoadd.detectTemplate.background.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# Ignore NaNs when estimating the background
config.assembleCoadd.detectTemplate.background.isNanSafe=False

# Use Approximate (Chebyshev) to model background.
config.assembleCoadd.detectTemplate.background.useApprox=True

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.assembleCoadd.detectTemplate.background.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.assembleCoadd.detectTemplate.background.approxOrderY=-1

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.assembleCoadd.detectTemplate.background.weighting=True

# type of statistic to use for grid points
config.assembleCoadd.detectTemplate.tempLocalBackground.statisticsProperty='MEANCLIP'

# behaviour if there are too few points in grid for requested interpolation style
config.assembleCoadd.detectTemplate.tempLocalBackground.undersampleStyle='REDUCE_INTERP_ORDER'

# how large a region of the sky should be used for each background point
config.assembleCoadd.detectTemplate.tempLocalBackground.binSize=64

# Sky region size to be used for each background point in X direction. If 0, the binSize config is used.
config.assembleCoadd.detectTemplate.tempLocalBackground.binSizeX=0

# Sky region size to be used for each background point in Y direction. If 0, the binSize config is used.
config.assembleCoadd.detectTemplate.tempLocalBackground.binSizeY=0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
config.assembleCoadd.detectTemplate.tempLocalBackground.algorithm='AKIMA_SPLINE'

# Names of mask planes to ignore while estimating the background
config.assembleCoadd.detectTemplate.tempLocalBackground.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# Ignore NaNs when estimating the background
config.assembleCoadd.detectTemplate.tempLocalBackground.isNanSafe=False

# Use Approximate (Chebyshev) to model background.
config.assembleCoadd.detectTemplate.tempLocalBackground.useApprox=False

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.assembleCoadd.detectTemplate.tempLocalBackground.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.assembleCoadd.detectTemplate.tempLocalBackground.approxOrderY=-1

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.assembleCoadd.detectTemplate.tempLocalBackground.weighting=True

# Enable temporary local background subtraction? (see tempLocalBackground)
config.assembleCoadd.detectTemplate.doTempLocalBackground=False

# type of statistic to use for grid points
config.assembleCoadd.detectTemplate.tempWideBackground.statisticsProperty='MEANCLIP'

# behaviour if there are too few points in grid for requested interpolation style
config.assembleCoadd.detectTemplate.tempWideBackground.undersampleStyle='REDUCE_INTERP_ORDER'

# how large a region of the sky should be used for each background point
config.assembleCoadd.detectTemplate.tempWideBackground.binSize=512

# Sky region size to be used for each background point in X direction. If 0, the binSize config is used.
config.assembleCoadd.detectTemplate.tempWideBackground.binSizeX=0

# Sky region size to be used for each background point in Y direction. If 0, the binSize config is used.
config.assembleCoadd.detectTemplate.tempWideBackground.binSizeY=0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
config.assembleCoadd.detectTemplate.tempWideBackground.algorithm='AKIMA_SPLINE'

# Names of mask planes to ignore while estimating the background
config.assembleCoadd.detectTemplate.tempWideBackground.ignoredPixelMask=['BAD', 'EDGE', 'NO_DATA']

# Ignore NaNs when estimating the background
config.assembleCoadd.detectTemplate.tempWideBackground.isNanSafe=False

# Use Approximate (Chebyshev) to model background.
config.assembleCoadd.detectTemplate.tempWideBackground.useApprox=False

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.assembleCoadd.detectTemplate.tempWideBackground.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.assembleCoadd.detectTemplate.tempWideBackground.approxOrderY=-1

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.assembleCoadd.detectTemplate.tempWideBackground.weighting=True

# Do temporary wide (large-scale) background subtraction before footprint detection?
config.assembleCoadd.detectTemplate.doTempWideBackground=False

# The maximum number of peaks in a Footprint before trying to replace its peaks using the temporary local background
config.assembleCoadd.detectTemplate.nPeaksMaxSimple=1

# Multiple of PSF RMS size to use for convolution kernel bounding box size; note that this is not a half-size. The size will be rounded up to the nearest odd integer
config.assembleCoadd.detectTemplate.nSigmaForKernel=7.0

# Mask planes to ignore when calculating statistics of image (for thresholdType=stdev)
config.assembleCoadd.detectTemplate.statsMask=['BAD', 'SAT', 'EDGE', 'NO_DATA']

# Charactistic maximum local number of epochs/visits in which an artifact candidate can appear  and still be masked.  The effective maxNumEpochs is a broken linear function of local number of epochs (N): min(maxFractionEpochsLow*N, maxNumEpochs + maxFractionEpochsHigh*N). For each footprint detected on the image difference between the psfMatched warp and static sky model, if a significant fraction of pixels (defined by spatialThreshold) are residuals in more than the computed effective maxNumEpochs, the artifact candidate is deemed persistant rather than transient and not masked.
config.assembleCoadd.maxNumEpochs=2

# Fraction of local number of epochs (N) to use as effective maxNumEpochs for low N. Effective maxNumEpochs = min(maxFractionEpochsLow * N, maxNumEpochs + maxFractionEpochsHigh * N)
config.assembleCoadd.maxFractionEpochsLow=0.4

# Fraction of local number of epochs (N) to use as effective maxNumEpochs for high N. Effective maxNumEpochs = min(maxFractionEpochsLow * N, maxNumEpochs + maxFractionEpochsHigh * N)
config.assembleCoadd.maxFractionEpochsHigh=0.03

# Unitless fraction of pixels defining how much of the outlier region has to meet the temporal criteria. If 0, clip all. If 1, clip none.
config.assembleCoadd.spatialThreshold=0.5

# Rescale Warp variance plane using empirical noise?
config.assembleCoadd.doScaleWarpVariance=True

# type of statistic to use for grid points
config.assembleCoadd.scaleWarpVariance.background.statisticsProperty='MEANCLIP'

# behaviour if there are too few points in grid for requested interpolation style
config.assembleCoadd.scaleWarpVariance.background.undersampleStyle='REDUCE_INTERP_ORDER'

# how large a region of the sky should be used for each background point
config.assembleCoadd.scaleWarpVariance.background.binSize=32

# Sky region size to be used for each background point in X direction. If 0, the binSize config is used.
config.assembleCoadd.scaleWarpVariance.background.binSizeX=0

# Sky region size to be used for each background point in Y direction. If 0, the binSize config is used.
config.assembleCoadd.scaleWarpVariance.background.binSizeY=0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
config.assembleCoadd.scaleWarpVariance.background.algorithm='AKIMA_SPLINE'

# Names of mask planes to ignore while estimating the background
config.assembleCoadd.scaleWarpVariance.background.ignoredPixelMask=['DETECTED', 'DETECTED_NEGATIVE', 'BAD', 'SAT', 'NO_DATA', 'INTRP']

# Ignore NaNs when estimating the background
config.assembleCoadd.scaleWarpVariance.background.isNanSafe=False

# Use Approximate (Chebyshev) to model background.
config.assembleCoadd.scaleWarpVariance.background.useApprox=False

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.assembleCoadd.scaleWarpVariance.background.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.assembleCoadd.scaleWarpVariance.background.approxOrderY=-1

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.assembleCoadd.scaleWarpVariance.background.weighting=True

# Mask planes for pixels to ignore when scaling variance
config.assembleCoadd.scaleWarpVariance.maskPlanes=['DETECTED', 'DETECTED_NEGATIVE', 'BAD', 'SAT', 'NO_DATA', 'INTRP']

# Maximum variance scaling value to permit
config.assembleCoadd.scaleWarpVariance.limit=10.0

# Rescue artifacts from clipping that completely lie within a footprint detectedon the PsfMatched Template Coadd. Replicates a behavior of SafeClip.
config.assembleCoadd.doPreserveContainedBySource=True

# Ignore artifact candidates that are mostly covered by the bad pixel mask, because they will be excluded anyway. This prevents them from contributing to the outlier epoch count image and potentially being labeled as persistant.'Mostly' is defined by the config 'prefilterArtifactsRatio'.
config.assembleCoadd.doPrefilterArtifacts=True

# Prefilter artifact candidates that are mostly covered by these bad mask planes.
config.assembleCoadd.prefilterArtifactsMaskPlanes=['NO_DATA', 'BAD', 'SAT', 'SUSPECT']

# Prefilter artifact candidates with less than this fraction overlapping good pixels
config.assembleCoadd.prefilterArtifactsRatio=0.05

# Run detection on the coaddition product
config.doDetection=True

# Scale variance plane using empirical noise?
config.detectCoaddSources.doScaleVariance=True

# type of statistic to use for grid points
config.detectCoaddSources.scaleVariance.background.statisticsProperty='MEANCLIP'

# behaviour if there are too few points in grid for requested interpolation style
config.detectCoaddSources.scaleVariance.background.undersampleStyle='REDUCE_INTERP_ORDER'

# how large a region of the sky should be used for each background point
config.detectCoaddSources.scaleVariance.background.binSize=32

# Sky region size to be used for each background point in X direction. If 0, the binSize config is used.
config.detectCoaddSources.scaleVariance.background.binSizeX=0

# Sky region size to be used for each background point in Y direction. If 0, the binSize config is used.
config.detectCoaddSources.scaleVariance.background.binSizeY=0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
config.detectCoaddSources.scaleVariance.background.algorithm='AKIMA_SPLINE'

# Names of mask planes to ignore while estimating the background
config.detectCoaddSources.scaleVariance.background.ignoredPixelMask=['DETECTED', 'DETECTED_NEGATIVE', 'BAD', 'SAT', 'NO_DATA', 'INTRP']

# Ignore NaNs when estimating the background
config.detectCoaddSources.scaleVariance.background.isNanSafe=False

# Use Approximate (Chebyshev) to model background.
config.detectCoaddSources.scaleVariance.background.useApprox=False

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.detectCoaddSources.scaleVariance.background.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.detectCoaddSources.scaleVariance.background.approxOrderY=-1

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.detectCoaddSources.scaleVariance.background.weighting=True

# Mask planes for pixels to ignore when scaling variance
config.detectCoaddSources.scaleVariance.maskPlanes=['DETECTED', 'DETECTED_NEGATIVE', 'BAD', 'SAT', 'NO_DATA', 'INTRP']

# Maximum variance scaling value to permit
config.detectCoaddSources.scaleVariance.limit=10.0

# detected sources with fewer than the specified number of pixels will be ignored
config.detectCoaddSources.detection.minPixels=1

# Pixels should be grown as isotropically as possible (slower)
config.detectCoaddSources.detection.isotropicGrow=True

# Grow all footprints at the same time? This allows disconnected footprints to merge.
config.detectCoaddSources.detection.combinedGrow=True

# Grow detections by nSigmaToGrow * [PSF RMS width]; if 0 then do not grow
config.detectCoaddSources.detection.nSigmaToGrow=2.4

# Grow detections to set the image mask bits, but return the original (not-grown) footprints
config.detectCoaddSources.detection.returnOriginalFootprints=False

# Threshold for footprints
config.detectCoaddSources.detection.thresholdValue=5.0

# Include threshold relative to thresholdValue
config.detectCoaddSources.detection.includeThresholdMultiplier=1.0

# specifies the desired flavor of Threshold
config.detectCoaddSources.detection.thresholdType='pixel_stdev'

# specifies whether to detect positive, or negative sources, or both
config.detectCoaddSources.detection.thresholdPolarity='positive'

# Fiddle factor to add to the background; debugging only
config.detectCoaddSources.detection.adjustBackground=0.0

# Estimate the background again after final source detection?
config.detectCoaddSources.detection.reEstimateBackground=False

# type of statistic to use for grid points
config.detectCoaddSources.detection.background.statisticsProperty='MEANCLIP'

# behaviour if there are too few points in grid for requested interpolation style
config.detectCoaddSources.detection.background.undersampleStyle='REDUCE_INTERP_ORDER'

# how large a region of the sky should be used for each background point
config.detectCoaddSources.detection.background.binSize=4096

# Sky region size to be used for each background point in X direction. If 0, the binSize config is used.
config.detectCoaddSources.detection.background.binSizeX=0

# Sky region size to be used for each background point in Y direction. If 0, the binSize config is used.
config.detectCoaddSources.detection.background.binSizeY=0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
config.detectCoaddSources.detection.background.algorithm='AKIMA_SPLINE'

# Names of mask planes to ignore while estimating the background
config.detectCoaddSources.detection.background.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# Ignore NaNs when estimating the background
config.detectCoaddSources.detection.background.isNanSafe=False

# Use Approximate (Chebyshev) to model background.
config.detectCoaddSources.detection.background.useApprox=False

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.detectCoaddSources.detection.background.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.detectCoaddSources.detection.background.approxOrderY=-1

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.detectCoaddSources.detection.background.weighting=True

# type of statistic to use for grid points
config.detectCoaddSources.detection.tempLocalBackground.statisticsProperty='MEANCLIP'

# behaviour if there are too few points in grid for requested interpolation style
config.detectCoaddSources.detection.tempLocalBackground.undersampleStyle='REDUCE_INTERP_ORDER'

# how large a region of the sky should be used for each background point
config.detectCoaddSources.detection.tempLocalBackground.binSize=64

# Sky region size to be used for each background point in X direction. If 0, the binSize config is used.
config.detectCoaddSources.detection.tempLocalBackground.binSizeX=0

# Sky region size to be used for each background point in Y direction. If 0, the binSize config is used.
config.detectCoaddSources.detection.tempLocalBackground.binSizeY=0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
config.detectCoaddSources.detection.tempLocalBackground.algorithm='AKIMA_SPLINE'

# Names of mask planes to ignore while estimating the background
config.detectCoaddSources.detection.tempLocalBackground.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# Ignore NaNs when estimating the background
config.detectCoaddSources.detection.tempLocalBackground.isNanSafe=False

# Use Approximate (Chebyshev) to model background.
config.detectCoaddSources.detection.tempLocalBackground.useApprox=False

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.detectCoaddSources.detection.tempLocalBackground.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.detectCoaddSources.detection.tempLocalBackground.approxOrderY=-1

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.detectCoaddSources.detection.tempLocalBackground.weighting=True

# Enable temporary local background subtraction? (see tempLocalBackground)
config.detectCoaddSources.detection.doTempLocalBackground=True

# type of statistic to use for grid points
config.detectCoaddSources.detection.tempWideBackground.statisticsProperty='MEANCLIP'

# behaviour if there are too few points in grid for requested interpolation style
config.detectCoaddSources.detection.tempWideBackground.undersampleStyle='REDUCE_INTERP_ORDER'

# how large a region of the sky should be used for each background point
config.detectCoaddSources.detection.tempWideBackground.binSize=512

# Sky region size to be used for each background point in X direction. If 0, the binSize config is used.
config.detectCoaddSources.detection.tempWideBackground.binSizeX=0

# Sky region size to be used for each background point in Y direction. If 0, the binSize config is used.
config.detectCoaddSources.detection.tempWideBackground.binSizeY=0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
config.detectCoaddSources.detection.tempWideBackground.algorithm='AKIMA_SPLINE'

# Names of mask planes to ignore while estimating the background
config.detectCoaddSources.detection.tempWideBackground.ignoredPixelMask=['BAD', 'EDGE', 'NO_DATA']

# Ignore NaNs when estimating the background
config.detectCoaddSources.detection.tempWideBackground.isNanSafe=False

# Use Approximate (Chebyshev) to model background.
config.detectCoaddSources.detection.tempWideBackground.useApprox=False

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.detectCoaddSources.detection.tempWideBackground.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.detectCoaddSources.detection.tempWideBackground.approxOrderY=-1

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.detectCoaddSources.detection.tempWideBackground.weighting=True

# Do temporary wide (large-scale) background subtraction before footprint detection?
config.detectCoaddSources.detection.doTempWideBackground=True

# The maximum number of peaks in a Footprint before trying to replace its peaks using the temporary local background
config.detectCoaddSources.detection.nPeaksMaxSimple=1

# Multiple of PSF RMS size to use for convolution kernel bounding box size; note that this is not a half-size. The size will be rounded up to the nearest odd integer
config.detectCoaddSources.detection.nSigmaForKernel=7.0

# Mask planes to ignore when calculating statistics of image (for thresholdType=stdev)
config.detectCoaddSources.detection.statsMask=['BAD', 'SAT', 'EDGE', 'NO_DATA']

# Fraction of the threshold to use for first pass (to find sky objects)
config.detectCoaddSources.detection.prelimThresholdFactor=0.5

# Avoid pixels masked with these mask planes
config.detectCoaddSources.detection.skyObjects.avoidMask=['DETECTED', 'DETECTED_NEGATIVE', 'BAD', 'NO_DATA']

# Number of pixels to grow the masked pixels when adding sky objects
config.detectCoaddSources.detection.skyObjects.growMask=0

# Radius, in pixels, of sky objects
config.detectCoaddSources.detection.skyObjects.sourceRadius=8.0

# Try to add this many sky objects
config.detectCoaddSources.detection.skyObjects.nSources=1000

# Maximum number of trial sky object positions
# (default: nSkySources*nTrialSkySourcesMultiplier)
config.detectCoaddSources.detection.skyObjects.nTrialSources=None

# Set nTrialSkySources to
#     nSkySources*nTrialSkySourcesMultiplier
# if nTrialSkySources is None
config.detectCoaddSources.detection.skyObjects.nTrialSourcesMultiplier=5

# Tweak background level so median PSF flux of sky objects is zero?
config.detectCoaddSources.detection.doBackgroundTweak=True

# Minimum number of sky sources in statistical sample; if below this number, we refuse to modify the threshold.
config.detectCoaddSources.detection.minNumSources=10

# Name of coadd
config.detectCoaddSources.coaddName='deep'

# Run fake sources injection task
config.detectCoaddSources.doInsertFakes=False

# Mask plane to set on pixels affected by fakes.  Will be added if not already present.
config.detectCoaddSources.insertFakes.maskPlaneName='FAKE'

