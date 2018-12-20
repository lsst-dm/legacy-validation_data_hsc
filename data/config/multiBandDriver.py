import lsst.pipe.drivers.multiBandDriver
assert type(config)==lsst.pipe.drivers.multiBandDriver.MultiBandDriverConfig, 'config is of type %s.%s instead of lsst.pipe.drivers.multiBandDriver.MultiBandDriverConfig' % (type(config).__module__, type(config).__name__)
import lsst.shapelet.radialProfile
import lsst.shapelet.radialProfile.radialProfile
import lsst.meas.modelfit.pixelFitRegion.pixelFitRegionContinued
import lsst.shapelet.multiShapeletFunction
import lsst.shapelet.multiShapeletBasis
import lsst.shapelet.functorKeys
import lsst.meas.modelfit.psf
import lsst.meas.modelfit.cmodel.cmodelContinued
import lsst.meas.extensions.photometryKron.version
import lsst.shapelet.tractor
import lsst.meas.modelfit.priors.priors
import lsst.meas.modelfit.priors.priorsContinued
import lsst.meas.modelfit.mixture
import lsst.shapelet
import lsst.meas.modelfit.optimizer.optimizer
import lsst.meas.extensions.convolved.convolved
import lsst.shapelet.version
import lsst.meas.modelfit.cmodel
import lsst.meas.modelfit.integrals
import lsst.meas.extensions.shapeHSM.version
import lsst.shapelet.constants.constants
import lsst.meas.modelfit.optimizer
import lsst.meas.modelfit.version
import lsst.shapelet.constants
import lsst.shapelet.constants.constantsContinued
import lsst.meas.modelfit.common
import lsst.shapelet.gaussHermiteProjection
import lsst.shapelet.basisEvaluator
import lsst.meas.extensions.shapeHSM
import lsst.meas.extensions.shapeHSM.hsmShapeControl
import lsst.meas.modelfit.truncatedGaussian
import lsst.meas.extensions.photometryKron.photometryKron
import lsst.meas.extensions
import lsst.meas.modelfit.multiModel
import lsst.meas.modelfit.pixelFitRegion
import lsst.meas.modelfit.pixelFitRegion.pixelFitRegion
import lsst.meas.extensions.shapeHSM.hsmMomentsControl
import lsst.meas.modelfit.psf.psf
import lsst.meas.modelfit.sampler
import lsst.meas.extensions.convolved
import lsst.shapelet.shapeletFunction
import lsst.shapelet.generator
import lsst.meas.modelfit.likelihood
import lsst.meas.modelfit.optimizer.optimizerContinued
import lsst.shapelet.multiShapeletFunction.multiShapeletFunctionContinued
import lsst.meas.modelfit.priors
import lsst.meas.modelfit.adaptiveImportanceSampler
import lsst.meas.modelfit
import lsst.shapelet.multiShapeletFunction.multiShapeletFunction
import lsst.meas.modelfit.cmodel.cmodel
import lsst.meas.modelfit.model
import lsst.shapelet.radialProfile.radialProfileContinued
import lsst.shapelet.shapeletFunction.shapeletFunctionContinued
import lsst.meas.modelfit.unitSystem
import lsst.shapelet.hermiteTransformMatrix
import lsst.meas.modelfit.unitTransformedLikelihood
import lsst.shapelet.gaussHermiteConvolution
import lsst.meas.modelfit.psf.psfContinued
import lsst.meas.extensions.photometryKron
import lsst.shapelet.shapeletFunction.shapeletFunction
import lsst.shapelet.matrixBuilder
import lsst.meas.extensions.convolved.version
# Name of coadd
config.coaddName='deep'

# Re-run detection? (requires *Coadd dataset to have been written)
config.doDetection=False

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

# Priority-ordered list of bands for the merge.
config.mergeCoaddDetections.priorityList=['HSC-I2', 'HSC-I', 'HSC-R2', 'HSC-R', 'HSC-Z', 'HSC-Y', 'HSC-G', 'NB0921', 'NB0816', 'NB1010', 'NB0387', 'NB0515']

# Name of coadd
config.mergeCoaddDetections.coaddName='deep'

# Minimum distance from closest peak to create a new one (in arcsec).
config.mergeCoaddDetections.minNewPeak=1.0

# When adding new catalogs to the merge, all peaks less than this distance  (in arcsec) to an existing peak will be flagged as detected in that catalog.
config.mergeCoaddDetections.maxSamePeak=0.3

# Always keep peaks detected in this many bands
config.mergeCoaddDetections.cullPeaks.nBandsSufficient=2

# Always keep this many peaks in each family
config.mergeCoaddDetections.cullPeaks.rankSufficient=20

# Keep peaks with less than this rank that also match the rankNormalizedConsidered condition.
config.mergeCoaddDetections.cullPeaks.rankConsidered=30

# Keep peaks with less than this normalized rank that also match the rankConsidered condition.
config.mergeCoaddDetections.cullPeaks.rankNormalizedConsidered=0.7

# Name of `filter' used to label sky objects (e.g. flag merge_peak_sky is set)
# (N.b. should be in MergeMeasurementsConfig.pseudoFilterList)
config.mergeCoaddDetections.skyFilterName='sky'

# Avoid pixels masked with these mask planes
config.mergeCoaddDetections.skyObjects.avoidMask=['DETECTED']

# Number of pixels to grow the masked pixels when adding sky objects
config.mergeCoaddDetections.skyObjects.growMask=0

# Radius, in pixels, of sky objects
config.mergeCoaddDetections.skyObjects.sourceRadius=8.0

# Try to add this many sky objects
config.mergeCoaddDetections.skyObjects.nSources=100

# Maximum number of trial sky object positions
# (default: nSkySources*nTrialSkySourcesMultiplier)
config.mergeCoaddDetections.skyObjects.nTrialSources=None

# Set nTrialSkySources to
#     nSkySources*nTrialSkySourcesMultiplier
# if nTrialSkySources is None
config.mergeCoaddDetections.skyObjects.nTrialSourcesMultiplier=5

# What to do when a peak to be deblended is close to the edge of the image
config.deblendCoaddSources.singleBandDeblend.edgeHandling='ramp'

# When the deblender should attribute stray flux to point sources
config.deblendCoaddSources.singleBandDeblend.strayFluxToPointSources='necessary'

# Assign stray flux (not claimed by any child in the deblender) to deblend children.
config.deblendCoaddSources.singleBandDeblend.assignStrayFlux=True

# How to split flux among peaks
config.deblendCoaddSources.singleBandDeblend.strayFluxRule='trim'

# When splitting stray flux, clip fractions below this value to zero.
config.deblendCoaddSources.singleBandDeblend.clipStrayFluxFraction=0.001

# Chi-squared per DOF cut for deciding a source is a PSF during deblending (un-shifted PSF model)
config.deblendCoaddSources.singleBandDeblend.psfChisq1=1.5

# Chi-squared per DOF cut for deciding a source is PSF during deblending (shifted PSF model)
config.deblendCoaddSources.singleBandDeblend.psfChisq2=1.5

# Chi-squared per DOF cut for deciding a source is a PSF during deblending (shifted PSF model #2)
config.deblendCoaddSources.singleBandDeblend.psfChisq2b=1.5

# Only deblend the brightest maxNumberOfPeaks peaks in the parent (<= 0: unlimited)
config.deblendCoaddSources.singleBandDeblend.maxNumberOfPeaks=0

# Maximum area for footprints before they are ignored as large; non-positive means no threshold applied
config.deblendCoaddSources.singleBandDeblend.maxFootprintArea=1000000

# Maximum linear dimension for footprints before they are ignored as large; non-positive means no threshold applied
config.deblendCoaddSources.singleBandDeblend.maxFootprintSize=0

# Minimum axis ratio for footprints before they are ignored as large; non-positive means no threshold applied
config.deblendCoaddSources.singleBandDeblend.minFootprintAxisRatio=0.0

# Mask name for footprints not deblended, or None
config.deblendCoaddSources.singleBandDeblend.notDeblendedMask='NOT_DEBLENDED'

# Footprints smaller in width or height than this value will be ignored; minimum of 2 due to PSF gradient calculation.
config.deblendCoaddSources.singleBandDeblend.tinyFootprintSize=2

# Guarantee that all peaks produce a child source.
config.deblendCoaddSources.singleBandDeblend.propagateAllPeaks=True

# If True, catch exceptions thrown by the deblender, log them, and set a flag on the parent, instead of letting them propagate up
config.deblendCoaddSources.singleBandDeblend.catchFailures=False

# Mask planes to ignore when performing statistics
config.deblendCoaddSources.singleBandDeblend.maskPlanes=['SAT', 'INTRP', 'NO_DATA']

# Mask planes with the corresponding limit on the fraction of masked pixels. Sources violating this limit will not be deblended.
config.deblendCoaddSources.singleBandDeblend.maskLimits={}

# If true, a least-squares fit of the templates will be done to the full image. The templates will be re-weighted based on this fit.
config.deblendCoaddSources.singleBandDeblend.weightTemplates=False

# Try to remove similar templates?
config.deblendCoaddSources.singleBandDeblend.removeDegenerateTemplates=False

# If the dot product between two templates is larger than this value, we consider them to be describing the same object (i.e. they are degenerate).  If one of the objects has been labeled as a PSF it will be removed, otherwise the template with the lowest value will be removed.
config.deblendCoaddSources.singleBandDeblend.maxTempDotProd=0.5

# Apply a smoothing filter to all of the template images
config.deblendCoaddSources.singleBandDeblend.medianSmoothTemplate=True

# Maximum number of iterations to deblend a single parent
config.deblendCoaddSources.multiBandDeblend.maxIter=200

# Relative error to use when determining stopping criteria
config.deblendCoaddSources.multiBandDeblend.relativeError=0.001

# A peak must be updated by at least 'minTranslation' (pixels)or no update is performed.This field is ignored if fitPositions is False.
config.deblendCoaddSources.multiBandDeblend.minTranslation=0.001

# If fitPositions is True, the positions and box sizes areupdated on every 'refinementSkip' iterations.
config.deblendCoaddSources.multiBandDeblend.refinementSkip=10

# Method to use for fitting translations.Currently 'default' is the only available option,which performs a linear fit, but it is possible that wewill use galsim or some other method as a future option
config.deblendCoaddSources.multiBandDeblend.translationMethod='default'

# Boxes are resized when the flux at an edge is > edgeFluxThresh * background RMS
config.deblendCoaddSources.multiBandDeblend.edgeFluxThresh=1.0

# Calculate exact Lipschitz constant in every step(True) or only calculate the approximateLipschitz constant with significant changes in A,S(False)
config.deblendCoaddSources.multiBandDeblend.exactLipschitz=False

# A fractional measure of how much a value (like the exactLipschitz)can change before it needs to be recalculated.This must be between 0 and 1.
config.deblendCoaddSources.multiBandDeblend.stepSlack=0.2

# List of constraints to use for each object(order does not matter)Current options are all used by default:
# S: symmetry
# M: monotonicity
# 1: normalized SED to unity+: non-negative morphology
config.deblendCoaddSources.multiBandDeblend.constraints='1,+,S,M'

# Strictness of symmetry, from0 (no symmetry enforced) to1 (perfect symmetry required).If 'S' is not in `constraints`, this argument is ignored
config.deblendCoaddSources.multiBandDeblend.symmetryThresh=1.0

# L0 threshold. NaN results in no L0 penalty.
config.deblendCoaddSources.multiBandDeblend.l0Thresh=float('nan')

# L1 threshold. NaN results in no L1 penalty.
config.deblendCoaddSources.multiBandDeblend.l1Thresh=float('nan')

# Threshold for TV (total variation) constraint in the x-direction.NaN results in no TVx penalty.
config.deblendCoaddSources.multiBandDeblend.tvxThresh=float('nan')

# Threshold for TV (total variation) constraint in the y-direction.NaN results in no TVy penalty.
config.deblendCoaddSources.multiBandDeblend.tvyThresh=float('nan')

# Use inverse variance as deblender weights
config.deblendCoaddSources.multiBandDeblend.useWeights=False

# Fraction of background RMS level to use as acutoff for defining the background of the imageThis is used to initialize the model for each sourceand to set the size of the bounding box for each sourceevery `refinementSkip` iteration.
config.deblendCoaddSources.multiBandDeblend.bgScale=0.5

# Whether or not to convolve the morphology with thePSF in each band or use the same morphology in all bands
config.deblendCoaddSources.multiBandDeblend.usePsfConvolution=True

# Whether or not to save the SEDs and templates
config.deblendCoaddSources.multiBandDeblend.saveTemplates=True

# Whether or not to process isolated sources in the deblender
config.deblendCoaddSources.multiBandDeblend.processSingles=False

# Whether or not to process isolated sources in the deblender
config.deblendCoaddSources.multiBandDeblend.badMask='BAD,CR,NO_DATA,SAT,SUSPECT'

# Only deblend the brightest maxNumberOfPeaks peaks in the parent (<= 0: unlimited)
config.deblendCoaddSources.multiBandDeblend.maxNumberOfPeaks=0

# Maximum area for footprints before they are ignored as large; non-positive means no threshold applied
config.deblendCoaddSources.multiBandDeblend.maxFootprintArea=1000000

# Maximum linear dimension for footprints before they are ignored as large; non-positive means no threshold applied
config.deblendCoaddSources.multiBandDeblend.maxFootprintSize=0

# Minimum axis ratio for footprints before they are ignored as large; non-positive means no threshold applied
config.deblendCoaddSources.multiBandDeblend.minFootprintAxisRatio=0.0

# Mask name for footprints not deblended, or None
config.deblendCoaddSources.multiBandDeblend.notDeblendedMask='NOT_DEBLENDED'

# Footprints smaller in width or height than this value will be ignored; minimum of 2 due to PSF gradient calculation.
config.deblendCoaddSources.multiBandDeblend.tinyFootprintSize=2

# If True, catch exceptions thrown by the deblender, log them, and set a flag on the parent, instead of letting them propagate up
config.deblendCoaddSources.multiBandDeblend.catchFailures=False

# Guarantee that all peaks produce a child source.
config.deblendCoaddSources.multiBandDeblend.propagateAllPeaks=False

# Mask planes to ignore when performing statistics
config.deblendCoaddSources.multiBandDeblend.maskPlanes=['SAT', 'INTRP', 'NO_DATA']

# Mask planes with the corresponding limit on the fraction of masked pixels. Sources violating this limit will not be deblended.
config.deblendCoaddSources.multiBandDeblend.maskLimits={}

# What to do when a peak to be deblended is close to the edge of the image
config.deblendCoaddSources.multiBandDeblend.edgeHandling='noclip'

# Apply a smoothing filter to all of the template images
config.deblendCoaddSources.multiBandDeblend.medianSmoothTemplate=False

# Half size of the median smoothing filter
config.deblendCoaddSources.multiBandDeblend.medianFilterHalfsize=2.0

# Clip non-zero spans in the footprints
config.deblendCoaddSources.multiBandDeblend.clipFootprintToNonzero=False

# Reapportion flux to the footprints so that flux is conserved
config.deblendCoaddSources.multiBandDeblend.conserveFlux=True

# If true, a least-squares fit of the templates will be done to the full image. The templates will be re-weighted based on this fit.
config.deblendCoaddSources.multiBandDeblend.weightTemplates=False

# When the deblender should attribute stray flux to point sources
config.deblendCoaddSources.multiBandDeblend.strayFluxToPointSources='necessary'

# Assign stray flux (not claimed by any child in the deblender) to deblend children.
config.deblendCoaddSources.multiBandDeblend.assignStrayFlux=True

# How to split flux among peaks
config.deblendCoaddSources.multiBandDeblend.strayFluxRule='trim'

# When splitting stray flux, clip fractions below this value to zero.
config.deblendCoaddSources.multiBandDeblend.clipStrayFluxFraction=0.001

# As part of the flux calculation, the sum of the templates iscalculated. If 'getTemplateSum==True' then the sum of thetemplates is stored in the result (a 'PerFootprint').
config.deblendCoaddSources.multiBandDeblend.getTemplateSum=False

# Simultaneously deblend all bands?
config.deblendCoaddSources.simultaneous=False

# Name of coadd
config.deblendCoaddSources.coaddName='deep'

import lsst.meas.extensions.convolved.version
import lsst.meas.modelfit.pixelFitRegion.pixelFitRegionContinued
import lsst.shapelet.multiShapeletFunction
import lsst.meas.modelfit.psf
import lsst.meas.extensions.photometryKron.version
import lsst.meas.modelfit.priors.priors
import lsst.meas.modelfit.mixture
import lsst.shapelet
import lsst.meas.modelfit.optimizer.optimizer
import lsst.meas.extensions.convolved.convolved
import lsst.shapelet.constants.constants
import lsst.meas.modelfit.optimizer
import lsst.shapelet.constants.constantsContinued
import lsst.meas.extensions.shapeHSM.hsmShapeControl
import lsst.meas.extensions.photometryKron.photometryKron
import lsst.meas.extensions
import lsst.meas.modelfit.pixelFitRegion
import lsst.meas.modelfit.sampler
import lsst.shapelet.generator
import lsst.meas.extensions.convolved
import lsst.meas.modelfit.optimizer.optimizerContinued
import lsst.meas.modelfit.priors
import lsst.meas.modelfit.adaptiveImportanceSampler
import lsst.meas.modelfit.cmodel.cmodel
import lsst.meas.modelfit.model
import lsst.shapelet.shapeletFunction.shapeletFunctionContinued
import lsst.shapelet.hermiteTransformMatrix
import lsst.shapelet.gaussHermiteConvolution
import lsst.meas.extensions.photometryKron
import lsst.shapelet.shapeletFunction.shapeletFunction
import lsst.shapelet.radialProfile
import lsst.shapelet.radialProfile.radialProfile
import lsst.shapelet.multiShapeletBasis
import lsst.shapelet.functorKeys
import lsst.meas.modelfit.cmodel.cmodelContinued
import lsst.shapelet.tractor
import lsst.meas.modelfit.priors.priorsContinued
import lsst.meas.modelfit.cmodel
import lsst.meas.modelfit.integrals
import lsst.shapelet.version
import lsst.meas.extensions.shapeHSM.version
import lsst.meas.modelfit.version
import lsst.shapelet.constants
import lsst.meas.modelfit.common
import lsst.shapelet.gaussHermiteProjection
import lsst.meas.extensions.shapeHSM
import lsst.meas.modelfit.truncatedGaussian
import lsst.meas.modelfit.multiModel
import lsst.meas.modelfit.pixelFitRegion.pixelFitRegion
import lsst.meas.modelfit.psf.psf
import lsst.meas.extensions.shapeHSM.hsmMomentsControl
import lsst.shapelet.shapeletFunction
import lsst.meas.modelfit.likelihood
import lsst.shapelet.multiShapeletFunction.multiShapeletFunctionContinued
import lsst.meas.modelfit
import lsst.shapelet.multiShapeletFunction.multiShapeletFunction
import lsst.shapelet.radialProfile.radialProfileContinued
import lsst.meas.modelfit.unitSystem
import lsst.meas.modelfit.unitTransformedLikelihood
import lsst.meas.modelfit.psf.psfContinued
import lsst.shapelet.matrixBuilder
import lsst.shapelet.basisEvaluator
# Name of the input catalog to use.If the single band deblender was used this should be 'deblendedFlux.If the multi-band deblender was used this should be 'deblendedModel.If no deblending was performed this should be 'mergeDet'
config.measureCoaddSources.inputCatalog='deblendedFlux'

import lsst.meas.extensions.shapeHSM.hsmShapeControl
import lsst.meas.extensions.photometryKron.photometryKron
import lsst.meas.extensions
import lsst.meas.extensions.convolved.convolved
import lsst.meas.extensions.shapeHSM.version
import lsst.meas.extensions.shapeHSM.hsmMomentsControl
import lsst.meas.extensions.convolved
import lsst.meas.extensions.photometryKron
import lsst.meas.extensions.photometryKron.version
import lsst.meas.extensions.shapeHSM
import lsst.meas.extensions.convolved.version
# the name of the centroiding algorithm used to set source x,y
config.measureCoaddSources.measurement.slots.centroid='base_SdssCentroid'

# the name of the algorithm used to set source moments parameters
config.measureCoaddSources.measurement.slots.shape='ext_shapeHSM_HsmSourceMoments'

# the name of the algorithm used to set PSF moments parameters
config.measureCoaddSources.measurement.slots.psfShape='ext_shapeHSM_HsmPsfMoments'

# the name of the algorithm used to set the source aperture flux slot
config.measureCoaddSources.measurement.slots.apFlux='base_CircularApertureFlux_12_0'

# the name of the algorithm used to set the source model flux slot
config.measureCoaddSources.measurement.slots.modelFlux='modelfit_CModel'

# the name of the algorithm used to set the source psf flux slot
config.measureCoaddSources.measurement.slots.psfFlux='base_PsfFlux'

# the name of the algorithm used to set the source inst flux slot
config.measureCoaddSources.measurement.slots.instFlux='base_GaussianFlux'

# the name of the flux measurement algorithm used for calibration
config.measureCoaddSources.measurement.slots.calibFlux='base_CircularApertureFlux_12_0'

# When measuring, replace other detected footprints with noise?
config.measureCoaddSources.measurement.doReplaceWithNoise=True

# How to choose mean and variance of the Gaussian noise we generate?
config.measureCoaddSources.measurement.noiseReplacer.noiseSource='measure'

# Add ann offset to the generated noise.
config.measureCoaddSources.measurement.noiseReplacer.noiseOffset=0.0

# The seed multiplier value to use for random number generation
#    >= 1: set the seed deterministically based on exposureId
#       0: fall back to the afw.math.Random default constructor (which uses a seed value of 1)
config.measureCoaddSources.measurement.noiseReplacer.noiseSeedMultiplier=1

# Prefix to give undeblended plugins
config.measureCoaddSources.measurement.undeblendedPrefix='undeblended_'

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.plugins['base_PsfFlux'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.measureCoaddSources.measurement.plugins['base_PsfFlux'].badMaskPlanes=[]

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.plugins['base_PeakLikelihoodFlux'].doMeasure=True

# Name of warping kernel (e.g. "lanczos4") used to compute the peak
config.measureCoaddSources.measurement.plugins['base_PeakLikelihoodFlux'].warpingKernelName='lanczos4'

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.plugins['base_GaussianFlux'].doMeasure=True

# FIXME! NEVER DOCUMENTED!
config.measureCoaddSources.measurement.plugins['base_GaussianFlux'].background=0.0

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.plugins['base_NaiveCentroid'].doMeasure=True

# Value to subtract from the image pixel values
config.measureCoaddSources.measurement.plugins['base_NaiveCentroid'].background=0.0

# Do check that the centroid is contained in footprint.
config.measureCoaddSources.measurement.plugins['base_NaiveCentroid'].doFootprintCheck=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.measureCoaddSources.measurement.plugins['base_NaiveCentroid'].maxDistToPeak=-1.0

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.plugins['base_SdssCentroid'].doMeasure=True

# maximum allowed binning
config.measureCoaddSources.measurement.plugins['base_SdssCentroid'].binmax=16

# Do check that the centroid is contained in footprint.
config.measureCoaddSources.measurement.plugins['base_SdssCentroid'].doFootprintCheck=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.measureCoaddSources.measurement.plugins['base_SdssCentroid'].maxDistToPeak=-1.0

# if the peak's less than this insist on binning at least once
config.measureCoaddSources.measurement.plugins['base_SdssCentroid'].peakMin=-1.0

# fiddle factor for adjusting the binning
config.measureCoaddSources.measurement.plugins['base_SdssCentroid'].wfac=1.5

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.plugins['base_PixelFlags'].doMeasure=True

# List of mask planes to be searched for which occur anywhere within a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.measureCoaddSources.measurement.plugins['base_PixelFlags'].masksFpAnywhere=['CLIPPED', 'SENSOR_EDGE', 'INEXACT_PSF', 'BRIGHT_OBJECT']

# List of mask planes to be searched for which occur in the center of a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.measureCoaddSources.measurement.plugins['base_PixelFlags'].masksFpCenter=['CLIPPED', 'SENSOR_EDGE', 'INEXACT_PSF', 'BRIGHT_OBJECT']

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.plugins['base_SdssShape'].doMeasure=True

# Additional value to add to background
config.measureCoaddSources.measurement.plugins['base_SdssShape'].background=0.0

# Whether to also compute the shape of the PSF model
config.measureCoaddSources.measurement.plugins['base_SdssShape'].doMeasurePsf=True

# Maximum number of iterations
config.measureCoaddSources.measurement.plugins['base_SdssShape'].maxIter=100

# Maximum centroid shift, limited to 2-10
config.measureCoaddSources.measurement.plugins['base_SdssShape'].maxShift=0.0

# Convergence tolerance for e1,e2
config.measureCoaddSources.measurement.plugins['base_SdssShape'].tol1=9.999999747378752e-06

# Convergence tolerance for FWHM
config.measureCoaddSources.measurement.plugins['base_SdssShape'].tol2=9.999999747378752e-05

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.plugins['base_ScaledApertureFlux'].doMeasure=True

# Scaling factor of PSF FWHM for aperture radius.
config.measureCoaddSources.measurement.plugins['base_ScaledApertureFlux'].scale=3.14

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.measureCoaddSources.measurement.plugins['base_ScaledApertureFlux'].shiftKernel='lanczos5'

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.plugins['base_CircularApertureFlux'].doMeasure=True

# Maximum radius (in pixels) for which the sinc algorithm should be used instead of the faster naive algorithm.  For elliptical apertures, this is the minor axis radius.
config.measureCoaddSources.measurement.plugins['base_CircularApertureFlux'].maxSincRadius=12.0

# Radius (in pixels) of apertures.
config.measureCoaddSources.measurement.plugins['base_CircularApertureFlux'].radii=[3.0, 4.5, 6.0, 9.0, 12.0, 17.0, 25.0, 35.0, 50.0, 70.0]

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.measureCoaddSources.measurement.plugins['base_CircularApertureFlux'].shiftKernel='lanczos5'

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.plugins['base_Blendedness'].doMeasure=True

# Whether to compute quantities related to the Gaussian-weighted flux
config.measureCoaddSources.measurement.plugins['base_Blendedness'].doFlux=True

# Whether to compute HeavyFootprint dot products (the old deblend.blendedness parameter)
config.measureCoaddSources.measurement.plugins['base_Blendedness'].doOld=True

# Whether to compute quantities related to the Gaussian-weighted shape
config.measureCoaddSources.measurement.plugins['base_Blendedness'].doShape=True

# Radius factor that sets the maximum extent of the weight function (and hence the flux measurements)
config.measureCoaddSources.measurement.plugins['base_Blendedness'].nSigmaWeightMax=3.0

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.plugins['base_LocalBackground'].doMeasure=True

# Inner radius for background annulus as a multiple of the PSF sigma
config.measureCoaddSources.measurement.plugins['base_LocalBackground'].annulusInner=7.0

# Outer radius for background annulus as a multiple of the PSF sigma
config.measureCoaddSources.measurement.plugins['base_LocalBackground'].annulusOuter=15.0

# Mask planes that indicate pixels that should be excluded from the measurement
config.measureCoaddSources.measurement.plugins['base_LocalBackground'].badMaskPlanes=['BAD', 'SAT', 'NO_DATA']

# Number of sigma-clipping iterations for background measurement
config.measureCoaddSources.measurement.plugins['base_LocalBackground'].bgIter=3

# Rejection threshold (in standard deviations) for background measurement
config.measureCoaddSources.measurement.plugins['base_LocalBackground'].bgRej=3.0

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.plugins['base_FPPosition'].doMeasure=True

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.plugins['base_Jacobian'].doMeasure=True

# Nominal pixel size (arcsec)
config.measureCoaddSources.measurement.plugins['base_Jacobian'].pixelScale=0.5

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.plugins['base_Variance'].doMeasure=True

# Scale factor to apply to shape for aperture
config.measureCoaddSources.measurement.plugins['base_Variance'].scale=5.0

# Mask planes to ignore
config.measureCoaddSources.measurement.plugins['base_Variance'].mask=['DETECTED', 'DETECTED_NEGATIVE', 'BAD', 'SAT']

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.plugins['base_InputCount'].doMeasure=True

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.plugins['base_PeakCentroid'].doMeasure=True

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.plugins['base_SkyCoord'].doMeasure=True

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.plugins['subaru_FilterFraction'].doMeasure=True

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.plugins['ext_photometryKron_KronFlux'].doMeasure=True

# If true check that the Kron radius exceeds some minimum
config.measureCoaddSources.measurement.plugins['ext_photometryKron_KronFlux'].enforceMinimumRadius=True

# if true, use existing shape and centroid measurements instead of fitting
config.measureCoaddSources.measurement.plugins['ext_photometryKron_KronFlux'].fixed=False

# Largest aperture for which to use the slow, accurate, sinc aperture code
config.measureCoaddSources.measurement.plugins['ext_photometryKron_KronFlux'].maxSincRadius=10.0

# Minimum Kron radius (if == 0.0 use PSF's Kron radius) if enforceMinimumRadius. Also functions as fallback aperture radius if set.
config.measureCoaddSources.measurement.plugins['ext_photometryKron_KronFlux'].minimumRadius=0.0

# Number of times to iterate when setting the Kron radius
config.measureCoaddSources.measurement.plugins['ext_photometryKron_KronFlux'].nIterForRadius=1

# Number of Kron radii for Kron flux
config.measureCoaddSources.measurement.plugins['ext_photometryKron_KronFlux'].nRadiusForFlux=2.5

# Multiplier of rms size for aperture used to initially estimate the Kron radius
config.measureCoaddSources.measurement.plugins['ext_photometryKron_KronFlux'].nSigmaForRadius=6.0

# Name of field specifying reference Kron radius for forced measurement
config.measureCoaddSources.measurement.plugins['ext_photometryKron_KronFlux'].refRadiusName='ext_photometryKron_KronFlux_radius'

# Smooth image with N(0, smoothingSigma^2) Gaussian while estimating R_K
config.measureCoaddSources.measurement.plugins['ext_photometryKron_KronFlux'].smoothingSigma=-1.0

# Use the Footprint size as part of initial estimate of Kron radius
config.measureCoaddSources.measurement.plugins['ext_photometryKron_KronFlux'].useFootprintRadius=False

# list of target seeings (FWHM, pixels)
config.measureCoaddSources.measurement.plugins['ext_convolved_ConvolvedFlux'].seeing=[3.5, 5.0, 6.5, 8.0]

# scaling factor of kernel sigma for kernel size
config.measureCoaddSources.measurement.plugins['ext_convolved_ConvolvedFlux'].kernelScale=4.0

# Maximum radius (in pixels) for which the sinc algorithm should be used instead of the faster naive algorithm.  For elliptical apertures, this is the minor axis radius.
config.measureCoaddSources.measurement.plugins['ext_convolved_ConvolvedFlux'].aperture.maxSincRadius=10.0

# Radius (in pixels) of apertures.
config.measureCoaddSources.measurement.plugins['ext_convolved_ConvolvedFlux'].aperture.radii=[3.3, 4.5, 6.0]

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.measureCoaddSources.measurement.plugins['ext_convolved_ConvolvedFlux'].aperture.shiftKernel='lanczos5'

# name of Kron radius field in reference
config.measureCoaddSources.measurement.plugins['ext_convolved_ConvolvedFlux'].kronRadiusName='ext_photometryKron_KronFlux_radius'

# Largest aperture for which to use the sinc aperture code for Kron (pixels)
config.measureCoaddSources.measurement.plugins['ext_convolved_ConvolvedFlux'].maxSincRadius=10.0

# Number of Kron radii for Kron flux
config.measureCoaddSources.measurement.plugins['ext_convolved_ConvolvedFlux'].kronRadiusForFlux=2.5

# Register measurements for aperture correction?
# The aperture correction registration is done when the plugin is
# instantiated because the column names are derived from the configuration
# rather than being static. Sometimes you want to turn this off, e.g.,
# when you will use aperture corrections derived from somewhere else
# through the 'proxy' mechanism.
config.measureCoaddSources.measurement.plugins['ext_convolved_ConvolvedFlux'].registerForApCorr=True

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.plugins['ext_convolved_ConvolvedFlux'].doMeasure=True

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.plugins['ext_shapeHSM_HsmShapeBj'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.measureCoaddSources.measurement.plugins['ext_shapeHSM_HsmShapeBj'].badMaskPlanes=None

# Field name for number of deblend children
config.measureCoaddSources.measurement.plugins['ext_shapeHSM_HsmShapeBj'].deblendNChild=None

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.plugins['ext_shapeHSM_HsmShapeLinear'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.measureCoaddSources.measurement.plugins['ext_shapeHSM_HsmShapeLinear'].badMaskPlanes=None

# Field name for number of deblend children
config.measureCoaddSources.measurement.plugins['ext_shapeHSM_HsmShapeLinear'].deblendNChild=None

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.plugins['ext_shapeHSM_HsmShapeKsb'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.measureCoaddSources.measurement.plugins['ext_shapeHSM_HsmShapeKsb'].badMaskPlanes=None

# Field name for number of deblend children
config.measureCoaddSources.measurement.plugins['ext_shapeHSM_HsmShapeKsb'].deblendNChild=None

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.plugins['ext_shapeHSM_HsmShapeRegauss'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.measureCoaddSources.measurement.plugins['ext_shapeHSM_HsmShapeRegauss'].badMaskPlanes=None

# Field name for number of deblend children
config.measureCoaddSources.measurement.plugins['ext_shapeHSM_HsmShapeRegauss'].deblendNChild='deblend_nChild'

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.plugins['ext_shapeHSM_HsmSourceMoments'].doMeasure=True

# Store measured flux?
config.measureCoaddSources.measurement.plugins['ext_shapeHSM_HsmSourceMoments'].addFlux=None

# Mask planes used to reject bad pixels.
config.measureCoaddSources.measurement.plugins['ext_shapeHSM_HsmSourceMoments'].badMaskPlanes=None

# Use round weight function?
config.measureCoaddSources.measurement.plugins['ext_shapeHSM_HsmSourceMoments'].roundMoments=None

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.plugins['ext_shapeHSM_HsmSourceMomentsRound'].doMeasure=True

# Store measured flux?
config.measureCoaddSources.measurement.plugins['ext_shapeHSM_HsmSourceMomentsRound'].addFlux=None

# Mask planes used to reject bad pixels.
config.measureCoaddSources.measurement.plugins['ext_shapeHSM_HsmSourceMomentsRound'].badMaskPlanes=None

# Use round weight function?
config.measureCoaddSources.measurement.plugins['ext_shapeHSM_HsmSourceMomentsRound'].roundMoments=None

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.plugins['ext_shapeHSM_HsmPsfMoments'].doMeasure=True

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].doMeasure=True

# Shapelet order of inner expansion (0 == Gaussian)
config.measureCoaddSources.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].innerOrder=2

# Don't allow the semi-major radius of any component to go above this fraction of the PSF image width
config.measureCoaddSources.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].maxRadiusBoxFraction=0.4

# Don't allow the semi-minor radius of any component to drop below this value (pixels)
config.measureCoaddSources.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].minRadius=1.0

# Don't allow the determinant radii of the two components to differ by less than this (pixels)
config.measureCoaddSources.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].minRadiusDiff=0.5

# whether to save all iterations for debugging purposes
config.measureCoaddSources.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.measureCoaddSources.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.measureCoaddSources.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.maxInnerIterations=20

# maximum number of steps
config.measureCoaddSources.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.measureCoaddSources.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.measureCoaddSources.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.measureCoaddSources.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.measureCoaddSources.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.measureCoaddSources.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.measureCoaddSources.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.measureCoaddSources.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.measureCoaddSources.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.measureCoaddSources.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.measureCoaddSources.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.measureCoaddSources.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.measureCoaddSources.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.measureCoaddSources.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.measureCoaddSources.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionSolverTolerance=1e-08

# Shapelet order of outer expansion (0 == Gaussian)
config.measureCoaddSources.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].outerOrder=1

# Initial outer Gaussian peak height divided by inner Gaussian peak height
config.measureCoaddSources.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].peakRatio=0.1

# Initial outer radius divided by inner radius
config.measureCoaddSources.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].radiusRatio=2.0

config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models={}
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.maxInnerIterations=20

# maximum number of steps
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.radiusPriorSigma=0.5

config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.maxInnerIterations=20

# maximum number of steps
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.radiusPriorSigma=0.5

config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.maxInnerIterations=20

# maximum number of steps
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.order=2

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.order=1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.radiusPriorSigma=0.5

config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.maxInnerIterations=20

# maximum number of steps
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.order=4

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.order=4

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.radiusPriorSigma=0.5

# a sequence of model names indicating which models should be fit, and their order
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].sequence=['DoubleShapelet']

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].doMeasure=True

# Whether to record the steps the optimizer takes (or just the number, if running as a plugin)
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].dev.doRecordHistory=True

# Whether to record the time spent in this stage
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].dev.doRecordTime=True

# Softened core width for ellipticity distribution (conformal shear units).
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].dev.empiricalPriorConfig.ellipticityCore=0.001

# Width of exponential ellipticity distribution (conformal shear units).
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].dev.empiricalPriorConfig.ellipticitySigma=0.3

# ln(radius) at which the softened cutoff begins towards the minimum
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusMinInner=-6.0

# Minimum ln(radius).
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusMinOuter=-6.001

# Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T.
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusMu=-1.0

# Number of degrees of freedom for the Student's T distribution on ln(radius).
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusNu=50.0

# Width of the Student's T distribution in ln(radius).
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusSigma=0.45

# Ellipticity magnitude (conformal shear units) at which the softened cutoff begins
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].dev.linearPriorConfig.ellipticityMaxInner=2.0

# Maximum ellipticity magnitude (conformal shear units)
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].dev.linearPriorConfig.ellipticityMaxOuter=2.001

# ln(radius) at which the softened cutoff begins towards the maximum
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].dev.linearPriorConfig.logRadiusMaxInner=3.0

# Maximum ln(radius)
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].dev.linearPriorConfig.logRadiusMaxOuter=3.001

# ln(radius) at which the softened cutoff begins towards the minimum
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].dev.linearPriorConfig.logRadiusMinInner=-6.0

# The ratio P(logRadiusMinInner)/P(logRadiusMaxInner)
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].dev.linearPriorConfig.logRadiusMinMaxRatio=1.0

# Minimum ln(radius)
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].dev.linearPriorConfig.logRadiusMinOuter=-6.001

# Maximum radius used in approximating profile with Gaussians (0=default for this profile)
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].dev.maxRadius=0

# Number of Gaussian used to approximate the profile
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].dev.nComponents=8

# whether to save all iterations for debugging purposes
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].dev.optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].dev.optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].dev.optimizer.maxInnerIterations=20

# maximum number of steps
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].dev.optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].dev.optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].dev.optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].dev.optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].dev.optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].dev.optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].dev.optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].dev.optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].dev.optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].dev.optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].dev.optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].dev.optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].dev.optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].dev.optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].dev.optimizer.trustRegionSolverTolerance=1e-08

# Name of the Prior that defines the model to fit (a filename in $MEAS_MODELFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting.
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].dev.priorName=''

# One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].dev.priorSource='EMPIRICAL'

# Name of the shapelet.RadialProfile that defines the model to fit
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].dev.profileName='luv'

# Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances)
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].dev.usePixelWeights=False

# Scale the likelihood by this factor to artificially reweight it w.r.t. the prior.
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].dev.weightsMultiplier=1.0

# Whether to record the steps the optimizer takes (or just the number, if running as a plugin)
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].exp.doRecordHistory=True

# Whether to record the time spent in this stage
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].exp.doRecordTime=True

# Softened core width for ellipticity distribution (conformal shear units).
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].exp.empiricalPriorConfig.ellipticityCore=0.001

# Width of exponential ellipticity distribution (conformal shear units).
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].exp.empiricalPriorConfig.ellipticitySigma=0.3

# ln(radius) at which the softened cutoff begins towards the minimum
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusMinInner=-6.0

# Minimum ln(radius).
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusMinOuter=-6.001

# Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T.
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusMu=-1.0

# Number of degrees of freedom for the Student's T distribution on ln(radius).
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusNu=50.0

# Width of the Student's T distribution in ln(radius).
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusSigma=0.45

# Ellipticity magnitude (conformal shear units) at which the softened cutoff begins
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].exp.linearPriorConfig.ellipticityMaxInner=2.0

# Maximum ellipticity magnitude (conformal shear units)
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].exp.linearPriorConfig.ellipticityMaxOuter=2.001

# ln(radius) at which the softened cutoff begins towards the maximum
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].exp.linearPriorConfig.logRadiusMaxInner=3.0

# Maximum ln(radius)
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].exp.linearPriorConfig.logRadiusMaxOuter=3.001

# ln(radius) at which the softened cutoff begins towards the minimum
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].exp.linearPriorConfig.logRadiusMinInner=-6.0

# The ratio P(logRadiusMinInner)/P(logRadiusMaxInner)
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].exp.linearPriorConfig.logRadiusMinMaxRatio=1.0

# Minimum ln(radius)
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].exp.linearPriorConfig.logRadiusMinOuter=-6.001

# Maximum radius used in approximating profile with Gaussians (0=default for this profile)
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].exp.maxRadius=0

# Number of Gaussian used to approximate the profile
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].exp.nComponents=6

# whether to save all iterations for debugging purposes
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].exp.optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].exp.optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].exp.optimizer.maxInnerIterations=20

# maximum number of steps
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].exp.optimizer.maxOuterIterations=250

# If the trust radius falls below this threshold, consider the algorithm converged
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].exp.optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].exp.optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].exp.optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].exp.optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].exp.optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].exp.optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].exp.optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].exp.optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].exp.optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].exp.optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].exp.optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].exp.optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].exp.optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].exp.optimizer.trustRegionSolverTolerance=1e-08

# Name of the Prior that defines the model to fit (a filename in $MEAS_MODELFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting.
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].exp.priorName=''

# One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].exp.priorSource='EMPIRICAL'

# Name of the shapelet.RadialProfile that defines the model to fit
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].exp.profileName='lux'

# Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances)
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].exp.usePixelWeights=False

# Scale the likelihood by this factor to artificially reweight it w.r.t. the prior.
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].exp.weightsMultiplier=1.0

# If the 2nd-moments shape used to initialize the fit failed, use the PSF moments multiplied by this.  If <= 0.0, abort the fit early instead.
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].fallbackInitialMomentsPsfFactor=1.5

# Whether to record the steps the optimizer takes (or just the number, if running as a plugin)
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].initial.doRecordHistory=True

# Whether to record the time spent in this stage
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].initial.doRecordTime=True

# Softened core width for ellipticity distribution (conformal shear units).
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].initial.empiricalPriorConfig.ellipticityCore=0.001

# Width of exponential ellipticity distribution (conformal shear units).
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].initial.empiricalPriorConfig.ellipticitySigma=0.3

# ln(radius) at which the softened cutoff begins towards the minimum
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusMinInner=-6.0

# Minimum ln(radius).
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusMinOuter=-6.001

# Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T.
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusMu=-1.0

# Number of degrees of freedom for the Student's T distribution on ln(radius).
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusNu=50.0

# Width of the Student's T distribution in ln(radius).
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusSigma=0.45

# Ellipticity magnitude (conformal shear units) at which the softened cutoff begins
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].initial.linearPriorConfig.ellipticityMaxInner=2.0

# Maximum ellipticity magnitude (conformal shear units)
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].initial.linearPriorConfig.ellipticityMaxOuter=2.001

# ln(radius) at which the softened cutoff begins towards the maximum
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].initial.linearPriorConfig.logRadiusMaxInner=3.0

# Maximum ln(radius)
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].initial.linearPriorConfig.logRadiusMaxOuter=3.001

# ln(radius) at which the softened cutoff begins towards the minimum
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].initial.linearPriorConfig.logRadiusMinInner=-6.0

# The ratio P(logRadiusMinInner)/P(logRadiusMaxInner)
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].initial.linearPriorConfig.logRadiusMinMaxRatio=1.0

# Minimum ln(radius)
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].initial.linearPriorConfig.logRadiusMinOuter=-6.001

# Maximum radius used in approximating profile with Gaussians (0=default for this profile)
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].initial.maxRadius=0

# Number of Gaussian used to approximate the profile
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].initial.nComponents=3

# whether to save all iterations for debugging purposes
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].initial.optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].initial.optimizer.gradientThreshold=0.001

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].initial.optimizer.maxInnerIterations=20

# maximum number of steps
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].initial.optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].initial.optimizer.minTrustRadiusThreshold=0.01

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].initial.optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].initial.optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].initial.optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].initial.optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].initial.optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].initial.optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].initial.optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].initial.optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].initial.optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].initial.optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].initial.optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].initial.optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].initial.optimizer.trustRegionSolverTolerance=1e-08

# Name of the Prior that defines the model to fit (a filename in $MEAS_MODELFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting.
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].initial.priorName=''

# One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].initial.priorSource='EMPIRICAL'

# Name of the shapelet.RadialProfile that defines the model to fit
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].initial.profileName='lux'

# Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances)
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].initial.usePixelWeights=True

# Scale the likelihood by this factor to artificially reweight it w.r.t. the prior.
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].initial.weightsMultiplier=1.0

# Minimum initial radius in pixels (used to regularize initial moments-based PSF deconvolution)
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].minInitialRadius=0.1

# Field name prefix of the Shapelet PSF approximation used to convolve the galaxy model; must contain a set of fields matching the schema defined by shapelet.MultiShapeletFunctionKey.
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].psfName='modelfit_DoubleShapeletPsfApprox'

# Mask planes that indicate pixels that should be ignored in the fit.
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].region.badMaskPlanes=['EDGE', 'SAT', 'BAD', 'NO_DATA']

# Abort if the fit region grows beyond this many pixels.
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].region.maxArea=100000

# Maximum fraction of pixels that may be ignored due to masks; more than this and we don't even try.
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].region.maxBadPixelFraction=0.1

# Use this multiple of the initial fit ellipse then grow by the PSF width to determine the maximum final fit region size.
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].region.nFitRadiiMax=3.0

# Use this multiple of the initial fit ellipse then grow by the PSF width to determine the minimum final fit region size.
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].region.nFitRadiiMin=1.0

# Use this multiple of the Kron ellipse to set the fit region (for the final fit region, subject to the nFitRadiiMin and nFitRadiiMax constraints).
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].region.nKronRadii=1.5

# Grow the initial fit ellipses by this factor before comparing with the Kron/Footprint region
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].region.nPsfSigmaGrow=2.0

# If the Kron radius is less than this multiple of the PSF width, ignore it and fall back to a PSF-oriented ellipse scaled to match the area of the footprint or this radius (whichever is larger).
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].region.nPsfSigmaMin=4.0

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.plugins['modelfit_CModel'].doMeasure=True

config.measureCoaddSources.measurement.plugins.names=['base_LocalBackground', 'modelfit_DoubleShapeletPsfApprox', 'base_SkyCoord', 'base_PsfFlux', 'base_GaussianFlux', 'base_CircularApertureFlux', 'base_SdssShape', 'base_Variance', 'ext_shapeHSM_HsmSourceMoments', 'ext_photometryKron_KronFlux', 'base_PixelFlags', 'ext_convolved_ConvolvedFlux', 'base_InputCount', 'ext_shapeHSM_HsmSourceMomentsRound', 'base_SdssCentroid', 'ext_shapeHSM_HsmPsfMoments', 'base_NaiveCentroid', 'ext_shapeHSM_HsmShapeRegauss', 'subaru_FilterFraction', 'base_Blendedness', 'modelfit_CModel']
# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.undeblended['base_PsfFlux'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.measureCoaddSources.measurement.undeblended['base_PsfFlux'].badMaskPlanes=[]

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.undeblended['base_PeakLikelihoodFlux'].doMeasure=True

# Name of warping kernel (e.g. "lanczos4") used to compute the peak
config.measureCoaddSources.measurement.undeblended['base_PeakLikelihoodFlux'].warpingKernelName='lanczos4'

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.undeblended['base_GaussianFlux'].doMeasure=True

# FIXME! NEVER DOCUMENTED!
config.measureCoaddSources.measurement.undeblended['base_GaussianFlux'].background=0.0

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.undeblended['base_NaiveCentroid'].doMeasure=True

# Value to subtract from the image pixel values
config.measureCoaddSources.measurement.undeblended['base_NaiveCentroid'].background=0.0

# Do check that the centroid is contained in footprint.
config.measureCoaddSources.measurement.undeblended['base_NaiveCentroid'].doFootprintCheck=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.measureCoaddSources.measurement.undeblended['base_NaiveCentroid'].maxDistToPeak=-1.0

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.undeblended['base_SdssCentroid'].doMeasure=True

# maximum allowed binning
config.measureCoaddSources.measurement.undeblended['base_SdssCentroid'].binmax=16

# Do check that the centroid is contained in footprint.
config.measureCoaddSources.measurement.undeblended['base_SdssCentroid'].doFootprintCheck=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.measureCoaddSources.measurement.undeblended['base_SdssCentroid'].maxDistToPeak=-1.0

# if the peak's less than this insist on binning at least once
config.measureCoaddSources.measurement.undeblended['base_SdssCentroid'].peakMin=-1.0

# fiddle factor for adjusting the binning
config.measureCoaddSources.measurement.undeblended['base_SdssCentroid'].wfac=1.5

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.undeblended['base_PixelFlags'].doMeasure=True

# List of mask planes to be searched for which occur anywhere within a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.measureCoaddSources.measurement.undeblended['base_PixelFlags'].masksFpAnywhere=[]

# List of mask planes to be searched for which occur in the center of a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.measureCoaddSources.measurement.undeblended['base_PixelFlags'].masksFpCenter=[]

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.undeblended['base_SdssShape'].doMeasure=True

# Additional value to add to background
config.measureCoaddSources.measurement.undeblended['base_SdssShape'].background=0.0

# Whether to also compute the shape of the PSF model
config.measureCoaddSources.measurement.undeblended['base_SdssShape'].doMeasurePsf=True

# Maximum number of iterations
config.measureCoaddSources.measurement.undeblended['base_SdssShape'].maxIter=100

# Maximum centroid shift, limited to 2-10
config.measureCoaddSources.measurement.undeblended['base_SdssShape'].maxShift=0.0

# Convergence tolerance for e1,e2
config.measureCoaddSources.measurement.undeblended['base_SdssShape'].tol1=9.999999747378752e-06

# Convergence tolerance for FWHM
config.measureCoaddSources.measurement.undeblended['base_SdssShape'].tol2=9.999999747378752e-05

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.undeblended['base_ScaledApertureFlux'].doMeasure=True

# Scaling factor of PSF FWHM for aperture radius.
config.measureCoaddSources.measurement.undeblended['base_ScaledApertureFlux'].scale=3.14

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.measureCoaddSources.measurement.undeblended['base_ScaledApertureFlux'].shiftKernel='lanczos5'

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.undeblended['base_CircularApertureFlux'].doMeasure=True

# Maximum radius (in pixels) for which the sinc algorithm should be used instead of the faster naive algorithm.  For elliptical apertures, this is the minor axis radius.
config.measureCoaddSources.measurement.undeblended['base_CircularApertureFlux'].maxSincRadius=10.0

# Radius (in pixels) of apertures.
config.measureCoaddSources.measurement.undeblended['base_CircularApertureFlux'].radii=[3.0, 4.5, 6.0, 9.0, 12.0, 17.0, 25.0, 35.0, 50.0, 70.0]

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.measureCoaddSources.measurement.undeblended['base_CircularApertureFlux'].shiftKernel='lanczos5'

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.undeblended['base_Blendedness'].doMeasure=True

# Whether to compute quantities related to the Gaussian-weighted flux
config.measureCoaddSources.measurement.undeblended['base_Blendedness'].doFlux=True

# Whether to compute HeavyFootprint dot products (the old deblend.blendedness parameter)
config.measureCoaddSources.measurement.undeblended['base_Blendedness'].doOld=True

# Whether to compute quantities related to the Gaussian-weighted shape
config.measureCoaddSources.measurement.undeblended['base_Blendedness'].doShape=True

# Radius factor that sets the maximum extent of the weight function (and hence the flux measurements)
config.measureCoaddSources.measurement.undeblended['base_Blendedness'].nSigmaWeightMax=3.0

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.undeblended['base_LocalBackground'].doMeasure=True

# Inner radius for background annulus as a multiple of the PSF sigma
config.measureCoaddSources.measurement.undeblended['base_LocalBackground'].annulusInner=7.0

# Outer radius for background annulus as a multiple of the PSF sigma
config.measureCoaddSources.measurement.undeblended['base_LocalBackground'].annulusOuter=15.0

# Mask planes that indicate pixels that should be excluded from the measurement
config.measureCoaddSources.measurement.undeblended['base_LocalBackground'].badMaskPlanes=['BAD', 'SAT', 'NO_DATA']

# Number of sigma-clipping iterations for background measurement
config.measureCoaddSources.measurement.undeblended['base_LocalBackground'].bgIter=3

# Rejection threshold (in standard deviations) for background measurement
config.measureCoaddSources.measurement.undeblended['base_LocalBackground'].bgRej=3.0

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.undeblended['base_FPPosition'].doMeasure=True

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.undeblended['base_Jacobian'].doMeasure=True

# Nominal pixel size (arcsec)
config.measureCoaddSources.measurement.undeblended['base_Jacobian'].pixelScale=0.5

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.undeblended['base_Variance'].doMeasure=True

# Scale factor to apply to shape for aperture
config.measureCoaddSources.measurement.undeblended['base_Variance'].scale=5.0

# Mask planes to ignore
config.measureCoaddSources.measurement.undeblended['base_Variance'].mask=['DETECTED', 'DETECTED_NEGATIVE', 'BAD', 'SAT']

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.undeblended['base_InputCount'].doMeasure=True

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.undeblended['base_PeakCentroid'].doMeasure=True

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.undeblended['base_SkyCoord'].doMeasure=True

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.undeblended['subaru_FilterFraction'].doMeasure=True

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.undeblended['ext_photometryKron_KronFlux'].doMeasure=True

# If true check that the Kron radius exceeds some minimum
config.measureCoaddSources.measurement.undeblended['ext_photometryKron_KronFlux'].enforceMinimumRadius=True

# if true, use existing shape and centroid measurements instead of fitting
config.measureCoaddSources.measurement.undeblended['ext_photometryKron_KronFlux'].fixed=False

# Largest aperture for which to use the slow, accurate, sinc aperture code
config.measureCoaddSources.measurement.undeblended['ext_photometryKron_KronFlux'].maxSincRadius=10.0

# Minimum Kron radius (if == 0.0 use PSF's Kron radius) if enforceMinimumRadius. Also functions as fallback aperture radius if set.
config.measureCoaddSources.measurement.undeblended['ext_photometryKron_KronFlux'].minimumRadius=0.0

# Number of times to iterate when setting the Kron radius
config.measureCoaddSources.measurement.undeblended['ext_photometryKron_KronFlux'].nIterForRadius=1

# Number of Kron radii for Kron flux
config.measureCoaddSources.measurement.undeblended['ext_photometryKron_KronFlux'].nRadiusForFlux=2.5

# Multiplier of rms size for aperture used to initially estimate the Kron radius
config.measureCoaddSources.measurement.undeblended['ext_photometryKron_KronFlux'].nSigmaForRadius=6.0

# Name of field specifying reference Kron radius for forced measurement
config.measureCoaddSources.measurement.undeblended['ext_photometryKron_KronFlux'].refRadiusName='ext_photometryKron_KronFlux_radius'

# Smooth image with N(0, smoothingSigma^2) Gaussian while estimating R_K
config.measureCoaddSources.measurement.undeblended['ext_photometryKron_KronFlux'].smoothingSigma=-1.0

# Use the Footprint size as part of initial estimate of Kron radius
config.measureCoaddSources.measurement.undeblended['ext_photometryKron_KronFlux'].useFootprintRadius=False

# list of target seeings (FWHM, pixels)
config.measureCoaddSources.measurement.undeblended['ext_convolved_ConvolvedFlux'].seeing=[3.5, 5.0, 6.5]

# scaling factor of kernel sigma for kernel size
config.measureCoaddSources.measurement.undeblended['ext_convolved_ConvolvedFlux'].kernelScale=4.0

# Maximum radius (in pixels) for which the sinc algorithm should be used instead of the faster naive algorithm.  For elliptical apertures, this is the minor axis radius.
config.measureCoaddSources.measurement.undeblended['ext_convolved_ConvolvedFlux'].aperture.maxSincRadius=10.0

# Radius (in pixels) of apertures.
config.measureCoaddSources.measurement.undeblended['ext_convolved_ConvolvedFlux'].aperture.radii=[3.3, 4.5, 6.0]

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.measureCoaddSources.measurement.undeblended['ext_convolved_ConvolvedFlux'].aperture.shiftKernel='lanczos5'

# name of Kron radius field in reference
config.measureCoaddSources.measurement.undeblended['ext_convolved_ConvolvedFlux'].kronRadiusName='ext_photometryKron_KronFlux_radius'

# Largest aperture for which to use the sinc aperture code for Kron (pixels)
config.measureCoaddSources.measurement.undeblended['ext_convolved_ConvolvedFlux'].maxSincRadius=10.0

# Number of Kron radii for Kron flux
config.measureCoaddSources.measurement.undeblended['ext_convolved_ConvolvedFlux'].kronRadiusForFlux=2.5

# Register measurements for aperture correction?
# The aperture correction registration is done when the plugin is
# instantiated because the column names are derived from the configuration
# rather than being static. Sometimes you want to turn this off, e.g.,
# when you will use aperture corrections derived from somewhere else
# through the 'proxy' mechanism.
config.measureCoaddSources.measurement.undeblended['ext_convolved_ConvolvedFlux'].registerForApCorr=True

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.undeblended['ext_convolved_ConvolvedFlux'].doMeasure=True

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.undeblended['ext_shapeHSM_HsmShapeBj'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.measureCoaddSources.measurement.undeblended['ext_shapeHSM_HsmShapeBj'].badMaskPlanes=None

# Field name for number of deblend children
config.measureCoaddSources.measurement.undeblended['ext_shapeHSM_HsmShapeBj'].deblendNChild=None

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.undeblended['ext_shapeHSM_HsmShapeLinear'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.measureCoaddSources.measurement.undeblended['ext_shapeHSM_HsmShapeLinear'].badMaskPlanes=None

# Field name for number of deblend children
config.measureCoaddSources.measurement.undeblended['ext_shapeHSM_HsmShapeLinear'].deblendNChild=None

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.undeblended['ext_shapeHSM_HsmShapeKsb'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.measureCoaddSources.measurement.undeblended['ext_shapeHSM_HsmShapeKsb'].badMaskPlanes=None

# Field name for number of deblend children
config.measureCoaddSources.measurement.undeblended['ext_shapeHSM_HsmShapeKsb'].deblendNChild=None

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.undeblended['ext_shapeHSM_HsmShapeRegauss'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.measureCoaddSources.measurement.undeblended['ext_shapeHSM_HsmShapeRegauss'].badMaskPlanes=None

# Field name for number of deblend children
config.measureCoaddSources.measurement.undeblended['ext_shapeHSM_HsmShapeRegauss'].deblendNChild=None

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.undeblended['ext_shapeHSM_HsmSourceMoments'].doMeasure=True

# Store measured flux?
config.measureCoaddSources.measurement.undeblended['ext_shapeHSM_HsmSourceMoments'].addFlux=None

# Mask planes used to reject bad pixels.
config.measureCoaddSources.measurement.undeblended['ext_shapeHSM_HsmSourceMoments'].badMaskPlanes=None

# Use round weight function?
config.measureCoaddSources.measurement.undeblended['ext_shapeHSM_HsmSourceMoments'].roundMoments=None

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.undeblended['ext_shapeHSM_HsmSourceMomentsRound'].doMeasure=True

# Store measured flux?
config.measureCoaddSources.measurement.undeblended['ext_shapeHSM_HsmSourceMomentsRound'].addFlux=None

# Mask planes used to reject bad pixels.
config.measureCoaddSources.measurement.undeblended['ext_shapeHSM_HsmSourceMomentsRound'].badMaskPlanes=None

# Use round weight function?
config.measureCoaddSources.measurement.undeblended['ext_shapeHSM_HsmSourceMomentsRound'].roundMoments=None

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.undeblended['ext_shapeHSM_HsmPsfMoments'].doMeasure=True

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].doMeasure=True

# Shapelet order of inner expansion (0 == Gaussian)
config.measureCoaddSources.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].innerOrder=2

# Don't allow the semi-major radius of any component to go above this fraction of the PSF image width
config.measureCoaddSources.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].maxRadiusBoxFraction=0.4

# Don't allow the semi-minor radius of any component to drop below this value (pixels)
config.measureCoaddSources.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].minRadius=1.0

# Don't allow the determinant radii of the two components to differ by less than this (pixels)
config.measureCoaddSources.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].minRadiusDiff=0.5

# whether to save all iterations for debugging purposes
config.measureCoaddSources.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.measureCoaddSources.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.measureCoaddSources.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.maxInnerIterations=20

# maximum number of steps
config.measureCoaddSources.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.measureCoaddSources.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.measureCoaddSources.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.measureCoaddSources.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.measureCoaddSources.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.measureCoaddSources.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.measureCoaddSources.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.measureCoaddSources.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.measureCoaddSources.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.measureCoaddSources.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.measureCoaddSources.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.measureCoaddSources.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.measureCoaddSources.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.measureCoaddSources.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.measureCoaddSources.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionSolverTolerance=1e-08

# Shapelet order of outer expansion (0 == Gaussian)
config.measureCoaddSources.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].outerOrder=1

# Initial outer Gaussian peak height divided by inner Gaussian peak height
config.measureCoaddSources.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].peakRatio=0.1

# Initial outer radius divided by inner radius
config.measureCoaddSources.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].radiusRatio=2.0

config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models={}
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.maxInnerIterations=20

# maximum number of steps
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.radiusPriorSigma=0.5

config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.maxInnerIterations=20

# maximum number of steps
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.radiusPriorSigma=0.5

config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.maxInnerIterations=20

# maximum number of steps
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.order=2

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.order=1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.radiusPriorSigma=0.5

config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.maxInnerIterations=20

# maximum number of steps
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.order=4

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.order=4

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.radiusPriorSigma=0.5

# a sequence of model names indicating which models should be fit, and their order
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].sequence=['DoubleShapelet']

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].doMeasure=True

# Whether to record the steps the optimizer takes (or just the number, if running as a plugin)
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].dev.doRecordHistory=True

# Whether to record the time spent in this stage
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].dev.doRecordTime=True

# Softened core width for ellipticity distribution (conformal shear units).
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].dev.empiricalPriorConfig.ellipticityCore=0.001

# Width of exponential ellipticity distribution (conformal shear units).
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].dev.empiricalPriorConfig.ellipticitySigma=0.3

# ln(radius) at which the softened cutoff begins towards the minimum
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusMinInner=-6.0

# Minimum ln(radius).
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusMinOuter=-6.001

# Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T.
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusMu=-1.0

# Number of degrees of freedom for the Student's T distribution on ln(radius).
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusNu=50.0

# Width of the Student's T distribution in ln(radius).
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusSigma=0.45

# Ellipticity magnitude (conformal shear units) at which the softened cutoff begins
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].dev.linearPriorConfig.ellipticityMaxInner=2.0

# Maximum ellipticity magnitude (conformal shear units)
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].dev.linearPriorConfig.ellipticityMaxOuter=2.001

# ln(radius) at which the softened cutoff begins towards the maximum
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].dev.linearPriorConfig.logRadiusMaxInner=3.0

# Maximum ln(radius)
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].dev.linearPriorConfig.logRadiusMaxOuter=3.001

# ln(radius) at which the softened cutoff begins towards the minimum
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].dev.linearPriorConfig.logRadiusMinInner=-6.0

# The ratio P(logRadiusMinInner)/P(logRadiusMaxInner)
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].dev.linearPriorConfig.logRadiusMinMaxRatio=1.0

# Minimum ln(radius)
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].dev.linearPriorConfig.logRadiusMinOuter=-6.001

# Maximum radius used in approximating profile with Gaussians (0=default for this profile)
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].dev.maxRadius=0

# Number of Gaussian used to approximate the profile
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].dev.nComponents=8

# whether to save all iterations for debugging purposes
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].dev.optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].dev.optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].dev.optimizer.maxInnerIterations=20

# maximum number of steps
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].dev.optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].dev.optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].dev.optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].dev.optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].dev.optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].dev.optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].dev.optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].dev.optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].dev.optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].dev.optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].dev.optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].dev.optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].dev.optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].dev.optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].dev.optimizer.trustRegionSolverTolerance=1e-08

# Name of the Prior that defines the model to fit (a filename in $MEAS_MODELFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting.
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].dev.priorName=''

# One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].dev.priorSource='EMPIRICAL'

# Name of the shapelet.RadialProfile that defines the model to fit
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].dev.profileName='luv'

# Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances)
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].dev.usePixelWeights=False

# Scale the likelihood by this factor to artificially reweight it w.r.t. the prior.
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].dev.weightsMultiplier=1.0

# Whether to record the steps the optimizer takes (or just the number, if running as a plugin)
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].exp.doRecordHistory=True

# Whether to record the time spent in this stage
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].exp.doRecordTime=True

# Softened core width for ellipticity distribution (conformal shear units).
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].exp.empiricalPriorConfig.ellipticityCore=0.001

# Width of exponential ellipticity distribution (conformal shear units).
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].exp.empiricalPriorConfig.ellipticitySigma=0.3

# ln(radius) at which the softened cutoff begins towards the minimum
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusMinInner=-6.0

# Minimum ln(radius).
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusMinOuter=-6.001

# Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T.
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusMu=-1.0

# Number of degrees of freedom for the Student's T distribution on ln(radius).
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusNu=50.0

# Width of the Student's T distribution in ln(radius).
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusSigma=0.45

# Ellipticity magnitude (conformal shear units) at which the softened cutoff begins
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].exp.linearPriorConfig.ellipticityMaxInner=2.0

# Maximum ellipticity magnitude (conformal shear units)
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].exp.linearPriorConfig.ellipticityMaxOuter=2.001

# ln(radius) at which the softened cutoff begins towards the maximum
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].exp.linearPriorConfig.logRadiusMaxInner=3.0

# Maximum ln(radius)
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].exp.linearPriorConfig.logRadiusMaxOuter=3.001

# ln(radius) at which the softened cutoff begins towards the minimum
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].exp.linearPriorConfig.logRadiusMinInner=-6.0

# The ratio P(logRadiusMinInner)/P(logRadiusMaxInner)
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].exp.linearPriorConfig.logRadiusMinMaxRatio=1.0

# Minimum ln(radius)
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].exp.linearPriorConfig.logRadiusMinOuter=-6.001

# Maximum radius used in approximating profile with Gaussians (0=default for this profile)
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].exp.maxRadius=0

# Number of Gaussian used to approximate the profile
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].exp.nComponents=6

# whether to save all iterations for debugging purposes
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].exp.optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].exp.optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].exp.optimizer.maxInnerIterations=20

# maximum number of steps
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].exp.optimizer.maxOuterIterations=250

# If the trust radius falls below this threshold, consider the algorithm converged
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].exp.optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].exp.optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].exp.optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].exp.optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].exp.optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].exp.optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].exp.optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].exp.optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].exp.optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].exp.optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].exp.optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].exp.optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].exp.optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].exp.optimizer.trustRegionSolverTolerance=1e-08

# Name of the Prior that defines the model to fit (a filename in $MEAS_MODELFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting.
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].exp.priorName=''

# One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].exp.priorSource='EMPIRICAL'

# Name of the shapelet.RadialProfile that defines the model to fit
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].exp.profileName='lux'

# Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances)
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].exp.usePixelWeights=False

# Scale the likelihood by this factor to artificially reweight it w.r.t. the prior.
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].exp.weightsMultiplier=1.0

# If the 2nd-moments shape used to initialize the fit failed, use the PSF moments multiplied by this.  If <= 0.0, abort the fit early instead.
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].fallbackInitialMomentsPsfFactor=1.5

# Whether to record the steps the optimizer takes (or just the number, if running as a plugin)
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].initial.doRecordHistory=True

# Whether to record the time spent in this stage
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].initial.doRecordTime=True

# Softened core width for ellipticity distribution (conformal shear units).
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].initial.empiricalPriorConfig.ellipticityCore=0.001

# Width of exponential ellipticity distribution (conformal shear units).
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].initial.empiricalPriorConfig.ellipticitySigma=0.3

# ln(radius) at which the softened cutoff begins towards the minimum
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusMinInner=-6.0

# Minimum ln(radius).
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusMinOuter=-6.001

# Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T.
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusMu=-1.0

# Number of degrees of freedom for the Student's T distribution on ln(radius).
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusNu=50.0

# Width of the Student's T distribution in ln(radius).
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusSigma=0.45

# Ellipticity magnitude (conformal shear units) at which the softened cutoff begins
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].initial.linearPriorConfig.ellipticityMaxInner=2.0

# Maximum ellipticity magnitude (conformal shear units)
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].initial.linearPriorConfig.ellipticityMaxOuter=2.001

# ln(radius) at which the softened cutoff begins towards the maximum
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].initial.linearPriorConfig.logRadiusMaxInner=3.0

# Maximum ln(radius)
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].initial.linearPriorConfig.logRadiusMaxOuter=3.001

# ln(radius) at which the softened cutoff begins towards the minimum
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].initial.linearPriorConfig.logRadiusMinInner=-6.0

# The ratio P(logRadiusMinInner)/P(logRadiusMaxInner)
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].initial.linearPriorConfig.logRadiusMinMaxRatio=1.0

# Minimum ln(radius)
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].initial.linearPriorConfig.logRadiusMinOuter=-6.001

# Maximum radius used in approximating profile with Gaussians (0=default for this profile)
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].initial.maxRadius=0

# Number of Gaussian used to approximate the profile
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].initial.nComponents=3

# whether to save all iterations for debugging purposes
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].initial.optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].initial.optimizer.gradientThreshold=0.001

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].initial.optimizer.maxInnerIterations=20

# maximum number of steps
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].initial.optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].initial.optimizer.minTrustRadiusThreshold=0.01

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].initial.optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].initial.optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].initial.optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].initial.optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].initial.optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].initial.optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].initial.optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].initial.optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].initial.optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].initial.optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].initial.optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].initial.optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].initial.optimizer.trustRegionSolverTolerance=1e-08

# Name of the Prior that defines the model to fit (a filename in $MEAS_MODELFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting.
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].initial.priorName=''

# One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].initial.priorSource='EMPIRICAL'

# Name of the shapelet.RadialProfile that defines the model to fit
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].initial.profileName='lux'

# Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances)
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].initial.usePixelWeights=True

# Scale the likelihood by this factor to artificially reweight it w.r.t. the prior.
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].initial.weightsMultiplier=1.0

# Minimum initial radius in pixels (used to regularize initial moments-based PSF deconvolution)
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].minInitialRadius=0.1

# Field name prefix of the Shapelet PSF approximation used to convolve the galaxy model; must contain a set of fields matching the schema defined by shapelet.MultiShapeletFunctionKey.
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].psfName='modelfit_DoubleShapeletPsfApprox'

# Mask planes that indicate pixels that should be ignored in the fit.
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].region.badMaskPlanes=['EDGE', 'SAT', 'BAD', 'NO_DATA']

# Abort if the fit region grows beyond this many pixels.
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].region.maxArea=100000

# Maximum fraction of pixels that may be ignored due to masks; more than this and we don't even try.
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].region.maxBadPixelFraction=0.1

# Use this multiple of the initial fit ellipse then grow by the PSF width to determine the maximum final fit region size.
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].region.nFitRadiiMax=3.0

# Use this multiple of the initial fit ellipse then grow by the PSF width to determine the minimum final fit region size.
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].region.nFitRadiiMin=1.0

# Use this multiple of the Kron ellipse to set the fit region (for the final fit region, subject to the nFitRadiiMin and nFitRadiiMax constraints).
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].region.nKronRadii=1.5

# Grow the initial fit ellipses by this factor before comparing with the Kron/Footprint region
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].region.nPsfSigmaGrow=2.0

# If the Kron radius is less than this multiple of the PSF width, ignore it and fall back to a PSF-oriented ellipse scaled to match the area of the footprint or this radius (whichever is larger).
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].region.nPsfSigmaMin=4.0

# whether to run this plugin in single-object mode
config.measureCoaddSources.measurement.undeblended['modelfit_CModel'].doMeasure=True

config.measureCoaddSources.measurement.undeblended.names=[]
# Name of field in schema with number of deblended children
config.measureCoaddSources.setPrimaryFlags.nChildKeyName='deblend_nChild'

# Names of filters which should never be primary
config.measureCoaddSources.setPrimaryFlags.pseudoFilterList=['sky']

# Whether to match sources to CCD catalogs to propagate flags (to e.g. identify PSF stars)
config.measureCoaddSources.doPropagateFlags=True

# Source catalog flags to propagate, with the threshold of relative occurrence.
config.measureCoaddSources.propagateFlags.flags={'calib_psf_candidate': 0.2, 'calib_psf_used': 0.2}

# Source matching radius (arcsec)
config.measureCoaddSources.propagateFlags.matchRadius=0.2

# Match sources to reference catalog?
config.measureCoaddSources.doMatchSources=True

# Matching radius, arcsec
config.measureCoaddSources.match.matchRadius=0.25

# Apply flux limit?
config.measureCoaddSources.match.sourceSelection.doFluxLimit=False

# Apply flag limitation?
config.measureCoaddSources.match.sourceSelection.doFlags=False

# Apply unresolved limitation?
config.measureCoaddSources.match.sourceSelection.doUnresolved=False

# Apply signal-to-noise limit?
config.measureCoaddSources.match.sourceSelection.doSignalToNoise=False

# Apply isolated limitation?
config.measureCoaddSources.match.sourceSelection.doIsolated=False

# Select objects with value greater than this
config.measureCoaddSources.match.sourceSelection.fluxLimit.minimum=None

# Select objects with value less than this
config.measureCoaddSources.match.sourceSelection.fluxLimit.maximum=None

# Name of the source flux field to use.
config.measureCoaddSources.match.sourceSelection.fluxLimit.fluxField='slot_CalibFlux_flux'

# List of source flag fields that must be set for a source to be used.
config.measureCoaddSources.match.sourceSelection.flags.good=[]

# List of source flag fields that must NOT be set for a source to be used.
config.measureCoaddSources.match.sourceSelection.flags.bad=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_saturated', 'base_PsfFlux_flags']

# Select objects with value greater than this
config.measureCoaddSources.match.sourceSelection.unresolved.minimum=None

# Select objects with value less than this
config.measureCoaddSources.match.sourceSelection.unresolved.maximum=0.5

# Name of column for star/galaxy separation
config.measureCoaddSources.match.sourceSelection.unresolved.name='base_ClassificationExtendedness_value'

# Select objects with value greater than this
config.measureCoaddSources.match.sourceSelection.signalToNoise.minimum=None

# Select objects with value less than this
config.measureCoaddSources.match.sourceSelection.signalToNoise.maximum=None

# Name of the source flux field to use.
config.measureCoaddSources.match.sourceSelection.signalToNoise.fluxField='base_PsfFlux_flux'

# Name of the source flux error field to use.
config.measureCoaddSources.match.sourceSelection.signalToNoise.errField='base_PsfFlux_fluxErr'

# Name of column for parent
config.measureCoaddSources.match.sourceSelection.isolated.parentName='parent'

# Name of column for nChild
config.measureCoaddSources.match.sourceSelection.isolated.nChildName='deblend_nChild'

# Apply magnitude limit?
config.measureCoaddSources.match.referenceSelection.doMagLimit=False

# Apply flag limitation?
config.measureCoaddSources.match.referenceSelection.doFlags=False

# Apply signal-to-noise limit?
config.measureCoaddSources.match.referenceSelection.doSignalToNoise=False

# Apply magnitude error limit?
config.measureCoaddSources.match.referenceSelection.doMagError=False

# Select objects with value greater than this
config.measureCoaddSources.match.referenceSelection.magLimit.minimum=None

# Select objects with value less than this
config.measureCoaddSources.match.referenceSelection.magLimit.maximum=None

# Name of the source flux field to use.
config.measureCoaddSources.match.referenceSelection.magLimit.fluxField='flux'

# List of source flag fields that must be set for a source to be used.
config.measureCoaddSources.match.referenceSelection.flags.good=[]

# List of source flag fields that must NOT be set for a source to be used.
config.measureCoaddSources.match.referenceSelection.flags.bad=[]

# Select objects with value greater than this
config.measureCoaddSources.match.referenceSelection.signalToNoise.minimum=None

# Select objects with value less than this
config.measureCoaddSources.match.referenceSelection.signalToNoise.maximum=None

# Name of the source flux field to use.
config.measureCoaddSources.match.referenceSelection.signalToNoise.fluxField='flux'

# Name of the source flux error field to use.
config.measureCoaddSources.match.referenceSelection.signalToNoise.errField='flux_err'

# Select objects with value greater than this
config.measureCoaddSources.match.referenceSelection.magError.minimum=None

# Select objects with value less than this
config.measureCoaddSources.match.referenceSelection.magError.maximum=None

# Name of the source flux error field to use.
config.measureCoaddSources.match.referenceSelection.magError.magErrField='mag_err'

config.measureCoaddSources.match.referenceSelection.colorLimits={}
# Padding to add to 4 all edges of the bounding box (pixels)
config.measureCoaddSources.match.refObjLoader.pixelMargin=300

# Default reference catalog filter to use if filter not specified in exposure; if blank then filter must be specified in exposure
config.measureCoaddSources.match.refObjLoader.defaultFilter=''

# Mapping of camera filter name: reference catalog filter name; each reference filter must exist
config.measureCoaddSources.match.refObjLoader.filterMap={'B': 'g', 'V': 'r', 'R': 'r', 'I': 'i', 'r2': 'r', 'i2': 'i', 'N387': 'g', 'N468': 'g', 'N515': 'g', 'N527': 'g', 'N656': 'r', 'N718': 'i', 'N816': 'i', 'N921': 'z', 'N926': 'z', 'N973': 'y', 'N1010': 'y', 'I945': 'z', 'y': 'z'}

# Name of the ingested reference dataset
config.measureCoaddSources.match.refObjLoader.ref_dataset_name='ps1_pv3_3pi_20170110'

# Write reference matches in denormalized format? This format uses more disk space, but is more convenient to read.
config.measureCoaddSources.doWriteMatchesDenormalized=True

# Name of coadd
config.measureCoaddSources.coaddName='deep'

# Strictness of Astropy unit compatibility check, can be 'raise', 'warn' or 'silent'
config.measureCoaddSources.checkUnitsParseStrict='raise'

# Apply aperture corrections
config.measureCoaddSources.doApCorr=True

# flux measurement algorithms in getApCorrNameSet() to ignore; if a name is listed that does not appear in getApCorrNameSet() then a warning is logged
config.measureCoaddSources.applyApCorr.ignoreList=[]

# set the general failure flag for a flux when it cannot be aperture-corrected?
config.measureCoaddSources.applyApCorr.doFlagApCorrFailures=True

# flux measurement algorithms to be aperture-corrected by reference to another algorithm; this is a mapping alg1:alg2, where 'alg1' is the algorithm being corrected, and 'alg2' is the algorithm supplying the corrections
config.measureCoaddSources.applyApCorr.proxies={}

# Run catalogCalculation task
config.measureCoaddSources.doRunCatalogCalculation=True

# critical ratio of model to psf flux
config.measureCoaddSources.catalogCalculation.plugins['base_ClassificationExtendedness'].fluxRatio=0.985

# correction factor for modelFlux error
config.measureCoaddSources.catalogCalculation.plugins['base_ClassificationExtendedness'].modelErrFactor=0.0

# correction factor for psfFlux error
config.measureCoaddSources.catalogCalculation.plugins['base_ClassificationExtendedness'].psfErrFactor=0.0

config.measureCoaddSources.catalogCalculation.plugins.names=['base_ClassificationExtendedness', 'base_FootprintArea']
# Priority-ordered list of bands for the merge.
config.mergeCoaddMeasurements.priorityList=['HSC-I2', 'HSC-I', 'HSC-R2', 'HSC-R', 'HSC-Z', 'HSC-Y', 'HSC-G', 'NB0921', 'NB0816', 'NB1010', 'NB0387', 'NB0515']

# Name of coadd
config.mergeCoaddMeasurements.coaddName='deep'

# Names of filters which may have no associated detection
# (N.b. should include MergeDetectionsConfig.skyFilterName)
config.mergeCoaddMeasurements.pseudoFilterList=['sky']

# Name of flux measurement for calculating the S/N when choosing the reference band.
config.mergeCoaddMeasurements.snName='base_PsfFlux'

# If the S/N from the priority band is below this value (and the S/N is larger than minSNDiff compared to the priority band), use the band with the largest S/N as the reference band.
config.mergeCoaddMeasurements.minSN=10.0

# If the difference in S/N between another band and the priority band is larger than this value (and the S/N in the priority band is less than minSN) use the band with the largest S/N as the reference band
config.mergeCoaddMeasurements.minSNDiff=3.0

# Require that these flags, if available, are not set
config.mergeCoaddMeasurements.flags=['base_PixelFlags_flag_interpolatedCenter', 'base_PsfFlux_flag', 'ext_photometryKron_KronFlux_flag', 'modelfit_CModel_flag']

# Only include reference sources for each patch that lie within the patch's inner bbox
config.forcedPhotCoadd.references.removePatchOverlaps=False

# Bandpass for reference sources; None indicates chi-squared detections.
config.forcedPhotCoadd.references.filter=None

# Coadd name: typically one of deep or goodSeeing.
config.forcedPhotCoadd.references.coaddName='deep'

# Silently skip patches where the reference catalog does not exist.
config.forcedPhotCoadd.references.skipMissing=False

# the name of the centroiding algorithm used to set source x,y
config.forcedPhotCoadd.measurement.slots.centroid='base_TransformedCentroid'

# the name of the algorithm used to set source moments parameters
config.forcedPhotCoadd.measurement.slots.shape='base_TransformedShape'

# the name of the algorithm used to set PSF moments parameters
config.forcedPhotCoadd.measurement.slots.psfShape='base_SdssShape_psf'

# the name of the algorithm used to set the source aperture flux slot
config.forcedPhotCoadd.measurement.slots.apFlux=None

# the name of the algorithm used to set the source model flux slot
config.forcedPhotCoadd.measurement.slots.modelFlux='modelfit_CModel'

# the name of the algorithm used to set the source psf flux slot
config.forcedPhotCoadd.measurement.slots.psfFlux='base_PsfFlux'

# the name of the algorithm used to set the source inst flux slot
config.forcedPhotCoadd.measurement.slots.instFlux=None

# the name of the flux measurement algorithm used for calibration
config.forcedPhotCoadd.measurement.slots.calibFlux=None

# When measuring, replace other detected footprints with noise?
config.forcedPhotCoadd.measurement.doReplaceWithNoise=True

# How to choose mean and variance of the Gaussian noise we generate?
config.forcedPhotCoadd.measurement.noiseReplacer.noiseSource='measure'

# Add ann offset to the generated noise.
config.forcedPhotCoadd.measurement.noiseReplacer.noiseOffset=0.0

# The seed multiplier value to use for random number generation
#    >= 1: set the seed deterministically based on exposureId
#       0: fall back to the afw.math.Random default constructor (which uses a seed value of 1)
config.forcedPhotCoadd.measurement.noiseReplacer.noiseSeedMultiplier=1

# Prefix to give undeblended plugins
config.forcedPhotCoadd.measurement.undeblendedPrefix='undeblended_'

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.plugins['base_PsfFlux'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.forcedPhotCoadd.measurement.plugins['base_PsfFlux'].badMaskPlanes=[]

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.plugins['base_PeakLikelihoodFlux'].doMeasure=True

# Name of warping kernel (e.g. "lanczos4") used to compute the peak
config.forcedPhotCoadd.measurement.plugins['base_PeakLikelihoodFlux'].warpingKernelName='lanczos4'

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.plugins['base_GaussianFlux'].doMeasure=True

# FIXME! NEVER DOCUMENTED!
config.forcedPhotCoadd.measurement.plugins['base_GaussianFlux'].background=0.0

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.plugins['base_NaiveCentroid'].doMeasure=True

# Value to subtract from the image pixel values
config.forcedPhotCoadd.measurement.plugins['base_NaiveCentroid'].background=0.0

# Do check that the centroid is contained in footprint.
config.forcedPhotCoadd.measurement.plugins['base_NaiveCentroid'].doFootprintCheck=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.forcedPhotCoadd.measurement.plugins['base_NaiveCentroid'].maxDistToPeak=-1.0

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.plugins['base_SdssCentroid'].doMeasure=True

# maximum allowed binning
config.forcedPhotCoadd.measurement.plugins['base_SdssCentroid'].binmax=16

# Do check that the centroid is contained in footprint.
config.forcedPhotCoadd.measurement.plugins['base_SdssCentroid'].doFootprintCheck=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.forcedPhotCoadd.measurement.plugins['base_SdssCentroid'].maxDistToPeak=-1.0

# if the peak's less than this insist on binning at least once
config.forcedPhotCoadd.measurement.plugins['base_SdssCentroid'].peakMin=-1.0

# fiddle factor for adjusting the binning
config.forcedPhotCoadd.measurement.plugins['base_SdssCentroid'].wfac=1.5

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.plugins['base_PixelFlags'].doMeasure=True

# List of mask planes to be searched for which occur anywhere within a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.forcedPhotCoadd.measurement.plugins['base_PixelFlags'].masksFpAnywhere=['CLIPPED', 'SENSOR_EDGE', 'REJECTED', 'INEXACT_PSF', 'BRIGHT_OBJECT']

# List of mask planes to be searched for which occur in the center of a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.forcedPhotCoadd.measurement.plugins['base_PixelFlags'].masksFpCenter=['CLIPPED', 'SENSOR_EDGE', 'REJECTED', 'INEXACT_PSF', 'BRIGHT_OBJECT']

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.plugins['base_SdssShape'].doMeasure=True

# Additional value to add to background
config.forcedPhotCoadd.measurement.plugins['base_SdssShape'].background=0.0

# Whether to also compute the shape of the PSF model
config.forcedPhotCoadd.measurement.plugins['base_SdssShape'].doMeasurePsf=True

# Maximum number of iterations
config.forcedPhotCoadd.measurement.plugins['base_SdssShape'].maxIter=100

# Maximum centroid shift, limited to 2-10
config.forcedPhotCoadd.measurement.plugins['base_SdssShape'].maxShift=0.0

# Convergence tolerance for e1,e2
config.forcedPhotCoadd.measurement.plugins['base_SdssShape'].tol1=9.999999747378752e-06

# Convergence tolerance for FWHM
config.forcedPhotCoadd.measurement.plugins['base_SdssShape'].tol2=9.999999747378752e-05

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.plugins['base_ScaledApertureFlux'].doMeasure=True

# Scaling factor of PSF FWHM for aperture radius.
config.forcedPhotCoadd.measurement.plugins['base_ScaledApertureFlux'].scale=3.14

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.forcedPhotCoadd.measurement.plugins['base_ScaledApertureFlux'].shiftKernel='lanczos5'

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.plugins['base_CircularApertureFlux'].doMeasure=True

# Maximum radius (in pixels) for which the sinc algorithm should be used instead of the faster naive algorithm.  For elliptical apertures, this is the minor axis radius.
config.forcedPhotCoadd.measurement.plugins['base_CircularApertureFlux'].maxSincRadius=12.0

# Radius (in pixels) of apertures.
config.forcedPhotCoadd.measurement.plugins['base_CircularApertureFlux'].radii=[3.0, 4.5, 6.0, 9.0, 12.0, 17.0, 25.0, 35.0, 50.0, 70.0]

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.forcedPhotCoadd.measurement.plugins['base_CircularApertureFlux'].shiftKernel='lanczos5'

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.plugins['base_Blendedness'].doMeasure=True

# Whether to compute quantities related to the Gaussian-weighted flux
config.forcedPhotCoadd.measurement.plugins['base_Blendedness'].doFlux=True

# Whether to compute HeavyFootprint dot products (the old deblend.blendedness parameter)
config.forcedPhotCoadd.measurement.plugins['base_Blendedness'].doOld=True

# Whether to compute quantities related to the Gaussian-weighted shape
config.forcedPhotCoadd.measurement.plugins['base_Blendedness'].doShape=True

# Radius factor that sets the maximum extent of the weight function (and hence the flux measurements)
config.forcedPhotCoadd.measurement.plugins['base_Blendedness'].nSigmaWeightMax=3.0

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.plugins['base_LocalBackground'].doMeasure=True

# Inner radius for background annulus as a multiple of the PSF sigma
config.forcedPhotCoadd.measurement.plugins['base_LocalBackground'].annulusInner=7.0

# Outer radius for background annulus as a multiple of the PSF sigma
config.forcedPhotCoadd.measurement.plugins['base_LocalBackground'].annulusOuter=15.0

# Mask planes that indicate pixels that should be excluded from the measurement
config.forcedPhotCoadd.measurement.plugins['base_LocalBackground'].badMaskPlanes=['BAD', 'SAT', 'NO_DATA']

# Number of sigma-clipping iterations for background measurement
config.forcedPhotCoadd.measurement.plugins['base_LocalBackground'].bgIter=3

# Rejection threshold (in standard deviations) for background measurement
config.forcedPhotCoadd.measurement.plugins['base_LocalBackground'].bgRej=3.0

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.plugins['base_Variance'].doMeasure=True

# Scale factor to apply to shape for aperture
config.forcedPhotCoadd.measurement.plugins['base_Variance'].scale=5.0

# Mask planes to ignore
config.forcedPhotCoadd.measurement.plugins['base_Variance'].mask=['DETECTED', 'DETECTED_NEGATIVE', 'BAD', 'SAT']

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.plugins['base_InputCount'].doMeasure=True

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.plugins['base_PeakCentroid'].doMeasure=True

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.plugins['base_TransformedCentroid'].doMeasure=True

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.plugins['base_TransformedShape'].doMeasure=True

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.plugins['ext_photometryKron_KronFlux'].doMeasure=True

# If true check that the Kron radius exceeds some minimum
config.forcedPhotCoadd.measurement.plugins['ext_photometryKron_KronFlux'].enforceMinimumRadius=True

# if true, use existing shape and centroid measurements instead of fitting
config.forcedPhotCoadd.measurement.plugins['ext_photometryKron_KronFlux'].fixed=False

# Largest aperture for which to use the slow, accurate, sinc aperture code
config.forcedPhotCoadd.measurement.plugins['ext_photometryKron_KronFlux'].maxSincRadius=10.0

# Minimum Kron radius (if == 0.0 use PSF's Kron radius) if enforceMinimumRadius. Also functions as fallback aperture radius if set.
config.forcedPhotCoadd.measurement.plugins['ext_photometryKron_KronFlux'].minimumRadius=0.0

# Number of times to iterate when setting the Kron radius
config.forcedPhotCoadd.measurement.plugins['ext_photometryKron_KronFlux'].nIterForRadius=1

# Number of Kron radii for Kron flux
config.forcedPhotCoadd.measurement.plugins['ext_photometryKron_KronFlux'].nRadiusForFlux=2.5

# Multiplier of rms size for aperture used to initially estimate the Kron radius
config.forcedPhotCoadd.measurement.plugins['ext_photometryKron_KronFlux'].nSigmaForRadius=6.0

# Name of field specifying reference Kron radius for forced measurement
config.forcedPhotCoadd.measurement.plugins['ext_photometryKron_KronFlux'].refRadiusName='ext_photometryKron_KronFlux_radius'

# Smooth image with N(0, smoothingSigma^2) Gaussian while estimating R_K
config.forcedPhotCoadd.measurement.plugins['ext_photometryKron_KronFlux'].smoothingSigma=-1.0

# Use the Footprint size as part of initial estimate of Kron radius
config.forcedPhotCoadd.measurement.plugins['ext_photometryKron_KronFlux'].useFootprintRadius=False

# list of target seeings (FWHM, pixels)
config.forcedPhotCoadd.measurement.plugins['ext_convolved_ConvolvedFlux'].seeing=[3.5, 5.0, 6.5, 8.0]

# scaling factor of kernel sigma for kernel size
config.forcedPhotCoadd.measurement.plugins['ext_convolved_ConvolvedFlux'].kernelScale=4.0

# Maximum radius (in pixels) for which the sinc algorithm should be used instead of the faster naive algorithm.  For elliptical apertures, this is the minor axis radius.
config.forcedPhotCoadd.measurement.plugins['ext_convolved_ConvolvedFlux'].aperture.maxSincRadius=10.0

# Radius (in pixels) of apertures.
config.forcedPhotCoadd.measurement.plugins['ext_convolved_ConvolvedFlux'].aperture.radii=[3.3, 4.5, 6.0]

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.forcedPhotCoadd.measurement.plugins['ext_convolved_ConvolvedFlux'].aperture.shiftKernel='lanczos5'

# name of Kron radius field in reference
config.forcedPhotCoadd.measurement.plugins['ext_convolved_ConvolvedFlux'].kronRadiusName='ext_photometryKron_KronFlux_radius'

# Largest aperture for which to use the sinc aperture code for Kron (pixels)
config.forcedPhotCoadd.measurement.plugins['ext_convolved_ConvolvedFlux'].maxSincRadius=10.0

# Number of Kron radii for Kron flux
config.forcedPhotCoadd.measurement.plugins['ext_convolved_ConvolvedFlux'].kronRadiusForFlux=2.5

# Register measurements for aperture correction?
# The aperture correction registration is done when the plugin is
# instantiated because the column names are derived from the configuration
# rather than being static. Sometimes you want to turn this off, e.g.,
# when you will use aperture corrections derived from somewhere else
# through the 'proxy' mechanism.
config.forcedPhotCoadd.measurement.plugins['ext_convolved_ConvolvedFlux'].registerForApCorr=True

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.plugins['ext_convolved_ConvolvedFlux'].doMeasure=True

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.plugins['ext_shapeHSM_HsmShapeBj'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.forcedPhotCoadd.measurement.plugins['ext_shapeHSM_HsmShapeBj'].badMaskPlanes=None

# Field name for number of deblend children
config.forcedPhotCoadd.measurement.plugins['ext_shapeHSM_HsmShapeBj'].deblendNChild=None

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.plugins['ext_shapeHSM_HsmShapeLinear'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.forcedPhotCoadd.measurement.plugins['ext_shapeHSM_HsmShapeLinear'].badMaskPlanes=None

# Field name for number of deblend children
config.forcedPhotCoadd.measurement.plugins['ext_shapeHSM_HsmShapeLinear'].deblendNChild=None

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.plugins['ext_shapeHSM_HsmShapeKsb'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.forcedPhotCoadd.measurement.plugins['ext_shapeHSM_HsmShapeKsb'].badMaskPlanes=None

# Field name for number of deblend children
config.forcedPhotCoadd.measurement.plugins['ext_shapeHSM_HsmShapeKsb'].deblendNChild=None

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.plugins['ext_shapeHSM_HsmShapeRegauss'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.forcedPhotCoadd.measurement.plugins['ext_shapeHSM_HsmShapeRegauss'].badMaskPlanes=None

# Field name for number of deblend children
config.forcedPhotCoadd.measurement.plugins['ext_shapeHSM_HsmShapeRegauss'].deblendNChild=None

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.plugins['ext_shapeHSM_HsmSourceMoments'].doMeasure=True

# Store measured flux?
config.forcedPhotCoadd.measurement.plugins['ext_shapeHSM_HsmSourceMoments'].addFlux=None

# Mask planes used to reject bad pixels.
config.forcedPhotCoadd.measurement.plugins['ext_shapeHSM_HsmSourceMoments'].badMaskPlanes=None

# Use round weight function?
config.forcedPhotCoadd.measurement.plugins['ext_shapeHSM_HsmSourceMoments'].roundMoments=None

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.plugins['ext_shapeHSM_HsmSourceMomentsRound'].doMeasure=True

# Store measured flux?
config.forcedPhotCoadd.measurement.plugins['ext_shapeHSM_HsmSourceMomentsRound'].addFlux=None

# Mask planes used to reject bad pixels.
config.forcedPhotCoadd.measurement.plugins['ext_shapeHSM_HsmSourceMomentsRound'].badMaskPlanes=None

# Use round weight function?
config.forcedPhotCoadd.measurement.plugins['ext_shapeHSM_HsmSourceMomentsRound'].roundMoments=None

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.plugins['ext_shapeHSM_HsmPsfMoments'].doMeasure=True

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].doMeasure=True

# Shapelet order of inner expansion (0 == Gaussian)
config.forcedPhotCoadd.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].innerOrder=2

# Don't allow the semi-major radius of any component to go above this fraction of the PSF image width
config.forcedPhotCoadd.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].maxRadiusBoxFraction=0.4

# Don't allow the semi-minor radius of any component to drop below this value (pixels)
config.forcedPhotCoadd.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].minRadius=1.0

# Don't allow the determinant radii of the two components to differ by less than this (pixels)
config.forcedPhotCoadd.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].minRadiusDiff=0.5

# whether to save all iterations for debugging purposes
config.forcedPhotCoadd.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.forcedPhotCoadd.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.forcedPhotCoadd.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.maxInnerIterations=20

# maximum number of steps
config.forcedPhotCoadd.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.forcedPhotCoadd.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.forcedPhotCoadd.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.forcedPhotCoadd.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.forcedPhotCoadd.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.forcedPhotCoadd.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.forcedPhotCoadd.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.forcedPhotCoadd.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.forcedPhotCoadd.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.forcedPhotCoadd.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.forcedPhotCoadd.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.forcedPhotCoadd.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.forcedPhotCoadd.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.forcedPhotCoadd.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.forcedPhotCoadd.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionSolverTolerance=1e-08

# Shapelet order of outer expansion (0 == Gaussian)
config.forcedPhotCoadd.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].outerOrder=1

# Initial outer Gaussian peak height divided by inner Gaussian peak height
config.forcedPhotCoadd.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].peakRatio=0.1

# Initial outer radius divided by inner radius
config.forcedPhotCoadd.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].radiusRatio=2.0

config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models={}
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.maxInnerIterations=20

# maximum number of steps
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.radiusPriorSigma=0.5

config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.maxInnerIterations=20

# maximum number of steps
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.radiusPriorSigma=0.5

config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.maxInnerIterations=20

# maximum number of steps
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.order=2

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.order=1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.radiusPriorSigma=0.5

config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.maxInnerIterations=20

# maximum number of steps
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.order=4

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.order=4

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.radiusPriorSigma=0.5

# a sequence of model names indicating which models should be fit, and their order
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].sequence=['DoubleShapelet']

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].doMeasure=True

# Whether to record the steps the optimizer takes (or just the number, if running as a plugin)
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].dev.doRecordHistory=True

# Whether to record the time spent in this stage
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].dev.doRecordTime=True

# Softened core width for ellipticity distribution (conformal shear units).
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].dev.empiricalPriorConfig.ellipticityCore=0.001

# Width of exponential ellipticity distribution (conformal shear units).
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].dev.empiricalPriorConfig.ellipticitySigma=0.3

# ln(radius) at which the softened cutoff begins towards the minimum
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusMinInner=-6.0

# Minimum ln(radius).
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusMinOuter=-6.001

# Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T.
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusMu=-1.0

# Number of degrees of freedom for the Student's T distribution on ln(radius).
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusNu=50.0

# Width of the Student's T distribution in ln(radius).
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusSigma=0.45

# Ellipticity magnitude (conformal shear units) at which the softened cutoff begins
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].dev.linearPriorConfig.ellipticityMaxInner=2.0

# Maximum ellipticity magnitude (conformal shear units)
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].dev.linearPriorConfig.ellipticityMaxOuter=2.001

# ln(radius) at which the softened cutoff begins towards the maximum
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].dev.linearPriorConfig.logRadiusMaxInner=3.0

# Maximum ln(radius)
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].dev.linearPriorConfig.logRadiusMaxOuter=3.001

# ln(radius) at which the softened cutoff begins towards the minimum
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].dev.linearPriorConfig.logRadiusMinInner=-6.0

# The ratio P(logRadiusMinInner)/P(logRadiusMaxInner)
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].dev.linearPriorConfig.logRadiusMinMaxRatio=1.0

# Minimum ln(radius)
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].dev.linearPriorConfig.logRadiusMinOuter=-6.001

# Maximum radius used in approximating profile with Gaussians (0=default for this profile)
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].dev.maxRadius=0

# Number of Gaussian used to approximate the profile
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].dev.nComponents=8

# whether to save all iterations for debugging purposes
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].dev.optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].dev.optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].dev.optimizer.maxInnerIterations=20

# maximum number of steps
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].dev.optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].dev.optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].dev.optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].dev.optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].dev.optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].dev.optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].dev.optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].dev.optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].dev.optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].dev.optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].dev.optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].dev.optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].dev.optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].dev.optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].dev.optimizer.trustRegionSolverTolerance=1e-08

# Name of the Prior that defines the model to fit (a filename in $MEAS_MODELFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting.
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].dev.priorName=''

# One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].dev.priorSource='EMPIRICAL'

# Name of the shapelet.RadialProfile that defines the model to fit
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].dev.profileName='luv'

# Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances)
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].dev.usePixelWeights=False

# Scale the likelihood by this factor to artificially reweight it w.r.t. the prior.
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].dev.weightsMultiplier=1.0

# Whether to record the steps the optimizer takes (or just the number, if running as a plugin)
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].exp.doRecordHistory=True

# Whether to record the time spent in this stage
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].exp.doRecordTime=True

# Softened core width for ellipticity distribution (conformal shear units).
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].exp.empiricalPriorConfig.ellipticityCore=0.001

# Width of exponential ellipticity distribution (conformal shear units).
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].exp.empiricalPriorConfig.ellipticitySigma=0.3

# ln(radius) at which the softened cutoff begins towards the minimum
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusMinInner=-6.0

# Minimum ln(radius).
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusMinOuter=-6.001

# Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T.
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusMu=-1.0

# Number of degrees of freedom for the Student's T distribution on ln(radius).
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusNu=50.0

# Width of the Student's T distribution in ln(radius).
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusSigma=0.45

# Ellipticity magnitude (conformal shear units) at which the softened cutoff begins
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].exp.linearPriorConfig.ellipticityMaxInner=2.0

# Maximum ellipticity magnitude (conformal shear units)
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].exp.linearPriorConfig.ellipticityMaxOuter=2.001

# ln(radius) at which the softened cutoff begins towards the maximum
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].exp.linearPriorConfig.logRadiusMaxInner=3.0

# Maximum ln(radius)
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].exp.linearPriorConfig.logRadiusMaxOuter=3.001

# ln(radius) at which the softened cutoff begins towards the minimum
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].exp.linearPriorConfig.logRadiusMinInner=-6.0

# The ratio P(logRadiusMinInner)/P(logRadiusMaxInner)
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].exp.linearPriorConfig.logRadiusMinMaxRatio=1.0

# Minimum ln(radius)
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].exp.linearPriorConfig.logRadiusMinOuter=-6.001

# Maximum radius used in approximating profile with Gaussians (0=default for this profile)
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].exp.maxRadius=0

# Number of Gaussian used to approximate the profile
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].exp.nComponents=6

# whether to save all iterations for debugging purposes
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].exp.optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].exp.optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].exp.optimizer.maxInnerIterations=20

# maximum number of steps
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].exp.optimizer.maxOuterIterations=250

# If the trust radius falls below this threshold, consider the algorithm converged
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].exp.optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].exp.optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].exp.optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].exp.optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].exp.optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].exp.optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].exp.optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].exp.optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].exp.optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].exp.optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].exp.optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].exp.optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].exp.optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].exp.optimizer.trustRegionSolverTolerance=1e-08

# Name of the Prior that defines the model to fit (a filename in $MEAS_MODELFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting.
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].exp.priorName=''

# One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].exp.priorSource='EMPIRICAL'

# Name of the shapelet.RadialProfile that defines the model to fit
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].exp.profileName='lux'

# Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances)
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].exp.usePixelWeights=False

# Scale the likelihood by this factor to artificially reweight it w.r.t. the prior.
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].exp.weightsMultiplier=1.0

# If the 2nd-moments shape used to initialize the fit failed, use the PSF moments multiplied by this.  If <= 0.0, abort the fit early instead.
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].fallbackInitialMomentsPsfFactor=1.5

# Whether to record the steps the optimizer takes (or just the number, if running as a plugin)
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].initial.doRecordHistory=True

# Whether to record the time spent in this stage
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].initial.doRecordTime=True

# Softened core width for ellipticity distribution (conformal shear units).
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].initial.empiricalPriorConfig.ellipticityCore=0.001

# Width of exponential ellipticity distribution (conformal shear units).
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].initial.empiricalPriorConfig.ellipticitySigma=0.3

# ln(radius) at which the softened cutoff begins towards the minimum
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusMinInner=-6.0

# Minimum ln(radius).
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusMinOuter=-6.001

# Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T.
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusMu=-1.0

# Number of degrees of freedom for the Student's T distribution on ln(radius).
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusNu=50.0

# Width of the Student's T distribution in ln(radius).
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusSigma=0.45

# Ellipticity magnitude (conformal shear units) at which the softened cutoff begins
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].initial.linearPriorConfig.ellipticityMaxInner=2.0

# Maximum ellipticity magnitude (conformal shear units)
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].initial.linearPriorConfig.ellipticityMaxOuter=2.001

# ln(radius) at which the softened cutoff begins towards the maximum
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].initial.linearPriorConfig.logRadiusMaxInner=3.0

# Maximum ln(radius)
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].initial.linearPriorConfig.logRadiusMaxOuter=3.001

# ln(radius) at which the softened cutoff begins towards the minimum
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].initial.linearPriorConfig.logRadiusMinInner=-6.0

# The ratio P(logRadiusMinInner)/P(logRadiusMaxInner)
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].initial.linearPriorConfig.logRadiusMinMaxRatio=1.0

# Minimum ln(radius)
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].initial.linearPriorConfig.logRadiusMinOuter=-6.001

# Maximum radius used in approximating profile with Gaussians (0=default for this profile)
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].initial.maxRadius=0

# Number of Gaussian used to approximate the profile
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].initial.nComponents=3

# whether to save all iterations for debugging purposes
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].initial.optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].initial.optimizer.gradientThreshold=0.001

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].initial.optimizer.maxInnerIterations=20

# maximum number of steps
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].initial.optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].initial.optimizer.minTrustRadiusThreshold=0.01

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].initial.optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].initial.optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].initial.optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].initial.optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].initial.optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].initial.optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].initial.optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].initial.optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].initial.optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].initial.optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].initial.optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].initial.optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].initial.optimizer.trustRegionSolverTolerance=1e-08

# Name of the Prior that defines the model to fit (a filename in $MEAS_MODELFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting.
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].initial.priorName=''

# One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].initial.priorSource='EMPIRICAL'

# Name of the shapelet.RadialProfile that defines the model to fit
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].initial.profileName='lux'

# Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances)
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].initial.usePixelWeights=True

# Scale the likelihood by this factor to artificially reweight it w.r.t. the prior.
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].initial.weightsMultiplier=1.0

# Minimum initial radius in pixels (used to regularize initial moments-based PSF deconvolution)
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].minInitialRadius=0.1

# Field name prefix of the Shapelet PSF approximation used to convolve the galaxy model; must contain a set of fields matching the schema defined by shapelet.MultiShapeletFunctionKey.
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].psfName='modelfit_DoubleShapeletPsfApprox'

# Mask planes that indicate pixels that should be ignored in the fit.
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].region.badMaskPlanes=['EDGE', 'SAT', 'BAD', 'NO_DATA']

# Abort if the fit region grows beyond this many pixels.
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].region.maxArea=100000

# Maximum fraction of pixels that may be ignored due to masks; more than this and we don't even try.
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].region.maxBadPixelFraction=0.1

# Use this multiple of the initial fit ellipse then grow by the PSF width to determine the maximum final fit region size.
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].region.nFitRadiiMax=3.0

# Use this multiple of the initial fit ellipse then grow by the PSF width to determine the minimum final fit region size.
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].region.nFitRadiiMin=1.0

# Use this multiple of the Kron ellipse to set the fit region (for the final fit region, subject to the nFitRadiiMin and nFitRadiiMax constraints).
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].region.nKronRadii=1.5

# Grow the initial fit ellipses by this factor before comparing with the Kron/Footprint region
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].region.nPsfSigmaGrow=2.0

# If the Kron radius is less than this multiple of the PSF width, ignore it and fall back to a PSF-oriented ellipse scaled to match the area of the footprint or this radius (whichever is larger).
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].region.nPsfSigmaMin=4.0

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.plugins['modelfit_CModel'].doMeasure=True

config.forcedPhotCoadd.measurement.plugins.names=['ext_photometryKron_KronFlux', 'base_GaussianFlux', 'base_PixelFlags', 'ext_convolved_ConvolvedFlux', 'base_CircularApertureFlux', 'base_TransformedCentroid', 'base_LocalBackground', 'modelfit_DoubleShapeletPsfApprox', 'base_SdssShape', 'base_TransformedShape', 'base_InputCount', 'base_Variance', 'base_PsfFlux', 'base_SdssCentroid', 'modelfit_CModel']
# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.undeblended['base_PsfFlux'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.forcedPhotCoadd.measurement.undeblended['base_PsfFlux'].badMaskPlanes=[]

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.undeblended['base_PeakLikelihoodFlux'].doMeasure=True

# Name of warping kernel (e.g. "lanczos4") used to compute the peak
config.forcedPhotCoadd.measurement.undeblended['base_PeakLikelihoodFlux'].warpingKernelName='lanczos4'

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.undeblended['base_GaussianFlux'].doMeasure=True

# FIXME! NEVER DOCUMENTED!
config.forcedPhotCoadd.measurement.undeblended['base_GaussianFlux'].background=0.0

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.undeblended['base_NaiveCentroid'].doMeasure=True

# Value to subtract from the image pixel values
config.forcedPhotCoadd.measurement.undeblended['base_NaiveCentroid'].background=0.0

# Do check that the centroid is contained in footprint.
config.forcedPhotCoadd.measurement.undeblended['base_NaiveCentroid'].doFootprintCheck=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.forcedPhotCoadd.measurement.undeblended['base_NaiveCentroid'].maxDistToPeak=-1.0

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.undeblended['base_SdssCentroid'].doMeasure=True

# maximum allowed binning
config.forcedPhotCoadd.measurement.undeblended['base_SdssCentroid'].binmax=16

# Do check that the centroid is contained in footprint.
config.forcedPhotCoadd.measurement.undeblended['base_SdssCentroid'].doFootprintCheck=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.forcedPhotCoadd.measurement.undeblended['base_SdssCentroid'].maxDistToPeak=-1.0

# if the peak's less than this insist on binning at least once
config.forcedPhotCoadd.measurement.undeblended['base_SdssCentroid'].peakMin=-1.0

# fiddle factor for adjusting the binning
config.forcedPhotCoadd.measurement.undeblended['base_SdssCentroid'].wfac=1.5

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.undeblended['base_PixelFlags'].doMeasure=True

# List of mask planes to be searched for which occur anywhere within a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.forcedPhotCoadd.measurement.undeblended['base_PixelFlags'].masksFpAnywhere=[]

# List of mask planes to be searched for which occur in the center of a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.forcedPhotCoadd.measurement.undeblended['base_PixelFlags'].masksFpCenter=[]

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.undeblended['base_SdssShape'].doMeasure=True

# Additional value to add to background
config.forcedPhotCoadd.measurement.undeblended['base_SdssShape'].background=0.0

# Whether to also compute the shape of the PSF model
config.forcedPhotCoadd.measurement.undeblended['base_SdssShape'].doMeasurePsf=True

# Maximum number of iterations
config.forcedPhotCoadd.measurement.undeblended['base_SdssShape'].maxIter=100

# Maximum centroid shift, limited to 2-10
config.forcedPhotCoadd.measurement.undeblended['base_SdssShape'].maxShift=0.0

# Convergence tolerance for e1,e2
config.forcedPhotCoadd.measurement.undeblended['base_SdssShape'].tol1=9.999999747378752e-06

# Convergence tolerance for FWHM
config.forcedPhotCoadd.measurement.undeblended['base_SdssShape'].tol2=9.999999747378752e-05

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.undeblended['base_ScaledApertureFlux'].doMeasure=True

# Scaling factor of PSF FWHM for aperture radius.
config.forcedPhotCoadd.measurement.undeblended['base_ScaledApertureFlux'].scale=3.14

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.forcedPhotCoadd.measurement.undeblended['base_ScaledApertureFlux'].shiftKernel='lanczos5'

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.undeblended['base_CircularApertureFlux'].doMeasure=True

# Maximum radius (in pixels) for which the sinc algorithm should be used instead of the faster naive algorithm.  For elliptical apertures, this is the minor axis radius.
config.forcedPhotCoadd.measurement.undeblended['base_CircularApertureFlux'].maxSincRadius=12.0

# Radius (in pixels) of apertures.
config.forcedPhotCoadd.measurement.undeblended['base_CircularApertureFlux'].radii=[3.0, 4.5, 6.0, 9.0, 12.0, 17.0, 25.0, 35.0, 50.0, 70.0]

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.forcedPhotCoadd.measurement.undeblended['base_CircularApertureFlux'].shiftKernel='lanczos5'

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.undeblended['base_Blendedness'].doMeasure=True

# Whether to compute quantities related to the Gaussian-weighted flux
config.forcedPhotCoadd.measurement.undeblended['base_Blendedness'].doFlux=True

# Whether to compute HeavyFootprint dot products (the old deblend.blendedness parameter)
config.forcedPhotCoadd.measurement.undeblended['base_Blendedness'].doOld=True

# Whether to compute quantities related to the Gaussian-weighted shape
config.forcedPhotCoadd.measurement.undeblended['base_Blendedness'].doShape=True

# Radius factor that sets the maximum extent of the weight function (and hence the flux measurements)
config.forcedPhotCoadd.measurement.undeblended['base_Blendedness'].nSigmaWeightMax=3.0

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.undeblended['base_LocalBackground'].doMeasure=True

# Inner radius for background annulus as a multiple of the PSF sigma
config.forcedPhotCoadd.measurement.undeblended['base_LocalBackground'].annulusInner=7.0

# Outer radius for background annulus as a multiple of the PSF sigma
config.forcedPhotCoadd.measurement.undeblended['base_LocalBackground'].annulusOuter=15.0

# Mask planes that indicate pixels that should be excluded from the measurement
config.forcedPhotCoadd.measurement.undeblended['base_LocalBackground'].badMaskPlanes=['BAD', 'SAT', 'NO_DATA']

# Number of sigma-clipping iterations for background measurement
config.forcedPhotCoadd.measurement.undeblended['base_LocalBackground'].bgIter=3

# Rejection threshold (in standard deviations) for background measurement
config.forcedPhotCoadd.measurement.undeblended['base_LocalBackground'].bgRej=3.0

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.undeblended['base_Variance'].doMeasure=True

# Scale factor to apply to shape for aperture
config.forcedPhotCoadd.measurement.undeblended['base_Variance'].scale=5.0

# Mask planes to ignore
config.forcedPhotCoadd.measurement.undeblended['base_Variance'].mask=['DETECTED', 'DETECTED_NEGATIVE', 'BAD', 'SAT']

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.undeblended['base_InputCount'].doMeasure=True

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.undeblended['base_PeakCentroid'].doMeasure=True

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.undeblended['base_TransformedCentroid'].doMeasure=True

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.undeblended['base_TransformedShape'].doMeasure=True

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.undeblended['ext_photometryKron_KronFlux'].doMeasure=True

# If true check that the Kron radius exceeds some minimum
config.forcedPhotCoadd.measurement.undeblended['ext_photometryKron_KronFlux'].enforceMinimumRadius=True

# if true, use existing shape and centroid measurements instead of fitting
config.forcedPhotCoadd.measurement.undeblended['ext_photometryKron_KronFlux'].fixed=False

# Largest aperture for which to use the slow, accurate, sinc aperture code
config.forcedPhotCoadd.measurement.undeblended['ext_photometryKron_KronFlux'].maxSincRadius=10.0

# Minimum Kron radius (if == 0.0 use PSF's Kron radius) if enforceMinimumRadius. Also functions as fallback aperture radius if set.
config.forcedPhotCoadd.measurement.undeblended['ext_photometryKron_KronFlux'].minimumRadius=0.0

# Number of times to iterate when setting the Kron radius
config.forcedPhotCoadd.measurement.undeblended['ext_photometryKron_KronFlux'].nIterForRadius=1

# Number of Kron radii for Kron flux
config.forcedPhotCoadd.measurement.undeblended['ext_photometryKron_KronFlux'].nRadiusForFlux=2.5

# Multiplier of rms size for aperture used to initially estimate the Kron radius
config.forcedPhotCoadd.measurement.undeblended['ext_photometryKron_KronFlux'].nSigmaForRadius=6.0

# Name of field specifying reference Kron radius for forced measurement
config.forcedPhotCoadd.measurement.undeblended['ext_photometryKron_KronFlux'].refRadiusName='ext_photometryKron_KronFlux_radius'

# Smooth image with N(0, smoothingSigma^2) Gaussian while estimating R_K
config.forcedPhotCoadd.measurement.undeblended['ext_photometryKron_KronFlux'].smoothingSigma=-1.0

# Use the Footprint size as part of initial estimate of Kron radius
config.forcedPhotCoadd.measurement.undeblended['ext_photometryKron_KronFlux'].useFootprintRadius=False

# list of target seeings (FWHM, pixels)
config.forcedPhotCoadd.measurement.undeblended['ext_convolved_ConvolvedFlux'].seeing=[3.5, 5.0, 6.5, 8.0]

# scaling factor of kernel sigma for kernel size
config.forcedPhotCoadd.measurement.undeblended['ext_convolved_ConvolvedFlux'].kernelScale=4.0

# Maximum radius (in pixels) for which the sinc algorithm should be used instead of the faster naive algorithm.  For elliptical apertures, this is the minor axis radius.
config.forcedPhotCoadd.measurement.undeblended['ext_convolved_ConvolvedFlux'].aperture.maxSincRadius=10.0

# Radius (in pixels) of apertures.
config.forcedPhotCoadd.measurement.undeblended['ext_convolved_ConvolvedFlux'].aperture.radii=[3.3, 4.5, 6.0]

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.forcedPhotCoadd.measurement.undeblended['ext_convolved_ConvolvedFlux'].aperture.shiftKernel='lanczos5'

# name of Kron radius field in reference
config.forcedPhotCoadd.measurement.undeblended['ext_convolved_ConvolvedFlux'].kronRadiusName='ext_photometryKron_KronFlux_radius'

# Largest aperture for which to use the sinc aperture code for Kron (pixels)
config.forcedPhotCoadd.measurement.undeblended['ext_convolved_ConvolvedFlux'].maxSincRadius=10.0

# Number of Kron radii for Kron flux
config.forcedPhotCoadd.measurement.undeblended['ext_convolved_ConvolvedFlux'].kronRadiusForFlux=2.5

# Register measurements for aperture correction?
# The aperture correction registration is done when the plugin is
# instantiated because the column names are derived from the configuration
# rather than being static. Sometimes you want to turn this off, e.g.,
# when you will use aperture corrections derived from somewhere else
# through the 'proxy' mechanism.
config.forcedPhotCoadd.measurement.undeblended['ext_convolved_ConvolvedFlux'].registerForApCorr=False

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.undeblended['ext_convolved_ConvolvedFlux'].doMeasure=True

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.undeblended['ext_shapeHSM_HsmShapeBj'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.forcedPhotCoadd.measurement.undeblended['ext_shapeHSM_HsmShapeBj'].badMaskPlanes=None

# Field name for number of deblend children
config.forcedPhotCoadd.measurement.undeblended['ext_shapeHSM_HsmShapeBj'].deblendNChild=None

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.undeblended['ext_shapeHSM_HsmShapeLinear'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.forcedPhotCoadd.measurement.undeblended['ext_shapeHSM_HsmShapeLinear'].badMaskPlanes=None

# Field name for number of deblend children
config.forcedPhotCoadd.measurement.undeblended['ext_shapeHSM_HsmShapeLinear'].deblendNChild=None

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.undeblended['ext_shapeHSM_HsmShapeKsb'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.forcedPhotCoadd.measurement.undeblended['ext_shapeHSM_HsmShapeKsb'].badMaskPlanes=None

# Field name for number of deblend children
config.forcedPhotCoadd.measurement.undeblended['ext_shapeHSM_HsmShapeKsb'].deblendNChild=None

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.undeblended['ext_shapeHSM_HsmShapeRegauss'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.forcedPhotCoadd.measurement.undeblended['ext_shapeHSM_HsmShapeRegauss'].badMaskPlanes=None

# Field name for number of deblend children
config.forcedPhotCoadd.measurement.undeblended['ext_shapeHSM_HsmShapeRegauss'].deblendNChild=None

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.undeblended['ext_shapeHSM_HsmSourceMoments'].doMeasure=True

# Store measured flux?
config.forcedPhotCoadd.measurement.undeblended['ext_shapeHSM_HsmSourceMoments'].addFlux=None

# Mask planes used to reject bad pixels.
config.forcedPhotCoadd.measurement.undeblended['ext_shapeHSM_HsmSourceMoments'].badMaskPlanes=None

# Use round weight function?
config.forcedPhotCoadd.measurement.undeblended['ext_shapeHSM_HsmSourceMoments'].roundMoments=None

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.undeblended['ext_shapeHSM_HsmSourceMomentsRound'].doMeasure=True

# Store measured flux?
config.forcedPhotCoadd.measurement.undeblended['ext_shapeHSM_HsmSourceMomentsRound'].addFlux=None

# Mask planes used to reject bad pixels.
config.forcedPhotCoadd.measurement.undeblended['ext_shapeHSM_HsmSourceMomentsRound'].badMaskPlanes=None

# Use round weight function?
config.forcedPhotCoadd.measurement.undeblended['ext_shapeHSM_HsmSourceMomentsRound'].roundMoments=None

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.undeblended['ext_shapeHSM_HsmPsfMoments'].doMeasure=True

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].doMeasure=True

# Shapelet order of inner expansion (0 == Gaussian)
config.forcedPhotCoadd.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].innerOrder=2

# Don't allow the semi-major radius of any component to go above this fraction of the PSF image width
config.forcedPhotCoadd.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].maxRadiusBoxFraction=0.4

# Don't allow the semi-minor radius of any component to drop below this value (pixels)
config.forcedPhotCoadd.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].minRadius=1.0

# Don't allow the determinant radii of the two components to differ by less than this (pixels)
config.forcedPhotCoadd.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].minRadiusDiff=0.5

# whether to save all iterations for debugging purposes
config.forcedPhotCoadd.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.forcedPhotCoadd.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.forcedPhotCoadd.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.maxInnerIterations=20

# maximum number of steps
config.forcedPhotCoadd.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.forcedPhotCoadd.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.forcedPhotCoadd.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.forcedPhotCoadd.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.forcedPhotCoadd.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.forcedPhotCoadd.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.forcedPhotCoadd.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.forcedPhotCoadd.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.forcedPhotCoadd.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.forcedPhotCoadd.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.forcedPhotCoadd.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.forcedPhotCoadd.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.forcedPhotCoadd.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.forcedPhotCoadd.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.forcedPhotCoadd.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionSolverTolerance=1e-08

# Shapelet order of outer expansion (0 == Gaussian)
config.forcedPhotCoadd.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].outerOrder=1

# Initial outer Gaussian peak height divided by inner Gaussian peak height
config.forcedPhotCoadd.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].peakRatio=0.1

# Initial outer radius divided by inner radius
config.forcedPhotCoadd.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].radiusRatio=2.0

config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models={}
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.maxInnerIterations=20

# maximum number of steps
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.radiusPriorSigma=0.5

config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.maxInnerIterations=20

# maximum number of steps
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.radiusPriorSigma=0.5

config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.maxInnerIterations=20

# maximum number of steps
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.order=2

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.order=1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.radiusPriorSigma=0.5

config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.maxInnerIterations=20

# maximum number of steps
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.order=4

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.order=4

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.radiusPriorSigma=0.5

# a sequence of model names indicating which models should be fit, and their order
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].sequence=['DoubleShapelet']

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].doMeasure=True

# Whether to record the steps the optimizer takes (or just the number, if running as a plugin)
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].dev.doRecordHistory=True

# Whether to record the time spent in this stage
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].dev.doRecordTime=True

# Softened core width for ellipticity distribution (conformal shear units).
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].dev.empiricalPriorConfig.ellipticityCore=0.001

# Width of exponential ellipticity distribution (conformal shear units).
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].dev.empiricalPriorConfig.ellipticitySigma=0.3

# ln(radius) at which the softened cutoff begins towards the minimum
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusMinInner=-6.0

# Minimum ln(radius).
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusMinOuter=-6.001

# Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T.
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusMu=-1.0

# Number of degrees of freedom for the Student's T distribution on ln(radius).
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusNu=50.0

# Width of the Student's T distribution in ln(radius).
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusSigma=0.45

# Ellipticity magnitude (conformal shear units) at which the softened cutoff begins
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].dev.linearPriorConfig.ellipticityMaxInner=2.0

# Maximum ellipticity magnitude (conformal shear units)
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].dev.linearPriorConfig.ellipticityMaxOuter=2.001

# ln(radius) at which the softened cutoff begins towards the maximum
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].dev.linearPriorConfig.logRadiusMaxInner=3.0

# Maximum ln(radius)
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].dev.linearPriorConfig.logRadiusMaxOuter=3.001

# ln(radius) at which the softened cutoff begins towards the minimum
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].dev.linearPriorConfig.logRadiusMinInner=-6.0

# The ratio P(logRadiusMinInner)/P(logRadiusMaxInner)
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].dev.linearPriorConfig.logRadiusMinMaxRatio=1.0

# Minimum ln(radius)
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].dev.linearPriorConfig.logRadiusMinOuter=-6.001

# Maximum radius used in approximating profile with Gaussians (0=default for this profile)
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].dev.maxRadius=0

# Number of Gaussian used to approximate the profile
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].dev.nComponents=8

# whether to save all iterations for debugging purposes
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].dev.optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].dev.optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].dev.optimizer.maxInnerIterations=20

# maximum number of steps
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].dev.optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].dev.optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].dev.optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].dev.optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].dev.optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].dev.optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].dev.optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].dev.optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].dev.optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].dev.optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].dev.optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].dev.optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].dev.optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].dev.optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].dev.optimizer.trustRegionSolverTolerance=1e-08

# Name of the Prior that defines the model to fit (a filename in $MEAS_MODELFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting.
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].dev.priorName=''

# One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].dev.priorSource='EMPIRICAL'

# Name of the shapelet.RadialProfile that defines the model to fit
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].dev.profileName='luv'

# Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances)
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].dev.usePixelWeights=False

# Scale the likelihood by this factor to artificially reweight it w.r.t. the prior.
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].dev.weightsMultiplier=1.0

# Whether to record the steps the optimizer takes (or just the number, if running as a plugin)
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].exp.doRecordHistory=True

# Whether to record the time spent in this stage
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].exp.doRecordTime=True

# Softened core width for ellipticity distribution (conformal shear units).
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].exp.empiricalPriorConfig.ellipticityCore=0.001

# Width of exponential ellipticity distribution (conformal shear units).
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].exp.empiricalPriorConfig.ellipticitySigma=0.3

# ln(radius) at which the softened cutoff begins towards the minimum
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusMinInner=-6.0

# Minimum ln(radius).
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusMinOuter=-6.001

# Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T.
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusMu=-1.0

# Number of degrees of freedom for the Student's T distribution on ln(radius).
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusNu=50.0

# Width of the Student's T distribution in ln(radius).
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusSigma=0.45

# Ellipticity magnitude (conformal shear units) at which the softened cutoff begins
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].exp.linearPriorConfig.ellipticityMaxInner=2.0

# Maximum ellipticity magnitude (conformal shear units)
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].exp.linearPriorConfig.ellipticityMaxOuter=2.001

# ln(radius) at which the softened cutoff begins towards the maximum
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].exp.linearPriorConfig.logRadiusMaxInner=3.0

# Maximum ln(radius)
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].exp.linearPriorConfig.logRadiusMaxOuter=3.001

# ln(radius) at which the softened cutoff begins towards the minimum
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].exp.linearPriorConfig.logRadiusMinInner=-6.0

# The ratio P(logRadiusMinInner)/P(logRadiusMaxInner)
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].exp.linearPriorConfig.logRadiusMinMaxRatio=1.0

# Minimum ln(radius)
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].exp.linearPriorConfig.logRadiusMinOuter=-6.001

# Maximum radius used in approximating profile with Gaussians (0=default for this profile)
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].exp.maxRadius=0

# Number of Gaussian used to approximate the profile
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].exp.nComponents=6

# whether to save all iterations for debugging purposes
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].exp.optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].exp.optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].exp.optimizer.maxInnerIterations=20

# maximum number of steps
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].exp.optimizer.maxOuterIterations=250

# If the trust radius falls below this threshold, consider the algorithm converged
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].exp.optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].exp.optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].exp.optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].exp.optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].exp.optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].exp.optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].exp.optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].exp.optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].exp.optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].exp.optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].exp.optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].exp.optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].exp.optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].exp.optimizer.trustRegionSolverTolerance=1e-08

# Name of the Prior that defines the model to fit (a filename in $MEAS_MODELFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting.
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].exp.priorName=''

# One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].exp.priorSource='EMPIRICAL'

# Name of the shapelet.RadialProfile that defines the model to fit
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].exp.profileName='lux'

# Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances)
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].exp.usePixelWeights=False

# Scale the likelihood by this factor to artificially reweight it w.r.t. the prior.
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].exp.weightsMultiplier=1.0

# If the 2nd-moments shape used to initialize the fit failed, use the PSF moments multiplied by this.  If <= 0.0, abort the fit early instead.
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].fallbackInitialMomentsPsfFactor=1.5

# Whether to record the steps the optimizer takes (or just the number, if running as a plugin)
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].initial.doRecordHistory=True

# Whether to record the time spent in this stage
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].initial.doRecordTime=True

# Softened core width for ellipticity distribution (conformal shear units).
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].initial.empiricalPriorConfig.ellipticityCore=0.001

# Width of exponential ellipticity distribution (conformal shear units).
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].initial.empiricalPriorConfig.ellipticitySigma=0.3

# ln(radius) at which the softened cutoff begins towards the minimum
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusMinInner=-6.0

# Minimum ln(radius).
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusMinOuter=-6.001

# Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T.
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusMu=-1.0

# Number of degrees of freedom for the Student's T distribution on ln(radius).
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusNu=50.0

# Width of the Student's T distribution in ln(radius).
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusSigma=0.45

# Ellipticity magnitude (conformal shear units) at which the softened cutoff begins
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].initial.linearPriorConfig.ellipticityMaxInner=2.0

# Maximum ellipticity magnitude (conformal shear units)
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].initial.linearPriorConfig.ellipticityMaxOuter=2.001

# ln(radius) at which the softened cutoff begins towards the maximum
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].initial.linearPriorConfig.logRadiusMaxInner=3.0

# Maximum ln(radius)
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].initial.linearPriorConfig.logRadiusMaxOuter=3.001

# ln(radius) at which the softened cutoff begins towards the minimum
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].initial.linearPriorConfig.logRadiusMinInner=-6.0

# The ratio P(logRadiusMinInner)/P(logRadiusMaxInner)
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].initial.linearPriorConfig.logRadiusMinMaxRatio=1.0

# Minimum ln(radius)
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].initial.linearPriorConfig.logRadiusMinOuter=-6.001

# Maximum radius used in approximating profile with Gaussians (0=default for this profile)
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].initial.maxRadius=0

# Number of Gaussian used to approximate the profile
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].initial.nComponents=3

# whether to save all iterations for debugging purposes
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].initial.optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].initial.optimizer.gradientThreshold=0.001

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].initial.optimizer.maxInnerIterations=20

# maximum number of steps
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].initial.optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].initial.optimizer.minTrustRadiusThreshold=0.01

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].initial.optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].initial.optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].initial.optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].initial.optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].initial.optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].initial.optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].initial.optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].initial.optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].initial.optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].initial.optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].initial.optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].initial.optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].initial.optimizer.trustRegionSolverTolerance=1e-08

# Name of the Prior that defines the model to fit (a filename in $MEAS_MODELFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting.
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].initial.priorName=''

# One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].initial.priorSource='EMPIRICAL'

# Name of the shapelet.RadialProfile that defines the model to fit
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].initial.profileName='lux'

# Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances)
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].initial.usePixelWeights=True

# Scale the likelihood by this factor to artificially reweight it w.r.t. the prior.
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].initial.weightsMultiplier=1.0

# Minimum initial radius in pixels (used to regularize initial moments-based PSF deconvolution)
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].minInitialRadius=0.1

# Field name prefix of the Shapelet PSF approximation used to convolve the galaxy model; must contain a set of fields matching the schema defined by shapelet.MultiShapeletFunctionKey.
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].psfName='modelfit_DoubleShapeletPsfApprox'

# Mask planes that indicate pixels that should be ignored in the fit.
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].region.badMaskPlanes=['EDGE', 'SAT', 'BAD', 'NO_DATA']

# Abort if the fit region grows beyond this many pixels.
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].region.maxArea=100000

# Maximum fraction of pixels that may be ignored due to masks; more than this and we don't even try.
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].region.maxBadPixelFraction=0.1

# Use this multiple of the initial fit ellipse then grow by the PSF width to determine the maximum final fit region size.
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].region.nFitRadiiMax=3.0

# Use this multiple of the initial fit ellipse then grow by the PSF width to determine the minimum final fit region size.
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].region.nFitRadiiMin=1.0

# Use this multiple of the Kron ellipse to set the fit region (for the final fit region, subject to the nFitRadiiMin and nFitRadiiMax constraints).
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].region.nKronRadii=1.5

# Grow the initial fit ellipses by this factor before comparing with the Kron/Footprint region
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].region.nPsfSigmaGrow=2.0

# If the Kron radius is less than this multiple of the PSF width, ignore it and fall back to a PSF-oriented ellipse scaled to match the area of the footprint or this radius (whichever is larger).
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].region.nPsfSigmaMin=4.0

# whether to run this plugin in single-object mode
config.forcedPhotCoadd.measurement.undeblended['modelfit_CModel'].doMeasure=True

config.forcedPhotCoadd.measurement.undeblended.names=['ext_convolved_ConvolvedFlux', 'base_CircularApertureFlux', 'ext_photometryKron_KronFlux', 'base_PsfFlux']
# Mapping of reference columns to source columns
config.forcedPhotCoadd.measurement.copyColumns={'id': 'id', 'parent': 'parent', 'deblend_nChild': 'deblend_nChild', 'coord_ra': 'coord_ra', 'coord_dec': 'coord_dec'}

# Strictness of Astropy unit compatibility check, can be 'raise', 'warn' or 'silent'
config.forcedPhotCoadd.measurement.checkUnitsParseStrict='raise'

# coadd name: typically one of deep or goodSeeing
config.forcedPhotCoadd.coaddName='deep'

# Run subtask to apply aperture corrections
config.forcedPhotCoadd.doApCorr=True

# flux measurement algorithms in getApCorrNameSet() to ignore; if a name is listed that does not appear in getApCorrNameSet() then a warning is logged
config.forcedPhotCoadd.applyApCorr.ignoreList=[]

# set the general failure flag for a flux when it cannot be aperture-corrected?
config.forcedPhotCoadd.applyApCorr.doFlagApCorrFailures=True

# flux measurement algorithms to be aperture-corrected by reference to another algorithm; this is a mapping alg1:alg2, where 'alg1' is the algorithm being corrected, and 'alg2' is the algorithm supplying the corrections
config.forcedPhotCoadd.applyApCorr.proxies={'undeblended_base_PsfFlux_flux': 'base_PsfFlux_flux', 'undeblended_ext_photometryKron_KronFlux_flux': 'ext_photometryKron_KronFlux_flux', 'undeblended_ext_convolved_ConvolvedFlux_0_3_3': 'ext_convolved_ConvolvedFlux_0_3_3', 'undeblended_ext_convolved_ConvolvedFlux_0_4_5': 'ext_convolved_ConvolvedFlux_0_4_5', 'undeblended_ext_convolved_ConvolvedFlux_0_6_0': 'ext_convolved_ConvolvedFlux_0_6_0', 'undeblended_ext_convolved_ConvolvedFlux_1_3_3': 'ext_convolved_ConvolvedFlux_1_3_3', 'undeblended_ext_convolved_ConvolvedFlux_1_4_5': 'ext_convolved_ConvolvedFlux_1_4_5', 'undeblended_ext_convolved_ConvolvedFlux_1_6_0': 'ext_convolved_ConvolvedFlux_1_6_0', 'undeblended_ext_convolved_ConvolvedFlux_2_3_3': 'ext_convolved_ConvolvedFlux_2_3_3', 'undeblended_ext_convolved_ConvolvedFlux_2_4_5': 'ext_convolved_ConvolvedFlux_2_4_5', 'undeblended_ext_convolved_ConvolvedFlux_2_6_0': 'ext_convolved_ConvolvedFlux_2_6_0', 'undeblended_ext_convolved_ConvolvedFlux_3_3_3': 'ext_convolved_ConvolvedFlux_3_3_3', 'undeblended_ext_convolved_ConvolvedFlux_3_4_5': 'ext_convolved_ConvolvedFlux_3_4_5', 'undeblended_ext_convolved_ConvolvedFlux_3_6_0': 'ext_convolved_ConvolvedFlux_3_6_0', 'undeblended_ext_convolved_ConvolvedFlux_0_kron': 'ext_convolved_ConvolvedFlux_0_kron', 'undeblended_ext_convolved_ConvolvedFlux_1_kron': 'ext_convolved_ConvolvedFlux_1_kron', 'undeblended_ext_convolved_ConvolvedFlux_2_kron': 'ext_convolved_ConvolvedFlux_2_kron', 'undeblended_ext_convolved_ConvolvedFlux_3_kron': 'ext_convolved_ConvolvedFlux_3_kron'}

# critical ratio of model to psf flux
config.forcedPhotCoadd.catalogCalculation.plugins['base_ClassificationExtendedness'].fluxRatio=0.985

# correction factor for modelFlux error
config.forcedPhotCoadd.catalogCalculation.plugins['base_ClassificationExtendedness'].modelErrFactor=0.0

# correction factor for psfFlux error
config.forcedPhotCoadd.catalogCalculation.plugins['base_ClassificationExtendedness'].psfErrFactor=0.0

config.forcedPhotCoadd.catalogCalculation.plugins.names=['base_ClassificationExtendedness']
# Dataset (without coadd prefix) that should be used to obtain (Heavy)Footprints for sources. Must have IDs that match those of the reference catalog.If None, Footprints will be generated by transforming the reference Footprints.
config.forcedPhotCoadd.footprintDatasetName='meas'

# Are we reprocessing?
# 
# This exists as a workaround for large deblender footprints causing large memory use and/or very slow processing.  We refuse to deblend those footprints when running on a cluster and return to reprocess on a machine with larger memory or more time if we consider those footprints important to recover.
config.reprocessing=False

