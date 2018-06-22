import lsst.pipe.drivers.singleFrameDriver
assert type(config)==lsst.pipe.drivers.singleFrameDriver.SingleFrameDriverConfig, 'config is of type %s.%s instead of lsst.pipe.drivers.singleFrameDriver.SingleFrameDriverConfig' % (type(config).__module__, type(config).__name__)
import lsst.shapelet.gaussHermiteProjection
import lsst.meas.extensions.photometryKron
import lsst.meas.extensions.psfex.prefs
import lsst.shapelet.shapeletFunction
import lsst.meas.modelfit.cmodel.cmodelContinued
import lsst.meas.extensions.shapeHSM.hsmShapeControl
import lsst.meas.extensions.shapeHSM
import lsst.shapelet
import lsst.meas.modelfit.truncatedGaussian
import lsst.meas.extensions.psfex.version
import lsst.shapelet.shapeletFunction.shapeletFunctionContinued
import lsst.meas.modelfit.psf.psfContinued
import lsst.shapelet.version
import lsst.meas.extensions.shapeHSM.version
import eups.exceptions
import lsst.shapelet.multiShapeletFunction.multiShapeletFunctionContinued
import eups.distrib.builder
import lsst.shapelet.hermiteTransformMatrix
import lsst.pipe.tasks.setConfigFromEups
import lsst.obs.subaru.crosstalkYagi
import lsst.meas.modelfit.optimizer.optimizer
import eups.distrib.Repositories
import lsst.obs.subaru.crosstalk
import eups
import lsst.shapelet.multiShapeletBasis
import lsst.obs.subaru.strayLight.waveletCompression
import eups.VersionParser
import lsst.shapelet.radialProfile.radialProfile
import lsst.shapelet.tractor
import lsst.meas.extensions.psfex.field
import lsst.afw.display.rgb
import lsst.obs.subaru.strayLight.yStrayLight
import lsst.meas.modelfit.likelihood
import lsst.meas.extensions.convolved.convolved
import lsst.meas.modelfit.version
import lsst.shapelet.shapeletFunction.shapeletFunction
import lsst.meas.modelfit.mixture
import lsst.shapelet.matrixBuilder
import lsst.meas.modelfit.common
import lsst.meas.modelfit.pixelFitRegion.pixelFitRegion
import lsst.meas.extensions.photometryKron.version
import lsst.shapelet.multiShapeletFunction.multiShapeletFunction
import lsst.meas.modelfit.sampler
import lsst.shapelet.constants.constantsContinued
import lsst.shapelet.constants
import lsst.shapelet.basisEvaluator
import lsst.shapelet.radialProfile
import lsst.meas.modelfit.unitSystem
import lsst.meas.modelfit.multiModel
import eups.app
import eups.Product
import lsst.shapelet.radialProfile.radialProfileContinued
import lsst.meas.extensions.psfex
import lsst.obs.subaru.isr
import eups.table
import lsst.meas.modelfit.integrals
import lsst.meas.extensions.shapeHSM.hsmMomentsControl
import eups.distrib.eupspkg
import eups.tags
import eups.distrib.Distrib
import eups.distrib.tarball
import eups.distrib.pacman
import lsst.meas.modelfit.pixelFitRegion
import lsst.shapelet.constants.constants
import lsst.meas.extensions.psfex.psfexPsfDeterminer
import lsst.meas.extensions.convolved
import lsst.meas.extensions.convolved.version
import lsst.shapelet.generator
import eups.db
import eups.stack.ProductStack
import lsst.shapelet.functorKeys
import configparser
import eups.db.VersionFile
import lsst.meas.modelfit.pixelFitRegion.pixelFitRegionContinued
import lsst.obs.subaru.strayLight.rotatorAngle
import eups.stack.ProductFamily
import eups.db.Database
import lsst.meas.modelfit.psf.psf
import lsst.meas.modelfit.optimizer.optimizerContinued
import eups.db.ChainFile
import eups.hooks
import eups.Uses
import eups.cmd
import lsst.obs.hsc.vignette
import lsst.meas.extensions.psfex.psfexPsf
import lsst.meas.modelfit.cmodel
import eups.distrib.server
import lsst.meas.modelfit.model
import lsst.meas.modelfit.priors.priors
import lsst.shapelet.gaussHermiteConvolution
import eups.Eups
import lsst.meas.modelfit
import eups.stack
import lsst.shapelet.multiShapeletFunction
import eups.VersionCompare
import optparse
import lsst.meas.modelfit.cmodel.cmodel
import eups.lock
import lsst.meas.modelfit.psf
import lsst.meas.modelfit.priors.priorsContinued
import lsst.meas.extensions.photometryKron.photometryKron
import lsst.meas.modelfit.optimizer
import lsst.meas.modelfit.adaptiveImportanceSampler
import lsst.meas.modelfit.unitTransformedLikelihood
import lsst.obs.subaru.strayLight
import eups.distrib
import lsst.meas.extensions.psfex.psf
import eups.distrib.Repository
import eups.utils
import eups.distrib.DistribFactory
import lsst.meas.modelfit.priors
import lsst.shapelet.gaussHermiteProjection
import lsst.meas.extensions.photometryKron
import lsst.meas.extensions.psfex.prefs
import lsst.shapelet.shapeletFunction
import lsst.meas.modelfit.cmodel.cmodelContinued
import lsst.meas.extensions.shapeHSM.hsmShapeControl
import lsst.meas.extensions.shapeHSM
import lsst.shapelet
import lsst.meas.modelfit.truncatedGaussian
import lsst.meas.extensions.psfex.version
import lsst.shapelet.shapeletFunction.shapeletFunctionContinued
import lsst.meas.modelfit.psf.psfContinued
import lsst.shapelet.version
import lsst.meas.extensions.shapeHSM.version
import lsst.shapelet.multiShapeletFunction.multiShapeletFunctionContinued
import lsst.shapelet.hermiteTransformMatrix
import lsst.obs.subaru.crosstalkYagi
import lsst.meas.modelfit.optimizer.optimizer
import lsst.obs.subaru.crosstalk
import lsst.shapelet.multiShapeletBasis
import lsst.obs.subaru.strayLight.waveletCompression
import lsst.shapelet.radialProfile.radialProfile
import lsst.shapelet.tractor
import lsst.meas.extensions.psfex.field
import lsst.afw.display.rgb
import lsst.obs.subaru.strayLight.yStrayLight
import lsst.meas.modelfit.likelihood
import lsst.meas.extensions.convolved.convolved
import lsst.meas.modelfit.version
import lsst.shapelet.shapeletFunction.shapeletFunction
import lsst.meas.modelfit.mixture
import lsst.shapelet.matrixBuilder
import lsst.meas.modelfit.common
import lsst.meas.modelfit.pixelFitRegion.pixelFitRegion
import lsst.meas.extensions.photometryKron.version
import lsst.shapelet.multiShapeletFunction.multiShapeletFunction
import lsst.meas.modelfit.sampler
import lsst.shapelet.constants.constantsContinued
import lsst.shapelet.constants
import lsst.shapelet.basisEvaluator
import lsst.shapelet.radialProfile
import lsst.meas.modelfit.unitSystem
import lsst.meas.modelfit.multiModel
import lsst.shapelet.radialProfile.radialProfileContinued
import lsst.meas.extensions.psfex
import lsst.obs.subaru.isr
import lsst.meas.modelfit.integrals
import lsst.meas.extensions.shapeHSM.hsmMomentsControl
import lsst.meas.modelfit.pixelFitRegion
import lsst.shapelet.constants.constants
import lsst.meas.extensions.psfex.psfexPsfDeterminer
import lsst.meas.extensions.convolved
import lsst.meas.extensions.convolved.version
import lsst.shapelet.generator
import lsst.shapelet.functorKeys
import lsst.meas.modelfit.pixelFitRegion.pixelFitRegionContinued
import lsst.obs.subaru.strayLight.rotatorAngle
import lsst.meas.modelfit.psf.psf
import lsst.meas.modelfit.optimizer.optimizerContinued
import lsst.obs.hsc.vignette
import lsst.meas.extensions.psfex.psfexPsf
import lsst.meas.modelfit.cmodel
import lsst.meas.modelfit.model
import lsst.meas.modelfit.priors.priors
import lsst.shapelet.gaussHermiteConvolution
import lsst.meas.modelfit
import lsst.shapelet.multiShapeletFunction
import lsst.meas.modelfit.cmodel.cmodel
import lsst.meas.modelfit.psf
import lsst.meas.modelfit.priors.priorsContinued
import lsst.meas.extensions.photometryKron.photometryKron
import lsst.meas.modelfit.optimizer
import lsst.meas.modelfit.adaptiveImportanceSampler
import lsst.meas.modelfit.unitTransformedLikelihood
import lsst.obs.subaru.strayLight
import lsst.meas.extensions.psfex.psf
import lsst.meas.modelfit.priors
import lsst.obs.subaru.isr
config.processCcd.isr.retarget(target=lsst.obs.subaru.isr.SubaruIsrTask, ConfigClass=lsst.obs.subaru.isr.SubaruIsrConfig)

# Apply bias frame correction?
config.processCcd.isr.doBias=True

# Apply dark frame correction?
config.processCcd.isr.doDark=True

# Apply flat field correction?
config.processCcd.isr.doFlat=True

# Apply fringe correction?
config.processCcd.isr.doFringe=True

# Mask defect pixels?
config.processCcd.isr.doDefect=True

# Apply a distortion model based on camera geometry to the WCS?
config.processCcd.isr.doAddDistortionModel=True

# Persist postISRCCD?
config.processCcd.isr.doWrite=False

# Name of the bias data product
config.processCcd.isr.biasDataProductName='bias'

# Name of the dark data product
config.processCcd.isr.darkDataProductName='dark'

# Name of the flat data product
config.processCcd.isr.flatDataProductName='flat'

# trim out non-data regions?
config.processCcd.isr.assembleCcd.doTrim=True

# FITS headers to remove (in addition to DATASEC, BIASSEC, TRIMSEC and perhaps GAIN)
config.processCcd.isr.assembleCcd.keysToRemove=[]

# The gain to use if no Detector is present in the Exposure (ignored if NaN)
config.processCcd.isr.gain=float('nan')

# The read noise to use if no Detector is present in the Exposure
config.processCcd.isr.readNoise=0.0

# The saturation level to use if no Detector is present in the Exposure (ignored if NaN)
config.processCcd.isr.saturation=float('nan')

# Do fringe subtraction after flat-fielding?
config.processCcd.isr.fringeAfterFlat=True

# Only fringe-subtract these filters
config.processCcd.isr.fringe.filters=['y', 'N921', 'N926', 'N973', 'N1010']

# Number of fringe measurements
config.processCcd.isr.fringe.num=30000

# Half-size of small (fringe) measurements (pixels)
config.processCcd.isr.fringe.small=3

# Half-size of large (background) measurements (pixels)
config.processCcd.isr.fringe.large=30

# Number of fitting iterations
config.processCcd.isr.fringe.iterations=20

# Sigma clip threshold
config.processCcd.isr.fringe.clip=3.0

# Ignore pixels with these masks
config.processCcd.isr.fringe.stats.badMaskPlanes=['SAT', 'NO_DATA']

# Statistic to use
config.processCcd.isr.fringe.stats.stat=32

# Sigma clip threshold
config.processCcd.isr.fringe.stats.clip=3.0

# Number of fitting iterations
config.processCcd.isr.fringe.stats.iterations=3

# Offset to the random number generator seed (full seed includes exposure ID)
config.processCcd.isr.fringe.stats.rngSeedOffset=0

# Remove fringe pedestal?
config.processCcd.isr.fringe.pedestal=False

# FWHM of PSF (arcsec)
config.processCcd.isr.fwhm=1.0

# Name of mask plane to use in saturation detection and interpolation
config.processCcd.isr.saturatedMaskName='SAT'

# Name of mask plane to use for suspect pixels
config.processCcd.isr.suspectMaskName='SUSPECT'

# The method for scaling the flat on the fly.
config.processCcd.isr.flatScalingType='USER'

# If flatScalingType is 'USER' then scale flat by this amount; ignored otherwise
config.processCcd.isr.flatUserScale=1.0

# The method for fitting the overscan bias level.
config.processCcd.isr.overscanFitType='AKIMA_SPLINE'

# Order of polynomial or to fit if overscan fit type is a polynomial, or number of spline knots if overscan fit type is a spline.
config.processCcd.isr.overscanOrder=30

# Rejection threshold (sigma) for collapsing overscan before fit
config.processCcd.isr.overscanRej=3.0

# Number of pixels by which to grow the saturation footprints
config.processCcd.isr.growSaturationFootprintSize=1

# Perform interpolation over pixels masked as saturated?
config.processCcd.isr.doSaturationInterpolation=True

# If True, ensure we interpolate NaNs after flat-fielding, even if we also have to interpolate them before flat-fielding.
config.processCcd.isr.doNanInterpAfterFlat=False

# The approximate flux of a zero-magnitude object in a one-second exposure, per filter
config.processCcd.isr.fluxMag0T1={'g': 398107170553.49854, 'r': 398107170553.49854, 'i': 275422870333.81744, 'z': 120226443461.74132, 'y': 91201083935.59116, 'N515': 20892961308.54041, 'N816': 15848931924.611174, 'N921': 19054607179.632523}

# fields to remove from the metadata of the assembled ccd.
config.processCcd.isr.keysToRemoveFromAssembledCcd=[]

# Assemble amp-level calibration exposures into ccd-level exposure?
config.processCcd.isr.doAssembleIsrExposures=False

# Assemble amp-level exposures into a ccd-level exposure?
config.processCcd.isr.doAssembleCcd=True

# Expect input science images to have a WCS (set False for e.g. spectrographs)
config.processCcd.isr.expectWcs=True

# Correct for nonlinearity of the detector's response?
config.processCcd.isr.doLinearize=True

# Correct for crosstalk
config.processCcd.isr.doCrosstalk=True

import lsst.obs.subaru.crosstalk
config.processCcd.isr.crosstalk.retarget(target=lsst.obs.subaru.crosstalk.CrosstalkTask, ConfigClass=lsst.obs.subaru.crosstalk.CrosstalkConfig)

# Set crosstalk mask plane for pixels over this value
config.processCcd.isr.crosstalk.minPixelToMask=45000.0

# Name for crosstalk mask plane
config.processCcd.isr.crosstalk.crosstalkMaskPlane='CROSSTALK'

# Crosstalk coefficients
config.processCcd.isr.crosstalk.coeffs.values=[0.0, -0.000125, -0.000149, -0.000156, -0.000124, 0.0, -0.000132, -0.000157, -0.000171, -0.000134, 0.0, -0.000153, -0.000157, -0.000151, -0.000137, 0.0]

# Shape of coeffs array
config.processCcd.isr.crosstalk.coeffs.shape=[4, 4]

# Apply the brighter fatter correction
config.processCcd.isr.doBrighterFatter=True

# Kernel file used for the brighter fatter correction
config.processCcd.isr.brighterFatterKernelFile=''

# Maximum number of iterations for the brighter fatter correction
config.processCcd.isr.brighterFatterMaxIter=10

# Threshold used to stop iterating the brighter fatter correction.  It is the  absolute value of the difference between the current corrected image and the one from the previous iteration summed over all the pixels.
config.processCcd.isr.brighterFatterThreshold=1000.0

# Should the gain be applied when applying the brighter fatter correction?
config.processCcd.isr.brighterFatterApplyGain=True

# Dataset type for input data; users will typically leave this alone, but camera-specific ISR tasks will override it
config.processCcd.isr.datasetType='raw'

# Fallback default filter name for calibrations
config.processCcd.isr.fallbackFilterName=None

# Construct and attach a wavelength-dependent throughput curve for this CCD image?
config.processCcd.isr.doAttachTransmissionCurve=False

# Load and use transmission_optics (if doAttachTransmissionCurve is True)?
config.processCcd.isr.doUseOpticsTransmission=True

# Load and use transmission_filter (if doAttachTransmissionCurve is True)?
config.processCcd.isr.doUseFilterTransmission=True

# Load and use transmission_sensor (if doAttachTransmissionCurve is True)?
config.processCcd.isr.doUseSensorTransmission=True

# Load and use transmission_atmosphere (if doAttachTransmissionCurve is True)?
config.processCcd.isr.doUseAtmosphereTransmission=True

# Calculate empirical read noise instead of value from AmpInfo data?
config.processCcd.isr.doEmpiricalReadNoise=False

# Mesh size in X (pix) to calculate count statistics
config.processCcd.isr.qa.flatness.meshX=256

# Mesh size in Y (pix) to calculate count statistics
config.processCcd.isr.qa.flatness.meshY=256

# Do we clip outliers in calculate count statistics?
config.processCcd.isr.qa.flatness.doClip=True

# How many sigma is used to clip outliers in calculate count statistics?
config.processCcd.isr.qa.flatness.clipSigma=3.0

# How many times do we iterate clipping outliers in calculate count statistics?
config.processCcd.isr.qa.flatness.nIter=3

# Write OverScan-Subtracted image?
config.processCcd.isr.qa.doWriteOss=False

# Write OverScan-Subtracted thumbnail?
config.processCcd.isr.qa.doThumbnailOss=True

# Write flattened image?
config.processCcd.isr.qa.doWriteFlattened=False

# Write flattened thumbnail?
config.processCcd.isr.qa.doThumbnailFlattened=True

# Mask saturated pixels?
config.processCcd.isr.doSaturation=True

# Mask suspect pixels?
config.processCcd.isr.doSuspect=True

# Widen bleed trails based on their width?
config.processCcd.isr.doWidenSaturationTrails=True

# Do overscan subtraction?
config.processCcd.isr.doOverscan=True

# Calculate variance?
config.processCcd.isr.doVariance=True

# Trim guider shadow
config.processCcd.isr.doGuider=False

# Subtract stray light in the y-band (due to encoder LEDs)?
config.processCcd.isr.doStrayLight=False

# Calculate rotator start+end instead of using only the start value from the header?
config.processCcd.isr.strayLight.doRotatorAngleCorrection=True

# Correct the amplifiers for their gains
# 
# N.b. this is intended to be used *instead* of doFlat; it's useful if you're measuring system throughput
# 
config.processCcd.isr.doApplyGains=False

# Normalize all the amplifiers in each CCD to have the same gain
# 
# This does not measure the gains, it simply forces the median of each amplifier to be equal
# after applying the nominal gain
# 
config.processCcd.isr.normalizeGains=False

# Remove any PC cards in the header
config.processCcd.isr.removePcCards=True

# FWHM of PSF used when interpolating over bad columns (arcsec)
config.processCcd.isr.fwhmForBadColumnInterpolation=1.0

# Should we set the level of all BAD patches of the chip to the chip's average value?
config.processCcd.isr.doSetBadRegions=True

# How to estimate the average value for BAD regions.
config.processCcd.isr.badStatistic='MEANCLIP'

# Tweak flats to match observed amplifier ratios?
config.processCcd.isr.doTweakFlat=False

# Maximum deviation from the median for overscan
config.processCcd.isr.overscanMaxDev=1000.0

# Center of vignetting pattern, in x (focal plane coords)
config.processCcd.isr.vignette.xCenter=-100.0

# Center of vignetting pattern, in y (focal plane coords)
config.processCcd.isr.vignette.yCenter=100.0

# Radius of vignetting pattern, in focal plane coords
config.processCcd.isr.vignette.radius=17500.0

# Number of points to define the Vignette polygon
config.processCcd.isr.numPolygonPoints=100

# Persist Polygon used to define vignetted region?
config.processCcd.isr.doWriteVignettePolygon=True

# Default value for fluxMag0T1 (for an unrecognised filter)
config.processCcd.isr.defaultFluxMag0T1=158489319246.11172

# Binning factor for thumbnail
config.processCcd.isr.thumbnailBinning=4

# Number of stdev below the background to set thumbnail minimum
config.processCcd.isr.thumbnailStdev=3.0

# Range for thumbnail mapping
config.processCcd.isr.thumbnailRange=5.0

# Softening parameter for thumbnail mapping
config.processCcd.isr.thumbnailQ=20.0

# Border around saturated pixels for thumbnail
config.processCcd.isr.thumbnailSatBorder=2

import lsst.shapelet.gaussHermiteProjection
import lsst.shapelet.shapeletFunction
import lsst.meas.modelfit.cmodel.cmodelContinued
import lsst.shapelet.generator
import lsst.shapelet.functorKeys
import lsst.shapelet
import lsst.meas.modelfit.pixelFitRegion.pixelFitRegionContinued
import lsst.meas.modelfit.truncatedGaussian
import lsst.shapelet.shapeletFunction.shapeletFunctionContinued
import lsst.meas.modelfit.psf.psf
import lsst.shapelet.version
import lsst.meas.modelfit.optimizer.optimizerContinued
import lsst.meas.modelfit.psf.psfContinued
import lsst.shapelet.multiShapeletFunction.multiShapeletFunctionContinued
import lsst.shapelet.hermiteTransformMatrix
import lsst.meas.modelfit.optimizer.optimizer
import lsst.shapelet.multiShapeletBasis
import lsst.meas.modelfit.cmodel
import lsst.meas.modelfit.model
import lsst.meas.modelfit.priors.priors
import lsst.shapelet.gaussHermiteConvolution
import lsst.shapelet.radialProfile.radialProfile
import lsst.shapelet.tractor
import lsst.meas.modelfit
import lsst.shapelet.multiShapeletFunction
import lsst.meas.modelfit.cmodel.cmodel
import lsst.meas.modelfit.psf
import lsst.meas.modelfit.likelihood
import lsst.meas.modelfit.priors.priorsContinued
import lsst.meas.modelfit.version
import lsst.shapelet.shapeletFunction.shapeletFunction
import lsst.meas.modelfit.mixture
import lsst.shapelet.matrixBuilder
import lsst.meas.modelfit.common
import lsst.meas.modelfit.optimizer
import lsst.meas.modelfit.adaptiveImportanceSampler
import lsst.meas.modelfit.pixelFitRegion.pixelFitRegion
import lsst.shapelet.multiShapeletFunction.multiShapeletFunction
import lsst.meas.modelfit.unitTransformedLikelihood
import lsst.meas.modelfit.sampler
import lsst.shapelet.basisEvaluator
import lsst.shapelet.constants
import lsst.shapelet.constants.constantsContinued
import lsst.shapelet.radialProfile
import lsst.meas.modelfit.unitSystem
import lsst.meas.modelfit.multiModel
import lsst.shapelet.radialProfile.radialProfileContinued
import lsst.meas.modelfit.integrals
import lsst.meas.modelfit.pixelFitRegion
import lsst.meas.modelfit.priors
import lsst.shapelet.constants.constants
# Measure PSF? If False then for all subsequent operations use either existing PSF model when present, or install simple PSF model when not (see installSimplePsf config options)
config.processCcd.charImage.doMeasurePsf=True

# Persist results?
config.processCcd.charImage.doWrite=True

# Write icExp and icExpBackground in addition to icSrc? Ignored if doWrite False.
config.processCcd.charImage.doWriteExposure=False

# Number of iterations of detect sources, measure sources, estimate PSF. If useSimplePsf is True then 2 should be plenty; otherwise more may be wanted.
config.processCcd.charImage.psfIterations=2

# type of statistic to use for grid points
config.processCcd.charImage.background.statisticsProperty='MEANCLIP'

# behaviour if there are too few points in grid for requested interpolation style
config.processCcd.charImage.background.undersampleStyle='REDUCE_INTERP_ORDER'

# how large a region of the sky should be used for each background point
config.processCcd.charImage.background.binSize=128

# Sky region size to be used for each background point in X direction. If 0, the binSize config is used.
config.processCcd.charImage.background.binSizeX=0

# Sky region size to be used for each background point in Y direction. If 0, the binSize config is used.
config.processCcd.charImage.background.binSizeY=0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
config.processCcd.charImage.background.algorithm='AKIMA_SPLINE'

# Names of mask planes to ignore while estimating the background
config.processCcd.charImage.background.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# Ignore NaNs when estimating the background
config.processCcd.charImage.background.isNanSafe=False

# Use Approximate (Chebyshev) to model background.
config.processCcd.charImage.background.useApprox=True

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.processCcd.charImage.background.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.processCcd.charImage.background.approxOrderY=-1

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.processCcd.charImage.background.weighting=True

# detected sources with fewer than the specified number of pixels will be ignored
config.processCcd.charImage.detection.minPixels=1

# Pixels should be grown as isotropically as possible (slower)
config.processCcd.charImage.detection.isotropicGrow=True

# Grow all footprints at the same time? This allows disconnected footprints to merge.
config.processCcd.charImage.detection.combinedGrow=True

# Grow detections by nSigmaToGrow * [PSF RMS width]; if 0 then do not grow
config.processCcd.charImage.detection.nSigmaToGrow=2.4

# Grow detections to set the image mask bits, but return the original (not-grown) footprints
config.processCcd.charImage.detection.returnOriginalFootprints=False

# Threshold for footprints
config.processCcd.charImage.detection.thresholdValue=5.0

# Include threshold relative to thresholdValue
config.processCcd.charImage.detection.includeThresholdMultiplier=10.0

# specifies the desired flavor of Threshold
config.processCcd.charImage.detection.thresholdType='stdev'

# specifies whether to detect positive, or negative sources, or both
config.processCcd.charImage.detection.thresholdPolarity='positive'

# Fiddle factor to add to the background; debugging only
config.processCcd.charImage.detection.adjustBackground=0.0

# Estimate the background again after final source detection?
config.processCcd.charImage.detection.reEstimateBackground=True

# type of statistic to use for grid points
config.processCcd.charImage.detection.background.statisticsProperty='MEANCLIP'

# behaviour if there are too few points in grid for requested interpolation style
config.processCcd.charImage.detection.background.undersampleStyle='REDUCE_INTERP_ORDER'

# how large a region of the sky should be used for each background point
config.processCcd.charImage.detection.background.binSize=128

# Sky region size to be used for each background point in X direction. If 0, the binSize config is used.
config.processCcd.charImage.detection.background.binSizeX=0

# Sky region size to be used for each background point in Y direction. If 0, the binSize config is used.
config.processCcd.charImage.detection.background.binSizeY=0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
config.processCcd.charImage.detection.background.algorithm='AKIMA_SPLINE'

# Names of mask planes to ignore while estimating the background
config.processCcd.charImage.detection.background.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# Ignore NaNs when estimating the background
config.processCcd.charImage.detection.background.isNanSafe=False

# Use Approximate (Chebyshev) to model background.
config.processCcd.charImage.detection.background.useApprox=True

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.processCcd.charImage.detection.background.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.processCcd.charImage.detection.background.approxOrderY=-1

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.processCcd.charImage.detection.background.weighting=True

# type of statistic to use for grid points
config.processCcd.charImage.detection.tempLocalBackground.statisticsProperty='MEANCLIP'

# behaviour if there are too few points in grid for requested interpolation style
config.processCcd.charImage.detection.tempLocalBackground.undersampleStyle='REDUCE_INTERP_ORDER'

# how large a region of the sky should be used for each background point
config.processCcd.charImage.detection.tempLocalBackground.binSize=64

# Sky region size to be used for each background point in X direction. If 0, the binSize config is used.
config.processCcd.charImage.detection.tempLocalBackground.binSizeX=0

# Sky region size to be used for each background point in Y direction. If 0, the binSize config is used.
config.processCcd.charImage.detection.tempLocalBackground.binSizeY=0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
config.processCcd.charImage.detection.tempLocalBackground.algorithm='AKIMA_SPLINE'

# Names of mask planes to ignore while estimating the background
config.processCcd.charImage.detection.tempLocalBackground.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# Ignore NaNs when estimating the background
config.processCcd.charImage.detection.tempLocalBackground.isNanSafe=False

# Use Approximate (Chebyshev) to model background.
config.processCcd.charImage.detection.tempLocalBackground.useApprox=False

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.processCcd.charImage.detection.tempLocalBackground.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.processCcd.charImage.detection.tempLocalBackground.approxOrderY=-1

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.processCcd.charImage.detection.tempLocalBackground.weighting=True

# Enable temporary local background subtraction? (see tempLocalBackground)
config.processCcd.charImage.detection.doTempLocalBackground=False

# type of statistic to use for grid points
config.processCcd.charImage.detection.tempWideBackground.statisticsProperty='MEANCLIP'

# behaviour if there are too few points in grid for requested interpolation style
config.processCcd.charImage.detection.tempWideBackground.undersampleStyle='REDUCE_INTERP_ORDER'

# how large a region of the sky should be used for each background point
config.processCcd.charImage.detection.tempWideBackground.binSize=512

# Sky region size to be used for each background point in X direction. If 0, the binSize config is used.
config.processCcd.charImage.detection.tempWideBackground.binSizeX=0

# Sky region size to be used for each background point in Y direction. If 0, the binSize config is used.
config.processCcd.charImage.detection.tempWideBackground.binSizeY=0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
config.processCcd.charImage.detection.tempWideBackground.algorithm='AKIMA_SPLINE'

# Names of mask planes to ignore while estimating the background
config.processCcd.charImage.detection.tempWideBackground.ignoredPixelMask=['BAD', 'EDGE', 'NO_DATA']

# Ignore NaNs when estimating the background
config.processCcd.charImage.detection.tempWideBackground.isNanSafe=False

# Use Approximate (Chebyshev) to model background.
config.processCcd.charImage.detection.tempWideBackground.useApprox=False

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.processCcd.charImage.detection.tempWideBackground.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.processCcd.charImage.detection.tempWideBackground.approxOrderY=-1

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.processCcd.charImage.detection.tempWideBackground.weighting=True

# Do temporary wide (large-scale) background subtraction before footprint detection?
config.processCcd.charImage.detection.doTempWideBackground=False

# The maximum number of peaks in a Footprint before trying to replace its peaks using the temporary local background
config.processCcd.charImage.detection.nPeaksMaxSimple=1

# Multiple of PSF RMS size to use for convolution kernel bounding box size; note that this is not a half-size. The size will be rounded up to the nearest odd integer
config.processCcd.charImage.detection.nSigmaForKernel=7.0

# Mask planes to ignore when calculating statistics of image (for thresholdType=stdev)
config.processCcd.charImage.detection.statsMask=['BAD', 'SAT', 'EDGE', 'NO_DATA']

# Run deblender input exposure
config.processCcd.charImage.doDeblend=False

# What to do when a peak to be deblended is close to the edge of the image
config.processCcd.charImage.deblend.edgeHandling='ramp'

# When the deblender should attribute stray flux to point sources
config.processCcd.charImage.deblend.strayFluxToPointSources='necessary'

# Assign stray flux (not claimed by any child in the deblender) to deblend children.
config.processCcd.charImage.deblend.assignStrayFlux=True

# How to split flux among peaks
config.processCcd.charImage.deblend.strayFluxRule='trim'

# When splitting stray flux, clip fractions below this value to zero.
config.processCcd.charImage.deblend.clipStrayFluxFraction=0.001

# Chi-squared per DOF cut for deciding a source is a PSF during deblending (un-shifted PSF model)
config.processCcd.charImage.deblend.psfChisq1=1.5

# Chi-squared per DOF cut for deciding a source is PSF during deblending (shifted PSF model)
config.processCcd.charImage.deblend.psfChisq2=1.5

# Chi-squared per DOF cut for deciding a source is a PSF during deblending (shifted PSF model #2)
config.processCcd.charImage.deblend.psfChisq2b=1.5

# Only deblend the brightest maxNumberOfPeaks peaks in the parent (<= 0: unlimited)
config.processCcd.charImage.deblend.maxNumberOfPeaks=0

# Maximum area for footprints before they are ignored as large; non-positive means no threshold applied
config.processCcd.charImage.deblend.maxFootprintArea=10000

# Maximum linear dimension for footprints before they are ignored as large; non-positive means no threshold applied
config.processCcd.charImage.deblend.maxFootprintSize=0

# Minimum axis ratio for footprints before they are ignored as large; non-positive means no threshold applied
config.processCcd.charImage.deblend.minFootprintAxisRatio=0.0

# Mask name for footprints not deblended, or None
config.processCcd.charImage.deblend.notDeblendedMask='NOT_DEBLENDED'

# Footprints smaller in width or height than this value will be ignored; minimum of 2 due to PSF gradient calculation.
config.processCcd.charImage.deblend.tinyFootprintSize=2

# Guarantee that all peaks produce a child source.
config.processCcd.charImage.deblend.propagateAllPeaks=False

# If True, catch exceptions thrown by the deblender, log them, and set a flag on the parent, instead of letting them propagate up
config.processCcd.charImage.deblend.catchFailures=False

# Mask planes to ignore when performing statistics
config.processCcd.charImage.deblend.maskPlanes=['SAT', 'INTRP', 'NO_DATA']

# Mask planes with the corresponding limit on the fraction of masked pixels. Sources violating this limit will not be deblended.
config.processCcd.charImage.deblend.maskLimits={'NO_DATA': 0.25}

# If true, a least-squares fit of the templates will be done to the full image. The templates will be re-weighted based on this fit.
config.processCcd.charImage.deblend.weightTemplates=False

# Try to remove similar templates?
config.processCcd.charImage.deblend.removeDegenerateTemplates=False

# If the dot product between two templates is larger than this value, we consider them to be describing the same object (i.e. they are degenerate).  If one of the objects has been labeled as a PSF it will be removed, otherwise the template with the lowest value will be removed.
config.processCcd.charImage.deblend.maxTempDotProd=0.5

# Apply a smoothing filter to all of the template images
config.processCcd.charImage.deblend.medianSmoothTemplate=True

import lsst.meas.extensions.convolved
import lsst.meas.extensions.photometryKron
import lsst.meas.extensions.shapeHSM.version
import lsst.meas.extensions.convolved.version
import lsst.meas.extensions.shapeHSM.hsmMomentsControl
import lsst.meas.extensions.convolved.convolved
import lsst.meas.extensions.photometryKron.photometryKron
import lsst.meas.extensions.shapeHSM.hsmShapeControl
import lsst.meas.extensions.shapeHSM
import lsst.meas.extensions.photometryKron.version
# the name of the centroiding algorithm used to set source x,y
config.processCcd.charImage.measurement.slots.centroid='base_SdssCentroid'

# the name of the algorithm used to set source moments parameters
config.processCcd.charImage.measurement.slots.shape='ext_shapeHSM_HsmSourceMoments'

# the name of the algorithm used to set PSF moments parameters
config.processCcd.charImage.measurement.slots.psfShape='ext_shapeHSM_HsmPsfMoments'

# the name of the algorithm used to set the source aperture flux slot
config.processCcd.charImage.measurement.slots.apFlux='base_CircularApertureFlux_12_0'

# the name of the algorithm used to set the source model flux slot
config.processCcd.charImage.measurement.slots.modelFlux='modelfit_CModel'

# the name of the algorithm used to set the source psf flux slot
config.processCcd.charImage.measurement.slots.psfFlux='base_PsfFlux'

# the name of the algorithm used to set the source inst flux slot
config.processCcd.charImage.measurement.slots.instFlux='base_GaussianFlux'

# the name of the flux measurement algorithm used for calibration
config.processCcd.charImage.measurement.slots.calibFlux='base_CircularApertureFlux_12_0'

# When measuring, replace other detected footprints with noise?
config.processCcd.charImage.measurement.doReplaceWithNoise=True

# How to choose mean and variance of the Gaussian noise we generate?
config.processCcd.charImage.measurement.noiseReplacer.noiseSource='measure'

# Add ann offset to the generated noise.
config.processCcd.charImage.measurement.noiseReplacer.noiseOffset=0.0

# The seed multiplier value to use for random number generation
#    >= 1: set the seed deterministically based on exposureId
#       0: fall back to the afw.math.Random default constructor (which uses a seed value of 1)
config.processCcd.charImage.measurement.noiseReplacer.noiseSeedMultiplier=1

# Prefix to give undeblended plugins
config.processCcd.charImage.measurement.undeblendedPrefix='undeblended_'

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_PsfFlux'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.processCcd.charImage.measurement.plugins['base_PsfFlux'].badMaskPlanes=[]

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_PeakLikelihoodFlux'].doMeasure=True

# Name of warping kernel (e.g. "lanczos4") used to compute the peak
config.processCcd.charImage.measurement.plugins['base_PeakLikelihoodFlux'].warpingKernelName='lanczos4'

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_GaussianFlux'].doMeasure=True

# FIXME! NEVER DOCUMENTED!
config.processCcd.charImage.measurement.plugins['base_GaussianFlux'].background=0.0

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_NaiveCentroid'].doMeasure=True

# Value to subtract from the image pixel values
config.processCcd.charImage.measurement.plugins['base_NaiveCentroid'].background=0.0

# Do check that the centroid is contained in footprint.
config.processCcd.charImage.measurement.plugins['base_NaiveCentroid'].doFootprintCheck=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.processCcd.charImage.measurement.plugins['base_NaiveCentroid'].maxDistToPeak=-1.0

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_SdssCentroid'].doMeasure=True

# maximum allowed binning
config.processCcd.charImage.measurement.plugins['base_SdssCentroid'].binmax=16

# Do check that the centroid is contained in footprint.
config.processCcd.charImage.measurement.plugins['base_SdssCentroid'].doFootprintCheck=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.processCcd.charImage.measurement.plugins['base_SdssCentroid'].maxDistToPeak=-1.0

# if the peak's less than this insist on binning at least once
config.processCcd.charImage.measurement.plugins['base_SdssCentroid'].peakMin=-1.0

# fiddle factor for adjusting the binning
config.processCcd.charImage.measurement.plugins['base_SdssCentroid'].wfac=1.5

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_PixelFlags'].doMeasure=True

# List of mask planes to be searched for which occur anywhere within a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.processCcd.charImage.measurement.plugins['base_PixelFlags'].masksFpAnywhere=[]

# List of mask planes to be searched for which occur in the center of a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.processCcd.charImage.measurement.plugins['base_PixelFlags'].masksFpCenter=[]

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_SdssShape'].doMeasure=True

# Additional value to add to background
config.processCcd.charImage.measurement.plugins['base_SdssShape'].background=0.0

# Whether to also compute the shape of the PSF model
config.processCcd.charImage.measurement.plugins['base_SdssShape'].doMeasurePsf=True

# Maximum number of iterations
config.processCcd.charImage.measurement.plugins['base_SdssShape'].maxIter=100

# Maximum centroid shift, limited to 2-10
config.processCcd.charImage.measurement.plugins['base_SdssShape'].maxShift=0.0

# Convergence tolerance for e1,e2
config.processCcd.charImage.measurement.plugins['base_SdssShape'].tol1=9.999999747378752e-06

# Convergence tolerance for FWHM
config.processCcd.charImage.measurement.plugins['base_SdssShape'].tol2=9.999999747378752e-05

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_ScaledApertureFlux'].doMeasure=True

# Scaling factor of PSF FWHM for aperture radius.
config.processCcd.charImage.measurement.plugins['base_ScaledApertureFlux'].scale=3.14

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.processCcd.charImage.measurement.plugins['base_ScaledApertureFlux'].shiftKernel='lanczos5'

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_CircularApertureFlux'].doMeasure=True

# Maximum radius (in pixels) for which the sinc algorithm should be used instead of the faster naive algorithm.  For elliptical apertures, this is the minor axis radius.
config.processCcd.charImage.measurement.plugins['base_CircularApertureFlux'].maxSincRadius=12.0

# Radius (in pixels) of apertures.
config.processCcd.charImage.measurement.plugins['base_CircularApertureFlux'].radii=[3.0, 4.5, 6.0, 9.0, 12.0, 17.0, 25.0, 35.0, 50.0, 70.0]

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.processCcd.charImage.measurement.plugins['base_CircularApertureFlux'].shiftKernel='lanczos5'

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_Blendedness'].doMeasure=True

# Whether to compute quantities related to the Gaussian-weighted flux
config.processCcd.charImage.measurement.plugins['base_Blendedness'].doFlux=True

# Whether to compute HeavyFootprint dot products (the old deblend.blendedness parameter)
config.processCcd.charImage.measurement.plugins['base_Blendedness'].doOld=True

# Whether to compute quantities related to the Gaussian-weighted shape
config.processCcd.charImage.measurement.plugins['base_Blendedness'].doShape=True

# Radius factor that sets the maximum extent of the weight function (and hence the flux measurements)
config.processCcd.charImage.measurement.plugins['base_Blendedness'].nSigmaWeightMax=3.0

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_LocalBackground'].doMeasure=True

# Inner radius for background annulus as a multiple of the PSF sigma
config.processCcd.charImage.measurement.plugins['base_LocalBackground'].annulusInner=7.0

# Outer radius for background annulus as a multiple of the PSF sigma
config.processCcd.charImage.measurement.plugins['base_LocalBackground'].annulusOuter=15.0

# Mask planes that indicate pixels that should be excluded from the measurement
config.processCcd.charImage.measurement.plugins['base_LocalBackground'].badMaskPlanes=['BAD', 'SAT', 'NO_DATA']

# Number of sigma-clipping iterations for background measurement
config.processCcd.charImage.measurement.plugins['base_LocalBackground'].bgIter=3

# Rejection threshold (in standard deviations) for background measurement
config.processCcd.charImage.measurement.plugins['base_LocalBackground'].bgRej=3.0

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_FPPosition'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_Jacobian'].doMeasure=True

# Nominal pixel size (arcsec)
config.processCcd.charImage.measurement.plugins['base_Jacobian'].pixelScale=0.168

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_Variance'].doMeasure=True

# Scale factor to apply to shape for aperture
config.processCcd.charImage.measurement.plugins['base_Variance'].scale=5.0

# Mask planes to ignore
config.processCcd.charImage.measurement.plugins['base_Variance'].mask=['DETECTED', 'DETECTED_NEGATIVE', 'BAD', 'SAT']

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_InputCount'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_PeakCentroid'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_SkyCoord'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['subaru_FilterFraction'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].doMeasure=True

# Shapelet order of inner expansion (0 == Gaussian)
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].innerOrder=2

# Don't allow the semi-major radius of any component to go above this fraction of the PSF image width
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].maxRadiusBoxFraction=0.4

# Don't allow the semi-minor radius of any component to drop below this value (pixels)
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].minRadius=1.0

# Don't allow the determinant radii of the two components to differ by less than this (pixels)
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].minRadiusDiff=0.5

# whether to save all iterations for debugging purposes
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionSolverTolerance=1e-08

# Shapelet order of outer expansion (0 == Gaussian)
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].outerOrder=1

# Initial outer Gaussian peak height divided by inner Gaussian peak height
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].peakRatio=0.1

# Initial outer radius divided by inner radius
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].radiusRatio=2.0

config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models={}
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.radiusPriorSigma=0.5

config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.radiusPriorSigma=0.5

config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.order=2

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.order=1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.radiusPriorSigma=0.5

config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.order=4

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.order=4

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.radiusPriorSigma=0.5

# a sequence of model names indicating which models should be fit, and their order
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].sequence=['DoubleShapelet']

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].doMeasure=True

# Whether to record the steps the optimizer takes (or just the number, if running as a plugin)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.doRecordHistory=True

# Whether to record the time spent in this stage
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.doRecordTime=True

# Softened core width for ellipticity distribution (conformal shear units).
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.empiricalPriorConfig.ellipticityCore=0.001

# Width of exponential ellipticity distribution (conformal shear units).
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.empiricalPriorConfig.ellipticitySigma=0.3

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusMinInner=-6.0

# Minimum ln(radius).
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusMinOuter=-6.001

# Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T.
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusMu=-1.0

# Number of degrees of freedom for the Student's T distribution on ln(radius).
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusNu=50.0

# Width of the Student's T distribution in ln(radius).
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusSigma=0.45

# Ellipticity magnitude (conformal shear units) at which the softened cutoff begins
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.linearPriorConfig.ellipticityMaxInner=2.0

# Maximum ellipticity magnitude (conformal shear units)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.linearPriorConfig.ellipticityMaxOuter=2.001

# ln(radius) at which the softened cutoff begins towards the maximum
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.linearPriorConfig.logRadiusMaxInner=3.0

# Maximum ln(radius)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.linearPriorConfig.logRadiusMaxOuter=3.001

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.linearPriorConfig.logRadiusMinInner=-6.0

# The ratio P(logRadiusMinInner)/P(logRadiusMaxInner)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.linearPriorConfig.logRadiusMinMaxRatio=1.0

# Minimum ln(radius)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.linearPriorConfig.logRadiusMinOuter=-6.001

# Maximum radius used in approximating profile with Gaussians (0=default for this profile)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.maxRadius=0

# Number of Gaussian used to approximate the profile
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.nComponents=8

# whether to save all iterations for debugging purposes
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.optimizer.trustRegionSolverTolerance=1e-08

# Name of the Prior that defines the model to fit (a filename in $MEAS_MODELFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting.
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.priorName=''

# One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.priorSource='EMPIRICAL'

# Name of the shapelet.RadialProfile that defines the model to fit
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.profileName='luv'

# Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.usePixelWeights=False

# Scale the likelihood by this factor to artificially reweight it w.r.t. the prior.
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.weightsMultiplier=1.0

# Whether to record the steps the optimizer takes (or just the number, if running as a plugin)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.doRecordHistory=True

# Whether to record the time spent in this stage
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.doRecordTime=True

# Softened core width for ellipticity distribution (conformal shear units).
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.empiricalPriorConfig.ellipticityCore=0.001

# Width of exponential ellipticity distribution (conformal shear units).
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.empiricalPriorConfig.ellipticitySigma=0.3

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusMinInner=-6.0

# Minimum ln(radius).
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusMinOuter=-6.001

# Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T.
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusMu=-1.0

# Number of degrees of freedom for the Student's T distribution on ln(radius).
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusNu=50.0

# Width of the Student's T distribution in ln(radius).
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusSigma=0.45

# Ellipticity magnitude (conformal shear units) at which the softened cutoff begins
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.linearPriorConfig.ellipticityMaxInner=2.0

# Maximum ellipticity magnitude (conformal shear units)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.linearPriorConfig.ellipticityMaxOuter=2.001

# ln(radius) at which the softened cutoff begins towards the maximum
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.linearPriorConfig.logRadiusMaxInner=3.0

# Maximum ln(radius)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.linearPriorConfig.logRadiusMaxOuter=3.001

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.linearPriorConfig.logRadiusMinInner=-6.0

# The ratio P(logRadiusMinInner)/P(logRadiusMaxInner)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.linearPriorConfig.logRadiusMinMaxRatio=1.0

# Minimum ln(radius)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.linearPriorConfig.logRadiusMinOuter=-6.001

# Maximum radius used in approximating profile with Gaussians (0=default for this profile)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.maxRadius=0

# Number of Gaussian used to approximate the profile
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.nComponents=6

# whether to save all iterations for debugging purposes
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.optimizer.maxOuterIterations=250

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.optimizer.trustRegionSolverTolerance=1e-08

# Name of the Prior that defines the model to fit (a filename in $MEAS_MODELFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting.
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.priorName=''

# One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.priorSource='EMPIRICAL'

# Name of the shapelet.RadialProfile that defines the model to fit
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.profileName='lux'

# Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.usePixelWeights=False

# Scale the likelihood by this factor to artificially reweight it w.r.t. the prior.
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.weightsMultiplier=1.0

# If the 2nd-moments shape used to initialize the fit failed, use the PSF moments multiplied by this.  If <= 0.0, abort the fit early instead.
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].fallbackInitialMomentsPsfFactor=1.5

# Whether to record the steps the optimizer takes (or just the number, if running as a plugin)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.doRecordHistory=True

# Whether to record the time spent in this stage
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.doRecordTime=True

# Softened core width for ellipticity distribution (conformal shear units).
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.empiricalPriorConfig.ellipticityCore=0.001

# Width of exponential ellipticity distribution (conformal shear units).
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.empiricalPriorConfig.ellipticitySigma=0.3

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusMinInner=-6.0

# Minimum ln(radius).
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusMinOuter=-6.001

# Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T.
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusMu=-1.0

# Number of degrees of freedom for the Student's T distribution on ln(radius).
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusNu=50.0

# Width of the Student's T distribution in ln(radius).
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusSigma=0.45

# Ellipticity magnitude (conformal shear units) at which the softened cutoff begins
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.linearPriorConfig.ellipticityMaxInner=2.0

# Maximum ellipticity magnitude (conformal shear units)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.linearPriorConfig.ellipticityMaxOuter=2.001

# ln(radius) at which the softened cutoff begins towards the maximum
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.linearPriorConfig.logRadiusMaxInner=3.0

# Maximum ln(radius)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.linearPriorConfig.logRadiusMaxOuter=3.001

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.linearPriorConfig.logRadiusMinInner=-6.0

# The ratio P(logRadiusMinInner)/P(logRadiusMaxInner)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.linearPriorConfig.logRadiusMinMaxRatio=1.0

# Minimum ln(radius)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.linearPriorConfig.logRadiusMinOuter=-6.001

# Maximum radius used in approximating profile with Gaussians (0=default for this profile)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.maxRadius=0

# Number of Gaussian used to approximate the profile
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.nComponents=3

# whether to save all iterations for debugging purposes
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.optimizer.gradientThreshold=0.01

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.optimizer.minTrustRadiusThreshold=0.01

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.optimizer.trustRegionSolverTolerance=1e-08

# Name of the Prior that defines the model to fit (a filename in $MEAS_MODELFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting.
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.priorName=''

# One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.priorSource='EMPIRICAL'

# Name of the shapelet.RadialProfile that defines the model to fit
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.profileName='lux'

# Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.usePixelWeights=True

# Scale the likelihood by this factor to artificially reweight it w.r.t. the prior.
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.weightsMultiplier=1.0

# Minimum initial radius in pixels (used to regularize initial moments-based PSF deconvolution)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].minInitialRadius=0.1

# Field name prefix of the Shapelet PSF approximation used to convolve the galaxy model; must contain a set of fields matching the schema defined by shapelet.MultiShapeletFunctionKey.
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].psfName='modelfit_DoubleShapeletPsfApprox'

# Mask planes that indicate pixels that should be ignored in the fit.
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].region.badMaskPlanes=['EDGE', 'SAT', 'BAD', 'NO_DATA']

# Abort if the fit region grows beyond this many pixels.
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].region.maxArea=100000

# Maximum fraction of pixels that may be ignored due to masks; more than this and we don't even try.
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].region.maxBadPixelFraction=0.1

# Use this multiple of the initial fit ellipse then grow by the PSF width to determine the maximum final fit region size.
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].region.nFitRadiiMax=3.0

# Use this multiple of the initial fit ellipse then grow by the PSF width to determine the minimum final fit region size.
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].region.nFitRadiiMin=1.0

# Use this multiple of the Kron ellipse to set the fit region (for the final fit region, subject to the nFitRadiiMin and nFitRadiiMax constraints).
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].region.nKronRadii=1.5

# Grow the initial fit ellipses by this factor before comparing with the Kron/Footprint region
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].region.nPsfSigmaGrow=2.0

# If the Kron radius is less than this multiple of the PSF width, ignore it and fall back to a PSF-oriented ellipse scaled to match the area of the footprint or this radius (whichever is larger).
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].region.nPsfSigmaMin=4.0

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['ext_photometryKron_KronFlux'].doMeasure=True

# If true check that the Kron radius exceeds some minimum
config.processCcd.charImage.measurement.plugins['ext_photometryKron_KronFlux'].enforceMinimumRadius=True

# if true, use existing shape and centroid measurements instead of fitting
config.processCcd.charImage.measurement.plugins['ext_photometryKron_KronFlux'].fixed=False

# Largest aperture for which to use the slow, accurate, sinc aperture code
config.processCcd.charImage.measurement.plugins['ext_photometryKron_KronFlux'].maxSincRadius=10.0

# Minimum Kron radius (if == 0.0 use PSF's Kron radius) if enforceMinimumRadius. Also functions as fallback aperture radius if set.
config.processCcd.charImage.measurement.plugins['ext_photometryKron_KronFlux'].minimumRadius=0.0

# Number of times to iterate when setting the Kron radius
config.processCcd.charImage.measurement.plugins['ext_photometryKron_KronFlux'].nIterForRadius=1

# Number of Kron radii for Kron flux
config.processCcd.charImage.measurement.plugins['ext_photometryKron_KronFlux'].nRadiusForFlux=2.5

# Multiplier of rms size for aperture used to initially estimate the Kron radius
config.processCcd.charImage.measurement.plugins['ext_photometryKron_KronFlux'].nSigmaForRadius=6.0

# Name of field specifying reference Kron radius for forced measurement
config.processCcd.charImage.measurement.plugins['ext_photometryKron_KronFlux'].refRadiusName='ext_photometryKron_KronFlux_radius'

# Smooth image with N(0, smoothingSigma^2) Gaussian while estimating R_K
config.processCcd.charImage.measurement.plugins['ext_photometryKron_KronFlux'].smoothingSigma=-1.0

# Use the Footprint size as part of initial estimate of Kron radius
config.processCcd.charImage.measurement.plugins['ext_photometryKron_KronFlux'].useFootprintRadius=False

# list of target seeings (FWHM, pixels)
config.processCcd.charImage.measurement.plugins['ext_convolved_ConvolvedFlux'].seeing=[3.5, 5.0, 6.5, 8.0]

# scaling factor of kernel sigma for kernel size
config.processCcd.charImage.measurement.plugins['ext_convolved_ConvolvedFlux'].kernelScale=4.0

# Maximum radius (in pixels) for which the sinc algorithm should be used instead of the faster naive algorithm.  For elliptical apertures, this is the minor axis radius.
config.processCcd.charImage.measurement.plugins['ext_convolved_ConvolvedFlux'].aperture.maxSincRadius=10.0

# Radius (in pixels) of apertures.
config.processCcd.charImage.measurement.plugins['ext_convolved_ConvolvedFlux'].aperture.radii=[3.3, 4.5, 6.0]

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.processCcd.charImage.measurement.plugins['ext_convolved_ConvolvedFlux'].aperture.shiftKernel='lanczos5'

# name of Kron radius field in reference
config.processCcd.charImage.measurement.plugins['ext_convolved_ConvolvedFlux'].kronRadiusName='ext_photometryKron_KronFlux_radius'

# Largest aperture for which to use the sinc aperture code for Kron (pixels)
config.processCcd.charImage.measurement.plugins['ext_convolved_ConvolvedFlux'].maxSincRadius=10.0

# Number of Kron radii for Kron flux
config.processCcd.charImage.measurement.plugins['ext_convolved_ConvolvedFlux'].kronRadiusForFlux=2.5

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['ext_convolved_ConvolvedFlux'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['ext_shapeHSM_HsmShapeBj'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.processCcd.charImage.measurement.plugins['ext_shapeHSM_HsmShapeBj'].badMaskPlanes=None

# Field name for number of deblend children
config.processCcd.charImage.measurement.plugins['ext_shapeHSM_HsmShapeBj'].deblendNChild=None

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['ext_shapeHSM_HsmShapeLinear'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.processCcd.charImage.measurement.plugins['ext_shapeHSM_HsmShapeLinear'].badMaskPlanes=None

# Field name for number of deblend children
config.processCcd.charImage.measurement.plugins['ext_shapeHSM_HsmShapeLinear'].deblendNChild=None

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['ext_shapeHSM_HsmShapeKsb'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.processCcd.charImage.measurement.plugins['ext_shapeHSM_HsmShapeKsb'].badMaskPlanes=None

# Field name for number of deblend children
config.processCcd.charImage.measurement.plugins['ext_shapeHSM_HsmShapeKsb'].deblendNChild=None

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['ext_shapeHSM_HsmShapeRegauss'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.processCcd.charImage.measurement.plugins['ext_shapeHSM_HsmShapeRegauss'].badMaskPlanes=None

# Field name for number of deblend children
config.processCcd.charImage.measurement.plugins['ext_shapeHSM_HsmShapeRegauss'].deblendNChild=''

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['ext_shapeHSM_HsmSourceMoments'].doMeasure=True

# Store measured flux?
config.processCcd.charImage.measurement.plugins['ext_shapeHSM_HsmSourceMoments'].addFlux=None

# Mask planes used to reject bad pixels.
config.processCcd.charImage.measurement.plugins['ext_shapeHSM_HsmSourceMoments'].badMaskPlanes=None

# Use round weight function?
config.processCcd.charImage.measurement.plugins['ext_shapeHSM_HsmSourceMoments'].roundMoments=None

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['ext_shapeHSM_HsmSourceMomentsRound'].doMeasure=True

# Store measured flux?
config.processCcd.charImage.measurement.plugins['ext_shapeHSM_HsmSourceMomentsRound'].addFlux=None

# Mask planes used to reject bad pixels.
config.processCcd.charImage.measurement.plugins['ext_shapeHSM_HsmSourceMomentsRound'].badMaskPlanes=None

# Use round weight function?
config.processCcd.charImage.measurement.plugins['ext_shapeHSM_HsmSourceMomentsRound'].roundMoments=None

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['ext_shapeHSM_HsmPsfMoments'].doMeasure=True

config.processCcd.charImage.measurement.plugins.names=['ext_shapeHSM_HsmPsfMoments', 'ext_shapeHSM_HsmSourceMoments', 'modelfit_DoubleShapeletPsfApprox', 'base_SdssShape', 'base_PixelFlags', 'ext_shapeHSM_HsmShapeRegauss', 'base_CircularApertureFlux', 'ext_photometryKron_KronFlux', 'base_Jacobian', 'ext_shapeHSM_HsmSourceMomentsRound', 'base_PsfFlux', 'modelfit_CModel', 'ext_convolved_ConvolvedFlux', 'base_FPPosition', 'base_SdssCentroid', 'base_GaussianFlux']
# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['base_PsfFlux'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.processCcd.charImage.measurement.undeblended['base_PsfFlux'].badMaskPlanes=[]

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['base_PeakLikelihoodFlux'].doMeasure=True

# Name of warping kernel (e.g. "lanczos4") used to compute the peak
config.processCcd.charImage.measurement.undeblended['base_PeakLikelihoodFlux'].warpingKernelName='lanczos4'

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['base_GaussianFlux'].doMeasure=True

# FIXME! NEVER DOCUMENTED!
config.processCcd.charImage.measurement.undeblended['base_GaussianFlux'].background=0.0

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['base_NaiveCentroid'].doMeasure=True

# Value to subtract from the image pixel values
config.processCcd.charImage.measurement.undeblended['base_NaiveCentroid'].background=0.0

# Do check that the centroid is contained in footprint.
config.processCcd.charImage.measurement.undeblended['base_NaiveCentroid'].doFootprintCheck=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.processCcd.charImage.measurement.undeblended['base_NaiveCentroid'].maxDistToPeak=-1.0

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['base_SdssCentroid'].doMeasure=True

# maximum allowed binning
config.processCcd.charImage.measurement.undeblended['base_SdssCentroid'].binmax=16

# Do check that the centroid is contained in footprint.
config.processCcd.charImage.measurement.undeblended['base_SdssCentroid'].doFootprintCheck=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.processCcd.charImage.measurement.undeblended['base_SdssCentroid'].maxDistToPeak=-1.0

# if the peak's less than this insist on binning at least once
config.processCcd.charImage.measurement.undeblended['base_SdssCentroid'].peakMin=-1.0

# fiddle factor for adjusting the binning
config.processCcd.charImage.measurement.undeblended['base_SdssCentroid'].wfac=1.5

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['base_PixelFlags'].doMeasure=True

# List of mask planes to be searched for which occur anywhere within a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.processCcd.charImage.measurement.undeblended['base_PixelFlags'].masksFpAnywhere=[]

# List of mask planes to be searched for which occur in the center of a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.processCcd.charImage.measurement.undeblended['base_PixelFlags'].masksFpCenter=[]

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['base_SdssShape'].doMeasure=True

# Additional value to add to background
config.processCcd.charImage.measurement.undeblended['base_SdssShape'].background=0.0

# Whether to also compute the shape of the PSF model
config.processCcd.charImage.measurement.undeblended['base_SdssShape'].doMeasurePsf=True

# Maximum number of iterations
config.processCcd.charImage.measurement.undeblended['base_SdssShape'].maxIter=100

# Maximum centroid shift, limited to 2-10
config.processCcd.charImage.measurement.undeblended['base_SdssShape'].maxShift=0.0

# Convergence tolerance for e1,e2
config.processCcd.charImage.measurement.undeblended['base_SdssShape'].tol1=9.999999747378752e-06

# Convergence tolerance for FWHM
config.processCcd.charImage.measurement.undeblended['base_SdssShape'].tol2=9.999999747378752e-05

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['base_ScaledApertureFlux'].doMeasure=True

# Scaling factor of PSF FWHM for aperture radius.
config.processCcd.charImage.measurement.undeblended['base_ScaledApertureFlux'].scale=3.14

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.processCcd.charImage.measurement.undeblended['base_ScaledApertureFlux'].shiftKernel='lanczos5'

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['base_CircularApertureFlux'].doMeasure=True

# Maximum radius (in pixels) for which the sinc algorithm should be used instead of the faster naive algorithm.  For elliptical apertures, this is the minor axis radius.
config.processCcd.charImage.measurement.undeblended['base_CircularApertureFlux'].maxSincRadius=10.0

# Radius (in pixels) of apertures.
config.processCcd.charImage.measurement.undeblended['base_CircularApertureFlux'].radii=[3.0, 4.5, 6.0, 9.0, 12.0, 17.0, 25.0, 35.0, 50.0, 70.0]

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.processCcd.charImage.measurement.undeblended['base_CircularApertureFlux'].shiftKernel='lanczos5'

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['base_Blendedness'].doMeasure=True

# Whether to compute quantities related to the Gaussian-weighted flux
config.processCcd.charImage.measurement.undeblended['base_Blendedness'].doFlux=True

# Whether to compute HeavyFootprint dot products (the old deblend.blendedness parameter)
config.processCcd.charImage.measurement.undeblended['base_Blendedness'].doOld=True

# Whether to compute quantities related to the Gaussian-weighted shape
config.processCcd.charImage.measurement.undeblended['base_Blendedness'].doShape=True

# Radius factor that sets the maximum extent of the weight function (and hence the flux measurements)
config.processCcd.charImage.measurement.undeblended['base_Blendedness'].nSigmaWeightMax=3.0

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['base_LocalBackground'].doMeasure=True

# Inner radius for background annulus as a multiple of the PSF sigma
config.processCcd.charImage.measurement.undeblended['base_LocalBackground'].annulusInner=7.0

# Outer radius for background annulus as a multiple of the PSF sigma
config.processCcd.charImage.measurement.undeblended['base_LocalBackground'].annulusOuter=15.0

# Mask planes that indicate pixels that should be excluded from the measurement
config.processCcd.charImage.measurement.undeblended['base_LocalBackground'].badMaskPlanes=['BAD', 'SAT', 'NO_DATA']

# Number of sigma-clipping iterations for background measurement
config.processCcd.charImage.measurement.undeblended['base_LocalBackground'].bgIter=3

# Rejection threshold (in standard deviations) for background measurement
config.processCcd.charImage.measurement.undeblended['base_LocalBackground'].bgRej=3.0

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['base_FPPosition'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['base_Jacobian'].doMeasure=True

# Nominal pixel size (arcsec)
config.processCcd.charImage.measurement.undeblended['base_Jacobian'].pixelScale=0.5

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['base_Variance'].doMeasure=True

# Scale factor to apply to shape for aperture
config.processCcd.charImage.measurement.undeblended['base_Variance'].scale=5.0

# Mask planes to ignore
config.processCcd.charImage.measurement.undeblended['base_Variance'].mask=['DETECTED', 'DETECTED_NEGATIVE', 'BAD', 'SAT']

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['base_InputCount'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['base_PeakCentroid'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['base_SkyCoord'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['subaru_FilterFraction'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].doMeasure=True

# Shapelet order of inner expansion (0 == Gaussian)
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].innerOrder=2

# Don't allow the semi-major radius of any component to go above this fraction of the PSF image width
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].maxRadiusBoxFraction=0.4

# Don't allow the semi-minor radius of any component to drop below this value (pixels)
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].minRadius=1.0

# Don't allow the determinant radii of the two components to differ by less than this (pixels)
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].minRadiusDiff=0.5

# whether to save all iterations for debugging purposes
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionSolverTolerance=1e-08

# Shapelet order of outer expansion (0 == Gaussian)
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].outerOrder=1

# Initial outer Gaussian peak height divided by inner Gaussian peak height
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].peakRatio=0.1

# Initial outer radius divided by inner radius
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].radiusRatio=2.0

config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models={}
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.radiusPriorSigma=0.5

config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.radiusPriorSigma=0.5

config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.order=2

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.order=1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.radiusPriorSigma=0.5

config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.order=4

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.order=4

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.radiusPriorSigma=0.5

# a sequence of model names indicating which models should be fit, and their order
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].sequence=['DoubleShapelet']

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].doMeasure=True

# Whether to record the steps the optimizer takes (or just the number, if running as a plugin)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.doRecordHistory=True

# Whether to record the time spent in this stage
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.doRecordTime=True

# Softened core width for ellipticity distribution (conformal shear units).
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.empiricalPriorConfig.ellipticityCore=0.001

# Width of exponential ellipticity distribution (conformal shear units).
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.empiricalPriorConfig.ellipticitySigma=0.3

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusMinInner=-6.0

# Minimum ln(radius).
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusMinOuter=-6.001

# Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T.
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusMu=-1.0

# Number of degrees of freedom for the Student's T distribution on ln(radius).
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusNu=50.0

# Width of the Student's T distribution in ln(radius).
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusSigma=0.45

# Ellipticity magnitude (conformal shear units) at which the softened cutoff begins
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.linearPriorConfig.ellipticityMaxInner=2.0

# Maximum ellipticity magnitude (conformal shear units)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.linearPriorConfig.ellipticityMaxOuter=2.001

# ln(radius) at which the softened cutoff begins towards the maximum
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.linearPriorConfig.logRadiusMaxInner=3.0

# Maximum ln(radius)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.linearPriorConfig.logRadiusMaxOuter=3.001

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.linearPriorConfig.logRadiusMinInner=-6.0

# The ratio P(logRadiusMinInner)/P(logRadiusMaxInner)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.linearPriorConfig.logRadiusMinMaxRatio=1.0

# Minimum ln(radius)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.linearPriorConfig.logRadiusMinOuter=-6.001

# Maximum radius used in approximating profile with Gaussians (0=default for this profile)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.maxRadius=0

# Number of Gaussian used to approximate the profile
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.nComponents=8

# whether to save all iterations for debugging purposes
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.optimizer.trustRegionSolverTolerance=1e-08

# Name of the Prior that defines the model to fit (a filename in $MEAS_MODELFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting.
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.priorName=''

# One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.priorSource='EMPIRICAL'

# Name of the shapelet.RadialProfile that defines the model to fit
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.profileName='luv'

# Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.usePixelWeights=False

# Scale the likelihood by this factor to artificially reweight it w.r.t. the prior.
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.weightsMultiplier=1.0

# Whether to record the steps the optimizer takes (or just the number, if running as a plugin)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.doRecordHistory=True

# Whether to record the time spent in this stage
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.doRecordTime=True

# Softened core width for ellipticity distribution (conformal shear units).
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.empiricalPriorConfig.ellipticityCore=0.001

# Width of exponential ellipticity distribution (conformal shear units).
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.empiricalPriorConfig.ellipticitySigma=0.3

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusMinInner=-6.0

# Minimum ln(radius).
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusMinOuter=-6.001

# Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T.
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusMu=-1.0

# Number of degrees of freedom for the Student's T distribution on ln(radius).
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusNu=50.0

# Width of the Student's T distribution in ln(radius).
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusSigma=0.45

# Ellipticity magnitude (conformal shear units) at which the softened cutoff begins
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.linearPriorConfig.ellipticityMaxInner=2.0

# Maximum ellipticity magnitude (conformal shear units)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.linearPriorConfig.ellipticityMaxOuter=2.001

# ln(radius) at which the softened cutoff begins towards the maximum
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.linearPriorConfig.logRadiusMaxInner=3.0

# Maximum ln(radius)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.linearPriorConfig.logRadiusMaxOuter=3.001

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.linearPriorConfig.logRadiusMinInner=-6.0

# The ratio P(logRadiusMinInner)/P(logRadiusMaxInner)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.linearPriorConfig.logRadiusMinMaxRatio=1.0

# Minimum ln(radius)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.linearPriorConfig.logRadiusMinOuter=-6.001

# Maximum radius used in approximating profile with Gaussians (0=default for this profile)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.maxRadius=0

# Number of Gaussian used to approximate the profile
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.nComponents=6

# whether to save all iterations for debugging purposes
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.optimizer.maxOuterIterations=250

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.optimizer.trustRegionSolverTolerance=1e-08

# Name of the Prior that defines the model to fit (a filename in $MEAS_MODELFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting.
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.priorName=''

# One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.priorSource='EMPIRICAL'

# Name of the shapelet.RadialProfile that defines the model to fit
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.profileName='lux'

# Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.usePixelWeights=False

# Scale the likelihood by this factor to artificially reweight it w.r.t. the prior.
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.weightsMultiplier=1.0

# If the 2nd-moments shape used to initialize the fit failed, use the PSF moments multiplied by this.  If <= 0.0, abort the fit early instead.
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].fallbackInitialMomentsPsfFactor=1.5

# Whether to record the steps the optimizer takes (or just the number, if running as a plugin)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.doRecordHistory=True

# Whether to record the time spent in this stage
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.doRecordTime=True

# Softened core width for ellipticity distribution (conformal shear units).
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.empiricalPriorConfig.ellipticityCore=0.001

# Width of exponential ellipticity distribution (conformal shear units).
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.empiricalPriorConfig.ellipticitySigma=0.3

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusMinInner=-6.0

# Minimum ln(radius).
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusMinOuter=-6.001

# Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T.
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusMu=-1.0

# Number of degrees of freedom for the Student's T distribution on ln(radius).
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusNu=50.0

# Width of the Student's T distribution in ln(radius).
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusSigma=0.45

# Ellipticity magnitude (conformal shear units) at which the softened cutoff begins
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.linearPriorConfig.ellipticityMaxInner=2.0

# Maximum ellipticity magnitude (conformal shear units)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.linearPriorConfig.ellipticityMaxOuter=2.001

# ln(radius) at which the softened cutoff begins towards the maximum
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.linearPriorConfig.logRadiusMaxInner=3.0

# Maximum ln(radius)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.linearPriorConfig.logRadiusMaxOuter=3.001

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.linearPriorConfig.logRadiusMinInner=-6.0

# The ratio P(logRadiusMinInner)/P(logRadiusMaxInner)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.linearPriorConfig.logRadiusMinMaxRatio=1.0

# Minimum ln(radius)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.linearPriorConfig.logRadiusMinOuter=-6.001

# Maximum radius used in approximating profile with Gaussians (0=default for this profile)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.maxRadius=0

# Number of Gaussian used to approximate the profile
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.nComponents=3

# whether to save all iterations for debugging purposes
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.optimizer.gradientThreshold=0.01

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.optimizer.minTrustRadiusThreshold=0.01

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.optimizer.trustRegionSolverTolerance=1e-08

# Name of the Prior that defines the model to fit (a filename in $MEAS_MODELFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting.
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.priorName=''

# One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.priorSource='EMPIRICAL'

# Name of the shapelet.RadialProfile that defines the model to fit
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.profileName='lux'

# Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.usePixelWeights=True

# Scale the likelihood by this factor to artificially reweight it w.r.t. the prior.
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.weightsMultiplier=1.0

# Minimum initial radius in pixels (used to regularize initial moments-based PSF deconvolution)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].minInitialRadius=0.1

# Field name prefix of the Shapelet PSF approximation used to convolve the galaxy model; must contain a set of fields matching the schema defined by shapelet.MultiShapeletFunctionKey.
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].psfName='modelfit_DoubleShapeletPsfApprox'

# Mask planes that indicate pixels that should be ignored in the fit.
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].region.badMaskPlanes=['EDGE', 'SAT', 'BAD', 'NO_DATA']

# Abort if the fit region grows beyond this many pixels.
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].region.maxArea=100000

# Maximum fraction of pixels that may be ignored due to masks; more than this and we don't even try.
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].region.maxBadPixelFraction=0.1

# Use this multiple of the initial fit ellipse then grow by the PSF width to determine the maximum final fit region size.
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].region.nFitRadiiMax=3.0

# Use this multiple of the initial fit ellipse then grow by the PSF width to determine the minimum final fit region size.
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].region.nFitRadiiMin=1.0

# Use this multiple of the Kron ellipse to set the fit region (for the final fit region, subject to the nFitRadiiMin and nFitRadiiMax constraints).
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].region.nKronRadii=1.5

# Grow the initial fit ellipses by this factor before comparing with the Kron/Footprint region
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].region.nPsfSigmaGrow=2.0

# If the Kron radius is less than this multiple of the PSF width, ignore it and fall back to a PSF-oriented ellipse scaled to match the area of the footprint or this radius (whichever is larger).
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].region.nPsfSigmaMin=4.0

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['ext_photometryKron_KronFlux'].doMeasure=True

# If true check that the Kron radius exceeds some minimum
config.processCcd.charImage.measurement.undeblended['ext_photometryKron_KronFlux'].enforceMinimumRadius=True

# if true, use existing shape and centroid measurements instead of fitting
config.processCcd.charImage.measurement.undeblended['ext_photometryKron_KronFlux'].fixed=False

# Largest aperture for which to use the slow, accurate, sinc aperture code
config.processCcd.charImage.measurement.undeblended['ext_photometryKron_KronFlux'].maxSincRadius=10.0

# Minimum Kron radius (if == 0.0 use PSF's Kron radius) if enforceMinimumRadius. Also functions as fallback aperture radius if set.
config.processCcd.charImage.measurement.undeblended['ext_photometryKron_KronFlux'].minimumRadius=0.0

# Number of times to iterate when setting the Kron radius
config.processCcd.charImage.measurement.undeblended['ext_photometryKron_KronFlux'].nIterForRadius=1

# Number of Kron radii for Kron flux
config.processCcd.charImage.measurement.undeblended['ext_photometryKron_KronFlux'].nRadiusForFlux=2.5

# Multiplier of rms size for aperture used to initially estimate the Kron radius
config.processCcd.charImage.measurement.undeblended['ext_photometryKron_KronFlux'].nSigmaForRadius=6.0

# Name of field specifying reference Kron radius for forced measurement
config.processCcd.charImage.measurement.undeblended['ext_photometryKron_KronFlux'].refRadiusName='ext_photometryKron_KronFlux_radius'

# Smooth image with N(0, smoothingSigma^2) Gaussian while estimating R_K
config.processCcd.charImage.measurement.undeblended['ext_photometryKron_KronFlux'].smoothingSigma=-1.0

# Use the Footprint size as part of initial estimate of Kron radius
config.processCcd.charImage.measurement.undeblended['ext_photometryKron_KronFlux'].useFootprintRadius=False

# list of target seeings (FWHM, pixels)
config.processCcd.charImage.measurement.undeblended['ext_convolved_ConvolvedFlux'].seeing=[3.5, 5.0, 6.5]

# scaling factor of kernel sigma for kernel size
config.processCcd.charImage.measurement.undeblended['ext_convolved_ConvolvedFlux'].kernelScale=4.0

# Maximum radius (in pixels) for which the sinc algorithm should be used instead of the faster naive algorithm.  For elliptical apertures, this is the minor axis radius.
config.processCcd.charImage.measurement.undeblended['ext_convolved_ConvolvedFlux'].aperture.maxSincRadius=10.0

# Radius (in pixels) of apertures.
config.processCcd.charImage.measurement.undeblended['ext_convolved_ConvolvedFlux'].aperture.radii=[3.3, 4.5, 6.0]

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.processCcd.charImage.measurement.undeblended['ext_convolved_ConvolvedFlux'].aperture.shiftKernel='lanczos5'

# name of Kron radius field in reference
config.processCcd.charImage.measurement.undeblended['ext_convolved_ConvolvedFlux'].kronRadiusName='ext_photometryKron_KronFlux_radius'

# Largest aperture for which to use the sinc aperture code for Kron (pixels)
config.processCcd.charImage.measurement.undeblended['ext_convolved_ConvolvedFlux'].maxSincRadius=10.0

# Number of Kron radii for Kron flux
config.processCcd.charImage.measurement.undeblended['ext_convolved_ConvolvedFlux'].kronRadiusForFlux=2.5

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['ext_convolved_ConvolvedFlux'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['ext_shapeHSM_HsmShapeBj'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.processCcd.charImage.measurement.undeblended['ext_shapeHSM_HsmShapeBj'].badMaskPlanes=None

# Field name for number of deblend children
config.processCcd.charImage.measurement.undeblended['ext_shapeHSM_HsmShapeBj'].deblendNChild=None

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['ext_shapeHSM_HsmShapeLinear'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.processCcd.charImage.measurement.undeblended['ext_shapeHSM_HsmShapeLinear'].badMaskPlanes=None

# Field name for number of deblend children
config.processCcd.charImage.measurement.undeblended['ext_shapeHSM_HsmShapeLinear'].deblendNChild=None

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['ext_shapeHSM_HsmShapeKsb'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.processCcd.charImage.measurement.undeblended['ext_shapeHSM_HsmShapeKsb'].badMaskPlanes=None

# Field name for number of deblend children
config.processCcd.charImage.measurement.undeblended['ext_shapeHSM_HsmShapeKsb'].deblendNChild=None

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['ext_shapeHSM_HsmShapeRegauss'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.processCcd.charImage.measurement.undeblended['ext_shapeHSM_HsmShapeRegauss'].badMaskPlanes=None

# Field name for number of deblend children
config.processCcd.charImage.measurement.undeblended['ext_shapeHSM_HsmShapeRegauss'].deblendNChild=None

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['ext_shapeHSM_HsmSourceMoments'].doMeasure=True

# Store measured flux?
config.processCcd.charImage.measurement.undeblended['ext_shapeHSM_HsmSourceMoments'].addFlux=None

# Mask planes used to reject bad pixels.
config.processCcd.charImage.measurement.undeblended['ext_shapeHSM_HsmSourceMoments'].badMaskPlanes=None

# Use round weight function?
config.processCcd.charImage.measurement.undeblended['ext_shapeHSM_HsmSourceMoments'].roundMoments=None

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['ext_shapeHSM_HsmSourceMomentsRound'].doMeasure=True

# Store measured flux?
config.processCcd.charImage.measurement.undeblended['ext_shapeHSM_HsmSourceMomentsRound'].addFlux=None

# Mask planes used to reject bad pixels.
config.processCcd.charImage.measurement.undeblended['ext_shapeHSM_HsmSourceMomentsRound'].badMaskPlanes=None

# Use round weight function?
config.processCcd.charImage.measurement.undeblended['ext_shapeHSM_HsmSourceMomentsRound'].roundMoments=None

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['ext_shapeHSM_HsmPsfMoments'].doMeasure=True

config.processCcd.charImage.measurement.undeblended.names=[]
# Run subtasks to measure and apply aperture corrections
config.processCcd.charImage.doApCorr=True

# Field name prefix for the flux other measurements should be aperture corrected to match
config.processCcd.charImage.measureApCorr.refFluxName='slot_CalibFlux'

# Apply flux limit?
config.processCcd.charImage.measureApCorr.sourceSelector['science'].doFluxLimit=False

# Apply flag limitation?
config.processCcd.charImage.measureApCorr.sourceSelector['science'].doFlags=False

# Apply unresolved limitation?
config.processCcd.charImage.measureApCorr.sourceSelector['science'].doUnresolved=False

# Apply signal-to-noise limit?
config.processCcd.charImage.measureApCorr.sourceSelector['science'].doSignalToNoise=False

# Apply isolated limitation?
config.processCcd.charImage.measureApCorr.sourceSelector['science'].doIsolated=False

# Select objects with value greater than this
config.processCcd.charImage.measureApCorr.sourceSelector['science'].fluxLimit.minimum=None

# Select objects with value less than this
config.processCcd.charImage.measureApCorr.sourceSelector['science'].fluxLimit.maximum=None

# Name of the source flux field to use.
config.processCcd.charImage.measureApCorr.sourceSelector['science'].fluxLimit.fluxField='slot_CalibFlux_flux'

# List of source flag fields that must be set for a source to be used.
config.processCcd.charImage.measureApCorr.sourceSelector['science'].flags.good=[]

# List of source flag fields that must NOT be set for a source to be used.
config.processCcd.charImage.measureApCorr.sourceSelector['science'].flags.bad=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_saturated', 'base_PsfFlux_flags']

# Select objects with value greater than this
config.processCcd.charImage.measureApCorr.sourceSelector['science'].unresolved.minimum=None

# Select objects with value less than this
config.processCcd.charImage.measureApCorr.sourceSelector['science'].unresolved.maximum=None

# Name of column for star/galaxy separation
config.processCcd.charImage.measureApCorr.sourceSelector['science'].unresolved.name='base_ClassificationExtendedness_value'

# Select objects with value greater than this
config.processCcd.charImage.measureApCorr.sourceSelector['science'].signalToNoise.minimum=None

# Select objects with value less than this
config.processCcd.charImage.measureApCorr.sourceSelector['science'].signalToNoise.maximum=None

# Name of the source flux field to use.
config.processCcd.charImage.measureApCorr.sourceSelector['science'].signalToNoise.fluxField='base_PsfFlux_flux'

# Name of the source flux error field to use.
config.processCcd.charImage.measureApCorr.sourceSelector['science'].signalToNoise.errField='base_PsfFlux_fluxSigma'

# Name of column for parent
config.processCcd.charImage.measureApCorr.sourceSelector['science'].isolated.parentName='parent'

# Name of column for nChild
config.processCcd.charImage.measureApCorr.sourceSelector['science'].isolated.nChildName='deblend_nChild'

# Apply magnitude limit?
config.processCcd.charImage.measureApCorr.sourceSelector['references'].doMagLimit=False

# Apply flag limitation?
config.processCcd.charImage.measureApCorr.sourceSelector['references'].doFlags=False

# Apply signal-to-noise limit?
config.processCcd.charImage.measureApCorr.sourceSelector['references'].doSignalToNoise=False

# Apply magnitude error limit?
config.processCcd.charImage.measureApCorr.sourceSelector['references'].doMagError=False

# Select objects with value greater than this
config.processCcd.charImage.measureApCorr.sourceSelector['references'].magLimit.minimum=None

# Select objects with value less than this
config.processCcd.charImage.measureApCorr.sourceSelector['references'].magLimit.maximum=None

# Name of the source flux field to use.
config.processCcd.charImage.measureApCorr.sourceSelector['references'].magLimit.fluxField='flux'

# List of source flag fields that must be set for a source to be used.
config.processCcd.charImage.measureApCorr.sourceSelector['references'].flags.good=[]

# List of source flag fields that must NOT be set for a source to be used.
config.processCcd.charImage.measureApCorr.sourceSelector['references'].flags.bad=[]

# Select objects with value greater than this
config.processCcd.charImage.measureApCorr.sourceSelector['references'].signalToNoise.minimum=None

# Select objects with value less than this
config.processCcd.charImage.measureApCorr.sourceSelector['references'].signalToNoise.maximum=None

# Name of the source flux field to use.
config.processCcd.charImage.measureApCorr.sourceSelector['references'].signalToNoise.fluxField='flux'

# Name of the source flux error field to use.
config.processCcd.charImage.measureApCorr.sourceSelector['references'].signalToNoise.errField='flux_err'

# Select objects with value greater than this
config.processCcd.charImage.measureApCorr.sourceSelector['references'].magError.minimum=None

# Select objects with value less than this
config.processCcd.charImage.measureApCorr.sourceSelector['references'].magError.maximum=None

# Name of the source flux error field to use.
config.processCcd.charImage.measureApCorr.sourceSelector['references'].magError.magErrField='mag_err'

config.processCcd.charImage.measureApCorr.sourceSelector['references'].colorLimits={}
# specify the minimum psfFlux for good Psf Candidates
config.processCcd.charImage.measureApCorr.sourceSelector['objectSize'].fluxMin=12500.0

# specify the maximum psfFlux for good Psf Candidates (ignored if == 0)
config.processCcd.charImage.measureApCorr.sourceSelector['objectSize'].fluxMax=0.0

# minimum width to include in histogram
config.processCcd.charImage.measureApCorr.sourceSelector['objectSize'].widthMin=0.0

# maximum width to include in histogram
config.processCcd.charImage.measureApCorr.sourceSelector['objectSize'].widthMax=10.0

# Name of field in Source to use for flux measurement
config.processCcd.charImage.measureApCorr.sourceSelector['objectSize'].sourceFluxField='base_GaussianFlux_flux'

# Standard deviation of width allowed to be interpreted as good stars
config.processCcd.charImage.measureApCorr.sourceSelector['objectSize'].widthStdAllowed=0.15

# Keep objects within this many sigma of cluster 0's median
config.processCcd.charImage.measureApCorr.sourceSelector['objectSize'].nSigmaClip=2.0

# List of flags which cause a source to be rejected as bad
config.processCcd.charImage.measureApCorr.sourceSelector['objectSize'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter', 'base_PixelFlags_flag_bad', 'base_PixelFlags_flag_interpolated']

# Name of a flag field that is True for Sources that should be used.
config.processCcd.charImage.measureApCorr.sourceSelector['flagged'].field='calib_psfUsed'

# List of flags which cause a source to be rejected as bad
config.processCcd.charImage.measureApCorr.sourceSelector['astrometry'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter', 'base_PixelFlags_flag_bad']

# Type of source flux; typically one of Ap or Psf
config.processCcd.charImage.measureApCorr.sourceSelector['astrometry'].sourceFluxType='Ap'

# Minimum allowed signal-to-noise ratio for sources used for matching (in the flux specified by sourceFluxType); <= 0 for no limit
config.processCcd.charImage.measureApCorr.sourceSelector['astrometry'].minSnr=10.0

# Type of source flux; typically one of Ap or Psf
config.processCcd.charImage.measureApCorr.sourceSelector['matcher'].sourceFluxType='Ap'

# Minimum allowed signal-to-noise ratio for sources used for matching (in the flux specified by sourceFluxType); <= 0 for no limit
config.processCcd.charImage.measureApCorr.sourceSelector['matcher'].minSnr=40.0

# Type of source flux; typically one of Ap or Psf
config.processCcd.charImage.measureApCorr.sourceSelector['matcherPessimistic'].sourceFluxType='Ap'

# Minimum allowed signal-to-noise ratio for sources used for matching (in the flux specified by sourceFluxType); <= 0 for no limit
config.processCcd.charImage.measureApCorr.sourceSelector['matcherPessimistic'].minSnr=40.0

# specify the minimum psfFlux for good Psf Candidates
config.processCcd.charImage.measureApCorr.sourceSelector['catalog'].fluxLim=0.0

# specify the maximum psfFlux for good Psf Candidates (ignored if == 0)
config.processCcd.charImage.measureApCorr.sourceSelector['catalog'].fluxMax=0.0

# List of flags which cause a source to be rejected as bad
config.processCcd.charImage.measureApCorr.sourceSelector['catalog'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter']

config.processCcd.charImage.measureApCorr.sourceSelector.name='flagged'
# Minimum number of degrees of freedom (# of valid data points - # of parameters); if this is exceeded, the order of the fit is decreased (in both dimensions), and if we can't decrease it enough, we'll raise ValueError.
config.processCcd.charImage.measureApCorr.minDegreesOfFreedom=1

# maximum Chebyshev function order in x
config.processCcd.charImage.measureApCorr.fitConfig.orderX=2

# maximum Chebyshev function order in y
config.processCcd.charImage.measureApCorr.fitConfig.orderY=2

# if true, only include terms where the sum of the x and y order is less than or equal to max(orderX, orderY)
config.processCcd.charImage.measureApCorr.fitConfig.triangular=True

# Number of iterations for sigma clipping
config.processCcd.charImage.measureApCorr.numIter=4

# Number of standard devisations to clip at
config.processCcd.charImage.measureApCorr.numSigmaClip=3.0

# Allow these measurement algorithms to fail without an exception
config.processCcd.charImage.measureApCorr.allowFailure=['ext_convolved_ConvolvedFlux_0_3_3', 'ext_convolved_ConvolvedFlux_0_4_5', 'ext_convolved_ConvolvedFlux_0_6_0', 'ext_convolved_ConvolvedFlux_1_3_3', 'ext_convolved_ConvolvedFlux_1_4_5', 'ext_convolved_ConvolvedFlux_1_6_0', 'ext_convolved_ConvolvedFlux_2_3_3', 'ext_convolved_ConvolvedFlux_2_4_5', 'ext_convolved_ConvolvedFlux_2_6_0', 'ext_convolved_ConvolvedFlux_3_3_3', 'ext_convolved_ConvolvedFlux_3_4_5', 'ext_convolved_ConvolvedFlux_3_6_0', 'ext_convolved_ConvolvedFlux_0_kron', 'ext_convolved_ConvolvedFlux_1_kron', 'ext_convolved_ConvolvedFlux_2_kron', 'ext_convolved_ConvolvedFlux_3_kron']

# flux measurement algorithms in getApCorrNameSet() to ignore; if a name is listed that does not appear in getApCorrNameSet() then a warning is logged
config.processCcd.charImage.applyApCorr.ignoreList=[]

# set the general failure flag for a flux when it cannot be aperture-corrected?
config.processCcd.charImage.applyApCorr.doFlagApCorrFailures=True

# flux measurement algorithms to be aperture-corrected by reference to another algorithm; this is a mapping alg1:alg2, where 'alg1' is the algorithm being corrected, and 'alg2' is the algorithm supplying the corrections
config.processCcd.charImage.applyApCorr.proxies={}

# critical ratio of model to psf flux
config.processCcd.charImage.catalogCalculation.plugins['base_ClassificationExtendedness'].fluxRatio=0.985

# correction factor for modelFlux error
config.processCcd.charImage.catalogCalculation.plugins['base_ClassificationExtendedness'].modelErrFactor=0.0

# correction factor for psfFlux error
config.processCcd.charImage.catalogCalculation.plugins['base_ClassificationExtendedness'].psfErrFactor=0.0

config.processCcd.charImage.catalogCalculation.plugins.names=['base_ClassificationExtendedness', 'base_FootprintArea']
# Replace the existing PSF model with a simplified version that has the same sigma at the start of each PSF determination iteration? Doing so makes PSF determination converge more robustly and quickly.
config.processCcd.charImage.useSimplePsf=True

# Estimated FWHM of simple Gaussian PSF model, in pixels. Ignored if input exposure has a PSF model.
config.processCcd.charImage.installSimplePsf.fwhm=3.5322300675464238

# Width and height of PSF model, in pixels. Must be odd.
config.processCcd.charImage.installSimplePsf.width=11

import lsst.meas.algorithms.loadIndexedReferenceObjects
config.processCcd.charImage.refObjLoader.retarget(target=lsst.meas.algorithms.loadIndexedReferenceObjects.LoadIndexedReferenceObjectsTask, ConfigClass=lsst.meas.algorithms.loadIndexedReferenceObjects.LoadIndexedReferenceObjectsConfig)

# Padding to add to 4 all edges of the bounding box (pixels)
config.processCcd.charImage.refObjLoader.pixelMargin=300

# Default reference catalog filter to use if filter not specified in exposure; if blank then filter must be specified in exposure
config.processCcd.charImage.refObjLoader.defaultFilter=''

# Mapping of camera filter name: reference catalog filter name; each reference filter must exist
config.processCcd.charImage.refObjLoader.filterMap={'B': 'g', 'V': 'r', 'R': 'r', 'I': 'i', 'r2': 'r', 'i2': 'i', 'N387': 'g', 'N468': 'g', 'N515': 'g', 'N527': 'g', 'N656': 'r', 'N718': 'i', 'N816': 'i', 'N921': 'z', 'N926': 'z', 'N973': 'y', 'N1010': 'y', 'I945': 'z'}

# Name of the ingested reference dataset
config.processCcd.charImage.refObjLoader.ref_dataset_name='ps1_pv3_3pi_20170110'

# Maximum separation between reference objects and sources beyond which they will not be considered a match (arcsec)
config.processCcd.charImage.ref_match.matcher.maxMatchDistArcSec=2.0

# Number of bright stars to use
config.processCcd.charImage.ref_match.matcher.numBrightStars=50

# Minimum number of matched pairs; see also minFracMatchedPairs
config.processCcd.charImage.ref_match.matcher.minMatchedPairs=30

# Minimum number of matched pairs as a fraction of the smaller of the number of reference stars or the number of good sources; the actual minimum is the smaller of this value or minMatchedPairs
config.processCcd.charImage.ref_match.matcher.minFracMatchedPairs=0.3

# Maximum allowed shift of WCS, due to matching (pixel). When changing this value, the LoadReferenceObjectsConfig.pixelMargin should also be updated.
config.processCcd.charImage.ref_match.matcher.maxOffsetPix=250

# Rotation angle allowed between sources and position reference objects (degrees)
config.processCcd.charImage.ref_match.matcher.maxRotationDeg=1.145916

# Allowed non-perpendicularity of x and y (degree)
config.processCcd.charImage.ref_match.matcher.allowedNonperpDeg=0.2

# number of points to define a shape for matching
config.processCcd.charImage.ref_match.matcher.numPointsForShape=6

# maximum determinant of linear transformation matrix for a usable solution
config.processCcd.charImage.ref_match.matcher.maxDeterminant=0.02

# Apply flux limit?
config.processCcd.charImage.ref_match.matcher.sourceSelector['science'].doFluxLimit=False

# Apply flag limitation?
config.processCcd.charImage.ref_match.matcher.sourceSelector['science'].doFlags=False

# Apply unresolved limitation?
config.processCcd.charImage.ref_match.matcher.sourceSelector['science'].doUnresolved=False

# Apply signal-to-noise limit?
config.processCcd.charImage.ref_match.matcher.sourceSelector['science'].doSignalToNoise=False

# Apply isolated limitation?
config.processCcd.charImage.ref_match.matcher.sourceSelector['science'].doIsolated=False

# Select objects with value greater than this
config.processCcd.charImage.ref_match.matcher.sourceSelector['science'].fluxLimit.minimum=None

# Select objects with value less than this
config.processCcd.charImage.ref_match.matcher.sourceSelector['science'].fluxLimit.maximum=None

# Name of the source flux field to use.
config.processCcd.charImage.ref_match.matcher.sourceSelector['science'].fluxLimit.fluxField='slot_CalibFlux_flux'

# List of source flag fields that must be set for a source to be used.
config.processCcd.charImage.ref_match.matcher.sourceSelector['science'].flags.good=[]

# List of source flag fields that must NOT be set for a source to be used.
config.processCcd.charImage.ref_match.matcher.sourceSelector['science'].flags.bad=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_saturated', 'base_PsfFlux_flags']

# Select objects with value greater than this
config.processCcd.charImage.ref_match.matcher.sourceSelector['science'].unresolved.minimum=None

# Select objects with value less than this
config.processCcd.charImage.ref_match.matcher.sourceSelector['science'].unresolved.maximum=None

# Name of column for star/galaxy separation
config.processCcd.charImage.ref_match.matcher.sourceSelector['science'].unresolved.name='base_ClassificationExtendedness_value'

# Select objects with value greater than this
config.processCcd.charImage.ref_match.matcher.sourceSelector['science'].signalToNoise.minimum=None

# Select objects with value less than this
config.processCcd.charImage.ref_match.matcher.sourceSelector['science'].signalToNoise.maximum=None

# Name of the source flux field to use.
config.processCcd.charImage.ref_match.matcher.sourceSelector['science'].signalToNoise.fluxField='base_PsfFlux_flux'

# Name of the source flux error field to use.
config.processCcd.charImage.ref_match.matcher.sourceSelector['science'].signalToNoise.errField='base_PsfFlux_fluxSigma'

# Name of column for parent
config.processCcd.charImage.ref_match.matcher.sourceSelector['science'].isolated.parentName='parent'

# Name of column for nChild
config.processCcd.charImage.ref_match.matcher.sourceSelector['science'].isolated.nChildName='deblend_nChild'

# Apply magnitude limit?
config.processCcd.charImage.ref_match.matcher.sourceSelector['references'].doMagLimit=False

# Apply flag limitation?
config.processCcd.charImage.ref_match.matcher.sourceSelector['references'].doFlags=False

# Apply signal-to-noise limit?
config.processCcd.charImage.ref_match.matcher.sourceSelector['references'].doSignalToNoise=False

# Apply magnitude error limit?
config.processCcd.charImage.ref_match.matcher.sourceSelector['references'].doMagError=False

# Select objects with value greater than this
config.processCcd.charImage.ref_match.matcher.sourceSelector['references'].magLimit.minimum=None

# Select objects with value less than this
config.processCcd.charImage.ref_match.matcher.sourceSelector['references'].magLimit.maximum=None

# Name of the source flux field to use.
config.processCcd.charImage.ref_match.matcher.sourceSelector['references'].magLimit.fluxField='flux'

# List of source flag fields that must be set for a source to be used.
config.processCcd.charImage.ref_match.matcher.sourceSelector['references'].flags.good=[]

# List of source flag fields that must NOT be set for a source to be used.
config.processCcd.charImage.ref_match.matcher.sourceSelector['references'].flags.bad=[]

# Select objects with value greater than this
config.processCcd.charImage.ref_match.matcher.sourceSelector['references'].signalToNoise.minimum=None

# Select objects with value less than this
config.processCcd.charImage.ref_match.matcher.sourceSelector['references'].signalToNoise.maximum=None

# Name of the source flux field to use.
config.processCcd.charImage.ref_match.matcher.sourceSelector['references'].signalToNoise.fluxField='flux'

# Name of the source flux error field to use.
config.processCcd.charImage.ref_match.matcher.sourceSelector['references'].signalToNoise.errField='flux_err'

# Select objects with value greater than this
config.processCcd.charImage.ref_match.matcher.sourceSelector['references'].magError.minimum=None

# Select objects with value less than this
config.processCcd.charImage.ref_match.matcher.sourceSelector['references'].magError.maximum=None

# Name of the source flux error field to use.
config.processCcd.charImage.ref_match.matcher.sourceSelector['references'].magError.magErrField='mag_err'

config.processCcd.charImage.ref_match.matcher.sourceSelector['references'].colorLimits={}
# specify the minimum psfFlux for good Psf Candidates
config.processCcd.charImage.ref_match.matcher.sourceSelector['objectSize'].fluxMin=12500.0

# specify the maximum psfFlux for good Psf Candidates (ignored if == 0)
config.processCcd.charImage.ref_match.matcher.sourceSelector['objectSize'].fluxMax=0.0

# minimum width to include in histogram
config.processCcd.charImage.ref_match.matcher.sourceSelector['objectSize'].widthMin=0.0

# maximum width to include in histogram
config.processCcd.charImage.ref_match.matcher.sourceSelector['objectSize'].widthMax=10.0

# Name of field in Source to use for flux measurement
config.processCcd.charImage.ref_match.matcher.sourceSelector['objectSize'].sourceFluxField='base_GaussianFlux_flux'

# Standard deviation of width allowed to be interpreted as good stars
config.processCcd.charImage.ref_match.matcher.sourceSelector['objectSize'].widthStdAllowed=0.15

# Keep objects within this many sigma of cluster 0's median
config.processCcd.charImage.ref_match.matcher.sourceSelector['objectSize'].nSigmaClip=2.0

# List of flags which cause a source to be rejected as bad
config.processCcd.charImage.ref_match.matcher.sourceSelector['objectSize'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter', 'base_PixelFlags_flag_bad', 'base_PixelFlags_flag_interpolated']

# Name of a flag field that is True for Sources that should be used.
config.processCcd.charImage.ref_match.matcher.sourceSelector['flagged'].field='calib_psfUsed'

# List of flags which cause a source to be rejected as bad
config.processCcd.charImage.ref_match.matcher.sourceSelector['astrometry'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter', 'base_PixelFlags_flag_bad']

# Type of source flux; typically one of Ap or Psf
config.processCcd.charImage.ref_match.matcher.sourceSelector['astrometry'].sourceFluxType='Ap'

# Minimum allowed signal-to-noise ratio for sources used for matching (in the flux specified by sourceFluxType); <= 0 for no limit
config.processCcd.charImage.ref_match.matcher.sourceSelector['astrometry'].minSnr=10.0

# Type of source flux; typically one of Ap or Psf
config.processCcd.charImage.ref_match.matcher.sourceSelector['matcher'].sourceFluxType='Psf'

# Minimum allowed signal-to-noise ratio for sources used for matching (in the flux specified by sourceFluxType); <= 0 for no limit
config.processCcd.charImage.ref_match.matcher.sourceSelector['matcher'].minSnr=40.0

# Type of source flux; typically one of Ap or Psf
config.processCcd.charImage.ref_match.matcher.sourceSelector['matcherPessimistic'].sourceFluxType='Ap'

# Minimum allowed signal-to-noise ratio for sources used for matching (in the flux specified by sourceFluxType); <= 0 for no limit
config.processCcd.charImage.ref_match.matcher.sourceSelector['matcherPessimistic'].minSnr=40.0

# specify the minimum psfFlux for good Psf Candidates
config.processCcd.charImage.ref_match.matcher.sourceSelector['catalog'].fluxLim=0.0

# specify the maximum psfFlux for good Psf Candidates (ignored if == 0)
config.processCcd.charImage.ref_match.matcher.sourceSelector['catalog'].fluxMax=0.0

# List of flags which cause a source to be rejected as bad
config.processCcd.charImage.ref_match.matcher.sourceSelector['catalog'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter']

config.processCcd.charImage.ref_match.matcher.sourceSelector.name='matcher'
# the maximum match distance is set to  mean_match_distance + matchDistanceSigma*std_dev_match_distance; ignored if not fitting a WCS
config.processCcd.charImage.ref_match.matchDistanceSigma=2.0

# Apply flux limit?
config.processCcd.charImage.ref_match.sourceSelection.doFluxLimit=False

# Apply flag limitation?
config.processCcd.charImage.ref_match.sourceSelection.doFlags=False

# Apply unresolved limitation?
config.processCcd.charImage.ref_match.sourceSelection.doUnresolved=False

# Apply signal-to-noise limit?
config.processCcd.charImage.ref_match.sourceSelection.doSignalToNoise=False

# Apply isolated limitation?
config.processCcd.charImage.ref_match.sourceSelection.doIsolated=False

# Select objects with value greater than this
config.processCcd.charImage.ref_match.sourceSelection.fluxLimit.minimum=None

# Select objects with value less than this
config.processCcd.charImage.ref_match.sourceSelection.fluxLimit.maximum=None

# Name of the source flux field to use.
config.processCcd.charImage.ref_match.sourceSelection.fluxLimit.fluxField='slot_CalibFlux_flux'

# List of source flag fields that must be set for a source to be used.
config.processCcd.charImage.ref_match.sourceSelection.flags.good=[]

# List of source flag fields that must NOT be set for a source to be used.
config.processCcd.charImage.ref_match.sourceSelection.flags.bad=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_saturated', 'base_PsfFlux_flags']

# Select objects with value greater than this
config.processCcd.charImage.ref_match.sourceSelection.unresolved.minimum=None

# Select objects with value less than this
config.processCcd.charImage.ref_match.sourceSelection.unresolved.maximum=None

# Name of column for star/galaxy separation
config.processCcd.charImage.ref_match.sourceSelection.unresolved.name='base_ClassificationExtendedness_value'

# Select objects with value greater than this
config.processCcd.charImage.ref_match.sourceSelection.signalToNoise.minimum=None

# Select objects with value less than this
config.processCcd.charImage.ref_match.sourceSelection.signalToNoise.maximum=None

# Name of the source flux field to use.
config.processCcd.charImage.ref_match.sourceSelection.signalToNoise.fluxField='base_PsfFlux_flux'

# Name of the source flux error field to use.
config.processCcd.charImage.ref_match.sourceSelection.signalToNoise.errField='base_PsfFlux_fluxSigma'

# Name of column for parent
config.processCcd.charImage.ref_match.sourceSelection.isolated.parentName='parent'

# Name of column for nChild
config.processCcd.charImage.ref_match.sourceSelection.isolated.nChildName='deblend_nChild'

# Apply magnitude limit?
config.processCcd.charImage.ref_match.referenceSelection.doMagLimit=False

# Apply flag limitation?
config.processCcd.charImage.ref_match.referenceSelection.doFlags=False

# Apply signal-to-noise limit?
config.processCcd.charImage.ref_match.referenceSelection.doSignalToNoise=False

# Apply magnitude error limit?
config.processCcd.charImage.ref_match.referenceSelection.doMagError=False

# Select objects with value greater than this
config.processCcd.charImage.ref_match.referenceSelection.magLimit.minimum=None

# Select objects with value less than this
config.processCcd.charImage.ref_match.referenceSelection.magLimit.maximum=None

# Name of the source flux field to use.
config.processCcd.charImage.ref_match.referenceSelection.magLimit.fluxField='flux'

# List of source flag fields that must be set for a source to be used.
config.processCcd.charImage.ref_match.referenceSelection.flags.good=[]

# List of source flag fields that must NOT be set for a source to be used.
config.processCcd.charImage.ref_match.referenceSelection.flags.bad=[]

# Select objects with value greater than this
config.processCcd.charImage.ref_match.referenceSelection.signalToNoise.minimum=None

# Select objects with value less than this
config.processCcd.charImage.ref_match.referenceSelection.signalToNoise.maximum=None

# Name of the source flux field to use.
config.processCcd.charImage.ref_match.referenceSelection.signalToNoise.fluxField='flux'

# Name of the source flux error field to use.
config.processCcd.charImage.ref_match.referenceSelection.signalToNoise.errField='flux_err'

# Select objects with value greater than this
config.processCcd.charImage.ref_match.referenceSelection.magError.minimum=None

# Select objects with value less than this
config.processCcd.charImage.ref_match.referenceSelection.magError.maximum=None

# Name of the source flux error field to use.
config.processCcd.charImage.ref_match.referenceSelection.magError.magErrField='mag_err'

config.processCcd.charImage.ref_match.referenceSelection.colorLimits={}
# Apply flux limit?
config.processCcd.charImage.measurePsf.starSelector['science'].doFluxLimit=False

# Apply flag limitation?
config.processCcd.charImage.measurePsf.starSelector['science'].doFlags=False

# Apply unresolved limitation?
config.processCcd.charImage.measurePsf.starSelector['science'].doUnresolved=False

# Apply signal-to-noise limit?
config.processCcd.charImage.measurePsf.starSelector['science'].doSignalToNoise=False

# Apply isolated limitation?
config.processCcd.charImage.measurePsf.starSelector['science'].doIsolated=False

# Select objects with value greater than this
config.processCcd.charImage.measurePsf.starSelector['science'].fluxLimit.minimum=None

# Select objects with value less than this
config.processCcd.charImage.measurePsf.starSelector['science'].fluxLimit.maximum=None

# Name of the source flux field to use.
config.processCcd.charImage.measurePsf.starSelector['science'].fluxLimit.fluxField='slot_CalibFlux_flux'

# List of source flag fields that must be set for a source to be used.
config.processCcd.charImage.measurePsf.starSelector['science'].flags.good=[]

# List of source flag fields that must NOT be set for a source to be used.
config.processCcd.charImage.measurePsf.starSelector['science'].flags.bad=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_saturated', 'base_PsfFlux_flags']

# Select objects with value greater than this
config.processCcd.charImage.measurePsf.starSelector['science'].unresolved.minimum=None

# Select objects with value less than this
config.processCcd.charImage.measurePsf.starSelector['science'].unresolved.maximum=None

# Name of column for star/galaxy separation
config.processCcd.charImage.measurePsf.starSelector['science'].unresolved.name='base_ClassificationExtendedness_value'

# Select objects with value greater than this
config.processCcd.charImage.measurePsf.starSelector['science'].signalToNoise.minimum=None

# Select objects with value less than this
config.processCcd.charImage.measurePsf.starSelector['science'].signalToNoise.maximum=None

# Name of the source flux field to use.
config.processCcd.charImage.measurePsf.starSelector['science'].signalToNoise.fluxField='base_PsfFlux_flux'

# Name of the source flux error field to use.
config.processCcd.charImage.measurePsf.starSelector['science'].signalToNoise.errField='base_PsfFlux_fluxSigma'

# Name of column for parent
config.processCcd.charImage.measurePsf.starSelector['science'].isolated.parentName='parent'

# Name of column for nChild
config.processCcd.charImage.measurePsf.starSelector['science'].isolated.nChildName='deblend_nChild'

# Apply magnitude limit?
config.processCcd.charImage.measurePsf.starSelector['references'].doMagLimit=False

# Apply flag limitation?
config.processCcd.charImage.measurePsf.starSelector['references'].doFlags=False

# Apply signal-to-noise limit?
config.processCcd.charImage.measurePsf.starSelector['references'].doSignalToNoise=False

# Apply magnitude error limit?
config.processCcd.charImage.measurePsf.starSelector['references'].doMagError=False

# Select objects with value greater than this
config.processCcd.charImage.measurePsf.starSelector['references'].magLimit.minimum=None

# Select objects with value less than this
config.processCcd.charImage.measurePsf.starSelector['references'].magLimit.maximum=None

# Name of the source flux field to use.
config.processCcd.charImage.measurePsf.starSelector['references'].magLimit.fluxField='flux'

# List of source flag fields that must be set for a source to be used.
config.processCcd.charImage.measurePsf.starSelector['references'].flags.good=[]

# List of source flag fields that must NOT be set for a source to be used.
config.processCcd.charImage.measurePsf.starSelector['references'].flags.bad=[]

# Select objects with value greater than this
config.processCcd.charImage.measurePsf.starSelector['references'].signalToNoise.minimum=None

# Select objects with value less than this
config.processCcd.charImage.measurePsf.starSelector['references'].signalToNoise.maximum=None

# Name of the source flux field to use.
config.processCcd.charImage.measurePsf.starSelector['references'].signalToNoise.fluxField='flux'

# Name of the source flux error field to use.
config.processCcd.charImage.measurePsf.starSelector['references'].signalToNoise.errField='flux_err'

# Select objects with value greater than this
config.processCcd.charImage.measurePsf.starSelector['references'].magError.minimum=None

# Select objects with value less than this
config.processCcd.charImage.measurePsf.starSelector['references'].magError.maximum=None

# Name of the source flux error field to use.
config.processCcd.charImage.measurePsf.starSelector['references'].magError.magErrField='mag_err'

config.processCcd.charImage.measurePsf.starSelector['references'].colorLimits={}
# specify the minimum psfFlux for good Psf Candidates
config.processCcd.charImage.measurePsf.starSelector['objectSize'].fluxMin=4000.0

# specify the maximum psfFlux for good Psf Candidates (ignored if == 0)
config.processCcd.charImage.measurePsf.starSelector['objectSize'].fluxMax=0.0

# minimum width to include in histogram
config.processCcd.charImage.measurePsf.starSelector['objectSize'].widthMin=0.9

# maximum width to include in histogram
config.processCcd.charImage.measurePsf.starSelector['objectSize'].widthMax=10.0

# Name of field in Source to use for flux measurement
config.processCcd.charImage.measurePsf.starSelector['objectSize'].sourceFluxField='base_PsfFlux_flux'

# Standard deviation of width allowed to be interpreted as good stars
config.processCcd.charImage.measurePsf.starSelector['objectSize'].widthStdAllowed=0.15

# Keep objects within this many sigma of cluster 0's median
config.processCcd.charImage.measurePsf.starSelector['objectSize'].nSigmaClip=2.0

# List of flags which cause a source to be rejected as bad
config.processCcd.charImage.measurePsf.starSelector['objectSize'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter', 'base_PixelFlags_flag_bad', 'base_PixelFlags_flag_interpolated']

# Name of a flag field that is True for Sources that should be used.
config.processCcd.charImage.measurePsf.starSelector['flagged'].field='calib_psfUsed'

# List of flags which cause a source to be rejected as bad
config.processCcd.charImage.measurePsf.starSelector['astrometry'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter', 'base_PixelFlags_flag_bad']

# Type of source flux; typically one of Ap or Psf
config.processCcd.charImage.measurePsf.starSelector['astrometry'].sourceFluxType='Ap'

# Minimum allowed signal-to-noise ratio for sources used for matching (in the flux specified by sourceFluxType); <= 0 for no limit
config.processCcd.charImage.measurePsf.starSelector['astrometry'].minSnr=10.0

# Type of source flux; typically one of Ap or Psf
config.processCcd.charImage.measurePsf.starSelector['matcher'].sourceFluxType='Ap'

# Minimum allowed signal-to-noise ratio for sources used for matching (in the flux specified by sourceFluxType); <= 0 for no limit
config.processCcd.charImage.measurePsf.starSelector['matcher'].minSnr=40.0

# Type of source flux; typically one of Ap or Psf
config.processCcd.charImage.measurePsf.starSelector['matcherPessimistic'].sourceFluxType='Ap'

# Minimum allowed signal-to-noise ratio for sources used for matching (in the flux specified by sourceFluxType); <= 0 for no limit
config.processCcd.charImage.measurePsf.starSelector['matcherPessimistic'].minSnr=40.0

# specify the minimum psfFlux for good Psf Candidates
config.processCcd.charImage.measurePsf.starSelector['catalog'].fluxLim=0.0

# specify the maximum psfFlux for good Psf Candidates (ignored if == 0)
config.processCcd.charImage.measurePsf.starSelector['catalog'].fluxMax=0.0

# List of flags which cause a source to be rejected as bad
config.processCcd.charImage.measurePsf.starSelector['catalog'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter']

config.processCcd.charImage.measurePsf.starSelector.name='objectSize'
# size of the kernel to create
config.processCcd.charImage.measurePsf.makePsfCandidates.kernelSize=21

# number of pixels to ignore around the edge of PSF candidate postage stamps
config.processCcd.charImage.measurePsf.makePsfCandidates.borderWidth=0

# radius of the kernel to create, relative to the square root of the stellar quadrupole moments
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].kernelSize=10.0

# Minimum radius of the kernel
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].kernelSizeMin=25

# Maximum radius of the kernel
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].kernelSizeMax=45

# Use non-linear fitter for spatial variation of Kernel
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].nonLinearSpatialFit=False

# number of eigen components for PSF kernel creation
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].nEigenComponents=4

# specify spatial order for PSF kernel creation
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].spatialOrder=2

# size of cell used to determine PSF (pixels, column direction)
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].sizeCellX=256

# size of cell used to determine PSF (pixels, row direction)
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].sizeCellY=256

# number of stars per psf cell for PSF kernel creation
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].nStarPerCell=3

# Number of pixels to ignore around the edge of PSF candidate postage stamps
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].borderWidth=0

# number of stars per psf Cell for spatial fitting
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].nStarPerCellSpatialFit=5

# Should each PSF candidate be given the same weight, independent of magnitude?
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].constantWeight=True

# number of iterations of PSF candidate star list
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].nIterForPsf=3

# tolerance of spatial fitting
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].tolerance=0.01

# floor for variance is lam*data
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].lam=0.05

# for psf candidate evaluation
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].reducedChi2ForPsfCandidates=2.0

# Rejection threshold (stdev) for candidates based on spatial fit
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].spatialReject=3.0

# Threshold (stdev) for rejecting extraneous pixels around candidate; applied if positive
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].pixelThreshold=0.0

# Reject candidates that are blended?
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].doRejectBlends=False

# Mask blends in image?
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].doMaskBlends=True

# radius of the kernel to create, relative to the square root of the stellar quadrupole moments
config.processCcd.charImage.measurePsf.psfDeterminer['psfex'].kernelSize=81.0

# Minimum radius of the kernel
config.processCcd.charImage.measurePsf.psfDeterminer['psfex'].kernelSizeMin=25

# Maximum radius of the kernel
config.processCcd.charImage.measurePsf.psfDeterminer['psfex'].kernelSizeMax=45

# number of eigen components for PSF kernel creation
config.processCcd.charImage.measurePsf.psfDeterminer['psfex']._PsfexPsfDeterminerConfig__nEigenComponents=4

# specify spatial order for PSF kernel creation
config.processCcd.charImage.measurePsf.psfDeterminer['psfex'].spatialOrder=2

# size of cell used to determine PSF (pixels, column direction)
config.processCcd.charImage.measurePsf.psfDeterminer['psfex'].sizeCellX=256

# size of cell used to determine PSF (pixels, row direction)
config.processCcd.charImage.measurePsf.psfDeterminer['psfex'].sizeCellY=256

# number of stars per psf cell for PSF kernel creation
config.processCcd.charImage.measurePsf.psfDeterminer['psfex']._PsfexPsfDeterminerConfig__nStarPerCell=3

# Resolution of the internal PSF model relative to the pixel size; e.g. 0.5 is equal to 2x oversampling
config.processCcd.charImage.measurePsf.psfDeterminer['psfex'].samplingSize=0.5

# List of mask bits which cause a source to be rejected as bad N.b. INTRP is used specially in PsfCandidateSet; it means "Contaminated by neighbour"
config.processCcd.charImage.measurePsf.psfDeterminer['psfex'].badMaskBits=['INTRP', 'SAT']

# BASIS value given to psfex.  PIXEL_AUTO will use the requested samplingSize only if the FWHM < 3 pixels.  Otherwise, it will use samplingSize=1.  PIXEL will always use the requested samplingSize
config.processCcd.charImage.measurePsf.psfDeterminer['psfex'].psfexBasis='PIXEL_AUTO'

# Number of pixels to ignore around the edge of PSF candidate postage stamps
config.processCcd.charImage.measurePsf.psfDeterminer['psfex']._PsfexPsfDeterminerConfig__borderWidth=0

# number of stars per psf Cell for spatial fitting
config.processCcd.charImage.measurePsf.psfDeterminer['psfex']._PsfexPsfDeterminerConfig__nStarPerCellSpatialFit=5

# Should each PSF candidate be given the same weight, independent of magnitude?
config.processCcd.charImage.measurePsf.psfDeterminer['psfex']._PsfexPsfDeterminerConfig__constantWeight=True

# number of iterations of PSF candidate star list
config.processCcd.charImage.measurePsf.psfDeterminer['psfex']._PsfexPsfDeterminerConfig__nIterForPsf=3

# tolerance of spatial fitting
config.processCcd.charImage.measurePsf.psfDeterminer['psfex'].tolerance=0.01

# floor for variance is lam*data
config.processCcd.charImage.measurePsf.psfDeterminer['psfex'].lam=0.05

# for psf candidate evaluation
config.processCcd.charImage.measurePsf.psfDeterminer['psfex'].reducedChi2ForPsfCandidates=2.0

# Rejection threshold (stdev) for candidates based on spatial fit
config.processCcd.charImage.measurePsf.psfDeterminer['psfex'].spatialReject=3.0

# Should PSFEX be permitted to recentroid PSF candidates?
config.processCcd.charImage.measurePsf.psfDeterminer['psfex'].recentroid=False

config.processCcd.charImage.measurePsf.psfDeterminer.name='psfex'
# Fraction of candidates to reserve from fitting; none if <= 0
config.processCcd.charImage.measurePsf.reserve.fraction=0.2

# This number will be added to the exposure ID to set the random seed for reserving candidates
config.processCcd.charImage.measurePsf.reserve.seed=1

# Interpolate over defects? (ignored unless you provide a list of defects)
config.processCcd.charImage.repair.doInterpolate=True

# Find and mask out cosmic rays?
config.processCcd.charImage.repair.doCosmicRay=True

# maximum number of contaminated pixels
config.processCcd.charImage.repair.cosmicray.nCrPixelMax=1000000

# CRs must be > this many sky-sig above sky
config.processCcd.charImage.repair.cosmicray.minSigma=6.0

# CRs must have > this many DN (== electrons/gain) in initial detection
config.processCcd.charImage.repair.cosmicray.min_DN=150.0

# used in condition 3 for CR; see CR.cc code
config.processCcd.charImage.repair.cosmicray.cond3_fac=2.5

# used in condition 3 for CR; see CR.cc code
config.processCcd.charImage.repair.cosmicray.cond3_fac2=0.4

# number of times to look for contaminated pixels near known CR pixels
config.processCcd.charImage.repair.cosmicray.niteration=3

# Don't interpolate over CR pixels
config.processCcd.charImage.repair.cosmicray.keepCRs=False

# type of statistic to use for grid points
config.processCcd.charImage.repair.cosmicray.background.statisticsProperty='MEDIAN'

# behaviour if there are too few points in grid for requested interpolation style
config.processCcd.charImage.repair.cosmicray.background.undersampleStyle='REDUCE_INTERP_ORDER'

# how large a region of the sky should be used for each background point
config.processCcd.charImage.repair.cosmicray.background.binSize=100000

# Sky region size to be used for each background point in X direction. If 0, the binSize config is used.
config.processCcd.charImage.repair.cosmicray.background.binSizeX=0

# Sky region size to be used for each background point in Y direction. If 0, the binSize config is used.
config.processCcd.charImage.repair.cosmicray.background.binSizeY=0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
config.processCcd.charImage.repair.cosmicray.background.algorithm='AKIMA_SPLINE'

# Names of mask planes to ignore while estimating the background
config.processCcd.charImage.repair.cosmicray.background.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# Ignore NaNs when estimating the background
config.processCcd.charImage.repair.cosmicray.background.isNanSafe=False

# Use Approximate (Chebyshev) to model background.
config.processCcd.charImage.repair.cosmicray.background.useApprox=False

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.processCcd.charImage.repair.cosmicray.background.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.processCcd.charImage.repair.cosmicray.background.approxOrderY=-1

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.processCcd.charImage.repair.cosmicray.background.weighting=True

# Kernel size (width and height) (pixels); if None then sizeFactor is used
config.processCcd.charImage.repair.interp.modelPsf.size=None

# Kernel size as a factor of fwhm (dimensionless); size = sizeFactor * fwhm; ignored if size is not None
config.processCcd.charImage.repair.interp.modelPsf.sizeFactor=3.0

# Minimum kernel size if using sizeFactor (pixels); ignored if size is not None
config.processCcd.charImage.repair.interp.modelPsf.minSize=5

# Maximum kernel size if using sizeFactor (pixels); ignored if size is not None
config.processCcd.charImage.repair.interp.modelPsf.maxSize=None

# Default FWHM of Gaussian model of core of star (pixels)
config.processCcd.charImage.repair.interp.modelPsf.defaultFwhm=3.0

# Add a Gaussian to represent wings?
config.processCcd.charImage.repair.interp.modelPsf.addWing=True

# wing width, as a multiple of core width (dimensionless); ignored if addWing false
config.processCcd.charImage.repair.interp.modelPsf.wingFwhmFactor=2.5

# wing amplitude, as a multiple of core amplitude (dimensionless); ignored if addWing false
config.processCcd.charImage.repair.interp.modelPsf.wingAmplitude=0.1

# Smoothly taper to the fallback value at the edge of the image?
config.processCcd.charImage.repair.interp.useFallbackValueAtEdge=True

# Type of statistic to calculate edge fallbackValue for interpolation
config.processCcd.charImage.repair.interp.fallbackValueType='MEANCLIP'

# If fallbackValueType is 'USER' then use this as the fallbackValue; ignored otherwise
config.processCcd.charImage.repair.interp.fallbackUserValue=0.0

# Allow negative values for egde interpolation fallbackValue?  If False, set fallbackValue to max(fallbackValue, 0.0)
config.processCcd.charImage.repair.interp.negativeFallbackAllowed=True

# Strictness of Astropy unit compatibility check, can be 'raise', 'warn' or 'silent'
config.processCcd.charImage.checkUnitsParseStrict='raise'

# Perform calibration?
config.processCcd.doCalibrate=True

# Save calibration results?
config.processCcd.calibrate.doWrite=True

# Include HeavyFootprint data in source table? If false then heavy footprints are saved as normal footprints, which saves some space
config.processCcd.calibrate.doWriteHeavyFootprintsInSources=True

# Write reference matches (ignored if doWrite false)?
config.processCcd.calibrate.doWriteMatches=True

# Write reference matches in denormalized format? This format uses more disk space, but is more convenient to read. Ignored if doWriteMatches=False or doWrite=False.
config.processCcd.calibrate.doWriteMatchesDenormalized=True

# Perform astrometric calibration?
config.processCcd.calibrate.doAstrometry=True

# Padding to add to 4 all edges of the bounding box (pixels)
config.processCcd.calibrate.astromRefObjLoader.pixelMargin=300

# Default reference catalog filter to use if filter not specified in exposure; if blank then filter must be specified in exposure
config.processCcd.calibrate.astromRefObjLoader.defaultFilter=''

# Mapping of camera filter name: reference catalog filter name; each reference filter must exist
config.processCcd.calibrate.astromRefObjLoader.filterMap={'y': 'z'}

# Padding to add to 4 all edges of the bounding box (pixels)
config.processCcd.calibrate.photoRefObjLoader.pixelMargin=300

# Default reference catalog filter to use if filter not specified in exposure; if blank then filter must be specified in exposure
config.processCcd.calibrate.photoRefObjLoader.defaultFilter=''

# Mapping of camera filter name: reference catalog filter name; each reference filter must exist
config.processCcd.calibrate.photoRefObjLoader.filterMap={'y': 'z'}

# Maximum separation between reference objects and sources beyond which they will not be considered a match (arcsec)
config.processCcd.calibrate.astrometry.matcher.maxMatchDistArcSec=2.0

# Number of bright stars to use
config.processCcd.calibrate.astrometry.matcher.numBrightStars=150

# Minimum number of matched pairs; see also minFracMatchedPairs
config.processCcd.calibrate.astrometry.matcher.minMatchedPairs=30

# Minimum number of matched pairs as a fraction of the smaller of the number of reference stars or the number of good sources; the actual minimum is the smaller of this value or minMatchedPairs
config.processCcd.calibrate.astrometry.matcher.minFracMatchedPairs=0.3

# Maximum allowed shift of WCS, due to matching (pixel). When changing this value, the LoadReferenceObjectsConfig.pixelMargin should also be updated.
config.processCcd.calibrate.astrometry.matcher.maxOffsetPix=250

# Rotation angle allowed between sources and position reference objects (degrees)
config.processCcd.calibrate.astrometry.matcher.maxRotationDeg=1.145916

# Allowed non-perpendicularity of x and y (degree)
config.processCcd.calibrate.astrometry.matcher.allowedNonperpDeg=0.2

# number of points to define a shape for matching
config.processCcd.calibrate.astrometry.matcher.numPointsForShape=6

# maximum determinant of linear transformation matrix for a usable solution
config.processCcd.calibrate.astrometry.matcher.maxDeterminant=0.02

# Apply flux limit?
config.processCcd.calibrate.astrometry.matcher.sourceSelector['science'].doFluxLimit=False

# Apply flag limitation?
config.processCcd.calibrate.astrometry.matcher.sourceSelector['science'].doFlags=False

# Apply unresolved limitation?
config.processCcd.calibrate.astrometry.matcher.sourceSelector['science'].doUnresolved=False

# Apply signal-to-noise limit?
config.processCcd.calibrate.astrometry.matcher.sourceSelector['science'].doSignalToNoise=False

# Apply isolated limitation?
config.processCcd.calibrate.astrometry.matcher.sourceSelector['science'].doIsolated=False

# Select objects with value greater than this
config.processCcd.calibrate.astrometry.matcher.sourceSelector['science'].fluxLimit.minimum=None

# Select objects with value less than this
config.processCcd.calibrate.astrometry.matcher.sourceSelector['science'].fluxLimit.maximum=None

# Name of the source flux field to use.
config.processCcd.calibrate.astrometry.matcher.sourceSelector['science'].fluxLimit.fluxField='slot_CalibFlux_flux'

# List of source flag fields that must be set for a source to be used.
config.processCcd.calibrate.astrometry.matcher.sourceSelector['science'].flags.good=[]

# List of source flag fields that must NOT be set for a source to be used.
config.processCcd.calibrate.astrometry.matcher.sourceSelector['science'].flags.bad=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_saturated', 'base_PsfFlux_flags']

# Select objects with value greater than this
config.processCcd.calibrate.astrometry.matcher.sourceSelector['science'].unresolved.minimum=None

# Select objects with value less than this
config.processCcd.calibrate.astrometry.matcher.sourceSelector['science'].unresolved.maximum=None

# Name of column for star/galaxy separation
config.processCcd.calibrate.astrometry.matcher.sourceSelector['science'].unresolved.name='base_ClassificationExtendedness_value'

# Select objects with value greater than this
config.processCcd.calibrate.astrometry.matcher.sourceSelector['science'].signalToNoise.minimum=None

# Select objects with value less than this
config.processCcd.calibrate.astrometry.matcher.sourceSelector['science'].signalToNoise.maximum=None

# Name of the source flux field to use.
config.processCcd.calibrate.astrometry.matcher.sourceSelector['science'].signalToNoise.fluxField='base_PsfFlux_flux'

# Name of the source flux error field to use.
config.processCcd.calibrate.astrometry.matcher.sourceSelector['science'].signalToNoise.errField='base_PsfFlux_fluxSigma'

# Name of column for parent
config.processCcd.calibrate.astrometry.matcher.sourceSelector['science'].isolated.parentName='parent'

# Name of column for nChild
config.processCcd.calibrate.astrometry.matcher.sourceSelector['science'].isolated.nChildName='deblend_nChild'

# Apply magnitude limit?
config.processCcd.calibrate.astrometry.matcher.sourceSelector['references'].doMagLimit=False

# Apply flag limitation?
config.processCcd.calibrate.astrometry.matcher.sourceSelector['references'].doFlags=False

# Apply signal-to-noise limit?
config.processCcd.calibrate.astrometry.matcher.sourceSelector['references'].doSignalToNoise=False

# Apply magnitude error limit?
config.processCcd.calibrate.astrometry.matcher.sourceSelector['references'].doMagError=False

# Select objects with value greater than this
config.processCcd.calibrate.astrometry.matcher.sourceSelector['references'].magLimit.minimum=None

# Select objects with value less than this
config.processCcd.calibrate.astrometry.matcher.sourceSelector['references'].magLimit.maximum=None

# Name of the source flux field to use.
config.processCcd.calibrate.astrometry.matcher.sourceSelector['references'].magLimit.fluxField='flux'

# List of source flag fields that must be set for a source to be used.
config.processCcd.calibrate.astrometry.matcher.sourceSelector['references'].flags.good=[]

# List of source flag fields that must NOT be set for a source to be used.
config.processCcd.calibrate.astrometry.matcher.sourceSelector['references'].flags.bad=[]

# Select objects with value greater than this
config.processCcd.calibrate.astrometry.matcher.sourceSelector['references'].signalToNoise.minimum=None

# Select objects with value less than this
config.processCcd.calibrate.astrometry.matcher.sourceSelector['references'].signalToNoise.maximum=None

# Name of the source flux field to use.
config.processCcd.calibrate.astrometry.matcher.sourceSelector['references'].signalToNoise.fluxField='flux'

# Name of the source flux error field to use.
config.processCcd.calibrate.astrometry.matcher.sourceSelector['references'].signalToNoise.errField='flux_err'

# Select objects with value greater than this
config.processCcd.calibrate.astrometry.matcher.sourceSelector['references'].magError.minimum=None

# Select objects with value less than this
config.processCcd.calibrate.astrometry.matcher.sourceSelector['references'].magError.maximum=None

# Name of the source flux error field to use.
config.processCcd.calibrate.astrometry.matcher.sourceSelector['references'].magError.magErrField='mag_err'

config.processCcd.calibrate.astrometry.matcher.sourceSelector['references'].colorLimits={}
# specify the minimum psfFlux for good Psf Candidates
config.processCcd.calibrate.astrometry.matcher.sourceSelector['objectSize'].fluxMin=12500.0

# specify the maximum psfFlux for good Psf Candidates (ignored if == 0)
config.processCcd.calibrate.astrometry.matcher.sourceSelector['objectSize'].fluxMax=0.0

# minimum width to include in histogram
config.processCcd.calibrate.astrometry.matcher.sourceSelector['objectSize'].widthMin=0.0

# maximum width to include in histogram
config.processCcd.calibrate.astrometry.matcher.sourceSelector['objectSize'].widthMax=10.0

# Name of field in Source to use for flux measurement
config.processCcd.calibrate.astrometry.matcher.sourceSelector['objectSize'].sourceFluxField='base_GaussianFlux_flux'

# Standard deviation of width allowed to be interpreted as good stars
config.processCcd.calibrate.astrometry.matcher.sourceSelector['objectSize'].widthStdAllowed=0.15

# Keep objects within this many sigma of cluster 0's median
config.processCcd.calibrate.astrometry.matcher.sourceSelector['objectSize'].nSigmaClip=2.0

# List of flags which cause a source to be rejected as bad
config.processCcd.calibrate.astrometry.matcher.sourceSelector['objectSize'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter', 'base_PixelFlags_flag_bad', 'base_PixelFlags_flag_interpolated']

# Name of a flag field that is True for Sources that should be used.
config.processCcd.calibrate.astrometry.matcher.sourceSelector['flagged'].field='calib_psfUsed'

# List of flags which cause a source to be rejected as bad
config.processCcd.calibrate.astrometry.matcher.sourceSelector['astrometry'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter', 'base_PixelFlags_flag_bad']

# Type of source flux; typically one of Ap or Psf
config.processCcd.calibrate.astrometry.matcher.sourceSelector['astrometry'].sourceFluxType='Ap'

# Minimum allowed signal-to-noise ratio for sources used for matching (in the flux specified by sourceFluxType); <= 0 for no limit
config.processCcd.calibrate.astrometry.matcher.sourceSelector['astrometry'].minSnr=10.0

# Type of source flux; typically one of Ap or Psf
config.processCcd.calibrate.astrometry.matcher.sourceSelector['matcher'].sourceFluxType='Psf'

# Minimum allowed signal-to-noise ratio for sources used for matching (in the flux specified by sourceFluxType); <= 0 for no limit
config.processCcd.calibrate.astrometry.matcher.sourceSelector['matcher'].minSnr=40.0

# Type of source flux; typically one of Ap or Psf
config.processCcd.calibrate.astrometry.matcher.sourceSelector['matcherPessimistic'].sourceFluxType='Ap'

# Minimum allowed signal-to-noise ratio for sources used for matching (in the flux specified by sourceFluxType); <= 0 for no limit
config.processCcd.calibrate.astrometry.matcher.sourceSelector['matcherPessimistic'].minSnr=40.0

# specify the minimum psfFlux for good Psf Candidates
config.processCcd.calibrate.astrometry.matcher.sourceSelector['catalog'].fluxLim=0.0

# specify the maximum psfFlux for good Psf Candidates (ignored if == 0)
config.processCcd.calibrate.astrometry.matcher.sourceSelector['catalog'].fluxMax=0.0

# List of flags which cause a source to be rejected as bad
config.processCcd.calibrate.astrometry.matcher.sourceSelector['catalog'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter']

config.processCcd.calibrate.astrometry.matcher.sourceSelector.name='matcher'
# the maximum match distance is set to  mean_match_distance + matchDistanceSigma*std_dev_match_distance; ignored if not fitting a WCS
config.processCcd.calibrate.astrometry.matchDistanceSigma=2.0

# Apply flux limit?
config.processCcd.calibrate.astrometry.sourceSelection.doFluxLimit=False

# Apply flag limitation?
config.processCcd.calibrate.astrometry.sourceSelection.doFlags=False

# Apply unresolved limitation?
config.processCcd.calibrate.astrometry.sourceSelection.doUnresolved=False

# Apply signal-to-noise limit?
config.processCcd.calibrate.astrometry.sourceSelection.doSignalToNoise=False

# Apply isolated limitation?
config.processCcd.calibrate.astrometry.sourceSelection.doIsolated=False

# Select objects with value greater than this
config.processCcd.calibrate.astrometry.sourceSelection.fluxLimit.minimum=None

# Select objects with value less than this
config.processCcd.calibrate.astrometry.sourceSelection.fluxLimit.maximum=None

# Name of the source flux field to use.
config.processCcd.calibrate.astrometry.sourceSelection.fluxLimit.fluxField='slot_CalibFlux_flux'

# List of source flag fields that must be set for a source to be used.
config.processCcd.calibrate.astrometry.sourceSelection.flags.good=[]

# List of source flag fields that must NOT be set for a source to be used.
config.processCcd.calibrate.astrometry.sourceSelection.flags.bad=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_saturated', 'base_PsfFlux_flags']

# Select objects with value greater than this
config.processCcd.calibrate.astrometry.sourceSelection.unresolved.minimum=None

# Select objects with value less than this
config.processCcd.calibrate.astrometry.sourceSelection.unresolved.maximum=None

# Name of column for star/galaxy separation
config.processCcd.calibrate.astrometry.sourceSelection.unresolved.name='base_ClassificationExtendedness_value'

# Select objects with value greater than this
config.processCcd.calibrate.astrometry.sourceSelection.signalToNoise.minimum=None

# Select objects with value less than this
config.processCcd.calibrate.astrometry.sourceSelection.signalToNoise.maximum=None

# Name of the source flux field to use.
config.processCcd.calibrate.astrometry.sourceSelection.signalToNoise.fluxField='base_PsfFlux_flux'

# Name of the source flux error field to use.
config.processCcd.calibrate.astrometry.sourceSelection.signalToNoise.errField='base_PsfFlux_fluxSigma'

# Name of column for parent
config.processCcd.calibrate.astrometry.sourceSelection.isolated.parentName='parent'

# Name of column for nChild
config.processCcd.calibrate.astrometry.sourceSelection.isolated.nChildName='deblend_nChild'

# Apply magnitude limit?
config.processCcd.calibrate.astrometry.referenceSelection.doMagLimit=False

# Apply flag limitation?
config.processCcd.calibrate.astrometry.referenceSelection.doFlags=False

# Apply signal-to-noise limit?
config.processCcd.calibrate.astrometry.referenceSelection.doSignalToNoise=False

# Apply magnitude error limit?
config.processCcd.calibrate.astrometry.referenceSelection.doMagError=False

# Select objects with value greater than this
config.processCcd.calibrate.astrometry.referenceSelection.magLimit.minimum=None

# Select objects with value less than this
config.processCcd.calibrate.astrometry.referenceSelection.magLimit.maximum=None

# Name of the source flux field to use.
config.processCcd.calibrate.astrometry.referenceSelection.magLimit.fluxField='flux'

# List of source flag fields that must be set for a source to be used.
config.processCcd.calibrate.astrometry.referenceSelection.flags.good=[]

# List of source flag fields that must NOT be set for a source to be used.
config.processCcd.calibrate.astrometry.referenceSelection.flags.bad=[]

# Select objects with value greater than this
config.processCcd.calibrate.astrometry.referenceSelection.signalToNoise.minimum=None

# Select objects with value less than this
config.processCcd.calibrate.astrometry.referenceSelection.signalToNoise.maximum=None

# Name of the source flux field to use.
config.processCcd.calibrate.astrometry.referenceSelection.signalToNoise.fluxField='flux'

# Name of the source flux error field to use.
config.processCcd.calibrate.astrometry.referenceSelection.signalToNoise.errField='flux_err'

# Select objects with value greater than this
config.processCcd.calibrate.astrometry.referenceSelection.magError.minimum=None

# Select objects with value less than this
config.processCcd.calibrate.astrometry.referenceSelection.magError.maximum=None

# Name of the source flux error field to use.
config.processCcd.calibrate.astrometry.referenceSelection.magError.magErrField='mag_err'

config.processCcd.calibrate.astrometry.referenceSelection.colorLimits={}
# order of SIP polynomial
config.processCcd.calibrate.astrometry.wcsFitter.order=3

# number of iterations of fitter (which fits X and Y separately, and so benefits from a few iterations
config.processCcd.calibrate.astrometry.wcsFitter.numIter=3

# number of rejection iterations
config.processCcd.calibrate.astrometry.wcsFitter.numRejIter=3

# Number of standard deviations for clipping level
config.processCcd.calibrate.astrometry.wcsFitter.rejSigma=3.0

# maximum median scatter of a WCS fit beyond which the fit fails (arcsec); be generous, as this is only intended to catch catastrophic failures
config.processCcd.calibrate.astrometry.wcsFitter.maxScatterArcsec=10.0

# If True then load reference objects and match sources but do not fit a WCS;  this simply controls whether 'run' calls 'solve' or 'loadAndMatch'
config.processCcd.calibrate.astrometry.forceKnownWcs=False

# maximum number of iterations of match sources and fit WCSignored if not fitting a WCS
config.processCcd.calibrate.astrometry.maxIter=3

# the match distance below which further iteration is pointless (arcsec); ignored if not fitting a WCS
config.processCcd.calibrate.astrometry.minMatchDistanceArcSec=0.001

# Raise an exception if astrometry fails? Ignored if doAstrometry false.
config.processCcd.calibrate.requireAstrometry=True

# Perform phometric calibration?
config.processCcd.calibrate.doPhotoCal=True

# Raise an exception if photoCal fails? Ignored if doPhotoCal false.
config.processCcd.calibrate.requirePhotoCal=True

# Matching radius, arcsec
config.processCcd.calibrate.photoCal.match.matchRadius=0.25

# Apply flux limit?
config.processCcd.calibrate.photoCal.match.sourceSelection.doFluxLimit=False

# Apply flag limitation?
config.processCcd.calibrate.photoCal.match.sourceSelection.doFlags=True

# Apply unresolved limitation?
config.processCcd.calibrate.photoCal.match.sourceSelection.doUnresolved=True

# Apply signal-to-noise limit?
config.processCcd.calibrate.photoCal.match.sourceSelection.doSignalToNoise=False

# Apply isolated limitation?
config.processCcd.calibrate.photoCal.match.sourceSelection.doIsolated=False

# Select objects with value greater than this
config.processCcd.calibrate.photoCal.match.sourceSelection.fluxLimit.minimum=None

# Select objects with value less than this
config.processCcd.calibrate.photoCal.match.sourceSelection.fluxLimit.maximum=None

# Name of the source flux field to use.
config.processCcd.calibrate.photoCal.match.sourceSelection.fluxLimit.fluxField='slot_CalibFlux_flux'

# List of source flag fields that must be set for a source to be used.
config.processCcd.calibrate.photoCal.match.sourceSelection.flags.good=[]

# List of source flag fields that must NOT be set for a source to be used.
config.processCcd.calibrate.photoCal.match.sourceSelection.flags.bad=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolated', 'base_PixelFlags_flag_saturated']

# Select objects with value greater than this
config.processCcd.calibrate.photoCal.match.sourceSelection.unresolved.minimum=None

# Select objects with value less than this
config.processCcd.calibrate.photoCal.match.sourceSelection.unresolved.maximum=None

# Name of column for star/galaxy separation
config.processCcd.calibrate.photoCal.match.sourceSelection.unresolved.name='base_ClassificationExtendedness_value'

# Select objects with value greater than this
config.processCcd.calibrate.photoCal.match.sourceSelection.signalToNoise.minimum=None

# Select objects with value less than this
config.processCcd.calibrate.photoCal.match.sourceSelection.signalToNoise.maximum=None

# Name of the source flux field to use.
config.processCcd.calibrate.photoCal.match.sourceSelection.signalToNoise.fluxField='base_PsfFlux_flux'

# Name of the source flux error field to use.
config.processCcd.calibrate.photoCal.match.sourceSelection.signalToNoise.errField='base_PsfFlux_fluxSigma'

# Name of column for parent
config.processCcd.calibrate.photoCal.match.sourceSelection.isolated.parentName='parent'

# Name of column for nChild
config.processCcd.calibrate.photoCal.match.sourceSelection.isolated.nChildName='deblend_nChild'

# Apply magnitude limit?
config.processCcd.calibrate.photoCal.match.referenceSelection.doMagLimit=True

# Apply flag limitation?
config.processCcd.calibrate.photoCal.match.referenceSelection.doFlags=False

# Apply signal-to-noise limit?
config.processCcd.calibrate.photoCal.match.referenceSelection.doSignalToNoise=False

# Apply magnitude error limit?
config.processCcd.calibrate.photoCal.match.referenceSelection.doMagError=False

# Select objects with value greater than this
config.processCcd.calibrate.photoCal.match.referenceSelection.magLimit.minimum=None

# Select objects with value less than this
config.processCcd.calibrate.photoCal.match.referenceSelection.magLimit.maximum=22.0

# Name of the source flux field to use.
config.processCcd.calibrate.photoCal.match.referenceSelection.magLimit.fluxField='i_flux'

# List of source flag fields that must be set for a source to be used.
config.processCcd.calibrate.photoCal.match.referenceSelection.flags.good=[]

# List of source flag fields that must NOT be set for a source to be used.
config.processCcd.calibrate.photoCal.match.referenceSelection.flags.bad=[]

# Select objects with value greater than this
config.processCcd.calibrate.photoCal.match.referenceSelection.signalToNoise.minimum=None

# Select objects with value less than this
config.processCcd.calibrate.photoCal.match.referenceSelection.signalToNoise.maximum=None

# Name of the source flux field to use.
config.processCcd.calibrate.photoCal.match.referenceSelection.signalToNoise.fluxField='flux'

# Name of the source flux error field to use.
config.processCcd.calibrate.photoCal.match.referenceSelection.signalToNoise.errField='flux_err'

# Select objects with value greater than this
config.processCcd.calibrate.photoCal.match.referenceSelection.magError.minimum=None

# Select objects with value less than this
config.processCcd.calibrate.photoCal.match.referenceSelection.magError.maximum=None

# Name of the source flux error field to use.
config.processCcd.calibrate.photoCal.match.referenceSelection.magError.magErrField='mag_err'

config.processCcd.calibrate.photoCal.match.referenceSelection.colorLimits={}
config.processCcd.calibrate.photoCal.match.referenceSelection.colorLimits['g-r']=lsst.meas.algorithms.sourceSelector.ColorLimit()
# Select objects with value greater than this
config.processCcd.calibrate.photoCal.match.referenceSelection.colorLimits['g-r'].minimum=0.0

# Select objects with value less than this
config.processCcd.calibrate.photoCal.match.referenceSelection.colorLimits['g-r'].maximum=None

# Name of column with primary flux measurement
config.processCcd.calibrate.photoCal.match.referenceSelection.colorLimits['g-r'].primary='g_flux'

# Name of column with secondary flux measurement
config.processCcd.calibrate.photoCal.match.referenceSelection.colorLimits['g-r'].secondary='r_flux'

config.processCcd.calibrate.photoCal.match.referenceSelection.colorLimits['r-i']=lsst.meas.algorithms.sourceSelector.ColorLimit()
# Select objects with value greater than this
config.processCcd.calibrate.photoCal.match.referenceSelection.colorLimits['r-i'].minimum=None

# Select objects with value less than this
config.processCcd.calibrate.photoCal.match.referenceSelection.colorLimits['r-i'].maximum=0.5

# Name of column with primary flux measurement
config.processCcd.calibrate.photoCal.match.referenceSelection.colorLimits['r-i'].primary='r_flux'

# Name of column with secondary flux measurement
config.processCcd.calibrate.photoCal.match.referenceSelection.colorLimits['r-i'].secondary='i_flux'

# Fraction of candidates to reserve from fitting; none if <= 0
config.processCcd.calibrate.photoCal.reserve.fraction=0.0

# This number will be added to the exposure ID to set the random seed for reserving candidates
config.processCcd.calibrate.photoCal.reserve.seed=1

# Name of the source flux field to use.  The associated flag field
# ('<name>_flags') will be implicitly included in badFlags.
config.processCcd.calibrate.photoCal.fluxField='slot_CalibFlux_flux'

# Apply photometric color terms to reference stars? One of:
# None: apply if colorterms and photoCatName are not None;
#       fail if color term data is not available for the specified ref catalog and filter.
# True: always apply colorterms; fail if color term data is not available for the
#       specified reference catalog and filter.
# False: do not apply.
config.processCcd.calibrate.photoCal.applyColorTerms=True

# maximum sigma to use when clipping
config.processCcd.calibrate.photoCal.sigmaMax=0.25

# clip at nSigma
config.processCcd.calibrate.photoCal.nSigma=3.0

# use median instead of mean to compute zeropoint
config.processCcd.calibrate.photoCal.useMedian=True

# number of iterations
config.processCcd.calibrate.photoCal.nIter=20

config.processCcd.calibrate.photoCal.colorterms.data={}
config.processCcd.calibrate.photoCal.colorterms.data['hsc*']=lsst.pipe.tasks.colorterms.ColortermDict()
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data={}
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['g']=lsst.pipe.tasks.colorterms.Colorterm()
# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['g'].primary='g'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['g'].secondary='g'

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['g'].c0=0.0

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['g'].c1=0.0

# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['g'].c2=0.0

config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['r']=lsst.pipe.tasks.colorterms.Colorterm()
# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['r'].primary='r'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['r'].secondary='r'

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['r'].c0=0.0

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['r'].c1=0.0

# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['r'].c2=0.0

config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['i']=lsst.pipe.tasks.colorterms.Colorterm()
# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['i'].primary='i'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['i'].secondary='i'

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['i'].c0=0.0

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['i'].c1=0.0

# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['i'].c2=0.0

config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['z']=lsst.pipe.tasks.colorterms.Colorterm()
# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['z'].primary='z'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['z'].secondary='z'

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['z'].c0=0.0

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['z'].c1=0.0

# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['z'].c2=0.0

config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['y']=lsst.pipe.tasks.colorterms.Colorterm()
# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['y'].primary='y'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['y'].secondary='y'

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['y'].c0=0.0

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['y'].c1=0.0

# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['y'].c2=0.0

config.processCcd.calibrate.photoCal.colorterms.data['sdss*']=lsst.pipe.tasks.colorterms.ColortermDict()
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data={}
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['g']=lsst.pipe.tasks.colorterms.Colorterm()
# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['g'].primary='g'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['g'].secondary='r'

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['g'].c0=-0.00816446

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['g'].c1=-0.08366937

# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['g'].c2=-0.00726883

config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['r']=lsst.pipe.tasks.colorterms.Colorterm()
# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['r'].primary='r'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['r'].secondary='i'

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['r'].c0=0.0023181

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['r'].c1=0.01284177

# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['r'].c2=-0.03068248

config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['r2']=lsst.pipe.tasks.colorterms.Colorterm()
# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['r2'].primary='r'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['r2'].secondary='i'

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['r2'].c0=0.00074087

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['r2'].c1=-0.00830543

# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['r2'].c2=-0.0284842

config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['i']=lsst.pipe.tasks.colorterms.Colorterm()
# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['i'].primary='i'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['i'].secondary='z'

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['i'].c0=0.00130204

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['i'].c1=-0.16922042

# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['i'].c2=-0.01374245

config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['i2']=lsst.pipe.tasks.colorterms.Colorterm()
# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['i2'].primary='i'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['i2'].secondary='z'

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['i2'].c0=0.00124676

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['i2'].c1=-0.20739606

# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['i2'].c2=-0.01067212

config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['z']=lsst.pipe.tasks.colorterms.Colorterm()
# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['z'].primary='z'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['z'].secondary='i'

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['z'].c0=-0.0068062

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['z'].c1=0.01353969

# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['z'].c2=0.01479369

config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['y']=lsst.pipe.tasks.colorterms.Colorterm()
# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['y'].primary='z'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['y'].secondary='i'

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['y'].c0=0.01739708

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['y'].c1=0.35652971

# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['y'].c2=0.00574408

config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['N816']=lsst.pipe.tasks.colorterms.Colorterm()
# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['N816'].primary='i'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['N816'].secondary='z'

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['N816'].c0=0.00927133

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['N816'].c1=-0.63558358

# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['N816'].c2=-0.05474862

config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['N921']=lsst.pipe.tasks.colorterms.Colorterm()
# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['N921'].primary='z'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['N921'].secondary='i'

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['N921'].c0=0.00752972

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['N921'].c1=0.0986353

# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['N921'].c2=-0.05451118

config.processCcd.calibrate.photoCal.colorterms.data['ps1*']=lsst.pipe.tasks.colorterms.ColortermDict()
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data={}
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['g']=lsst.pipe.tasks.colorterms.Colorterm()
# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['g'].primary='g'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['g'].secondary='r'

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['g'].c0=0.005905

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['g'].c1=0.063651

# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['g'].c2=-0.000716

config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['r']=lsst.pipe.tasks.colorterms.Colorterm()
# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['r'].primary='r'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['r'].secondary='i'

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['r'].c0=-0.000294

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['r'].c1=-0.005458

# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['r'].c2=-0.009451

config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['r2']=lsst.pipe.tasks.colorterms.Colorterm()
# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['r2'].primary='r'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['r2'].secondary='i'

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['r2'].c0=0.000118

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['r2'].c1=-0.00279

# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['r2'].c2=-0.014363

config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['i']=lsst.pipe.tasks.colorterms.Colorterm()
# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['i'].primary='i'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['i'].secondary='z'

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['i'].c0=0.000979

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['i'].c1=-0.154608

# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['i'].c2=-0.010429

config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['i2']=lsst.pipe.tasks.colorterms.Colorterm()
# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['i2'].primary='i'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['i2'].secondary='z'

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['i2'].c0=0.001653

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['i2'].c1=-0.206313

# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['i2'].c2=-0.016085

config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['z']=lsst.pipe.tasks.colorterms.Colorterm()
# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['z'].primary='z'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['z'].secondary='y'

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['z'].c0=-0.005585

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['z'].c1=-0.220704

# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['z'].c2=-0.298211

config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['y']=lsst.pipe.tasks.colorterms.Colorterm()
# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['y'].primary='y'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['y'].secondary='z'

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['y'].c0=-0.001952

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['y'].c1=0.19957

# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['y'].c2=0.216821

config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['I945']=lsst.pipe.tasks.colorterms.Colorterm()
# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['I945'].primary='y'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['I945'].secondary='z'

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['I945'].c0=0.005275

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['I945'].c1=-0.194285

# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['I945'].c2=-0.125424

config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N387']=lsst.pipe.tasks.colorterms.Colorterm()
# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N387'].primary='g'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N387'].secondary='r'

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N387'].c0=0.427879

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N387'].c1=1.869068

# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N387'].c2=0.54058

config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N468']=lsst.pipe.tasks.colorterms.Colorterm()
# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N468'].primary='g'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N468'].secondary='r'

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N468'].c0=-0.04224

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N468'].c1=0.121756

# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N468'].c2=0.027599

config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N515']=lsst.pipe.tasks.colorterms.Colorterm()
# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N515'].primary='g'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N515'].secondary='r'

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N515'].c0=-0.021913

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N515'].c1=-0.253159

# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N515'].c2=0.151553

config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N527']=lsst.pipe.tasks.colorterms.Colorterm()
# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N527'].primary='g'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N527'].secondary='r'

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N527'].c0=-0.020641

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N527'].c1=-0.366167

# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N527'].c2=0.038497

config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N656']=lsst.pipe.tasks.colorterms.Colorterm()
# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N656'].primary='r'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N656'].secondary='i'

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N656'].c0=0.035658

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N656'].c1=-0.512071

# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N656'].c2=0.042824

config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N718']=lsst.pipe.tasks.colorterms.Colorterm()
# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N718'].primary='i'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N718'].secondary='r'

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N718'].c0=-0.016294

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N718'].c1=-0.233139

# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N718'].c2=0.252505

config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N816']=lsst.pipe.tasks.colorterms.Colorterm()
# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N816'].primary='i'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N816'].secondary='z'

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N816'].c0=0.013806

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N816'].c1=-0.717681

# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N816'].c2=0.049289

config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N921']=lsst.pipe.tasks.colorterms.Colorterm()
# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N921'].primary='z'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N921'].secondary='y'

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N921'].c0=0.002039

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N921'].c1=-0.477412

# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N921'].c2=-0.492151

config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N973']=lsst.pipe.tasks.colorterms.Colorterm()
# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N973'].primary='y'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N973'].secondary='z'

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N973'].c0=-0.007775

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N973'].c1=-0.050972

# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N973'].c2=-0.197278

# Name of photometric reference catalog; used to select a color term dict in colorterms. see also applyColorTerms
config.processCcd.calibrate.photoCal.photoCatName='sdss-dr9-fink-v5b'

# Additional magnitude uncertainty to be added in quadrature with measurement errors.
config.processCcd.calibrate.photoCal.magErrFloor=0.0

# Fields to copy from the icSource catalog to the output catalog for matching sources Any missing fields will trigger a RuntimeError exception. Ignored if icSourceCat is not provided.
config.processCcd.calibrate.icSourceFieldsToCopy=['calib_psfCandidate', 'calib_psfUsed', 'calib_psf_reserved']

# Match radius for matching icSourceCat objects to sourceCat objects (pixels)
config.processCcd.calibrate.matchRadiusPix=3.0

# Strictness of Astropy unit compatibility check, can be 'raise', 'warn' or 'silent'
config.processCcd.calibrate.checkUnitsParseStrict='raise'

# detected sources with fewer than the specified number of pixels will be ignored
config.processCcd.calibrate.detection.minPixels=1

# Pixels should be grown as isotropically as possible (slower)
config.processCcd.calibrate.detection.isotropicGrow=True

# Grow all footprints at the same time? This allows disconnected footprints to merge.
config.processCcd.calibrate.detection.combinedGrow=True

# Grow detections by nSigmaToGrow * [PSF RMS width]; if 0 then do not grow
config.processCcd.calibrate.detection.nSigmaToGrow=2.4

# Grow detections to set the image mask bits, but return the original (not-grown) footprints
config.processCcd.calibrate.detection.returnOriginalFootprints=False

# Threshold for footprints
config.processCcd.calibrate.detection.thresholdValue=5.0

# Include threshold relative to thresholdValue
config.processCcd.calibrate.detection.includeThresholdMultiplier=1.0

# specifies the desired flavor of Threshold
config.processCcd.calibrate.detection.thresholdType='stdev'

# specifies whether to detect positive, or negative sources, or both
config.processCcd.calibrate.detection.thresholdPolarity='positive'

# Fiddle factor to add to the background; debugging only
config.processCcd.calibrate.detection.adjustBackground=0.0

# Estimate the background again after final source detection?
config.processCcd.calibrate.detection.reEstimateBackground=True

# type of statistic to use for grid points
config.processCcd.calibrate.detection.background.statisticsProperty='MEANCLIP'

# behaviour if there are too few points in grid for requested interpolation style
config.processCcd.calibrate.detection.background.undersampleStyle='REDUCE_INTERP_ORDER'

# how large a region of the sky should be used for each background point
config.processCcd.calibrate.detection.background.binSize=128

# Sky region size to be used for each background point in X direction. If 0, the binSize config is used.
config.processCcd.calibrate.detection.background.binSizeX=0

# Sky region size to be used for each background point in Y direction. If 0, the binSize config is used.
config.processCcd.calibrate.detection.background.binSizeY=0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
config.processCcd.calibrate.detection.background.algorithm='AKIMA_SPLINE'

# Names of mask planes to ignore while estimating the background
config.processCcd.calibrate.detection.background.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# Ignore NaNs when estimating the background
config.processCcd.calibrate.detection.background.isNanSafe=False

# Use Approximate (Chebyshev) to model background.
config.processCcd.calibrate.detection.background.useApprox=True

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.processCcd.calibrate.detection.background.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.processCcd.calibrate.detection.background.approxOrderY=-1

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.processCcd.calibrate.detection.background.weighting=True

# type of statistic to use for grid points
config.processCcd.calibrate.detection.tempLocalBackground.statisticsProperty='MEANCLIP'

# behaviour if there are too few points in grid for requested interpolation style
config.processCcd.calibrate.detection.tempLocalBackground.undersampleStyle='REDUCE_INTERP_ORDER'

# how large a region of the sky should be used for each background point
config.processCcd.calibrate.detection.tempLocalBackground.binSize=64

# Sky region size to be used for each background point in X direction. If 0, the binSize config is used.
config.processCcd.calibrate.detection.tempLocalBackground.binSizeX=0

# Sky region size to be used for each background point in Y direction. If 0, the binSize config is used.
config.processCcd.calibrate.detection.tempLocalBackground.binSizeY=0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
config.processCcd.calibrate.detection.tempLocalBackground.algorithm='AKIMA_SPLINE'

# Names of mask planes to ignore while estimating the background
config.processCcd.calibrate.detection.tempLocalBackground.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# Ignore NaNs when estimating the background
config.processCcd.calibrate.detection.tempLocalBackground.isNanSafe=False

# Use Approximate (Chebyshev) to model background.
config.processCcd.calibrate.detection.tempLocalBackground.useApprox=False

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.processCcd.calibrate.detection.tempLocalBackground.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.processCcd.calibrate.detection.tempLocalBackground.approxOrderY=-1

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.processCcd.calibrate.detection.tempLocalBackground.weighting=True

# Enable temporary local background subtraction? (see tempLocalBackground)
config.processCcd.calibrate.detection.doTempLocalBackground=False

# type of statistic to use for grid points
config.processCcd.calibrate.detection.tempWideBackground.statisticsProperty='MEANCLIP'

# behaviour if there are too few points in grid for requested interpolation style
config.processCcd.calibrate.detection.tempWideBackground.undersampleStyle='REDUCE_INTERP_ORDER'

# how large a region of the sky should be used for each background point
config.processCcd.calibrate.detection.tempWideBackground.binSize=512

# Sky region size to be used for each background point in X direction. If 0, the binSize config is used.
config.processCcd.calibrate.detection.tempWideBackground.binSizeX=0

# Sky region size to be used for each background point in Y direction. If 0, the binSize config is used.
config.processCcd.calibrate.detection.tempWideBackground.binSizeY=0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
config.processCcd.calibrate.detection.tempWideBackground.algorithm='AKIMA_SPLINE'

# Names of mask planes to ignore while estimating the background
config.processCcd.calibrate.detection.tempWideBackground.ignoredPixelMask=['BAD', 'EDGE', 'NO_DATA']

# Ignore NaNs when estimating the background
config.processCcd.calibrate.detection.tempWideBackground.isNanSafe=False

# Use Approximate (Chebyshev) to model background.
config.processCcd.calibrate.detection.tempWideBackground.useApprox=False

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.processCcd.calibrate.detection.tempWideBackground.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.processCcd.calibrate.detection.tempWideBackground.approxOrderY=-1

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.processCcd.calibrate.detection.tempWideBackground.weighting=True

# Do temporary wide (large-scale) background subtraction before footprint detection?
config.processCcd.calibrate.detection.doTempWideBackground=False

# The maximum number of peaks in a Footprint before trying to replace its peaks using the temporary local background
config.processCcd.calibrate.detection.nPeaksMaxSimple=1

# Multiple of PSF RMS size to use for convolution kernel bounding box size; note that this is not a half-size. The size will be rounded up to the nearest odd integer
config.processCcd.calibrate.detection.nSigmaForKernel=7.0

# Mask planes to ignore when calculating statistics of image (for thresholdType=stdev)
config.processCcd.calibrate.detection.statsMask=['BAD', 'SAT', 'EDGE', 'NO_DATA']

# Run deblender input exposure
config.processCcd.calibrate.doDeblend=True

# What to do when a peak to be deblended is close to the edge of the image
config.processCcd.calibrate.deblend.edgeHandling='ramp'

# When the deblender should attribute stray flux to point sources
config.processCcd.calibrate.deblend.strayFluxToPointSources='necessary'

# Assign stray flux (not claimed by any child in the deblender) to deblend children.
config.processCcd.calibrate.deblend.assignStrayFlux=True

# How to split flux among peaks
config.processCcd.calibrate.deblend.strayFluxRule='trim'

# When splitting stray flux, clip fractions below this value to zero.
config.processCcd.calibrate.deblend.clipStrayFluxFraction=0.001

# Chi-squared per DOF cut for deciding a source is a PSF during deblending (un-shifted PSF model)
config.processCcd.calibrate.deblend.psfChisq1=1.5

# Chi-squared per DOF cut for deciding a source is PSF during deblending (shifted PSF model)
config.processCcd.calibrate.deblend.psfChisq2=1.5

# Chi-squared per DOF cut for deciding a source is a PSF during deblending (shifted PSF model #2)
config.processCcd.calibrate.deblend.psfChisq2b=1.5

# Only deblend the brightest maxNumberOfPeaks peaks in the parent (<= 0: unlimited)
config.processCcd.calibrate.deblend.maxNumberOfPeaks=0

# Maximum area for footprints before they are ignored as large; non-positive means no threshold applied
config.processCcd.calibrate.deblend.maxFootprintArea=10000

# Maximum linear dimension for footprints before they are ignored as large; non-positive means no threshold applied
config.processCcd.calibrate.deblend.maxFootprintSize=0

# Minimum axis ratio for footprints before they are ignored as large; non-positive means no threshold applied
config.processCcd.calibrate.deblend.minFootprintAxisRatio=0.0

# Mask name for footprints not deblended, or None
config.processCcd.calibrate.deblend.notDeblendedMask='NOT_DEBLENDED'

# Footprints smaller in width or height than this value will be ignored; minimum of 2 due to PSF gradient calculation.
config.processCcd.calibrate.deblend.tinyFootprintSize=2

# Guarantee that all peaks produce a child source.
config.processCcd.calibrate.deblend.propagateAllPeaks=False

# If True, catch exceptions thrown by the deblender, log them, and set a flag on the parent, instead of letting them propagate up
config.processCcd.calibrate.deblend.catchFailures=False

# Mask planes to ignore when performing statistics
config.processCcd.calibrate.deblend.maskPlanes=['SAT', 'INTRP', 'NO_DATA']

# Mask planes with the corresponding limit on the fraction of masked pixels. Sources violating this limit will not be deblended.
config.processCcd.calibrate.deblend.maskLimits={'NO_DATA': 0.25}

# If true, a least-squares fit of the templates will be done to the full image. The templates will be re-weighted based on this fit.
config.processCcd.calibrate.deblend.weightTemplates=False

# Try to remove similar templates?
config.processCcd.calibrate.deblend.removeDegenerateTemplates=False

# If the dot product between two templates is larger than this value, we consider them to be describing the same object (i.e. they are degenerate).  If one of the objects has been labeled as a PSF it will be removed, otherwise the template with the lowest value will be removed.
config.processCcd.calibrate.deblend.maxTempDotProd=0.5

# Apply a smoothing filter to all of the template images
config.processCcd.calibrate.deblend.medianSmoothTemplate=True

# the name of the centroiding algorithm used to set source x,y
config.processCcd.calibrate.measurement.slots.centroid='base_SdssCentroid'

# the name of the algorithm used to set source moments parameters
config.processCcd.calibrate.measurement.slots.shape='ext_shapeHSM_HsmSourceMoments'

# the name of the algorithm used to set PSF moments parameters
config.processCcd.calibrate.measurement.slots.psfShape='ext_shapeHSM_HsmPsfMoments'

# the name of the algorithm used to set the source aperture flux slot
config.processCcd.calibrate.measurement.slots.apFlux='base_CircularApertureFlux_12_0'

# the name of the algorithm used to set the source model flux slot
config.processCcd.calibrate.measurement.slots.modelFlux='base_GaussianFlux'

# the name of the algorithm used to set the source psf flux slot
config.processCcd.calibrate.measurement.slots.psfFlux='base_PsfFlux'

# the name of the algorithm used to set the source inst flux slot
config.processCcd.calibrate.measurement.slots.instFlux='base_GaussianFlux'

# the name of the flux measurement algorithm used for calibration
config.processCcd.calibrate.measurement.slots.calibFlux='base_CircularApertureFlux_12_0'

# When measuring, replace other detected footprints with noise?
config.processCcd.calibrate.measurement.doReplaceWithNoise=True

# How to choose mean and variance of the Gaussian noise we generate?
config.processCcd.calibrate.measurement.noiseReplacer.noiseSource='measure'

# Add ann offset to the generated noise.
config.processCcd.calibrate.measurement.noiseReplacer.noiseOffset=0.0

# The seed multiplier value to use for random number generation
#    >= 1: set the seed deterministically based on exposureId
#       0: fall back to the afw.math.Random default constructor (which uses a seed value of 1)
config.processCcd.calibrate.measurement.noiseReplacer.noiseSeedMultiplier=1

# Prefix to give undeblended plugins
config.processCcd.calibrate.measurement.undeblendedPrefix='undeblended_'

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_PsfFlux'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.processCcd.calibrate.measurement.plugins['base_PsfFlux'].badMaskPlanes=[]

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_PeakLikelihoodFlux'].doMeasure=True

# Name of warping kernel (e.g. "lanczos4") used to compute the peak
config.processCcd.calibrate.measurement.plugins['base_PeakLikelihoodFlux'].warpingKernelName='lanczos4'

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_GaussianFlux'].doMeasure=True

# FIXME! NEVER DOCUMENTED!
config.processCcd.calibrate.measurement.plugins['base_GaussianFlux'].background=0.0

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_NaiveCentroid'].doMeasure=True

# Value to subtract from the image pixel values
config.processCcd.calibrate.measurement.plugins['base_NaiveCentroid'].background=0.0

# Do check that the centroid is contained in footprint.
config.processCcd.calibrate.measurement.plugins['base_NaiveCentroid'].doFootprintCheck=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.processCcd.calibrate.measurement.plugins['base_NaiveCentroid'].maxDistToPeak=-1.0

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_SdssCentroid'].doMeasure=True

# maximum allowed binning
config.processCcd.calibrate.measurement.plugins['base_SdssCentroid'].binmax=16

# Do check that the centroid is contained in footprint.
config.processCcd.calibrate.measurement.plugins['base_SdssCentroid'].doFootprintCheck=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.processCcd.calibrate.measurement.plugins['base_SdssCentroid'].maxDistToPeak=-1.0

# if the peak's less than this insist on binning at least once
config.processCcd.calibrate.measurement.plugins['base_SdssCentroid'].peakMin=-1.0

# fiddle factor for adjusting the binning
config.processCcd.calibrate.measurement.plugins['base_SdssCentroid'].wfac=1.5

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_PixelFlags'].doMeasure=True

# List of mask planes to be searched for which occur anywhere within a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.processCcd.calibrate.measurement.plugins['base_PixelFlags'].masksFpAnywhere=[]

# List of mask planes to be searched for which occur in the center of a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.processCcd.calibrate.measurement.plugins['base_PixelFlags'].masksFpCenter=[]

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_SdssShape'].doMeasure=True

# Additional value to add to background
config.processCcd.calibrate.measurement.plugins['base_SdssShape'].background=0.0

# Whether to also compute the shape of the PSF model
config.processCcd.calibrate.measurement.plugins['base_SdssShape'].doMeasurePsf=True

# Maximum number of iterations
config.processCcd.calibrate.measurement.plugins['base_SdssShape'].maxIter=100

# Maximum centroid shift, limited to 2-10
config.processCcd.calibrate.measurement.plugins['base_SdssShape'].maxShift=0.0

# Convergence tolerance for e1,e2
config.processCcd.calibrate.measurement.plugins['base_SdssShape'].tol1=9.999999747378752e-06

# Convergence tolerance for FWHM
config.processCcd.calibrate.measurement.plugins['base_SdssShape'].tol2=9.999999747378752e-05

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_ScaledApertureFlux'].doMeasure=True

# Scaling factor of PSF FWHM for aperture radius.
config.processCcd.calibrate.measurement.plugins['base_ScaledApertureFlux'].scale=3.14

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.processCcd.calibrate.measurement.plugins['base_ScaledApertureFlux'].shiftKernel='lanczos5'

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_CircularApertureFlux'].doMeasure=True

# Maximum radius (in pixels) for which the sinc algorithm should be used instead of the faster naive algorithm.  For elliptical apertures, this is the minor axis radius.
config.processCcd.calibrate.measurement.plugins['base_CircularApertureFlux'].maxSincRadius=12.0

# Radius (in pixels) of apertures.
config.processCcd.calibrate.measurement.plugins['base_CircularApertureFlux'].radii=[3.0, 4.5, 6.0, 9.0, 12.0, 17.0, 25.0, 35.0, 50.0, 70.0]

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.processCcd.calibrate.measurement.plugins['base_CircularApertureFlux'].shiftKernel='lanczos5'

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_Blendedness'].doMeasure=True

# Whether to compute quantities related to the Gaussian-weighted flux
config.processCcd.calibrate.measurement.plugins['base_Blendedness'].doFlux=True

# Whether to compute HeavyFootprint dot products (the old deblend.blendedness parameter)
config.processCcd.calibrate.measurement.plugins['base_Blendedness'].doOld=True

# Whether to compute quantities related to the Gaussian-weighted shape
config.processCcd.calibrate.measurement.plugins['base_Blendedness'].doShape=True

# Radius factor that sets the maximum extent of the weight function (and hence the flux measurements)
config.processCcd.calibrate.measurement.plugins['base_Blendedness'].nSigmaWeightMax=3.0

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_LocalBackground'].doMeasure=True

# Inner radius for background annulus as a multiple of the PSF sigma
config.processCcd.calibrate.measurement.plugins['base_LocalBackground'].annulusInner=7.0

# Outer radius for background annulus as a multiple of the PSF sigma
config.processCcd.calibrate.measurement.plugins['base_LocalBackground'].annulusOuter=15.0

# Mask planes that indicate pixels that should be excluded from the measurement
config.processCcd.calibrate.measurement.plugins['base_LocalBackground'].badMaskPlanes=['BAD', 'SAT', 'NO_DATA']

# Number of sigma-clipping iterations for background measurement
config.processCcd.calibrate.measurement.plugins['base_LocalBackground'].bgIter=3

# Rejection threshold (in standard deviations) for background measurement
config.processCcd.calibrate.measurement.plugins['base_LocalBackground'].bgRej=3.0

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_FPPosition'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_Jacobian'].doMeasure=True

# Nominal pixel size (arcsec)
config.processCcd.calibrate.measurement.plugins['base_Jacobian'].pixelScale=0.168

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_Variance'].doMeasure=True

# Scale factor to apply to shape for aperture
config.processCcd.calibrate.measurement.plugins['base_Variance'].scale=5.0

# Mask planes to ignore
config.processCcd.calibrate.measurement.plugins['base_Variance'].mask=['DETECTED', 'DETECTED_NEGATIVE', 'BAD', 'SAT']

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_InputCount'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_PeakCentroid'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_SkyCoord'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['subaru_FilterFraction'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].doMeasure=True

# Shapelet order of inner expansion (0 == Gaussian)
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].innerOrder=2

# Don't allow the semi-major radius of any component to go above this fraction of the PSF image width
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].maxRadiusBoxFraction=0.4

# Don't allow the semi-minor radius of any component to drop below this value (pixels)
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].minRadius=1.0

# Don't allow the determinant radii of the two components to differ by less than this (pixels)
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].minRadiusDiff=0.5

# whether to save all iterations for debugging purposes
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionSolverTolerance=1e-08

# Shapelet order of outer expansion (0 == Gaussian)
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].outerOrder=1

# Initial outer Gaussian peak height divided by inner Gaussian peak height
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].peakRatio=0.1

# Initial outer radius divided by inner radius
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].radiusRatio=2.0

config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models={}
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.radiusPriorSigma=0.5

config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.radiusPriorSigma=0.5

config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.order=2

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.order=1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.radiusPriorSigma=0.5

config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.order=4

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.order=4

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.radiusPriorSigma=0.5

# a sequence of model names indicating which models should be fit, and their order
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].sequence=['DoubleShapelet']

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].doMeasure=True

# Whether to record the steps the optimizer takes (or just the number, if running as a plugin)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.doRecordHistory=True

# Whether to record the time spent in this stage
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.doRecordTime=True

# Softened core width for ellipticity distribution (conformal shear units).
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.empiricalPriorConfig.ellipticityCore=0.001

# Width of exponential ellipticity distribution (conformal shear units).
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.empiricalPriorConfig.ellipticitySigma=0.3

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusMinInner=-6.0

# Minimum ln(radius).
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusMinOuter=-6.001

# Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T.
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusMu=-1.0

# Number of degrees of freedom for the Student's T distribution on ln(radius).
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusNu=50.0

# Width of the Student's T distribution in ln(radius).
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusSigma=0.45

# Ellipticity magnitude (conformal shear units) at which the softened cutoff begins
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.linearPriorConfig.ellipticityMaxInner=2.0

# Maximum ellipticity magnitude (conformal shear units)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.linearPriorConfig.ellipticityMaxOuter=2.001

# ln(radius) at which the softened cutoff begins towards the maximum
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.linearPriorConfig.logRadiusMaxInner=3.0

# Maximum ln(radius)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.linearPriorConfig.logRadiusMaxOuter=3.001

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.linearPriorConfig.logRadiusMinInner=-6.0

# The ratio P(logRadiusMinInner)/P(logRadiusMaxInner)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.linearPriorConfig.logRadiusMinMaxRatio=1.0

# Minimum ln(radius)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.linearPriorConfig.logRadiusMinOuter=-6.001

# Maximum radius used in approximating profile with Gaussians (0=default for this profile)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.maxRadius=0

# Number of Gaussian used to approximate the profile
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.nComponents=8

# whether to save all iterations for debugging purposes
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.optimizer.trustRegionSolverTolerance=1e-08

# Name of the Prior that defines the model to fit (a filename in $MEAS_MODELFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting.
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.priorName=''

# One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.priorSource='EMPIRICAL'

# Name of the shapelet.RadialProfile that defines the model to fit
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.profileName='luv'

# Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.usePixelWeights=False

# Scale the likelihood by this factor to artificially reweight it w.r.t. the prior.
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.weightsMultiplier=1.0

# Whether to record the steps the optimizer takes (or just the number, if running as a plugin)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.doRecordHistory=True

# Whether to record the time spent in this stage
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.doRecordTime=True

# Softened core width for ellipticity distribution (conformal shear units).
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.empiricalPriorConfig.ellipticityCore=0.001

# Width of exponential ellipticity distribution (conformal shear units).
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.empiricalPriorConfig.ellipticitySigma=0.3

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusMinInner=-6.0

# Minimum ln(radius).
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusMinOuter=-6.001

# Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T.
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusMu=-1.0

# Number of degrees of freedom for the Student's T distribution on ln(radius).
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusNu=50.0

# Width of the Student's T distribution in ln(radius).
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusSigma=0.45

# Ellipticity magnitude (conformal shear units) at which the softened cutoff begins
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.linearPriorConfig.ellipticityMaxInner=2.0

# Maximum ellipticity magnitude (conformal shear units)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.linearPriorConfig.ellipticityMaxOuter=2.001

# ln(radius) at which the softened cutoff begins towards the maximum
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.linearPriorConfig.logRadiusMaxInner=3.0

# Maximum ln(radius)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.linearPriorConfig.logRadiusMaxOuter=3.001

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.linearPriorConfig.logRadiusMinInner=-6.0

# The ratio P(logRadiusMinInner)/P(logRadiusMaxInner)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.linearPriorConfig.logRadiusMinMaxRatio=1.0

# Minimum ln(radius)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.linearPriorConfig.logRadiusMinOuter=-6.001

# Maximum radius used in approximating profile with Gaussians (0=default for this profile)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.maxRadius=0

# Number of Gaussian used to approximate the profile
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.nComponents=6

# whether to save all iterations for debugging purposes
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.optimizer.maxOuterIterations=250

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.optimizer.trustRegionSolverTolerance=1e-08

# Name of the Prior that defines the model to fit (a filename in $MEAS_MODELFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting.
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.priorName=''

# One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.priorSource='EMPIRICAL'

# Name of the shapelet.RadialProfile that defines the model to fit
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.profileName='lux'

# Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.usePixelWeights=False

# Scale the likelihood by this factor to artificially reweight it w.r.t. the prior.
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.weightsMultiplier=1.0

# If the 2nd-moments shape used to initialize the fit failed, use the PSF moments multiplied by this.  If <= 0.0, abort the fit early instead.
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].fallbackInitialMomentsPsfFactor=1.5

# Whether to record the steps the optimizer takes (or just the number, if running as a plugin)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.doRecordHistory=True

# Whether to record the time spent in this stage
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.doRecordTime=True

# Softened core width for ellipticity distribution (conformal shear units).
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.empiricalPriorConfig.ellipticityCore=0.001

# Width of exponential ellipticity distribution (conformal shear units).
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.empiricalPriorConfig.ellipticitySigma=0.3

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusMinInner=-6.0

# Minimum ln(radius).
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusMinOuter=-6.001

# Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T.
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusMu=-1.0

# Number of degrees of freedom for the Student's T distribution on ln(radius).
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusNu=50.0

# Width of the Student's T distribution in ln(radius).
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusSigma=0.45

# Ellipticity magnitude (conformal shear units) at which the softened cutoff begins
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.linearPriorConfig.ellipticityMaxInner=2.0

# Maximum ellipticity magnitude (conformal shear units)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.linearPriorConfig.ellipticityMaxOuter=2.001

# ln(radius) at which the softened cutoff begins towards the maximum
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.linearPriorConfig.logRadiusMaxInner=3.0

# Maximum ln(radius)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.linearPriorConfig.logRadiusMaxOuter=3.001

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.linearPriorConfig.logRadiusMinInner=-6.0

# The ratio P(logRadiusMinInner)/P(logRadiusMaxInner)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.linearPriorConfig.logRadiusMinMaxRatio=1.0

# Minimum ln(radius)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.linearPriorConfig.logRadiusMinOuter=-6.001

# Maximum radius used in approximating profile with Gaussians (0=default for this profile)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.maxRadius=0

# Number of Gaussian used to approximate the profile
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.nComponents=3

# whether to save all iterations for debugging purposes
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.optimizer.gradientThreshold=0.01

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.optimizer.minTrustRadiusThreshold=0.01

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.optimizer.trustRegionSolverTolerance=1e-08

# Name of the Prior that defines the model to fit (a filename in $MEAS_MODELFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting.
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.priorName=''

# One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.priorSource='EMPIRICAL'

# Name of the shapelet.RadialProfile that defines the model to fit
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.profileName='lux'

# Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.usePixelWeights=True

# Scale the likelihood by this factor to artificially reweight it w.r.t. the prior.
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.weightsMultiplier=1.0

# Minimum initial radius in pixels (used to regularize initial moments-based PSF deconvolution)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].minInitialRadius=0.1

# Field name prefix of the Shapelet PSF approximation used to convolve the galaxy model; must contain a set of fields matching the schema defined by shapelet.MultiShapeletFunctionKey.
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].psfName='modelfit_DoubleShapeletPsfApprox'

# Mask planes that indicate pixels that should be ignored in the fit.
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].region.badMaskPlanes=['EDGE', 'SAT', 'BAD', 'NO_DATA']

# Abort if the fit region grows beyond this many pixels.
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].region.maxArea=100000

# Maximum fraction of pixels that may be ignored due to masks; more than this and we don't even try.
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].region.maxBadPixelFraction=0.1

# Use this multiple of the initial fit ellipse then grow by the PSF width to determine the maximum final fit region size.
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].region.nFitRadiiMax=3.0

# Use this multiple of the initial fit ellipse then grow by the PSF width to determine the minimum final fit region size.
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].region.nFitRadiiMin=1.0

# Use this multiple of the Kron ellipse to set the fit region (for the final fit region, subject to the nFitRadiiMin and nFitRadiiMax constraints).
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].region.nKronRadii=1.5

# Grow the initial fit ellipses by this factor before comparing with the Kron/Footprint region
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].region.nPsfSigmaGrow=2.0

# If the Kron radius is less than this multiple of the PSF width, ignore it and fall back to a PSF-oriented ellipse scaled to match the area of the footprint or this radius (whichever is larger).
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].region.nPsfSigmaMin=4.0

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['ext_photometryKron_KronFlux'].doMeasure=True

# If true check that the Kron radius exceeds some minimum
config.processCcd.calibrate.measurement.plugins['ext_photometryKron_KronFlux'].enforceMinimumRadius=True

# if true, use existing shape and centroid measurements instead of fitting
config.processCcd.calibrate.measurement.plugins['ext_photometryKron_KronFlux'].fixed=False

# Largest aperture for which to use the slow, accurate, sinc aperture code
config.processCcd.calibrate.measurement.plugins['ext_photometryKron_KronFlux'].maxSincRadius=10.0

# Minimum Kron radius (if == 0.0 use PSF's Kron radius) if enforceMinimumRadius. Also functions as fallback aperture radius if set.
config.processCcd.calibrate.measurement.plugins['ext_photometryKron_KronFlux'].minimumRadius=0.0

# Number of times to iterate when setting the Kron radius
config.processCcd.calibrate.measurement.plugins['ext_photometryKron_KronFlux'].nIterForRadius=1

# Number of Kron radii for Kron flux
config.processCcd.calibrate.measurement.plugins['ext_photometryKron_KronFlux'].nRadiusForFlux=2.5

# Multiplier of rms size for aperture used to initially estimate the Kron radius
config.processCcd.calibrate.measurement.plugins['ext_photometryKron_KronFlux'].nSigmaForRadius=6.0

# Name of field specifying reference Kron radius for forced measurement
config.processCcd.calibrate.measurement.plugins['ext_photometryKron_KronFlux'].refRadiusName='ext_photometryKron_KronFlux_radius'

# Smooth image with N(0, smoothingSigma^2) Gaussian while estimating R_K
config.processCcd.calibrate.measurement.plugins['ext_photometryKron_KronFlux'].smoothingSigma=-1.0

# Use the Footprint size as part of initial estimate of Kron radius
config.processCcd.calibrate.measurement.plugins['ext_photometryKron_KronFlux'].useFootprintRadius=False

# list of target seeings (FWHM, pixels)
config.processCcd.calibrate.measurement.plugins['ext_convolved_ConvolvedFlux'].seeing=[3.5, 5.0, 6.5]

# scaling factor of kernel sigma for kernel size
config.processCcd.calibrate.measurement.plugins['ext_convolved_ConvolvedFlux'].kernelScale=4.0

# Maximum radius (in pixels) for which the sinc algorithm should be used instead of the faster naive algorithm.  For elliptical apertures, this is the minor axis radius.
config.processCcd.calibrate.measurement.plugins['ext_convolved_ConvolvedFlux'].aperture.maxSincRadius=10.0

# Radius (in pixels) of apertures.
config.processCcd.calibrate.measurement.plugins['ext_convolved_ConvolvedFlux'].aperture.radii=[3.3, 4.5, 6.0]

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.processCcd.calibrate.measurement.plugins['ext_convolved_ConvolvedFlux'].aperture.shiftKernel='lanczos5'

# name of Kron radius field in reference
config.processCcd.calibrate.measurement.plugins['ext_convolved_ConvolvedFlux'].kronRadiusName='ext_photometryKron_KronFlux_radius'

# Largest aperture for which to use the sinc aperture code for Kron (pixels)
config.processCcd.calibrate.measurement.plugins['ext_convolved_ConvolvedFlux'].maxSincRadius=10.0

# Number of Kron radii for Kron flux
config.processCcd.calibrate.measurement.plugins['ext_convolved_ConvolvedFlux'].kronRadiusForFlux=2.5

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['ext_convolved_ConvolvedFlux'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['ext_shapeHSM_HsmShapeBj'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.processCcd.calibrate.measurement.plugins['ext_shapeHSM_HsmShapeBj'].badMaskPlanes=None

# Field name for number of deblend children
config.processCcd.calibrate.measurement.plugins['ext_shapeHSM_HsmShapeBj'].deblendNChild=None

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['ext_shapeHSM_HsmShapeLinear'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.processCcd.calibrate.measurement.plugins['ext_shapeHSM_HsmShapeLinear'].badMaskPlanes=None

# Field name for number of deblend children
config.processCcd.calibrate.measurement.plugins['ext_shapeHSM_HsmShapeLinear'].deblendNChild=None

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['ext_shapeHSM_HsmShapeKsb'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.processCcd.calibrate.measurement.plugins['ext_shapeHSM_HsmShapeKsb'].badMaskPlanes=None

# Field name for number of deblend children
config.processCcd.calibrate.measurement.plugins['ext_shapeHSM_HsmShapeKsb'].deblendNChild=None

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['ext_shapeHSM_HsmShapeRegauss'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.processCcd.calibrate.measurement.plugins['ext_shapeHSM_HsmShapeRegauss'].badMaskPlanes=None

# Field name for number of deblend children
config.processCcd.calibrate.measurement.plugins['ext_shapeHSM_HsmShapeRegauss'].deblendNChild='deblend_nChild'

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['ext_shapeHSM_HsmSourceMoments'].doMeasure=True

# Store measured flux?
config.processCcd.calibrate.measurement.plugins['ext_shapeHSM_HsmSourceMoments'].addFlux=None

# Mask planes used to reject bad pixels.
config.processCcd.calibrate.measurement.plugins['ext_shapeHSM_HsmSourceMoments'].badMaskPlanes=None

# Use round weight function?
config.processCcd.calibrate.measurement.plugins['ext_shapeHSM_HsmSourceMoments'].roundMoments=None

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['ext_shapeHSM_HsmSourceMomentsRound'].doMeasure=True

# Store measured flux?
config.processCcd.calibrate.measurement.plugins['ext_shapeHSM_HsmSourceMomentsRound'].addFlux=None

# Mask planes used to reject bad pixels.
config.processCcd.calibrate.measurement.plugins['ext_shapeHSM_HsmSourceMomentsRound'].badMaskPlanes=None

# Use round weight function?
config.processCcd.calibrate.measurement.plugins['ext_shapeHSM_HsmSourceMomentsRound'].roundMoments=None

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['ext_shapeHSM_HsmPsfMoments'].doMeasure=True

config.processCcd.calibrate.measurement.plugins.names=['ext_shapeHSM_HsmPsfMoments', 'ext_shapeHSM_HsmSourceMoments', 'base_LocalBackground', 'base_SdssShape', 'base_PixelFlags', 'base_Blendedness', 'ext_shapeHSM_HsmShapeRegauss', 'base_SdssCentroid', 'base_CircularApertureFlux', 'ext_photometryKron_KronFlux', 'base_Jacobian', 'ext_shapeHSM_HsmSourceMomentsRound', 'base_PsfFlux', 'base_FPPosition', 'base_Variance', 'base_SkyCoord', 'base_NaiveCentroid', 'base_GaussianFlux']
# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['base_PsfFlux'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.processCcd.calibrate.measurement.undeblended['base_PsfFlux'].badMaskPlanes=[]

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['base_PeakLikelihoodFlux'].doMeasure=True

# Name of warping kernel (e.g. "lanczos4") used to compute the peak
config.processCcd.calibrate.measurement.undeblended['base_PeakLikelihoodFlux'].warpingKernelName='lanczos4'

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['base_GaussianFlux'].doMeasure=True

# FIXME! NEVER DOCUMENTED!
config.processCcd.calibrate.measurement.undeblended['base_GaussianFlux'].background=0.0

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['base_NaiveCentroid'].doMeasure=True

# Value to subtract from the image pixel values
config.processCcd.calibrate.measurement.undeblended['base_NaiveCentroid'].background=0.0

# Do check that the centroid is contained in footprint.
config.processCcd.calibrate.measurement.undeblended['base_NaiveCentroid'].doFootprintCheck=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.processCcd.calibrate.measurement.undeblended['base_NaiveCentroid'].maxDistToPeak=-1.0

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['base_SdssCentroid'].doMeasure=True

# maximum allowed binning
config.processCcd.calibrate.measurement.undeblended['base_SdssCentroid'].binmax=16

# Do check that the centroid is contained in footprint.
config.processCcd.calibrate.measurement.undeblended['base_SdssCentroid'].doFootprintCheck=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.processCcd.calibrate.measurement.undeblended['base_SdssCentroid'].maxDistToPeak=-1.0

# if the peak's less than this insist on binning at least once
config.processCcd.calibrate.measurement.undeblended['base_SdssCentroid'].peakMin=-1.0

# fiddle factor for adjusting the binning
config.processCcd.calibrate.measurement.undeblended['base_SdssCentroid'].wfac=1.5

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['base_PixelFlags'].doMeasure=True

# List of mask planes to be searched for which occur anywhere within a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.processCcd.calibrate.measurement.undeblended['base_PixelFlags'].masksFpAnywhere=[]

# List of mask planes to be searched for which occur in the center of a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.processCcd.calibrate.measurement.undeblended['base_PixelFlags'].masksFpCenter=[]

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['base_SdssShape'].doMeasure=True

# Additional value to add to background
config.processCcd.calibrate.measurement.undeblended['base_SdssShape'].background=0.0

# Whether to also compute the shape of the PSF model
config.processCcd.calibrate.measurement.undeblended['base_SdssShape'].doMeasurePsf=True

# Maximum number of iterations
config.processCcd.calibrate.measurement.undeblended['base_SdssShape'].maxIter=100

# Maximum centroid shift, limited to 2-10
config.processCcd.calibrate.measurement.undeblended['base_SdssShape'].maxShift=0.0

# Convergence tolerance for e1,e2
config.processCcd.calibrate.measurement.undeblended['base_SdssShape'].tol1=9.999999747378752e-06

# Convergence tolerance for FWHM
config.processCcd.calibrate.measurement.undeblended['base_SdssShape'].tol2=9.999999747378752e-05

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['base_ScaledApertureFlux'].doMeasure=True

# Scaling factor of PSF FWHM for aperture radius.
config.processCcd.calibrate.measurement.undeblended['base_ScaledApertureFlux'].scale=3.14

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.processCcd.calibrate.measurement.undeblended['base_ScaledApertureFlux'].shiftKernel='lanczos5'

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['base_CircularApertureFlux'].doMeasure=True

# Maximum radius (in pixels) for which the sinc algorithm should be used instead of the faster naive algorithm.  For elliptical apertures, this is the minor axis radius.
config.processCcd.calibrate.measurement.undeblended['base_CircularApertureFlux'].maxSincRadius=10.0

# Radius (in pixels) of apertures.
config.processCcd.calibrate.measurement.undeblended['base_CircularApertureFlux'].radii=[3.0, 4.5, 6.0, 9.0, 12.0, 17.0, 25.0, 35.0, 50.0, 70.0]

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.processCcd.calibrate.measurement.undeblended['base_CircularApertureFlux'].shiftKernel='lanczos5'

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['base_Blendedness'].doMeasure=True

# Whether to compute quantities related to the Gaussian-weighted flux
config.processCcd.calibrate.measurement.undeblended['base_Blendedness'].doFlux=True

# Whether to compute HeavyFootprint dot products (the old deblend.blendedness parameter)
config.processCcd.calibrate.measurement.undeblended['base_Blendedness'].doOld=True

# Whether to compute quantities related to the Gaussian-weighted shape
config.processCcd.calibrate.measurement.undeblended['base_Blendedness'].doShape=True

# Radius factor that sets the maximum extent of the weight function (and hence the flux measurements)
config.processCcd.calibrate.measurement.undeblended['base_Blendedness'].nSigmaWeightMax=3.0

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['base_LocalBackground'].doMeasure=True

# Inner radius for background annulus as a multiple of the PSF sigma
config.processCcd.calibrate.measurement.undeblended['base_LocalBackground'].annulusInner=7.0

# Outer radius for background annulus as a multiple of the PSF sigma
config.processCcd.calibrate.measurement.undeblended['base_LocalBackground'].annulusOuter=15.0

# Mask planes that indicate pixels that should be excluded from the measurement
config.processCcd.calibrate.measurement.undeblended['base_LocalBackground'].badMaskPlanes=['BAD', 'SAT', 'NO_DATA']

# Number of sigma-clipping iterations for background measurement
config.processCcd.calibrate.measurement.undeblended['base_LocalBackground'].bgIter=3

# Rejection threshold (in standard deviations) for background measurement
config.processCcd.calibrate.measurement.undeblended['base_LocalBackground'].bgRej=3.0

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['base_FPPosition'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['base_Jacobian'].doMeasure=True

# Nominal pixel size (arcsec)
config.processCcd.calibrate.measurement.undeblended['base_Jacobian'].pixelScale=0.5

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['base_Variance'].doMeasure=True

# Scale factor to apply to shape for aperture
config.processCcd.calibrate.measurement.undeblended['base_Variance'].scale=5.0

# Mask planes to ignore
config.processCcd.calibrate.measurement.undeblended['base_Variance'].mask=['DETECTED', 'DETECTED_NEGATIVE', 'BAD', 'SAT']

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['base_InputCount'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['base_PeakCentroid'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['base_SkyCoord'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['subaru_FilterFraction'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].doMeasure=True

# Shapelet order of inner expansion (0 == Gaussian)
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].innerOrder=2

# Don't allow the semi-major radius of any component to go above this fraction of the PSF image width
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].maxRadiusBoxFraction=0.4

# Don't allow the semi-minor radius of any component to drop below this value (pixels)
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].minRadius=1.0

# Don't allow the determinant radii of the two components to differ by less than this (pixels)
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].minRadiusDiff=0.5

# whether to save all iterations for debugging purposes
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionSolverTolerance=1e-08

# Shapelet order of outer expansion (0 == Gaussian)
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].outerOrder=1

# Initial outer Gaussian peak height divided by inner Gaussian peak height
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].peakRatio=0.1

# Initial outer radius divided by inner radius
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].radiusRatio=2.0

config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models={}
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.radiusPriorSigma=0.5

config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.radiusPriorSigma=0.5

config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.order=2

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.order=1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.radiusPriorSigma=0.5

config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.order=4

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.order=4

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.radiusPriorSigma=0.5

# a sequence of model names indicating which models should be fit, and their order
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].sequence=['DoubleShapelet']

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].doMeasure=True

# Whether to record the steps the optimizer takes (or just the number, if running as a plugin)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.doRecordHistory=True

# Whether to record the time spent in this stage
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.doRecordTime=True

# Softened core width for ellipticity distribution (conformal shear units).
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.empiricalPriorConfig.ellipticityCore=0.001

# Width of exponential ellipticity distribution (conformal shear units).
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.empiricalPriorConfig.ellipticitySigma=0.3

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusMinInner=-6.0

# Minimum ln(radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusMinOuter=-6.001

# Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T.
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusMu=-1.0

# Number of degrees of freedom for the Student's T distribution on ln(radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusNu=50.0

# Width of the Student's T distribution in ln(radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusSigma=0.45

# Ellipticity magnitude (conformal shear units) at which the softened cutoff begins
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.linearPriorConfig.ellipticityMaxInner=2.0

# Maximum ellipticity magnitude (conformal shear units)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.linearPriorConfig.ellipticityMaxOuter=2.001

# ln(radius) at which the softened cutoff begins towards the maximum
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.linearPriorConfig.logRadiusMaxInner=3.0

# Maximum ln(radius)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.linearPriorConfig.logRadiusMaxOuter=3.001

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.linearPriorConfig.logRadiusMinInner=-6.0

# The ratio P(logRadiusMinInner)/P(logRadiusMaxInner)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.linearPriorConfig.logRadiusMinMaxRatio=1.0

# Minimum ln(radius)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.linearPriorConfig.logRadiusMinOuter=-6.001

# Maximum radius used in approximating profile with Gaussians (0=default for this profile)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.maxRadius=0

# Number of Gaussian used to approximate the profile
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.nComponents=8

# whether to save all iterations for debugging purposes
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.optimizer.trustRegionSolverTolerance=1e-08

# Name of the Prior that defines the model to fit (a filename in $MEAS_MODELFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting.
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.priorName=''

# One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.priorSource='EMPIRICAL'

# Name of the shapelet.RadialProfile that defines the model to fit
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.profileName='luv'

# Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.usePixelWeights=False

# Scale the likelihood by this factor to artificially reweight it w.r.t. the prior.
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.weightsMultiplier=1.0

# Whether to record the steps the optimizer takes (or just the number, if running as a plugin)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.doRecordHistory=True

# Whether to record the time spent in this stage
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.doRecordTime=True

# Softened core width for ellipticity distribution (conformal shear units).
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.empiricalPriorConfig.ellipticityCore=0.001

# Width of exponential ellipticity distribution (conformal shear units).
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.empiricalPriorConfig.ellipticitySigma=0.3

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusMinInner=-6.0

# Minimum ln(radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusMinOuter=-6.001

# Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T.
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusMu=-1.0

# Number of degrees of freedom for the Student's T distribution on ln(radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusNu=50.0

# Width of the Student's T distribution in ln(radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusSigma=0.45

# Ellipticity magnitude (conformal shear units) at which the softened cutoff begins
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.linearPriorConfig.ellipticityMaxInner=2.0

# Maximum ellipticity magnitude (conformal shear units)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.linearPriorConfig.ellipticityMaxOuter=2.001

# ln(radius) at which the softened cutoff begins towards the maximum
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.linearPriorConfig.logRadiusMaxInner=3.0

# Maximum ln(radius)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.linearPriorConfig.logRadiusMaxOuter=3.001

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.linearPriorConfig.logRadiusMinInner=-6.0

# The ratio P(logRadiusMinInner)/P(logRadiusMaxInner)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.linearPriorConfig.logRadiusMinMaxRatio=1.0

# Minimum ln(radius)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.linearPriorConfig.logRadiusMinOuter=-6.001

# Maximum radius used in approximating profile with Gaussians (0=default for this profile)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.maxRadius=0

# Number of Gaussian used to approximate the profile
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.nComponents=6

# whether to save all iterations for debugging purposes
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.optimizer.maxOuterIterations=250

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.optimizer.trustRegionSolverTolerance=1e-08

# Name of the Prior that defines the model to fit (a filename in $MEAS_MODELFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting.
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.priorName=''

# One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.priorSource='EMPIRICAL'

# Name of the shapelet.RadialProfile that defines the model to fit
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.profileName='lux'

# Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.usePixelWeights=False

# Scale the likelihood by this factor to artificially reweight it w.r.t. the prior.
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.weightsMultiplier=1.0

# If the 2nd-moments shape used to initialize the fit failed, use the PSF moments multiplied by this.  If <= 0.0, abort the fit early instead.
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].fallbackInitialMomentsPsfFactor=1.5

# Whether to record the steps the optimizer takes (or just the number, if running as a plugin)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.doRecordHistory=True

# Whether to record the time spent in this stage
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.doRecordTime=True

# Softened core width for ellipticity distribution (conformal shear units).
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.empiricalPriorConfig.ellipticityCore=0.001

# Width of exponential ellipticity distribution (conformal shear units).
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.empiricalPriorConfig.ellipticitySigma=0.3

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusMinInner=-6.0

# Minimum ln(radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusMinOuter=-6.001

# Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T.
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusMu=-1.0

# Number of degrees of freedom for the Student's T distribution on ln(radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusNu=50.0

# Width of the Student's T distribution in ln(radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusSigma=0.45

# Ellipticity magnitude (conformal shear units) at which the softened cutoff begins
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.linearPriorConfig.ellipticityMaxInner=2.0

# Maximum ellipticity magnitude (conformal shear units)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.linearPriorConfig.ellipticityMaxOuter=2.001

# ln(radius) at which the softened cutoff begins towards the maximum
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.linearPriorConfig.logRadiusMaxInner=3.0

# Maximum ln(radius)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.linearPriorConfig.logRadiusMaxOuter=3.001

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.linearPriorConfig.logRadiusMinInner=-6.0

# The ratio P(logRadiusMinInner)/P(logRadiusMaxInner)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.linearPriorConfig.logRadiusMinMaxRatio=1.0

# Minimum ln(radius)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.linearPriorConfig.logRadiusMinOuter=-6.001

# Maximum radius used in approximating profile with Gaussians (0=default for this profile)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.maxRadius=0

# Number of Gaussian used to approximate the profile
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.nComponents=3

# whether to save all iterations for debugging purposes
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.optimizer.gradientThreshold=0.01

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.optimizer.minTrustRadiusThreshold=0.01

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.optimizer.trustRegionSolverTolerance=1e-08

# Name of the Prior that defines the model to fit (a filename in $MEAS_MODELFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting.
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.priorName=''

# One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.priorSource='EMPIRICAL'

# Name of the shapelet.RadialProfile that defines the model to fit
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.profileName='lux'

# Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.usePixelWeights=True

# Scale the likelihood by this factor to artificially reweight it w.r.t. the prior.
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.weightsMultiplier=1.0

# Minimum initial radius in pixels (used to regularize initial moments-based PSF deconvolution)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].minInitialRadius=0.1

# Field name prefix of the Shapelet PSF approximation used to convolve the galaxy model; must contain a set of fields matching the schema defined by shapelet.MultiShapeletFunctionKey.
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].psfName='modelfit_DoubleShapeletPsfApprox'

# Mask planes that indicate pixels that should be ignored in the fit.
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].region.badMaskPlanes=['EDGE', 'SAT', 'BAD', 'NO_DATA']

# Abort if the fit region grows beyond this many pixels.
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].region.maxArea=100000

# Maximum fraction of pixels that may be ignored due to masks; more than this and we don't even try.
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].region.maxBadPixelFraction=0.1

# Use this multiple of the initial fit ellipse then grow by the PSF width to determine the maximum final fit region size.
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].region.nFitRadiiMax=3.0

# Use this multiple of the initial fit ellipse then grow by the PSF width to determine the minimum final fit region size.
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].region.nFitRadiiMin=1.0

# Use this multiple of the Kron ellipse to set the fit region (for the final fit region, subject to the nFitRadiiMin and nFitRadiiMax constraints).
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].region.nKronRadii=1.5

# Grow the initial fit ellipses by this factor before comparing with the Kron/Footprint region
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].region.nPsfSigmaGrow=2.0

# If the Kron radius is less than this multiple of the PSF width, ignore it and fall back to a PSF-oriented ellipse scaled to match the area of the footprint or this radius (whichever is larger).
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].region.nPsfSigmaMin=4.0

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['ext_photometryKron_KronFlux'].doMeasure=True

# If true check that the Kron radius exceeds some minimum
config.processCcd.calibrate.measurement.undeblended['ext_photometryKron_KronFlux'].enforceMinimumRadius=True

# if true, use existing shape and centroid measurements instead of fitting
config.processCcd.calibrate.measurement.undeblended['ext_photometryKron_KronFlux'].fixed=False

# Largest aperture for which to use the slow, accurate, sinc aperture code
config.processCcd.calibrate.measurement.undeblended['ext_photometryKron_KronFlux'].maxSincRadius=10.0

# Minimum Kron radius (if == 0.0 use PSF's Kron radius) if enforceMinimumRadius. Also functions as fallback aperture radius if set.
config.processCcd.calibrate.measurement.undeblended['ext_photometryKron_KronFlux'].minimumRadius=0.0

# Number of times to iterate when setting the Kron radius
config.processCcd.calibrate.measurement.undeblended['ext_photometryKron_KronFlux'].nIterForRadius=1

# Number of Kron radii for Kron flux
config.processCcd.calibrate.measurement.undeblended['ext_photometryKron_KronFlux'].nRadiusForFlux=2.5

# Multiplier of rms size for aperture used to initially estimate the Kron radius
config.processCcd.calibrate.measurement.undeblended['ext_photometryKron_KronFlux'].nSigmaForRadius=6.0

# Name of field specifying reference Kron radius for forced measurement
config.processCcd.calibrate.measurement.undeblended['ext_photometryKron_KronFlux'].refRadiusName='ext_photometryKron_KronFlux_radius'

# Smooth image with N(0, smoothingSigma^2) Gaussian while estimating R_K
config.processCcd.calibrate.measurement.undeblended['ext_photometryKron_KronFlux'].smoothingSigma=-1.0

# Use the Footprint size as part of initial estimate of Kron radius
config.processCcd.calibrate.measurement.undeblended['ext_photometryKron_KronFlux'].useFootprintRadius=False

# list of target seeings (FWHM, pixels)
config.processCcd.calibrate.measurement.undeblended['ext_convolved_ConvolvedFlux'].seeing=[3.5, 5.0, 6.5]

# scaling factor of kernel sigma for kernel size
config.processCcd.calibrate.measurement.undeblended['ext_convolved_ConvolvedFlux'].kernelScale=4.0

# Maximum radius (in pixels) for which the sinc algorithm should be used instead of the faster naive algorithm.  For elliptical apertures, this is the minor axis radius.
config.processCcd.calibrate.measurement.undeblended['ext_convolved_ConvolvedFlux'].aperture.maxSincRadius=10.0

# Radius (in pixels) of apertures.
config.processCcd.calibrate.measurement.undeblended['ext_convolved_ConvolvedFlux'].aperture.radii=[3.3, 4.5, 6.0]

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.processCcd.calibrate.measurement.undeblended['ext_convolved_ConvolvedFlux'].aperture.shiftKernel='lanczos5'

# name of Kron radius field in reference
config.processCcd.calibrate.measurement.undeblended['ext_convolved_ConvolvedFlux'].kronRadiusName='ext_photometryKron_KronFlux_radius'

# Largest aperture for which to use the sinc aperture code for Kron (pixels)
config.processCcd.calibrate.measurement.undeblended['ext_convolved_ConvolvedFlux'].maxSincRadius=10.0

# Number of Kron radii for Kron flux
config.processCcd.calibrate.measurement.undeblended['ext_convolved_ConvolvedFlux'].kronRadiusForFlux=2.5

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['ext_convolved_ConvolvedFlux'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['ext_shapeHSM_HsmShapeBj'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.processCcd.calibrate.measurement.undeblended['ext_shapeHSM_HsmShapeBj'].badMaskPlanes=None

# Field name for number of deblend children
config.processCcd.calibrate.measurement.undeblended['ext_shapeHSM_HsmShapeBj'].deblendNChild=None

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['ext_shapeHSM_HsmShapeLinear'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.processCcd.calibrate.measurement.undeblended['ext_shapeHSM_HsmShapeLinear'].badMaskPlanes=None

# Field name for number of deblend children
config.processCcd.calibrate.measurement.undeblended['ext_shapeHSM_HsmShapeLinear'].deblendNChild=None

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['ext_shapeHSM_HsmShapeKsb'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.processCcd.calibrate.measurement.undeblended['ext_shapeHSM_HsmShapeKsb'].badMaskPlanes=None

# Field name for number of deblend children
config.processCcd.calibrate.measurement.undeblended['ext_shapeHSM_HsmShapeKsb'].deblendNChild=None

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['ext_shapeHSM_HsmShapeRegauss'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.processCcd.calibrate.measurement.undeblended['ext_shapeHSM_HsmShapeRegauss'].badMaskPlanes=None

# Field name for number of deblend children
config.processCcd.calibrate.measurement.undeblended['ext_shapeHSM_HsmShapeRegauss'].deblendNChild=None

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['ext_shapeHSM_HsmSourceMoments'].doMeasure=True

# Store measured flux?
config.processCcd.calibrate.measurement.undeblended['ext_shapeHSM_HsmSourceMoments'].addFlux=None

# Mask planes used to reject bad pixels.
config.processCcd.calibrate.measurement.undeblended['ext_shapeHSM_HsmSourceMoments'].badMaskPlanes=None

# Use round weight function?
config.processCcd.calibrate.measurement.undeblended['ext_shapeHSM_HsmSourceMoments'].roundMoments=None

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['ext_shapeHSM_HsmSourceMomentsRound'].doMeasure=True

# Store measured flux?
config.processCcd.calibrate.measurement.undeblended['ext_shapeHSM_HsmSourceMomentsRound'].addFlux=None

# Mask planes used to reject bad pixels.
config.processCcd.calibrate.measurement.undeblended['ext_shapeHSM_HsmSourceMomentsRound'].badMaskPlanes=None

# Use round weight function?
config.processCcd.calibrate.measurement.undeblended['ext_shapeHSM_HsmSourceMomentsRound'].roundMoments=None

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['ext_shapeHSM_HsmPsfMoments'].doMeasure=True

config.processCcd.calibrate.measurement.undeblended.names=[]
# Run subtask to apply aperture correction
config.processCcd.calibrate.doApCorr=True

# flux measurement algorithms in getApCorrNameSet() to ignore; if a name is listed that does not appear in getApCorrNameSet() then a warning is logged
config.processCcd.calibrate.applyApCorr.ignoreList=[]

# set the general failure flag for a flux when it cannot be aperture-corrected?
config.processCcd.calibrate.applyApCorr.doFlagApCorrFailures=True

# flux measurement algorithms to be aperture-corrected by reference to another algorithm; this is a mapping alg1:alg2, where 'alg1' is the algorithm being corrected, and 'alg2' is the algorithm supplying the corrections
config.processCcd.calibrate.applyApCorr.proxies={}

# critical ratio of model to psf flux
config.processCcd.calibrate.catalogCalculation.plugins['base_ClassificationExtendedness'].fluxRatio=0.95

# correction factor for modelFlux error
config.processCcd.calibrate.catalogCalculation.plugins['base_ClassificationExtendedness'].modelErrFactor=0.0

# correction factor for psfFlux error
config.processCcd.calibrate.catalogCalculation.plugins['base_ClassificationExtendedness'].psfErrFactor=0.0

config.processCcd.calibrate.catalogCalculation.plugins.names=['base_ClassificationExtendedness', 'base_FootprintArea']
# Run fake sources injection task
config.processCcd.calibrate.doInsertFakes=False

# Mask plane to set on pixels affected by fakes.  Will be added if not already present.
config.processCcd.calibrate.insertFakes.maskPlaneName='FAKE'

# List of CCDs to ignore when processing
config.ignoreCcdList=[104, 105, 106, 107, 108, 109, 110, 111]

# DataId key corresponding to a single sensor
config.ccdKey='ccd'

