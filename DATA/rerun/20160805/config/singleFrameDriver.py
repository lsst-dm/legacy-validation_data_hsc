import lsst.pipe.drivers.singleFrameDriver
assert type(config)==lsst.pipe.drivers.singleFrameDriver.SingleFrameDriverConfig, 'config is of type %s.%s instead of lsst.pipe.drivers.singleFrameDriver.SingleFrameDriverConfig' % (type(config).__module__, type(config).__name__)
import lsst.meas.extensions.photometryKron.version
import eups.distrib.server
import lsst.meas.extensions.psfex.psfexStarSelector
import eups.distrib.tarball
import lsst.meas.extensions.photometryKron.kronLib
import eups.cmd
import lsst.pipe.tasks.setConfigFromEups
import eups.distrib.DistribFactory
import eups.db.Database
import eups.tags
import lsst.meas.extensions.shapeHSM.version
import eups.Eups
import eups.distrib.Distrib
import lsst.meas.extensions.photometryKron
import eups
import eups.stack.ProductStack
import eups.exceptions
import eups.distrib.Repositories
import eups.Uses
import ConfigParser
import lsst.meas.extensions.psfex.psfex
import eups.VersionParser
import eups.VersionCompare
import lsst.meas.extensions.psfex
import eups.hooks
import eups.distrib.eupspkg
import eups.db
import lsst.obs.subaru.crosstalk
import lsst.meas.extensions.psfex.psfexLib
import eups.stack
import lsst.obs.hsc.vignette
import configparser
import eups.app
import lsst.obs.subaru.isr
import eups.distrib
import lsst.meas.extensions.psfex.psfexPsfDeterminer
import eups.distrib.Repository
import lsst.meas.extensions
import lsst.afw.display.rgb
import optparse
import eups.db.VersionFile
import eups.table
import eups.distrib.pacman
import eups.distrib.builder
import lsst.meas.extensions.shapeHSM
import eups.lock
import eups.utils
import eups.Product
import eups.db.ChainFile
import lsst.meas.extensions.shapeHSM.hsmLib
import lsst.obs.subaru.crosstalkYagi
import eups.stack.ProductFamily
# List of CCDs to ignore when processing
config.ignoreCcdList=[]

import lsst.meas.extensions.photometryKron.version
import eups.distrib.server
import lsst.meas.extensions.psfex.psfexStarSelector
import eups.distrib.tarball
import lsst.meas.extensions.photometryKron.kronLib
import eups.cmd
import lsst.pipe.tasks.setConfigFromEups
import eups.distrib.DistribFactory
import eups.db.Database
import eups.tags
import lsst.meas.extensions.shapeHSM.version
import eups.Eups
import eups.distrib.Distrib
import lsst.meas.extensions.photometryKron
import eups
import eups.stack.ProductStack
import eups.exceptions
import eups.distrib.Repositories
import eups.Uses
import ConfigParser
import lsst.meas.extensions.psfex.psfex
import eups.VersionParser
import eups.VersionCompare
import lsst.meas.extensions.psfex
import eups.hooks
import eups.distrib.eupspkg
import eups.db
import lsst.obs.subaru.crosstalk
import lsst.meas.extensions.psfex.psfexLib
import eups.stack
import lsst.obs.hsc.vignette
import configparser
import eups.app
import lsst.obs.subaru.isr
import eups.distrib
import lsst.meas.extensions.psfex.psfexPsfDeterminer
import eups.distrib.Repository
import lsst.meas.extensions
import lsst.afw.display.rgb
import optparse
import eups.db.VersionFile
import eups.table
import eups.distrib.pacman
import eups.distrib.builder
import lsst.meas.extensions.shapeHSM
import eups.lock
import eups.utils
import eups.Product
import eups.db.ChainFile
import lsst.meas.extensions.shapeHSM.hsmLib
import lsst.obs.subaru.crosstalkYagi
import eups.stack.ProductFamily
# Perform calibration?
config.processCcd.doCalibrate=True

# Interpolate over defects? (ignored unless you provide a list of defects)
config.processCcd.charImage.repair.doInterpolate=True

# Type of statistic to calculate edge fallbackValue for interpolation
# Allowed values:
# 	None	Field is optional
# 	MEAN	mean
# 	MEDIAN	median
# 	USER	user value set in fallbackUserValue config
# 	MEANCLIP	clipped mean
# 
config.processCcd.charImage.repair.interp.fallbackValueType='MEANCLIP'

# Add a Gaussian to represent wings?
config.processCcd.charImage.repair.interp.modelPsf.addWing=True

# wing amplitude, as a multiple of core amplitude (dimensionless); ignored if addWing false
config.processCcd.charImage.repair.interp.modelPsf.wingAmplitude=0.1

# Minimum kernel size if using sizeFactor (pixels); ignored if size is not None
config.processCcd.charImage.repair.interp.modelPsf.minSize=5

# Maximum kernel size if using sizeFactor (pixels); ignored if size is not None
config.processCcd.charImage.repair.interp.modelPsf.maxSize=None

# Kernel size as a factor of fwhm (dimensionless); size = sizeFactor * fwhm; ignored if size is not None
config.processCcd.charImage.repair.interp.modelPsf.sizeFactor=3.0

# Default FWHM of Gaussian model of core of star (pixels)
config.processCcd.charImage.repair.interp.modelPsf.defaultFwhm=3.0

# wing width, as a multiple of core width (dimensionless); ignored if addWing false
config.processCcd.charImage.repair.interp.modelPsf.wingFwhmFactor=2.5

# Kernel size (width and height) (pixels); if None then sizeFactor is used
config.processCcd.charImage.repair.interp.modelPsf.size=None

# If fallbackValueType is 'USER' then use this as the fallbackValue; ignored otherwise
config.processCcd.charImage.repair.interp.fallbackUserValue=0.0

# Smoothly taper to the fallback value at the edge of the image?
config.processCcd.charImage.repair.interp.useFallbackValueAtEdge=True

# Allow negative values for egde interpolation fallbackValue?  If False, set fallbackValue to max(fallbackValue, 0.0)
config.processCcd.charImage.repair.interp.negativeFallbackAllowed=True

# Find and mask out cosmic rays?
config.processCcd.charImage.repair.doCosmicRay=True

# Don't interpolate over CR pixels
config.processCcd.charImage.repair.cosmicray.keepCRs=False

# used in condition 3 for CR; see CR.cc code
config.processCcd.charImage.repair.cosmicray.cond3_fac=2.5

# used in condition 3 for CR; see CR.cc code
config.processCcd.charImage.repair.cosmicray.cond3_fac2=0.4

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.processCcd.charImage.repair.cosmicray.background.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.processCcd.charImage.repair.cosmicray.background.approxOrderY=-1

# Names of mask planes to ignore while estimating the background
config.processCcd.charImage.repair.cosmicray.background.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# behaviour if there are too few points in grid for requested interpolation style
# Allowed values:
# 	None	Field is optional
# 	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
# 	THROW_EXCEPTION	throw an exception if there are too few points
# 	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
# 
config.processCcd.charImage.repair.cosmicray.background.undersampleStyle='REDUCE_INTERP_ORDER'

# how to interpolate the background values. This maps to an enum; see afw::math::Background
# Allowed values:
# 	NONE	No background estimation is to be attempted
# 	CONSTANT	Use a single constant value
# 	LINEAR	Use linear interpolation
# 	None	Field is optional
# 	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
# 	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
# 
config.processCcd.charImage.repair.cosmicray.background.algorithm='AKIMA_SPLINE'

# Ignore NaNs when estimating the background
config.processCcd.charImage.repair.cosmicray.background.isNanSafe=False

# how large a region of the sky should be used for each background point
# 	Valid Range = [1,inf)
config.processCcd.charImage.repair.cosmicray.background.binSize=100000

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.processCcd.charImage.repair.cosmicray.background.weighting=True

# type of statistic to use for grid points
# Allowed values:
# 	None	Field is optional
# 	MEDIAN	median
# 	MEANCLIP	clipped mean
# 	MEAN	unclipped mean
# 
config.processCcd.charImage.repair.cosmicray.background.statisticsProperty='MEDIAN'

# Use Approximate (Chebyshev) to model background.
config.processCcd.charImage.repair.cosmicray.background.useApprox=False

# number of times to look for contaminated pixels near known CR pixels
config.processCcd.charImage.repair.cosmicray.niteration=3

# maximum number of contaminated pixels
config.processCcd.charImage.repair.cosmicray.nCrPixelMax=1000000

# CRs must be > this many sky-sig above sky
config.processCcd.charImage.repair.cosmicray.minSigma=6.0

# CRs must have > this many DN (== electrons/gain) in initial detection
config.processCcd.charImage.repair.cosmicray.min_DN=150.0

# Run subtasks to measure and apply aperture corrections
config.processCcd.charImage.doApCorr=True

# Width and height of PSF model, in pixels. Must be odd.
# 	Valid Range = [1,inf)
config.processCcd.charImage.installSimplePsf.width=11

# Estimated FWHM of simple Gaussian PSF model, in pixels. Ignored if input exposure has a PSF model.
config.processCcd.charImage.installSimplePsf.fwhm=3.5322300675464238

# Write icExp and icExpBackground in addition to icSrc? Ignored if doWrite False.
config.processCcd.charImage.doWriteExposure=True

# Maximum allowed shift of WCS, due to matching (pixel)
# 	Valid Range = [-inf,4000)
config.processCcd.charImage.astrometry.matcher.maxOffsetPix=300

# Type of source flux; typically one of Ap or Psf
config.processCcd.charImage.astrometry.matcher.sourceFluxType='Ap'

# Number of bright stars to use
# 	Valid Range = [2,inf)
config.processCcd.charImage.astrometry.matcher.numBrightStars=50

# number of points to define a shape for matching
config.processCcd.charImage.astrometry.matcher.numPointsForShape=6

# Minimum allowed signal-to-noise ratio for sources used for matching (in the flux specified by sourceFluxType); <=0 for no limit
config.processCcd.charImage.astrometry.matcher.minSnr=40.0

# Allowed non-perpendicularity of x and y (degree)
# 	Valid Range = [-inf,45.0)
config.processCcd.charImage.astrometry.matcher.allowedNonperpDeg=3.0

# Maximum separation between reference objects and sources beyond which they will not be considered a match (arcsec)
# 	Valid Range = [0,inf)
config.processCcd.charImage.astrometry.matcher.maxMatchDistArcSec=3.0

# Minimum number of matched pairs as a fraction of the smaller of the number of reference stars or the number of good sources; the actual minimum is the smaller of this value or minMatchedPairs
# 	Valid Range = [0,1)
config.processCcd.charImage.astrometry.matcher.minFracMatchedPairs=0.3

# maximum determinant of linear transformation matrix for a usable solution
config.processCcd.charImage.astrometry.matcher.maxDeterminant=0.02

# Rotation angle allowed between sources and position reference objects (degrees)
# 	Valid Range = [-inf,6.0)
config.processCcd.charImage.astrometry.matcher.maxRotationDeg=1.0

# Minimum number of matched pairs; see also minFracMatchedPairs
# 	Valid Range = [2,inf)
config.processCcd.charImage.astrometry.matcher.minMatchedPairs=30

# the match distance below which further iteration is pointless (arcsec); ignored if not fitting a WCS
# 	Valid Range = [0,inf)
config.processCcd.charImage.astrometry.minMatchDistanceArcSec=0.001

# If True then load reference objects and match sources but do not fit a WCS;  this simply controls whether 'run' calls 'solve' or 'loadAndMatch'
config.processCcd.charImage.astrometry.forceKnownWcs=False

# maximum number of iterations of match sources and fit WCSignored if not fitting a WCS
# 	Valid Range = [1,inf)
config.processCcd.charImage.astrometry.maxIter=3

# the maximum match distance is set to  mean_match_distance + matchDistanceSigma*std_dev_match_distance; ignored if not fitting a WCS
# 	Valid Range = [0,inf)
config.processCcd.charImage.astrometry.matchDistanceSigma=2.0

# Number of standard deviations for clipping level
# 	Valid Range = [0.0,inf)
config.processCcd.charImage.astrometry.wcsFitter.rejSigma=3.0

# maximum median scatter of a WCS fit beyond which the fit fails (arcsec); be generous, as this is only intended to catch catastrophic failures
# 	Valid Range = [0,inf)
config.processCcd.charImage.astrometry.wcsFitter.maxScatterArcsec=10.0

# number of rejection iterations
# 	Valid Range = [0,inf)
config.processCcd.charImage.astrometry.wcsFitter.numRejIter=1

# number of iterations of fitter (which fits X and Y separately, and so benefits from a few iterations
# 	Valid Range = [1,inf)
config.processCcd.charImage.astrometry.wcsFitter.numIter=3

# order of SIP polynomial
# 	Valid Range = [0,inf)
config.processCcd.charImage.astrometry.wcsFitter.order=4

# Number of iterations of detect sources, measure sources, estimate PSF. If useSimplePsf is True then 2 should be plenty; otherwise more may be wanted.
# 	Valid Range = [1,inf)
config.processCcd.charImage.psfIterations=2

import lsst.meas.extensions.photometryKron.version
import lsst.meas.extensions.photometryKron.kronLib
import lsst.meas.extensions.photometryKron
# The seed multiplier value to use for random number generation.  0 will not set seed.
config.processCcd.charImage.measurement.noiseReplacer.noiseSeedMultiplier=1

# Add ann offset to the generated noise.
config.processCcd.charImage.measurement.noiseReplacer.noiseOffset=0.0

# How to choose mean and variance of the Gaussian noise we generate?
# Allowed values:
# 	variance	Mean = 0, variance = the image's variance
# 	meta	Mean = 0, variance = the "BGMEAN" metadata entry
# 	measure	Measure clipped mean and variance from the whole image
# 
config.processCcd.charImage.measurement.noiseReplacer.noiseSource='measure'

# the name of the flux measurement algorithm used for calibration
config.processCcd.charImage.measurement.slots.calibFlux='base_CircularApertureFlux_12_0'

# the name of the algorithm used to set the source aperture flux slot
config.processCcd.charImage.measurement.slots.apFlux='base_CircularApertureFlux_3_0'

# the name of the algorithm used to set the source inst flux slot
config.processCcd.charImage.measurement.slots.instFlux='base_GaussianFlux'

# the name of the algorithm used to set source moments parameters
config.processCcd.charImage.measurement.slots.shape='base_SdssShape'

# the name of the centroiding algorithm used to set source x,y
config.processCcd.charImage.measurement.slots.centroid='base_SdssCentroid'

# the name of the algorithm used to set the source model flux slot
config.processCcd.charImage.measurement.slots.modelFlux='base_GaussianFlux'

# the name of the algorithm used to set the source psf flux slot
config.processCcd.charImage.measurement.slots.psfFlux='base_PsfFlux'

# When measuring, replace other detected footprints with noise?
config.processCcd.charImage.measurement.doReplaceWithNoise=True

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_FPPosition'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_NaiveCentroid'].doMeasure=True

# Value to subtract from the image pixel values
config.processCcd.charImage.measurement.plugins['base_NaiveCentroid'].background=0.0

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_InputCount'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_PeakLikelihoodFlux'].doMeasure=True

# Name of warping kernel (e.g. "lanczos4") used to compute the peak
config.processCcd.charImage.measurement.plugins['base_PeakLikelihoodFlux'].warpingKernelName='lanczos4'

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_SdssCentroid'].doMeasure=True

# if the peak's less than this insist on binning at least once
config.processCcd.charImage.measurement.plugins['base_SdssCentroid'].peakMin=-1.0

# fiddle factor for adjusting the binning
config.processCcd.charImage.measurement.plugins['base_SdssCentroid'].wfac=1.5

# maximum allowed binning
config.processCcd.charImage.measurement.plugins['base_SdssCentroid'].binmax=16

# Largest aperture for which to use the slow, accurate, sinc aperture code
config.processCcd.charImage.measurement.plugins['ext_photometryKron_KronFlux'].maxSincRadius=10.0

# Number of times to iterate when setting the Kron radius
config.processCcd.charImage.measurement.plugins['ext_photometryKron_KronFlux'].nIterForRadius=1

# Use the Footprint size as part of initial estimate of Kron radius
config.processCcd.charImage.measurement.plugins['ext_photometryKron_KronFlux'].useFootprintRadius=False

# Minimum Kron radius (if == 0.0 use PSF's Kron radius) if enforceMinimumRadius. Also functions as fallback aperture radius if set.
config.processCcd.charImage.measurement.plugins['ext_photometryKron_KronFlux'].minimumRadius=0.0

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['ext_photometryKron_KronFlux'].doMeasure=True

# Multiplier of rms size for aperture used to initially estimate the Kron radius
config.processCcd.charImage.measurement.plugins['ext_photometryKron_KronFlux'].nSigmaForRadius=6.0

# If true check that the Kron radius exceeds some minimum
config.processCcd.charImage.measurement.plugins['ext_photometryKron_KronFlux'].enforceMinimumRadius=True

# if true, use existing shape and centroid measurements instead of fitting
config.processCcd.charImage.measurement.plugins['ext_photometryKron_KronFlux'].fixed=False

# Smooth image with N(0, smoothingSigma^2) Gaussian while estimating R_K
config.processCcd.charImage.measurement.plugins['ext_photometryKron_KronFlux'].smoothingSigma=-1.0

# Number of Kron radii for Kron flux
config.processCcd.charImage.measurement.plugins['ext_photometryKron_KronFlux'].nRadiusForFlux=2.5

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['ext_shapeHSM_HsmSourceMoments'].doMeasure=True

# Mask planes used to reject bad pixels.
config.processCcd.charImage.measurement.plugins['ext_shapeHSM_HsmSourceMoments'].badMaskPlanes=['BAD', 'SAT', 'INTRP']

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_CircularApertureFlux'].doMeasure=True

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.processCcd.charImage.measurement.plugins['base_CircularApertureFlux'].shiftKernel='lanczos5'

# Radius (in pixels) of apertures.
config.processCcd.charImage.measurement.plugins['base_CircularApertureFlux'].radii=[3.0, 4.5, 6.0, 9.0, 12.0, 17.0, 25.0, 35.0, 50.0, 70.0]

# Maximum radius (in pixels) for which the sinc algorithm should be used instead of the faster naive algorithm.  For elliptical apertures, this is the minor axis radius.
config.processCcd.charImage.measurement.plugins['base_CircularApertureFlux'].maxSincRadius=12.0

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_Jacobian'].doMeasure=True

# Nominal pixel size (arcsec)
config.processCcd.charImage.measurement.plugins['base_Jacobian'].pixelScale=0.168

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_GaussianCentroid'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_GaussianFlux'].doMeasure=True

# FIXME! NEVER DOCUMENTED!
config.processCcd.charImage.measurement.plugins['base_GaussianFlux'].background=0.0

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['ext_shapeHSM_HsmPsfMoments'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['ext_shapeHSM_HsmShapeBj'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.processCcd.charImage.measurement.plugins['ext_shapeHSM_HsmShapeBj'].badMaskPlanes=['BAD', 'SAT', 'INTRP']

# Field name for number of deblend children
config.processCcd.charImage.measurement.plugins['ext_shapeHSM_HsmShapeBj'].deblendNChild=''

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_ScaledApertureFlux'].doMeasure=True

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.processCcd.charImage.measurement.plugins['base_ScaledApertureFlux'].shiftKernel='lanczos5'

# Scaling factor of PSF FWHM for aperture radius.
config.processCcd.charImage.measurement.plugins['base_ScaledApertureFlux'].scale=3.14

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_PsfFlux'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.processCcd.charImage.measurement.plugins['base_PsfFlux'].badMaskPlanes=[]

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['ext_shapeHSM_HsmShapeKsb'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.processCcd.charImage.measurement.plugins['ext_shapeHSM_HsmShapeKsb'].badMaskPlanes=['BAD', 'SAT', 'INTRP']

# Field name for number of deblend children
config.processCcd.charImage.measurement.plugins['ext_shapeHSM_HsmShapeKsb'].deblendNChild=''

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_Blendedness'].doMeasure=True

# Whether to compute quantities related to the Gaussian-weighted shape
config.processCcd.charImage.measurement.plugins['base_Blendedness'].doShape=True

# Whether to compute HeavyFootprint dot products (the old deblend.blendedness parameter)
config.processCcd.charImage.measurement.plugins['base_Blendedness'].doOld=True

# Whether to compute quantities related to the Gaussian-weighted flux
config.processCcd.charImage.measurement.plugins['base_Blendedness'].doFlux=True

# Radius factor that sets the maximum extent of the weight function (and hence the flux measurements)
config.processCcd.charImage.measurement.plugins['base_Blendedness'].nSigmaWeightMax=3.0

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['ext_shapeHSM_HsmShapeLinear'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.processCcd.charImage.measurement.plugins['ext_shapeHSM_HsmShapeLinear'].badMaskPlanes=['BAD', 'SAT', 'INTRP']

# Field name for number of deblend children
config.processCcd.charImage.measurement.plugins['ext_shapeHSM_HsmShapeLinear'].deblendNChild=''

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_PeakCentroid'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_PixelFlags'].doMeasure=True

# List of mask planes to be searched for which occur in the center of a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.processCcd.charImage.measurement.plugins['base_PixelFlags'].masksFpCenter=[]

# List of mask planes to be searched for which occur anywhere within a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.processCcd.charImage.measurement.plugins['base_PixelFlags'].masksFpAnywhere=[]

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_SkyCoord'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_Variance'].doMeasure=True

# Scale factor to apply to shape for aperture
config.processCcd.charImage.measurement.plugins['base_Variance'].scale=5.0

# Mask planes to ignore
config.processCcd.charImage.measurement.plugins['base_Variance'].mask=['DETECTED', 'DETECTED_NEGATIVE', 'BAD', 'SAT']

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['ext_shapeHSM_HsmShapeRegauss'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.processCcd.charImage.measurement.plugins['ext_shapeHSM_HsmShapeRegauss'].badMaskPlanes=['BAD', 'SAT', 'INTRP']

# Field name for number of deblend children
config.processCcd.charImage.measurement.plugins['ext_shapeHSM_HsmShapeRegauss'].deblendNChild=''

# Maximum centroid shift, limited to 2-10
config.processCcd.charImage.measurement.plugins['base_SdssShape'].maxShift=0.0

# Convergence tolerance for FWHM
config.processCcd.charImage.measurement.plugins['base_SdssShape'].tol2=9.999999747378752e-05

# Convergence tolerance for e1,e2
config.processCcd.charImage.measurement.plugins['base_SdssShape'].tol1=9.999999747378752e-06

# Whether to also compute the shape of the PSF model
config.processCcd.charImage.measurement.plugins['base_SdssShape'].doMeasurePsf=True

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_SdssShape'].doMeasure=True

# Additional value to add to background
config.processCcd.charImage.measurement.plugins['base_SdssShape'].background=0.0

# Maximum number of iterations
config.processCcd.charImage.measurement.plugins['base_SdssShape'].maxIter=100

config.processCcd.charImage.measurement.plugins.names=['base_CircularApertureFlux', 'base_FPPosition', 'base_PixelFlags', 'base_PsfFlux', 'ext_photometryKron_KronFlux', 'base_Jacobian', 'base_SdssCentroid', 'base_GaussianFlux', 'base_SdssShape']
# Measure PSF? If False then keep the existing PSF model (which must exist) and use that model for all operations.
config.processCcd.charImage.doMeasurePsf=True

# Run deblender input exposure
config.processCcd.charImage.doDeblend=False

# correction factor for psfFlux error
config.processCcd.charImage.afterburners.plugins['base_ClassificationExtendedness'].psfErrFactor=0.0

# correction factor for modelFlux error
config.processCcd.charImage.afterburners.plugins['base_ClassificationExtendedness'].modelErrFactor=0.0

# critical ratio of model to psf flux
config.processCcd.charImage.afterburners.plugins['base_ClassificationExtendedness'].fluxRatio=0.95

config.processCcd.charImage.afterburners.plugins.names=['base_ClassificationExtendedness']
# Estimate the background again after final source detection?
config.processCcd.charImage.detection.reEstimateBackground=True

# Grow detections by nSigmaToGrow * sigma; if 0 then do not grow
config.processCcd.charImage.detection.nSigmaToGrow=2.4

# detected sources with fewer than the specified number of pixels will be ignored
# 	Valid Range = [0,inf)
config.processCcd.charImage.detection.minPixels=1

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.processCcd.charImage.detection.tempLocalBackground.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.processCcd.charImage.detection.tempLocalBackground.approxOrderY=-1

# Names of mask planes to ignore while estimating the background
config.processCcd.charImage.detection.tempLocalBackground.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# behaviour if there are too few points in grid for requested interpolation style
# Allowed values:
# 	None	Field is optional
# 	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
# 	THROW_EXCEPTION	throw an exception if there are too few points
# 	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
# 
config.processCcd.charImage.detection.tempLocalBackground.undersampleStyle='REDUCE_INTERP_ORDER'

# how to interpolate the background values. This maps to an enum; see afw::math::Background
# Allowed values:
# 	NONE	No background estimation is to be attempted
# 	CONSTANT	Use a single constant value
# 	LINEAR	Use linear interpolation
# 	None	Field is optional
# 	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
# 	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
# 
config.processCcd.charImage.detection.tempLocalBackground.algorithm='AKIMA_SPLINE'

# Ignore NaNs when estimating the background
config.processCcd.charImage.detection.tempLocalBackground.isNanSafe=False

# how large a region of the sky should be used for each background point
# 	Valid Range = [1,inf)
config.processCcd.charImage.detection.tempLocalBackground.binSize=64

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.processCcd.charImage.detection.tempLocalBackground.weighting=True

# type of statistic to use for grid points
# Allowed values:
# 	None	Field is optional
# 	MEDIAN	median
# 	MEANCLIP	clipped mean
# 	MEAN	unclipped mean
# 
config.processCcd.charImage.detection.tempLocalBackground.statisticsProperty='MEANCLIP'

# Use Approximate (Chebyshev) to model background.
config.processCcd.charImage.detection.tempLocalBackground.useApprox=False

# Include threshold relative to thresholdValue
# 	Valid Range = [0.0,inf)
config.processCcd.charImage.detection.includeThresholdMultiplier=10.0

# Pixels should be grown as isotropically as possible (slower)
config.processCcd.charImage.detection.isotropicGrow=True

# Fiddle factor to add to the background; debugging only
config.processCcd.charImage.detection.adjustBackground=0.0

# Do temporary interpolated background subtraction before footprint detection?
config.processCcd.charImage.detection.doTempLocalBackground=False

# specifies whether to detect positive, or negative sources, or both
# Allowed values:
# 	positive	detect only positive sources
# 	negative	detect only negative sources
# 	both	detect both positive and negative sources
# 
config.processCcd.charImage.detection.thresholdPolarity='positive'

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.processCcd.charImage.detection.background.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.processCcd.charImage.detection.background.approxOrderY=-1

# Names of mask planes to ignore while estimating the background
config.processCcd.charImage.detection.background.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# behaviour if there are too few points in grid for requested interpolation style
# Allowed values:
# 	None	Field is optional
# 	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
# 	THROW_EXCEPTION	throw an exception if there are too few points
# 	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
# 
config.processCcd.charImage.detection.background.undersampleStyle='REDUCE_INTERP_ORDER'

# how to interpolate the background values. This maps to an enum; see afw::math::Background
# Allowed values:
# 	NONE	No background estimation is to be attempted
# 	CONSTANT	Use a single constant value
# 	LINEAR	Use linear interpolation
# 	None	Field is optional
# 	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
# 	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
# 
config.processCcd.charImage.detection.background.algorithm='NATURAL_SPLINE'

# Ignore NaNs when estimating the background
config.processCcd.charImage.detection.background.isNanSafe=False

# how large a region of the sky should be used for each background point
# 	Valid Range = [1,inf)
config.processCcd.charImage.detection.background.binSize=128

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.processCcd.charImage.detection.background.weighting=True

# type of statistic to use for grid points
# Allowed values:
# 	None	Field is optional
# 	MEDIAN	median
# 	MEANCLIP	clipped mean
# 	MEAN	unclipped mean
# 
config.processCcd.charImage.detection.background.statisticsProperty='MEANCLIP'

# Use Approximate (Chebyshev) to model background.
config.processCcd.charImage.detection.background.useApprox=True

# Grow detections to set the image mask bits, but return the original (not-grown) footprints
config.processCcd.charImage.detection.returnOriginalFootprints=False

# specifies the desired flavor of Threshold
# Allowed values:
# 	pixel_stdev	threshold applied to per-pixel std deviation
# 	variance	threshold applied to image variance
# 	value	threshold applied to image value
# 	stdev	threshold applied to image std deviation
# 
config.processCcd.charImage.detection.thresholdType='stdev'

# Threshold for footprints
# 	Valid Range = [0.0,inf)
config.processCcd.charImage.detection.thresholdValue=5.0

# Default reference catalog filter to use if filter not specified in exposure; if blank then filter must be specified in exposure
config.processCcd.charImage.refObjLoader.defaultFilter=''

# Padding to add to 4 all edges of the bounding box (pixels)
# 	Valid Range = [0,inf)
config.processCcd.charImage.refObjLoader.pixelMargin=50

# Mapping of camera filter name: reference catalog filter name; each reference filter must exist
config.processCcd.charImage.refObjLoader.filterMap={}

# Strictness of Astropy unit compatibility check, can be 'raise', 'warn' or 'silent'
config.processCcd.charImage.checkUnitsParseStrict='raise'

# flux measurement algorithms in getApCorrNameSet() to ignore; if a name is listed that does not appear in getApCorrNameSet() then a warning is logged
config.processCcd.charImage.applyApCorr.ignoreList=[]

# set the general failure flag for a flux when it cannot be aperture-corrected?
config.processCcd.charImage.applyApCorr.doFlagApCorrFailures=True

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.processCcd.charImage.background.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.processCcd.charImage.background.approxOrderY=-1

# Names of mask planes to ignore while estimating the background
config.processCcd.charImage.background.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# behaviour if there are too few points in grid for requested interpolation style
# Allowed values:
# 	None	Field is optional
# 	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
# 	THROW_EXCEPTION	throw an exception if there are too few points
# 	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
# 
config.processCcd.charImage.background.undersampleStyle='REDUCE_INTERP_ORDER'

# how to interpolate the background values. This maps to an enum; see afw::math::Background
# Allowed values:
# 	NONE	No background estimation is to be attempted
# 	CONSTANT	Use a single constant value
# 	LINEAR	Use linear interpolation
# 	None	Field is optional
# 	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
# 	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
# 
config.processCcd.charImage.background.algorithm='NATURAL_SPLINE'

# Ignore NaNs when estimating the background
config.processCcd.charImage.background.isNanSafe=False

# how large a region of the sky should be used for each background point
# 	Valid Range = [1,inf)
config.processCcd.charImage.background.binSize=128

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.processCcd.charImage.background.weighting=True

# type of statistic to use for grid points
# Allowed values:
# 	None	Field is optional
# 	MEDIAN	median
# 	MEANCLIP	clipped mean
# 	MEAN	unclipped mean
# 
config.processCcd.charImage.background.statisticsProperty='MEANCLIP'

# Use Approximate (Chebyshev) to model background.
config.processCcd.charImage.background.useApprox=True

# Persist results?
config.processCcd.charImage.doWrite=True

# What to do when a peak to be deblended is close to the edge of the image
# Allowed values:
# 	ramp	Ramp down flux at the image edge by the PSF
# 	noclip	Ignore the edge when building the symmetric template.
# 	clip	Clip the template at the edge AND the mirror of the edge.
# 	None	Field is optional
# 
config.processCcd.charImage.deblend.edgeHandling='ramp'

# Assign stray flux to deblend children.  Implies findStrayFlux.
config.processCcd.charImage.deblend.assignStrayFlux=True

# Mask planes with the corresponding limit on the fraction of masked pixels. Sources violating this limit will not be deblended.
config.processCcd.charImage.deblend.maskLimits={'NO_DATA': 0.25}

# How to split flux among peaks
# Allowed values:
# 	trim	Shrink the parent footprint to pixels that are not assigned to children
# 	r-to-peak	~ 1/(1+R^2) to the peak
# 	r-to-footprint	~ 1/(1+R^2) to the closest pixel in the footprint.  CAUTION: this can be computationally expensive on large footprints!
# 	nearest-footprint	Assign 100% to the nearest footprint (using L-1 norm aka Manhattan distance)
# 	None	Field is optional
# 
config.processCcd.charImage.deblend.strayFluxRule='trim'

# If True, catch exceptions thrown by the deblender, log them, and set a flag on the parent, instead of letting them propagate up
config.processCcd.charImage.deblend.catchFailures=False

# When the deblender should attribute stray flux to point sources
# Allowed values:
# 	always	Always
# 	never	Never; stray flux will not be attributed to any deblended child if the deblender thinks all peaks look like point sources
# 	necessary	When there is not an extended object in the footprint
# 	None	Field is optional
# 
config.processCcd.charImage.deblend.strayFluxToPointSources='necessary'

# Chi-squared per DOF cut for deciding a source is a PSF during deblending (shifted PSF model #2)
config.processCcd.charImage.deblend.psfChisq2b=1.5

# Mask planes to ignore when performing statistics
config.processCcd.charImage.deblend.maskPlanes=['SAT', 'INTRP', 'NO_DATA']

# Maximum area for footprints before they are ignored as large; non-positive means no threshold applied
config.processCcd.charImage.deblend.maxFootprintArea=10000

# Minimum axis ratio for footprints before they are ignored as large; non-positive means no threshold applied
config.processCcd.charImage.deblend.minFootprintAxisRatio=0.0

# Maximum linear dimension for footprints before they are ignored as large; non-positive means no threshold applied
config.processCcd.charImage.deblend.maxFootprintSize=0

# Chi-squared per DOF cut for deciding a source is a PSF during deblending (un-shifted PSF model)
config.processCcd.charImage.deblend.psfChisq1=1.5

# Find stray flux---flux not claimed by any child in the deblender.
config.processCcd.charImage.deblend.findStrayFlux=True

# Footprints smaller in width or height than this value will be ignored; minimum of 2 due to PSF gradient calculation.
# 	Valid Range = [2,inf)
config.processCcd.charImage.deblend.tinyFootprintSize=2

# Chi-squared per DOF cut for deciding a source is PSF during deblending (shifted PSF model)
config.processCcd.charImage.deblend.psfChisq2=1.5

# Guarantee that all peaks produce a child source.
config.processCcd.charImage.deblend.propagateAllPeaks=False

# Mask name for footprints not deblended, or None
config.processCcd.charImage.deblend.notDeblendedMask='NOT_DEBLENDED'

# Only deblend the brightest maxNumberOfPeaks peaks in the parent (<= 0: unlimited)
config.processCcd.charImage.deblend.maxNumberOfPeaks=0

# When splitting stray flux, clip fractions below this value to zero.
config.processCcd.charImage.deblend.clipStrayFluxFraction=0.001

# Replace the existing PSF model with a simplified version that has the same sigma at the start of each PSF determination iteration? Doing so makes PSF determination converge more robustly and quickly.
config.processCcd.charImage.useSimplePsf=True

# Field name prefix for the flux other measurements should be aperture corrected to match
config.processCcd.charImage.measureApCorr.refFluxName='slot_CalibFlux'

# Minimum number of degrees of freedom (# of valid data points - # of parameters); if this is exceeded, the order of the fit is decreased (in both dimensions), and if we can't decrease it enough, we'll raise ValueError.
# 	Valid Range = [1,inf)
config.processCcd.charImage.measureApCorr.minDegreesOfFreedom=1

# Standard deviation of width allowed to be interpreted as good stars
config.processCcd.charImage.measureApCorr.starSelector['objectSize'].widthStdAllowed=0.15

# maximum width to include in histogram
config.processCcd.charImage.measureApCorr.starSelector['objectSize'].widthMax=10.0

# minimum width to include in histogram
config.processCcd.charImage.measureApCorr.starSelector['objectSize'].widthMin=0.0

# size of the kernel to create
config.processCcd.charImage.measureApCorr.starSelector['objectSize'].kernelSize=21

# specify the minimum psfFlux for good Psf Candidates
config.processCcd.charImage.measureApCorr.starSelector['objectSize'].fluxMin=12500.0

# number of pixels to ignore around the edge of PSF candidate postage stamps
config.processCcd.charImage.measureApCorr.starSelector['objectSize'].borderWidth=0

# List of flags which cause a source to be rejected as bad
config.processCcd.charImage.measureApCorr.starSelector['objectSize'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter', 'base_PixelFlags_flag_bad', 'base_PixelFlags_flag_interpolated']

# specify the maximum psfFlux for good Psf Candidates (ignored if == 0)
config.processCcd.charImage.measureApCorr.starSelector['objectSize'].fluxMax=0.0

# Name of field in Source to use for flux measurement
config.processCcd.charImage.measureApCorr.starSelector['objectSize'].sourceFluxField='base_GaussianFlux_flux'

# Keep objects within this many sigma of cluster 0's median
config.processCcd.charImage.measureApCorr.starSelector['objectSize'].nSigmaClip=2.0

# specify the minimum psfFlux for good Psf Candidates
# 	Valid Range = [0.0,inf)
config.processCcd.charImage.measureApCorr.starSelector['catalog'].fluxLim=0.0

# List of flags which cause a source to be rejected as bad
config.processCcd.charImage.measureApCorr.starSelector['catalog'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter']

# specify the maximum psfFlux for good Psf Candidates (ignored if == 0)
# 	Valid Range = [0.0,inf)
config.processCcd.charImage.measureApCorr.starSelector['catalog'].fluxMax=0.0

# size of the kernel to create
config.processCcd.charImage.measureApCorr.starSelector['catalog'].kernelSize=21

# number of pixels to ignore around the edge of PSF candidate postage stamps
config.processCcd.charImage.measureApCorr.starSelector['catalog'].borderWidth=0

# size of the kernel to create
config.processCcd.charImage.measureApCorr.starSelector['secondMoment'].kernelSize=21

# Number of bins in moment histogram
config.processCcd.charImage.measureApCorr.starSelector['secondMoment'].histSize=64

# Clipping threshold for moments histogram range
config.processCcd.charImage.measureApCorr.starSelector['secondMoment'].histMomentClip=5.0

# Multiplier of mean for minimum moments histogram range
config.processCcd.charImage.measureApCorr.starSelector['secondMoment'].histMomentMinMultiplier=2.0

# List of flags which cause a source to be rejected as bad
config.processCcd.charImage.measureApCorr.starSelector['secondMoment'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter']

# number of pixels to ignore around the edge of PSF candidate postage stamps
config.processCcd.charImage.measureApCorr.starSelector['secondMoment'].borderWidth=0

# Multiplier of mean for maximum moments histogram range
config.processCcd.charImage.measureApCorr.starSelector['secondMoment'].histMomentMaxMultiplier=5.0

# specify the maximum psfFlux for good Psf Candidates (ignored if == 0)
config.processCcd.charImage.measureApCorr.starSelector['secondMoment'].fluxMax=0.0

# candidate PSF's shapes must lie within this many sigma of the average shape
config.processCcd.charImage.measureApCorr.starSelector['secondMoment'].clumpNSigma=2.0

# specify the minimum psfFlux for good Psf Candidates
config.processCcd.charImage.measureApCorr.starSelector['secondMoment'].fluxLim=12500.0

# Maximum moment to consider
config.processCcd.charImage.measureApCorr.starSelector['secondMoment'].histMomentMax=100.0

# Name of a flag field that is True for stars that should be used.
config.processCcd.charImage.measureApCorr.starSelector['flagged'].field='calib_psfUsed'

# List of flags which cause a source to be rejected as bad
config.processCcd.charImage.measureApCorr.starSelector['flagged'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter', 'base_PixelFlags_flag_bad', 'base_PixelFlags_flag_interpolated']

# number of pixels to ignore around the edge of PSF candidate postage stamps
config.processCcd.charImage.measureApCorr.starSelector['flagged'].borderWidth=0

# size of the kernel to create
config.processCcd.charImage.measureApCorr.starSelector['flagged'].kernelSize=21

# Filter bad pixels? 
config.processCcd.charImage.measureApCorr.starSelector['psfex'].maxbadflag=True

# List of flags which cause a source to be rejected as bad
config.processCcd.charImage.measureApCorr.starSelector['psfex'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter', 'base_PixelFlags_flag_bad', 'base_PixelFlags_flag_suspectCenter', 'base_PsfFlux_flag']

# Allowed FWHM variability (1.0 = 100%)
config.processCcd.charImage.measureApCorr.starSelector['psfex'].maxFwhmVariability=0.2

# size of the kernel to create
config.processCcd.charImage.measureApCorr.starSelector['psfex'].kernelSize=21

# Maximum (A-B)/(A+B) 
config.processCcd.charImage.measureApCorr.starSelector['psfex'].maxellip=0.3

# Maximum allowed FWHM 
config.processCcd.charImage.measureApCorr.starSelector['psfex'].minFwhm=2.0

# number of pixels to ignore around the edge of PSF candidate postage stamps
config.processCcd.charImage.measureApCorr.starSelector['psfex'].borderWidth=0

# Minimum S/N for candidates
config.processCcd.charImage.measureApCorr.starSelector['psfex'].minsn=100.0

# Name of phot. flux err. key
config.processCcd.charImage.measureApCorr.starSelector['psfex'].fluxErrName=''

# Name of photometric flux key 
config.processCcd.charImage.measureApCorr.starSelector['psfex'].fluxName='base_PsfFlux'

# Max number of bad pixels 
config.processCcd.charImage.measureApCorr.starSelector['psfex'].maxbad=0

# Minimum allowed FWHM 
config.processCcd.charImage.measureApCorr.starSelector['psfex'].maxFwhm=10.0

config.processCcd.charImage.measureApCorr.starSelector.name='flagged'
# Number of standard devisations to clip at
config.processCcd.charImage.measureApCorr.numSigmaClip=3.0

# if true, only include terms where the sum of the x and y order is less than or equal to max(orderX, orderY)
config.processCcd.charImage.measureApCorr.fitConfig.triangular=True

# maximum Chebyshev function order in x
config.processCcd.charImage.measureApCorr.fitConfig.orderX=2

# maximum Chebyshev function order in y
config.processCcd.charImage.measureApCorr.fitConfig.orderY=2

# Number of iterations for sigma clipping
config.processCcd.charImage.measureApCorr.numIter=4

# Maximum radius of the kernel
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].kernelSizeMax=45

# floor for variance is lam*data
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].lam=0.05

# for psf candidate evaluation
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].reducedChi2ForPsfCandidates=2.0

# Use non-linear fitter for spatial variation of Kernel
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].nonLinearSpatialFit=False

# Mask blends in image?
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].doMaskBlends=True

# size of cell used to determine PSF (pixels, column direction)
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].sizeCellX=256

# size of cell used to determine PSF (pixels, row direction)
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].sizeCellY=256

# Rejection threshold (stdev) for candidates based on spatial fit
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].spatialReject=3.0

# radius of the kernel to create, relative to the square root of the stellar quadrupole moments
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].kernelSize=10.0

# Reject candidates that are blended?
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].doRejectBlends=False

# Should each PSF candidate be given the same weight, independent of magnitude?
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].constantWeight=True

# number of iterations of PSF candidate star list
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].nIterForPsf=3

# number of stars per psf Cell for spatial fitting
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].nStarPerCellSpatialFit=5

# Number of pixels to ignore around the edge of PSF candidate postage stamps
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].borderWidth=0

# number of eigen components for PSF kernel creation
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].nEigenComponents=4

# Threshold (stdev) for rejecting extraneous pixels around candidate; applied if positive
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].pixelThreshold=0.0

# specify spatial order for PSF kernel creation
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].spatialOrder=2

# tolerance of spatial fitting
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].tolerance=0.01

# Minimum radius of the kernel
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].kernelSizeMin=25

# number of stars per psf cell for PSF kernel creation
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].nStarPerCell=3

# floor for variance is lam*data
config.processCcd.charImage.measurePsf.psfDeterminer['psfex'].lam=0.05

# for psf candidate evaluation
config.processCcd.charImage.measurePsf.psfDeterminer['psfex'].reducedChi2ForPsfCandidates=2.0

# specify spatial order for PSF kernel creation
config.processCcd.charImage.measurePsf.psfDeterminer['psfex'].spatialOrder=2

# number of eigen components for PSF kernel creation
config.processCcd.charImage.measurePsf.psfDeterminer['psfex']._PsfexPsfDeterminerConfig__nEigenComponents=4

# Should PSFEX be permitted to recentroid PSF candidates?
config.processCcd.charImage.measurePsf.psfDeterminer['psfex'].recentroid=False

# Resolution of the internal PSF model relative to the pixel size; e.g. 0.5 is equal to 2x oversampling
config.processCcd.charImage.measurePsf.psfDeterminer['psfex'].samplingSize=1.0

# size of cell used to determine PSF (pixels, row direction)
config.processCcd.charImage.measurePsf.psfDeterminer['psfex'].sizeCellY=256

# List of mask bits which cause a source to be rejected as bad N.b. INTRP is used specially in PsfCandidateSet; it means "Contaminated by neighbour"
config.processCcd.charImage.measurePsf.psfDeterminer['psfex'].badMaskBits=['INTRP', 'SAT']

# Rejection threshold (stdev) for candidates based on spatial fit
config.processCcd.charImage.measurePsf.psfDeterminer['psfex'].spatialReject=3.0

# Number of pixels to ignore around the edge of PSF candidate postage stamps
config.processCcd.charImage.measurePsf.psfDeterminer['psfex']._PsfexPsfDeterminerConfig__borderWidth=0

# radius of the kernel to create, relative to the square root of the stellar quadrupole moments
config.processCcd.charImage.measurePsf.psfDeterminer['psfex'].kernelSize=41.0

# size of cell used to determine PSF (pixels, column direction)
config.processCcd.charImage.measurePsf.psfDeterminer['psfex'].sizeCellX=256

# Should each PSF candidate be given the same weight, independent of magnitude?
config.processCcd.charImage.measurePsf.psfDeterminer['psfex']._PsfexPsfDeterminerConfig__constantWeight=True

# Maximum radius of the kernel
config.processCcd.charImage.measurePsf.psfDeterminer['psfex'].kernelSizeMax=45

# number of stars per psf Cell for spatial fitting
config.processCcd.charImage.measurePsf.psfDeterminer['psfex']._PsfexPsfDeterminerConfig__nStarPerCellSpatialFit=5

# number of stars per psf cell for PSF kernel creation
config.processCcd.charImage.measurePsf.psfDeterminer['psfex']._PsfexPsfDeterminerConfig__nStarPerCell=3

# number of iterations of PSF candidate star list
config.processCcd.charImage.measurePsf.psfDeterminer['psfex']._PsfexPsfDeterminerConfig__nIterForPsf=3

# tolerance of spatial fitting
config.processCcd.charImage.measurePsf.psfDeterminer['psfex'].tolerance=0.01

# Minimum radius of the kernel
config.processCcd.charImage.measurePsf.psfDeterminer['psfex'].kernelSizeMin=25

config.processCcd.charImage.measurePsf.psfDeterminer.name='psfex'
# Standard deviation of width allowed to be interpreted as good stars
config.processCcd.charImage.measurePsf.starSelector['objectSize'].widthStdAllowed=0.15

# maximum width to include in histogram
config.processCcd.charImage.measurePsf.starSelector['objectSize'].widthMax=10.0

# minimum width to include in histogram
config.processCcd.charImage.measurePsf.starSelector['objectSize'].widthMin=0.9

# size of the kernel to create
config.processCcd.charImage.measurePsf.starSelector['objectSize'].kernelSize=21

# specify the minimum psfFlux for good Psf Candidates
config.processCcd.charImage.measurePsf.starSelector['objectSize'].fluxMin=4000.0

# number of pixels to ignore around the edge of PSF candidate postage stamps
config.processCcd.charImage.measurePsf.starSelector['objectSize'].borderWidth=0

# List of flags which cause a source to be rejected as bad
config.processCcd.charImage.measurePsf.starSelector['objectSize'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter', 'base_PixelFlags_flag_bad', 'base_PixelFlags_flag_interpolated']

# specify the maximum psfFlux for good Psf Candidates (ignored if == 0)
config.processCcd.charImage.measurePsf.starSelector['objectSize'].fluxMax=0.0

# Name of field in Source to use for flux measurement
config.processCcd.charImage.measurePsf.starSelector['objectSize'].sourceFluxField='base_PsfFlux_flux'

# Keep objects within this many sigma of cluster 0's median
config.processCcd.charImage.measurePsf.starSelector['objectSize'].nSigmaClip=2.0

# specify the minimum psfFlux for good Psf Candidates
# 	Valid Range = [0.0,inf)
config.processCcd.charImage.measurePsf.starSelector['catalog'].fluxLim=0.0

# List of flags which cause a source to be rejected as bad
config.processCcd.charImage.measurePsf.starSelector['catalog'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter']

# specify the maximum psfFlux for good Psf Candidates (ignored if == 0)
# 	Valid Range = [0.0,inf)
config.processCcd.charImage.measurePsf.starSelector['catalog'].fluxMax=0.0

# size of the kernel to create
config.processCcd.charImage.measurePsf.starSelector['catalog'].kernelSize=21

# number of pixels to ignore around the edge of PSF candidate postage stamps
config.processCcd.charImage.measurePsf.starSelector['catalog'].borderWidth=0

# size of the kernel to create
config.processCcd.charImage.measurePsf.starSelector['secondMoment'].kernelSize=21

# Number of bins in moment histogram
config.processCcd.charImage.measurePsf.starSelector['secondMoment'].histSize=64

# Clipping threshold for moments histogram range
config.processCcd.charImage.measurePsf.starSelector['secondMoment'].histMomentClip=5.0

# Multiplier of mean for minimum moments histogram range
config.processCcd.charImage.measurePsf.starSelector['secondMoment'].histMomentMinMultiplier=2.0

# List of flags which cause a source to be rejected as bad
config.processCcd.charImage.measurePsf.starSelector['secondMoment'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter']

# number of pixels to ignore around the edge of PSF candidate postage stamps
config.processCcd.charImage.measurePsf.starSelector['secondMoment'].borderWidth=0

# Multiplier of mean for maximum moments histogram range
config.processCcd.charImage.measurePsf.starSelector['secondMoment'].histMomentMaxMultiplier=5.0

# specify the maximum psfFlux for good Psf Candidates (ignored if == 0)
config.processCcd.charImage.measurePsf.starSelector['secondMoment'].fluxMax=0.0

# candidate PSF's shapes must lie within this many sigma of the average shape
config.processCcd.charImage.measurePsf.starSelector['secondMoment'].clumpNSigma=2.0

# specify the minimum psfFlux for good Psf Candidates
config.processCcd.charImage.measurePsf.starSelector['secondMoment'].fluxLim=12500.0

# Maximum moment to consider
config.processCcd.charImage.measurePsf.starSelector['secondMoment'].histMomentMax=100.0

# Name of a flag field that is True for stars that should be used.
config.processCcd.charImage.measurePsf.starSelector['flagged'].field='calib_psfUsed'

# List of flags which cause a source to be rejected as bad
config.processCcd.charImage.measurePsf.starSelector['flagged'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter', 'base_PixelFlags_flag_bad', 'base_PixelFlags_flag_interpolated']

# number of pixels to ignore around the edge of PSF candidate postage stamps
config.processCcd.charImage.measurePsf.starSelector['flagged'].borderWidth=0

# size of the kernel to create
config.processCcd.charImage.measurePsf.starSelector['flagged'].kernelSize=21

# Filter bad pixels? 
config.processCcd.charImage.measurePsf.starSelector['psfex'].maxbadflag=True

# List of flags which cause a source to be rejected as bad
config.processCcd.charImage.measurePsf.starSelector['psfex'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter', 'base_PixelFlags_flag_bad', 'base_PixelFlags_flag_suspectCenter', 'base_PsfFlux_flag']

# Allowed FWHM variability (1.0 = 100%)
config.processCcd.charImage.measurePsf.starSelector['psfex'].maxFwhmVariability=0.2

# size of the kernel to create
config.processCcd.charImage.measurePsf.starSelector['psfex'].kernelSize=21

# Maximum (A-B)/(A+B) 
config.processCcd.charImage.measurePsf.starSelector['psfex'].maxellip=0.3

# Maximum allowed FWHM 
config.processCcd.charImage.measurePsf.starSelector['psfex'].minFwhm=2.0

# number of pixels to ignore around the edge of PSF candidate postage stamps
config.processCcd.charImage.measurePsf.starSelector['psfex'].borderWidth=0

# Minimum S/N for candidates
config.processCcd.charImage.measurePsf.starSelector['psfex'].minsn=100.0

# Name of phot. flux err. key
config.processCcd.charImage.measurePsf.starSelector['psfex'].fluxErrName=''

# Name of photometric flux key 
config.processCcd.charImage.measurePsf.starSelector['psfex'].fluxName='base_PsfFlux'

# Max number of bad pixels 
config.processCcd.charImage.measurePsf.starSelector['psfex'].maxbad=0

# Minimum allowed FWHM 
config.processCcd.charImage.measurePsf.starSelector['psfex'].maxFwhm=10.0

config.processCcd.charImage.measurePsf.starSelector.name='objectSize'
# This number will be multiplied by the exposure ID to set the random seed for reserving candidates
config.processCcd.charImage.measurePsf.reserveSeed=1

# Fraction of PSF candidates to reserve from fitting; none if <= 0
config.processCcd.charImage.measurePsf.reserveFraction=0.2

# Perform astrometric calibration?
config.processCcd.calibrate.doAstrometry=True

# Run subtask to apply aperture correction
config.processCcd.calibrate.doApCorr=True

# Include HeavyFootprint data in source table? If false then heavy footprints are saved as normal footprints, which saves some space
config.processCcd.calibrate.doWriteHeavyFootprintsInSources=True

# Perform phometric calibration?
config.processCcd.calibrate.doPhotoCal=True

# Raise an exception if astrometry fails? Ignored if doAstrometry false.
config.processCcd.calibrate.requireAstrometry=True

import lsst.meas.extensions.shapeHSM
import lsst.meas.extensions.shapeHSM.version
import lsst.meas.extensions.shapeHSM.hsmLib
# The seed multiplier value to use for random number generation.  0 will not set seed.
config.processCcd.calibrate.measurement.noiseReplacer.noiseSeedMultiplier=1

# Add ann offset to the generated noise.
config.processCcd.calibrate.measurement.noiseReplacer.noiseOffset=0.0

# How to choose mean and variance of the Gaussian noise we generate?
# Allowed values:
# 	variance	Mean = 0, variance = the image's variance
# 	meta	Mean = 0, variance = the "BGMEAN" metadata entry
# 	measure	Measure clipped mean and variance from the whole image
# 
config.processCcd.calibrate.measurement.noiseReplacer.noiseSource='measure'

# the name of the flux measurement algorithm used for calibration
config.processCcd.calibrate.measurement.slots.calibFlux='base_CircularApertureFlux_12_0'

# the name of the algorithm used to set the source aperture flux slot
config.processCcd.calibrate.measurement.slots.apFlux='base_CircularApertureFlux_3_0'

# the name of the algorithm used to set the source inst flux slot
config.processCcd.calibrate.measurement.slots.instFlux='base_GaussianFlux'

# the name of the algorithm used to set source moments parameters
config.processCcd.calibrate.measurement.slots.shape='ext_shapeHSM_HsmSourceMoments'

# the name of the centroiding algorithm used to set source x,y
config.processCcd.calibrate.measurement.slots.centroid='base_SdssCentroid'

# the name of the algorithm used to set the source model flux slot
config.processCcd.calibrate.measurement.slots.modelFlux='base_GaussianFlux'

# the name of the algorithm used to set the source psf flux slot
config.processCcd.calibrate.measurement.slots.psfFlux='base_PsfFlux'

# When measuring, replace other detected footprints with noise?
config.processCcd.calibrate.measurement.doReplaceWithNoise=True

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_FPPosition'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_NaiveCentroid'].doMeasure=True

# Value to subtract from the image pixel values
config.processCcd.calibrate.measurement.plugins['base_NaiveCentroid'].background=0.0

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_InputCount'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_PeakLikelihoodFlux'].doMeasure=True

# Name of warping kernel (e.g. "lanczos4") used to compute the peak
config.processCcd.calibrate.measurement.plugins['base_PeakLikelihoodFlux'].warpingKernelName='lanczos4'

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_SdssCentroid'].doMeasure=True

# if the peak's less than this insist on binning at least once
config.processCcd.calibrate.measurement.plugins['base_SdssCentroid'].peakMin=-1.0

# fiddle factor for adjusting the binning
config.processCcd.calibrate.measurement.plugins['base_SdssCentroid'].wfac=1.5

# maximum allowed binning
config.processCcd.calibrate.measurement.plugins['base_SdssCentroid'].binmax=16

# Largest aperture for which to use the slow, accurate, sinc aperture code
config.processCcd.calibrate.measurement.plugins['ext_photometryKron_KronFlux'].maxSincRadius=10.0

# Number of times to iterate when setting the Kron radius
config.processCcd.calibrate.measurement.plugins['ext_photometryKron_KronFlux'].nIterForRadius=1

# Use the Footprint size as part of initial estimate of Kron radius
config.processCcd.calibrate.measurement.plugins['ext_photometryKron_KronFlux'].useFootprintRadius=False

# Minimum Kron radius (if == 0.0 use PSF's Kron radius) if enforceMinimumRadius. Also functions as fallback aperture radius if set.
config.processCcd.calibrate.measurement.plugins['ext_photometryKron_KronFlux'].minimumRadius=0.0

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['ext_photometryKron_KronFlux'].doMeasure=True

# Multiplier of rms size for aperture used to initially estimate the Kron radius
config.processCcd.calibrate.measurement.plugins['ext_photometryKron_KronFlux'].nSigmaForRadius=6.0

# If true check that the Kron radius exceeds some minimum
config.processCcd.calibrate.measurement.plugins['ext_photometryKron_KronFlux'].enforceMinimumRadius=True

# if true, use existing shape and centroid measurements instead of fitting
config.processCcd.calibrate.measurement.plugins['ext_photometryKron_KronFlux'].fixed=False

# Smooth image with N(0, smoothingSigma^2) Gaussian while estimating R_K
config.processCcd.calibrate.measurement.plugins['ext_photometryKron_KronFlux'].smoothingSigma=-1.0

# Number of Kron radii for Kron flux
config.processCcd.calibrate.measurement.plugins['ext_photometryKron_KronFlux'].nRadiusForFlux=2.5

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['ext_shapeHSM_HsmSourceMoments'].doMeasure=True

# Mask planes used to reject bad pixels.
config.processCcd.calibrate.measurement.plugins['ext_shapeHSM_HsmSourceMoments'].badMaskPlanes=['BAD', 'SAT', 'INTRP']

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_CircularApertureFlux'].doMeasure=True

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.processCcd.calibrate.measurement.plugins['base_CircularApertureFlux'].shiftKernel='lanczos5'

# Radius (in pixels) of apertures.
config.processCcd.calibrate.measurement.plugins['base_CircularApertureFlux'].radii=[3.0, 4.5, 6.0, 9.0, 12.0, 17.0, 25.0, 35.0, 50.0, 70.0]

# Maximum radius (in pixels) for which the sinc algorithm should be used instead of the faster naive algorithm.  For elliptical apertures, this is the minor axis radius.
config.processCcd.calibrate.measurement.plugins['base_CircularApertureFlux'].maxSincRadius=12.0

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_Jacobian'].doMeasure=True

# Nominal pixel size (arcsec)
config.processCcd.calibrate.measurement.plugins['base_Jacobian'].pixelScale=0.168

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_GaussianCentroid'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_GaussianFlux'].doMeasure=True

# FIXME! NEVER DOCUMENTED!
config.processCcd.calibrate.measurement.plugins['base_GaussianFlux'].background=0.0

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['ext_shapeHSM_HsmPsfMoments'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['ext_shapeHSM_HsmShapeBj'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.processCcd.calibrate.measurement.plugins['ext_shapeHSM_HsmShapeBj'].badMaskPlanes=['BAD', 'SAT', 'INTRP']

# Field name for number of deblend children
config.processCcd.calibrate.measurement.plugins['ext_shapeHSM_HsmShapeBj'].deblendNChild=''

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_ScaledApertureFlux'].doMeasure=True

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.processCcd.calibrate.measurement.plugins['base_ScaledApertureFlux'].shiftKernel='lanczos5'

# Scaling factor of PSF FWHM for aperture radius.
config.processCcd.calibrate.measurement.plugins['base_ScaledApertureFlux'].scale=3.14

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_PsfFlux'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.processCcd.calibrate.measurement.plugins['base_PsfFlux'].badMaskPlanes=[]

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['ext_shapeHSM_HsmShapeKsb'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.processCcd.calibrate.measurement.plugins['ext_shapeHSM_HsmShapeKsb'].badMaskPlanes=['BAD', 'SAT', 'INTRP']

# Field name for number of deblend children
config.processCcd.calibrate.measurement.plugins['ext_shapeHSM_HsmShapeKsb'].deblendNChild=''

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_Blendedness'].doMeasure=True

# Whether to compute quantities related to the Gaussian-weighted shape
config.processCcd.calibrate.measurement.plugins['base_Blendedness'].doShape=True

# Whether to compute HeavyFootprint dot products (the old deblend.blendedness parameter)
config.processCcd.calibrate.measurement.plugins['base_Blendedness'].doOld=True

# Whether to compute quantities related to the Gaussian-weighted flux
config.processCcd.calibrate.measurement.plugins['base_Blendedness'].doFlux=True

# Radius factor that sets the maximum extent of the weight function (and hence the flux measurements)
config.processCcd.calibrate.measurement.plugins['base_Blendedness'].nSigmaWeightMax=3.0

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['ext_shapeHSM_HsmShapeLinear'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.processCcd.calibrate.measurement.plugins['ext_shapeHSM_HsmShapeLinear'].badMaskPlanes=['BAD', 'SAT', 'INTRP']

# Field name for number of deblend children
config.processCcd.calibrate.measurement.plugins['ext_shapeHSM_HsmShapeLinear'].deblendNChild=''

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_PeakCentroid'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_PixelFlags'].doMeasure=True

# List of mask planes to be searched for which occur in the center of a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.processCcd.calibrate.measurement.plugins['base_PixelFlags'].masksFpCenter=[]

# List of mask planes to be searched for which occur anywhere within a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.processCcd.calibrate.measurement.plugins['base_PixelFlags'].masksFpAnywhere=[]

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_SkyCoord'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_Variance'].doMeasure=True

# Scale factor to apply to shape for aperture
config.processCcd.calibrate.measurement.plugins['base_Variance'].scale=5.0

# Mask planes to ignore
config.processCcd.calibrate.measurement.plugins['base_Variance'].mask=['DETECTED', 'DETECTED_NEGATIVE', 'BAD', 'SAT']

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['ext_shapeHSM_HsmShapeRegauss'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.processCcd.calibrate.measurement.plugins['ext_shapeHSM_HsmShapeRegauss'].badMaskPlanes=['BAD', 'SAT', 'INTRP']

# Field name for number of deblend children
config.processCcd.calibrate.measurement.plugins['ext_shapeHSM_HsmShapeRegauss'].deblendNChild=''

# Maximum centroid shift, limited to 2-10
config.processCcd.calibrate.measurement.plugins['base_SdssShape'].maxShift=0.0

# Convergence tolerance for FWHM
config.processCcd.calibrate.measurement.plugins['base_SdssShape'].tol2=9.999999747378752e-05

# Convergence tolerance for e1,e2
config.processCcd.calibrate.measurement.plugins['base_SdssShape'].tol1=9.999999747378752e-06

# Whether to also compute the shape of the PSF model
config.processCcd.calibrate.measurement.plugins['base_SdssShape'].doMeasurePsf=True

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_SdssShape'].doMeasure=True

# Additional value to add to background
config.processCcd.calibrate.measurement.plugins['base_SdssShape'].background=0.0

# Maximum number of iterations
config.processCcd.calibrate.measurement.plugins['base_SdssShape'].maxIter=100

config.processCcd.calibrate.measurement.plugins.names=['base_CircularApertureFlux', 'base_FPPosition', 'base_NaiveCentroid', 'base_PixelFlags', 'base_SkyCoord', 'base_PsfFlux', 'base_Variance', 'base_Jacobian', 'base_GaussianCentroid', 'ext_shapeHSM_HsmShapeRegauss', 'base_SdssCentroid', 'base_GaussianFlux', 'ext_shapeHSM_HsmPsfMoments', 'base_SdssShape', 'ext_shapeHSM_HsmSourceMoments']
# Maximum allowed shift of WCS, due to matching (pixel)
# 	Valid Range = [-inf,4000)
config.processCcd.calibrate.astrometry.matcher.maxOffsetPix=750

# Type of source flux; typically one of Ap or Psf
config.processCcd.calibrate.astrometry.matcher.sourceFluxType='Psf'

# Number of bright stars to use
# 	Valid Range = [2,inf)
config.processCcd.calibrate.astrometry.matcher.numBrightStars=50

# number of points to define a shape for matching
config.processCcd.calibrate.astrometry.matcher.numPointsForShape=6

# Minimum allowed signal-to-noise ratio for sources used for matching (in the flux specified by sourceFluxType); <=0 for no limit
config.processCcd.calibrate.astrometry.matcher.minSnr=40.0

# Allowed non-perpendicularity of x and y (degree)
# 	Valid Range = [-inf,45.0)
config.processCcd.calibrate.astrometry.matcher.allowedNonperpDeg=0.2

# Maximum separation between reference objects and sources beyond which they will not be considered a match (arcsec)
# 	Valid Range = [0,inf)
config.processCcd.calibrate.astrometry.matcher.maxMatchDistArcSec=2.0

# Minimum number of matched pairs as a fraction of the smaller of the number of reference stars or the number of good sources; the actual minimum is the smaller of this value or minMatchedPairs
# 	Valid Range = [0,1)
config.processCcd.calibrate.astrometry.matcher.minFracMatchedPairs=0.3

# maximum determinant of linear transformation matrix for a usable solution
config.processCcd.calibrate.astrometry.matcher.maxDeterminant=0.02

# Rotation angle allowed between sources and position reference objects (degrees)
# 	Valid Range = [-inf,6.0)
config.processCcd.calibrate.astrometry.matcher.maxRotationDeg=1.145916

# Minimum number of matched pairs; see also minFracMatchedPairs
# 	Valid Range = [2,inf)
config.processCcd.calibrate.astrometry.matcher.minMatchedPairs=30

# the match distance below which further iteration is pointless (arcsec); ignored if not fitting a WCS
# 	Valid Range = [0,inf)
config.processCcd.calibrate.astrometry.minMatchDistanceArcSec=0.001

# If True then load reference objects and match sources but do not fit a WCS;  this simply controls whether 'run' calls 'solve' or 'loadAndMatch'
config.processCcd.calibrate.astrometry.forceKnownWcs=False

# maximum number of iterations of match sources and fit WCSignored if not fitting a WCS
# 	Valid Range = [1,inf)
config.processCcd.calibrate.astrometry.maxIter=3

# the maximum match distance is set to  mean_match_distance + matchDistanceSigma*std_dev_match_distance; ignored if not fitting a WCS
# 	Valid Range = [0,inf)
config.processCcd.calibrate.astrometry.matchDistanceSigma=2.0

# Number of standard deviations for clipping level
# 	Valid Range = [0.0,inf)
config.processCcd.calibrate.astrometry.wcsFitter.rejSigma=3.0

# maximum median scatter of a WCS fit beyond which the fit fails (arcsec); be generous, as this is only intended to catch catastrophic failures
# 	Valid Range = [0,inf)
config.processCcd.calibrate.astrometry.wcsFitter.maxScatterArcsec=10.0

# number of rejection iterations
# 	Valid Range = [0,inf)
config.processCcd.calibrate.astrometry.wcsFitter.numRejIter=3

# number of iterations of fitter (which fits X and Y separately, and so benefits from a few iterations
# 	Valid Range = [1,inf)
config.processCcd.calibrate.astrometry.wcsFitter.numIter=3

# order of SIP polynomial
# 	Valid Range = [0,inf)
config.processCcd.calibrate.astrometry.wcsFitter.order=3

# flux measurement algorithms in getApCorrNameSet() to ignore; if a name is listed that does not appear in getApCorrNameSet() then a warning is logged
config.processCcd.calibrate.applyApCorr.ignoreList=[]

# set the general failure flag for a flux when it cannot be aperture-corrected?
config.processCcd.calibrate.applyApCorr.doFlagApCorrFailures=True

# correction factor for psfFlux error
config.processCcd.calibrate.afterburners.plugins['base_ClassificationExtendedness'].psfErrFactor=0.0

# correction factor for modelFlux error
config.processCcd.calibrate.afterburners.plugins['base_ClassificationExtendedness'].modelErrFactor=0.0

# critical ratio of model to psf flux
config.processCcd.calibrate.afterburners.plugins['base_ClassificationExtendedness'].fluxRatio=0.95

config.processCcd.calibrate.afterburners.plugins.names=['base_ClassificationExtendedness']
# Estimate the background again after final source detection?
config.processCcd.calibrate.detection.reEstimateBackground=True

# Grow detections by nSigmaToGrow * sigma; if 0 then do not grow
config.processCcd.calibrate.detection.nSigmaToGrow=2.4

# detected sources with fewer than the specified number of pixels will be ignored
# 	Valid Range = [0,inf)
config.processCcd.calibrate.detection.minPixels=1

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.processCcd.calibrate.detection.tempLocalBackground.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.processCcd.calibrate.detection.tempLocalBackground.approxOrderY=-1

# Names of mask planes to ignore while estimating the background
config.processCcd.calibrate.detection.tempLocalBackground.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# behaviour if there are too few points in grid for requested interpolation style
# Allowed values:
# 	None	Field is optional
# 	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
# 	THROW_EXCEPTION	throw an exception if there are too few points
# 	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
# 
config.processCcd.calibrate.detection.tempLocalBackground.undersampleStyle='REDUCE_INTERP_ORDER'

# how to interpolate the background values. This maps to an enum; see afw::math::Background
# Allowed values:
# 	NONE	No background estimation is to be attempted
# 	CONSTANT	Use a single constant value
# 	LINEAR	Use linear interpolation
# 	None	Field is optional
# 	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
# 	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
# 
config.processCcd.calibrate.detection.tempLocalBackground.algorithm='AKIMA_SPLINE'

# Ignore NaNs when estimating the background
config.processCcd.calibrate.detection.tempLocalBackground.isNanSafe=False

# how large a region of the sky should be used for each background point
# 	Valid Range = [1,inf)
config.processCcd.calibrate.detection.tempLocalBackground.binSize=64

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.processCcd.calibrate.detection.tempLocalBackground.weighting=True

# type of statistic to use for grid points
# Allowed values:
# 	None	Field is optional
# 	MEDIAN	median
# 	MEANCLIP	clipped mean
# 	MEAN	unclipped mean
# 
config.processCcd.calibrate.detection.tempLocalBackground.statisticsProperty='MEANCLIP'

# Use Approximate (Chebyshev) to model background.
config.processCcd.calibrate.detection.tempLocalBackground.useApprox=False

# Include threshold relative to thresholdValue
# 	Valid Range = [0.0,inf)
config.processCcd.calibrate.detection.includeThresholdMultiplier=1.0

# Pixels should be grown as isotropically as possible (slower)
config.processCcd.calibrate.detection.isotropicGrow=True

# Fiddle factor to add to the background; debugging only
config.processCcd.calibrate.detection.adjustBackground=0.0

# Do temporary interpolated background subtraction before footprint detection?
config.processCcd.calibrate.detection.doTempLocalBackground=False

# specifies whether to detect positive, or negative sources, or both
# Allowed values:
# 	positive	detect only positive sources
# 	negative	detect only negative sources
# 	both	detect both positive and negative sources
# 
config.processCcd.calibrate.detection.thresholdPolarity='positive'

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.processCcd.calibrate.detection.background.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.processCcd.calibrate.detection.background.approxOrderY=-1

# Names of mask planes to ignore while estimating the background
config.processCcd.calibrate.detection.background.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# behaviour if there are too few points in grid for requested interpolation style
# Allowed values:
# 	None	Field is optional
# 	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
# 	THROW_EXCEPTION	throw an exception if there are too few points
# 	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
# 
config.processCcd.calibrate.detection.background.undersampleStyle='REDUCE_INTERP_ORDER'

# how to interpolate the background values. This maps to an enum; see afw::math::Background
# Allowed values:
# 	NONE	No background estimation is to be attempted
# 	CONSTANT	Use a single constant value
# 	LINEAR	Use linear interpolation
# 	None	Field is optional
# 	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
# 	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
# 
config.processCcd.calibrate.detection.background.algorithm='NATURAL_SPLINE'

# Ignore NaNs when estimating the background
config.processCcd.calibrate.detection.background.isNanSafe=False

# how large a region of the sky should be used for each background point
# 	Valid Range = [1,inf)
config.processCcd.calibrate.detection.background.binSize=128

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.processCcd.calibrate.detection.background.weighting=True

# type of statistic to use for grid points
# Allowed values:
# 	None	Field is optional
# 	MEDIAN	median
# 	MEANCLIP	clipped mean
# 	MEAN	unclipped mean
# 
config.processCcd.calibrate.detection.background.statisticsProperty='MEANCLIP'

# Use Approximate (Chebyshev) to model background.
config.processCcd.calibrate.detection.background.useApprox=True

# Grow detections to set the image mask bits, but return the original (not-grown) footprints
config.processCcd.calibrate.detection.returnOriginalFootprints=False

# specifies the desired flavor of Threshold
# Allowed values:
# 	pixel_stdev	threshold applied to per-pixel std deviation
# 	variance	threshold applied to image variance
# 	value	threshold applied to image value
# 	stdev	threshold applied to image std deviation
# 
config.processCcd.calibrate.detection.thresholdType='stdev'

# Threshold for footprints
# 	Valid Range = [0.0,inf)
config.processCcd.calibrate.detection.thresholdValue=5.0

# Write reference matches (ignored if doWrite false)?
config.processCcd.calibrate.doWriteMatches=True

# Default reference catalog filter to use if filter not specified in exposure; if blank then filter must be specified in exposure
config.processCcd.calibrate.refObjLoader.defaultFilter=''

# Padding to add to 4 all edges of the bounding box (pixels)
# 	Valid Range = [0,inf)
config.processCcd.calibrate.refObjLoader.pixelMargin=50

# Mapping of camera filter name: reference catalog filter name; each reference filter must exist
config.processCcd.calibrate.refObjLoader.filterMap={'N1010': 'z', 'N816': 'i', 'N387': 'g', 'i2': 'i', 'N921': 'z', 'N515': 'g', 'y': 'z'}

# Strictness of Astropy unit compatibility check, can be 'raise', 'warn' or 'silent'
config.processCcd.calibrate.checkUnitsParseStrict='raise'

# Run deblender input exposure
config.processCcd.calibrate.doDeblend=True

# Raise an exception if photoCal fails? Ignored if doPhotoCal false.
config.processCcd.calibrate.requirePhotoCal=True

# Save calibration results?
config.processCcd.calibrate.doWrite=True

# Use the extendedness parameter to select objects to use in photometric calibration?
# This applies only to the sources detected on the exposure, not the reference catalog
config.processCcd.calibrate.photoCal.doSelectUnresolved=True

# number of iterations
config.processCcd.calibrate.photoCal.nIter=20

# Name of the source flux field to use.  The associated flag field
# ('<name>_flags') will be implicitly included in badFlags.
config.processCcd.calibrate.photoCal.fluxField='slot_CalibFlux_flux'

config.processCcd.calibrate.photoCal.colorterms.data={}
config.processCcd.calibrate.photoCal.colorterms.data['hsc*']=lsst.pipe.tasks.colorterms.ColortermDict()
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data={}
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['i']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['i'].c2=0.0

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['i'].c1=0.0

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['i'].c0=0.0

# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['i'].primary='i'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['i'].secondary='i'

config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['y']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['y'].c2=0.0

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['y'].c1=0.0

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['y'].c0=0.0

# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['y'].primary='y'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['y'].secondary='y'

config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['r']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['r'].c2=0.0

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['r'].c1=0.0

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['r'].c0=0.0

# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['r'].primary='r'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['r'].secondary='r'

config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['z']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['z'].c2=0.0

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['z'].c1=0.0

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['z'].c0=0.0

# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['z'].primary='z'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['z'].secondary='z'

config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['g']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['g'].c2=0.0

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['g'].c1=0.0

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['g'].c0=0.0

# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['g'].primary='g'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['hsc*'].data['g'].secondary='g'

config.processCcd.calibrate.photoCal.colorterms.data['sdss*']=lsst.pipe.tasks.colorterms.ColortermDict()
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data={}
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['g']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['g'].c2=-0.00726883

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['g'].c1=-0.08366937

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['g'].c0=-0.00816446

# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['g'].primary='g'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['g'].secondary='r'

config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['N816']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['N816'].c2=-0.05474862

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['N816'].c1=-0.63558358

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['N816'].c0=0.00927133

# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['N816'].primary='i'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['N816'].secondary='z'

config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['i']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['i'].c2=-0.01374245

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['i'].c1=-0.16922042

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['i'].c0=0.00130204

# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['i'].primary='i'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['i'].secondary='z'

config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['i2']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['i2'].c2=-0.01067212

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['i2'].c1=-0.20739606

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['i2'].c0=0.00124676

# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['i2'].primary='i'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['i2'].secondary='z'

config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['r']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['r'].c2=-0.03068248

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['r'].c1=0.01284177

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['r'].c0=0.0023181

# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['r'].primary='r'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['r'].secondary='i'

config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['N921']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['N921'].c2=-0.05451118

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['N921'].c1=0.0986353

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['N921'].c0=0.00752972

# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['N921'].primary='z'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['N921'].secondary='i'

config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['y']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['y'].c2=0.00574408

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['y'].c1=0.35652971

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['y'].c0=0.01739708

# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['y'].primary='z'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['y'].secondary='i'

config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['z']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['z'].c2=0.01479369

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['z'].c1=0.01353969

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['z'].c0=-0.0068062

# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['z'].primary='z'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['sdss*'].data['z'].secondary='i'

config.processCcd.calibrate.photoCal.colorterms.data['ps1*']=lsst.pipe.tasks.colorterms.ColortermDict()
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data={}
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['g']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['g'].c2=-0.0151057

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['g'].c1=0.06508481

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['g'].c0=0.00730066

# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['g'].primary='g'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['g'].secondary='r'

config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N816']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N816'].c2=-0.10781564

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N816'].c1=-0.68757034

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N816'].c0=0.01191062

# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N816'].primary='i'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N816'].secondary='z'

config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['i']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['i'].c2=-0.03034094

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['i'].c1=-0.13944659

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['i'].c0=0.00166891

# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['i'].primary='i'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['i'].secondary='z'

config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['i2']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['i2'].c2=-0.02675511

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['i2'].c1=-0.18483562

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['i2'].c0=0.00180361

# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['i2'].primary='i'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['i2'].secondary='z'

config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['r']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['r'].c2=-0.01877566

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['r'].c1=0.02093734

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['r'].c0=0.00279757

# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['r'].primary='r'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['r'].secondary='i'

config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N921']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N921'].c2=-0.25059679

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N921'].c1=-0.59278367

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N921'].c0=0.00142051

# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N921'].primary='z'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['N921'].secondary='y'

config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['y']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['y'].c2=0.02880125

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['y'].c1=0.14747401

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['y'].c0=-0.00156858

# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['y'].primary='y'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['y'].secondary='z'

config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['z']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['z'].c2=-0.00316369

# First-order parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['z'].c1=-0.28840221

# Constant parameter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['z'].c0=-0.00907517

# name of primary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['z'].primary='z'

# name of secondary filter
config.processCcd.calibrate.photoCal.colorterms.data['ps1*'].data['z'].secondary='y'

# Name of photometric reference catalog; used to select a color term dict in colorterms. see also applyColorTerms
config.processCcd.calibrate.photoCal.photoCatName='sdss-dr9-fink-v5b'

# Don't use objects fainter than this magnitude
config.processCcd.calibrate.photoCal.magLimit=22.0

# maximum sigma to use when clipping
config.processCcd.calibrate.photoCal.sigmaMax=0.25

# List of source flag fields that must be set for a source to be used.
config.processCcd.calibrate.photoCal.goodFlags=[]

# Additional magnitude uncertainty to be added in quadrature with measurement errors.
# 	Valid Range = [0.0,inf)
config.processCcd.calibrate.photoCal.magErrFloor=0.0

# clip at nSigma
config.processCcd.calibrate.photoCal.nSigma=3.0

# List of source flag fields that will cause a source to be rejected when they are set.
config.processCcd.calibrate.photoCal.badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolated', 'base_PixelFlags_flag_saturated']

# use median instead of mean to compute zeropoint
config.processCcd.calibrate.photoCal.useMedian=True

# Write a field name astrom_usedByPhotoCal to the schema
config.processCcd.calibrate.photoCal.doWriteOutput=True

# Apply photometric color terms to reference stars? One of:
# None: apply if colorterms and photoCatName are not None;
#       fail if color term data is not available for the specified ref catalog and filter.
# True: always apply colorterms; fail if color term data is not available for the
#       specified reference catalog and filter.
# False: do not apply.
config.processCcd.calibrate.photoCal.applyColorTerms=True

# Fields to copy from the icSource catalog to the output catalog for matching sources Any missing fields will trigger a RuntimeError exception. Ignored if icSourceCat is not provided.
config.processCcd.calibrate.icSourceFieldsToCopy=['calib_psfCandidate', 'calib_psfUsed', 'calib_psfReserved']

# Match radius for matching icSourceCat objects to sourceCat objects (pixels)
config.processCcd.calibrate.matchRadiusPix=3.0

# What to do when a peak to be deblended is close to the edge of the image
# Allowed values:
# 	ramp	Ramp down flux at the image edge by the PSF
# 	noclip	Ignore the edge when building the symmetric template.
# 	clip	Clip the template at the edge AND the mirror of the edge.
# 	None	Field is optional
# 
config.processCcd.calibrate.deblend.edgeHandling='ramp'

# Assign stray flux to deblend children.  Implies findStrayFlux.
config.processCcd.calibrate.deblend.assignStrayFlux=True

# Mask planes with the corresponding limit on the fraction of masked pixels. Sources violating this limit will not be deblended.
config.processCcd.calibrate.deblend.maskLimits={'NO_DATA': 0.25}

# How to split flux among peaks
# Allowed values:
# 	trim	Shrink the parent footprint to pixels that are not assigned to children
# 	r-to-peak	~ 1/(1+R^2) to the peak
# 	r-to-footprint	~ 1/(1+R^2) to the closest pixel in the footprint.  CAUTION: this can be computationally expensive on large footprints!
# 	nearest-footprint	Assign 100% to the nearest footprint (using L-1 norm aka Manhattan distance)
# 	None	Field is optional
# 
config.processCcd.calibrate.deblend.strayFluxRule='trim'

# If True, catch exceptions thrown by the deblender, log them, and set a flag on the parent, instead of letting them propagate up
config.processCcd.calibrate.deblend.catchFailures=False

# When the deblender should attribute stray flux to point sources
# Allowed values:
# 	always	Always
# 	never	Never; stray flux will not be attributed to any deblended child if the deblender thinks all peaks look like point sources
# 	necessary	When there is not an extended object in the footprint
# 	None	Field is optional
# 
config.processCcd.calibrate.deblend.strayFluxToPointSources='necessary'

# Chi-squared per DOF cut for deciding a source is a PSF during deblending (shifted PSF model #2)
config.processCcd.calibrate.deblend.psfChisq2b=1.5

# Mask planes to ignore when performing statistics
config.processCcd.calibrate.deblend.maskPlanes=['SAT', 'INTRP', 'NO_DATA']

# Maximum area for footprints before they are ignored as large; non-positive means no threshold applied
config.processCcd.calibrate.deblend.maxFootprintArea=10000

# Minimum axis ratio for footprints before they are ignored as large; non-positive means no threshold applied
config.processCcd.calibrate.deblend.minFootprintAxisRatio=0.0

# Maximum linear dimension for footprints before they are ignored as large; non-positive means no threshold applied
config.processCcd.calibrate.deblend.maxFootprintSize=0

# Chi-squared per DOF cut for deciding a source is a PSF during deblending (un-shifted PSF model)
config.processCcd.calibrate.deblend.psfChisq1=1.5

# Find stray flux---flux not claimed by any child in the deblender.
config.processCcd.calibrate.deblend.findStrayFlux=True

# Footprints smaller in width or height than this value will be ignored; minimum of 2 due to PSF gradient calculation.
# 	Valid Range = [2,inf)
config.processCcd.calibrate.deblend.tinyFootprintSize=2

# Chi-squared per DOF cut for deciding a source is PSF during deblending (shifted PSF model)
config.processCcd.calibrate.deblend.psfChisq2=1.5

# Guarantee that all peaks produce a child source.
config.processCcd.calibrate.deblend.propagateAllPeaks=False

# Mask name for footprints not deblended, or None
config.processCcd.calibrate.deblend.notDeblendedMask='NOT_DEBLENDED'

# Only deblend the brightest maxNumberOfPeaks peaks in the parent (<= 0: unlimited)
config.processCcd.calibrate.deblend.maxNumberOfPeaks=0

# When splitting stray flux, clip fractions below this value to zero.
config.processCcd.calibrate.deblend.clipStrayFluxFraction=0.001

import lsst.obs.subaru.isr
config.processCcd.isr.retarget(target=lsst.obs.subaru.isr.SubaruIsrTask, ConfigClass=lsst.obs.subaru.isr.SubaruIsrConfig)
# Number of stdev below the background to set thumbnail minimum
config.processCcd.isr.thumbnailStdev=3.0

# How to estimate the average value for BAD regions.
# Allowed values:
# 	MEDIAN	Correct using the median of the good data
# 	MEANCLIP	Correct using the (clipped) mean of good data
# 	None	Field is optional
# 
config.processCcd.isr.badStatistic='MEANCLIP'

# Order of polynomial or to fit if overscan fit type is a polynomial, or number of spline knots if overscan fit type is a spline.
config.processCcd.isr.overscanOrder=30

# Widen bleed trails based on their width?
config.processCcd.isr.doWidenSaturationTrails=True

import lsst.obs.subaru.crosstalk
config.processCcd.isr.crosstalk.retarget(target=lsst.obs.subaru.crosstalk.CrosstalkTask, ConfigClass=lsst.obs.subaru.crosstalk.CrosstalkConfig)
# Shape of coeffs array
config.processCcd.isr.crosstalk.coeffs.shape=[4, 4]

# Crosstalk coefficients
config.processCcd.isr.crosstalk.coeffs.values=[0.0, -0.000125, -0.000149, -0.000156, -0.000124, 0.0, -0.000132, -0.000157, -0.000171, -0.000134, 0.0, -0.000153, -0.000157, -0.000151, -0.000137, 0.0]

# Name for crosstalk mask plane
config.processCcd.isr.crosstalk.crosstalkMaskPlane='CROSSTALK'

# Set crosstalk mask plane for pixels over this value
config.processCcd.isr.crosstalk.minPixelToMask=45000.0

# Apply dark frame correction?
config.processCcd.isr.doDark=True

# update exposure metadata in the assembled ccd to reflect the effective gain of the assembled chip
config.processCcd.isr.setGainAssembledCcd=True

# Fallback default filter name for calibrations
config.processCcd.isr.fallbackFilterName=None

# Correct for crosstalk
config.processCcd.isr.doCrosstalk=True

# FWHM of PSF used when interpolating over bad columns (arcsec)
config.processCcd.isr.fwhmForBadColumnInterpolation=1.0

# The read noise to use if no Detector is present in the Exposure
config.processCcd.isr.readNoise=0.0

# FWHM of PSF (arcsec)
config.processCcd.isr.fwhm=1.0

# Name of mask plane to use for suspect pixels
config.processCcd.isr.suspectMaskName='SUSPECT'

# Maximum number of iterations for the brighter fatter correction
config.processCcd.isr.brighterFatterMaxIter=10

# If flatScalingType is 'USER' then scale flat by this amount; ignored otherwise
config.processCcd.isr.flatUserScale=1.0

# Rejection threshold (sigma) for collapsing overscan before fit
config.processCcd.isr.overscanRej=3.0

# Dataset type for input data; users will typically leave this alone, but camera-specific ISR tasks will override it
config.processCcd.isr.datasetType='raw'

# Assemble amp-level exposures into a ccd-level exposure?
config.processCcd.isr.doAssembleCcd=True

# Name of mask plane to use in saturation detection and interpolation
config.processCcd.isr.saturatedMaskName='SAT'

# Should the gain be applied when applying the brighter fatter correction?
config.processCcd.isr.brighterFatterApplyGain=True

# Kernel file used for the brighter fatter correction
config.processCcd.isr.brighterFatterKernelFile=''

# Threshold used to stop iterating the brighter fatter correction.  It is the  absolute value of the difference between the current corrected image and the one from the previous iteration summed over all the pixels.
config.processCcd.isr.brighterFatterThreshold=1000.0

# Correct for nonlinearity of the detector's response?
config.processCcd.isr.doLinearize=True

# Should we set the level of all BAD patches of the chip to the chip's average value?
config.processCcd.isr.doSetBadRegions=True

# Apply the brighter fatter correction
config.processCcd.isr.doBrighterFatter=True

# Apply bias frame correction?
config.processCcd.isr.doBias=True

# Mask suspect pixels?
config.processCcd.isr.doSuspect=True

# Apply flat field correction?
config.processCcd.isr.doFlat=True

# Remove any PC cards in the header
config.processCcd.isr.removePcCards=True

# Apply fringe correction?
config.processCcd.isr.doFringe=True

# trim out non-data regions?
config.processCcd.isr.assembleCcd.doTrim=True

# FITS headers to remove (in addition to DATASEC, BIASSEC, TRIMSEC and perhaps GAIN)
config.processCcd.isr.assembleCcd.keysToRemove=[]

# renormalize to a gain of 1? (ignored if setGain false). Setting to True gives 1 ADU per electron. Setting to True is not recommended for mosaic cameras because it breaks normalization across the focal plane. However, if the CCDs are sufficiently flat then the resulting error may be acceptable.
config.processCcd.isr.assembleCcd.doRenorm=False

# set gain?
config.processCcd.isr.assembleCcd.setGain=True

# The saturation level to use if no Detector is present in the Exposure (ignored if NaN)
config.processCcd.isr.saturation=float('nan')

# Calculate variance?
config.processCcd.isr.doVariance=True

# Assemble amp-level calibration exposures into ccd-level exposure?
config.processCcd.isr.doAssembleIsrExposures=False

# Tweak flats to match observed amplifier ratios?
config.processCcd.isr.doTweakFlat=False

# fields to remove from the metadata of the assembled ccd.
config.processCcd.isr.keysToRemoveFromAssembledCcd=[]

# Center of vignetting pattern, in x (focal plane coords)
config.processCcd.isr.vignette.xCenter=-100.0

# Radius of vignetting pattern, in focal plane coords
config.processCcd.isr.vignette.radius=17500.0

# Center of vignetting pattern, in y (focal plane coords)
config.processCcd.isr.vignette.yCenter=100.0

# The approximate flux of a zero-magnitude object in a one-second exposure, per filter
config.processCcd.isr.fluxMag0T1={'g': 398107170553.49854, 'N816': 15848931924.611174, 'i': 275422870333.81744, 'r': 398107170553.49854, 'N921': 19054607179.632523, 'N515': 20892961308.54041, 'y': 91201083935.59116, 'z': 120226443461.74132}

# Do overscan subtraction?
config.processCcd.isr.doOverscan=True

# The gain to use if no Detector is present in the Exposure (ignored if NaN)
config.processCcd.isr.gain=float('nan')

# Do fringe subtraction after flat-fielding?
config.processCcd.isr.fringeAfterFlat=True

# Border around saturated pixels for thumbnail
config.processCcd.isr.thumbnailSatBorder=2

# Mask saturated pixels?
config.processCcd.isr.doSaturation=True

# Trim guider shadow
config.processCcd.isr.doGuider=False

# Softening parameter for thumbnail mapping
config.processCcd.isr.thumbnailQ=20.0

# Default value for fluxMag0T1 (for an unrecognised filter)
config.processCcd.isr.defaultFluxMag0T1=158489319246.11172

# Correct the amplifiers for their gains
# 
# N.b. this is intended to be used *instead* of doFlat; it's useful if you're measuring system throughput
# 
config.processCcd.isr.doApplyGains=False

# Number of points to define the Vignette polygon
config.processCcd.isr.numPolygonPoints=100

# Persist Polygon used to define vignetted region?
config.processCcd.isr.doWriteVignettePolygon=True

# Offset to the random number generator seed (full seed includes exposure ID)
config.processCcd.isr.fringe.stats.rngSeedOffset=0

# Ignore pixels with these masks
config.processCcd.isr.fringe.stats.badMaskPlanes=['SAT', 'NO_DATA']

# Statistic to use
config.processCcd.isr.fringe.stats.stat=32

# Number of fitting iterations
config.processCcd.isr.fringe.stats.iterations=3

# Sigma clip threshold
config.processCcd.isr.fringe.stats.clip=3.0

# Only fringe-subtract these filters
config.processCcd.isr.fringe.filters=['y', 'N921']

# Sigma clip threshold
config.processCcd.isr.fringe.clip=3.0

# Half-size of large (background) measurements (pixels)
config.processCcd.isr.fringe.large=30

# Number of fringe measurements
config.processCcd.isr.fringe.num=30000

# Number of fitting iterations
config.processCcd.isr.fringe.iterations=20

# Half-size of small (fringe) measurements (pixels)
config.processCcd.isr.fringe.small=3

# Remove fringe pedestal?
config.processCcd.isr.fringe.pedestal=False

# The method for scaling the flat on the fly.
# Allowed values:
# 	None	Field is optional
# 	MEDIAN	Scale by the inverse of the median
# 	USER	Scale by flatUserScale
# 	MEAN	Scale by the inverse of the mean
# 
config.processCcd.isr.flatScalingType='USER'

# Write OverScan-Subtracted thumbnail?
config.processCcd.isr.qa.doThumbnailOss=True

# Mesh size in X (pix) to calculate count statistics
config.processCcd.isr.qa.flatness.meshX=256

# How many times do we iterate clipping outliers in calculate count statistics?
config.processCcd.isr.qa.flatness.nIter=3

# Do we clip outliers in calculate count statistics?
config.processCcd.isr.qa.flatness.doClip=True

# How many sigma is used to clip outliers in calculate count statistics?
config.processCcd.isr.qa.flatness.clipSigma=3.0

# Mesh size in Y (pix) to calculate count statistics
config.processCcd.isr.qa.flatness.meshY=256

# Write flattened thumbnail?
config.processCcd.isr.qa.doThumbnailFlattened=True

# Write OverScan-Subtracted image?
config.processCcd.isr.qa.doWriteOss=False

# Write flattened image?
config.processCcd.isr.qa.doWriteFlattened=False

# Range for thumbnail mapping
config.processCcd.isr.thumbnailRange=5.0

# Persist postISRCCD?
config.processCcd.isr.doWrite=False

# Mask defect pixels?
config.processCcd.isr.doDefect=True

# Normalize all the amplifiers in each CCD to have the same gain
# 
# This does not measure the gains, it simply forces the median of each amplifier to be equal
# after applying the nominal gain
# 
config.processCcd.isr.normalizeGains=False

# Maximum deviation from the median for overscan
config.processCcd.isr.overscanMaxDev=1000.0

# The method for fitting the overscan bias level.
# Allowed values:
# 	LEG	Fit Legendre polynomial to the longest axis of the overscan region
# 	CUBIC_SPLINE	Fit cubic spline to the longest axis of the overscan region
# 	MEDIAN	Correct using the median of the overscan region
# 	None	Field is optional
# 	POLY	Fit ordinary polynomial to the longest axis of the overscan region
# 	CHEB	Fit Chebyshev polynomial to the longest axis of the overscan region
# 	AKIMA_SPLINE	Fit Akima spline to the longest axis of the overscan region
# 	NATURAL_SPLINE	Fit natural spline to the longest axis of the overscan region
# 	MEAN	Correct using the mean of the overscan region
# 
config.processCcd.isr.overscanFitType='AKIMA_SPLINE'

# Number of pixels by which to grow the saturation footprints
config.processCcd.isr.growSaturationFootprintSize=1

# Binning factor for thumbnail
config.processCcd.isr.thumbnailBinning=4

# DataId key corresponding to a single sensor
config.ccdKey='ccd'

