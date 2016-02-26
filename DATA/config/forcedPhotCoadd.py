import lsst.meas.base.forcedPhotCoadd
assert type(config)==lsst.meas.base.forcedPhotCoadd.ForcedPhotCoaddConfig, 'config is of type %s.%s instead of lsst.meas.base.forcedPhotCoadd.ForcedPhotCoaddConfig' % (type(config).__module__, type(config).__name__)
import lsst.meas.extensions
import lsst.meas.extensions
# Apply aperture corrections? Silently ignored if endOrder <= lsst.meas.base.APCORR_ORDER when calling run
# Allowed values:
# 	None	Field is optional
# 	yes	apply aperture corrections; fail if data not available
# 	noButWarn	do not apply aperture corrections, but warn if data available (since aperture corrections could have been applied)
# 	yesOrWarn	apply aperture corrections if data available, else warn
# 	no	do not apply aperture corrections
# 
config.measurement.doApplyApCorr='yes'

# Mapping of reference columns to source columns
config.measurement.copyColumns={'id': 'objectId', 'parent': 'parentObjectId'}

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
config.measurement.plugins['base_TransformedShape'].doMeasure=True

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

# whether to run this plugin in single-object mode
config.measurement.plugins['base_TransformedCentroid'].doMeasure=True

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

config.measurement.plugins.names=['base_CircularApertureFlux', 'base_GaussianFlux', 'base_PsfFlux', 'base_TransformedShape', 'base_TransformedCentroid']
# the name of the flux measurement algorithm used for calibration
config.measurement.slots.calibFlux=None

# the name of the algorithm used to set the source aperture flux slot
config.measurement.slots.apFlux=None

# the name of the algorithm used to set the source inst flux slot
config.measurement.slots.instFlux=None

# the name of the algorithm used to set source moments parameters
config.measurement.slots.shape='base_TransformedShape'

# the name of the centroiding algorithm used to set source x,y
config.measurement.slots.centroid='base_TransformedCentroid'

# the name of the algorithm used to set the source model flux slot
config.measurement.slots.modelFlux=None

# the name of the algorithm used to set the source psf flux slot
config.measurement.slots.psfFlux=None

# Dataset (without coadd prefix) that should be used to obtain (Heavy)Footprints for sources.  Must have IDs that match those of the reference catalog.  If None, Footprints will be generated by transforming the reference Footprints.
config.footprintDatasetName='meas'

# coadd name: typically one of deep or goodSeeing
config.coaddName='deep'

# Bandpass for reference sources; None indicates chi-squared detections.
config.references.filter=None

# Only include reference sources for each patch that lie within the patch's inner bbox
config.references.removePatchOverlaps=False

# Coadd name: typically one of deep or goodSeeing.
config.references.coaddName='deep'

# Mapping of reference columns to source columns
config.copyColumns={'id': 'objectId', 'parent': 'parentObjectId', 'deblend_nChild': 'deblend_nChild'}

