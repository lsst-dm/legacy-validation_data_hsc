import lsst.pipe.tasks.multiBand
assert type(config)==lsst.pipe.tasks.multiBand.MeasureMergedCoaddSourcesConfig, 'config is of type %s.%s instead of lsst.pipe.tasks.multiBand.MeasureMergedCoaddSourcesConfig' % (type(config).__module__, type(config).__name__)
import lsst.meas.extensions
# Source catalog flags to propagate, with the threshold of relative occurrence.
config.propagateFlags.flags={'calib_psfCandidate': 0.2, 'calib_psfUsed': 0.2}

# Source matching radius (arcsec)
config.propagateFlags.matchRadius=0.2

# Whether to match sources to CCD catalogs to propagate flags (to e.g. identify PSF stars)
config.doPropagateFlags=True

# Rejection iterations for Wcs fitting
# 	Valid Range = [0,inf)
config.astrometry.rejectIter=3

# Assume that the input image's WCS is correct, without comparing it to any external reality. (In contrast to using Astrometry.net).  NOTE, if you set this, you probably also want to un-set 'solver.calculateSip'; otherwise we'll still try to find a TAN-SIP WCS starting  from the existing WCS
config.astrometry.forceKnownWcs=False

# Rejection threshold for Wcs fitting
# 	Valid Range = (0.0,inf)
config.astrometry.rejectThresh=3.0

# Matching radius (arcsec) for matching sources to reference objects
# 	Valid Range = [0.0,inf)
config.astrometry.solver.catalogMatchDist=1.0

# Sigma-clipping parameter in sip/cleanBadPoints.py
# 	Valid Range = [0.0,inf)
config.astrometry.solver.cleaningParameter=3.0

# Range of pixel scales, around the value in the WCS header, to search. If the value of this field is X and the nominal scale is S, the range searched will be  S/X to S*X
# 	Valid Range = [1.001,inf)
config.astrometry.solver.pixelScaleUncertainty=1.1

# Default reference catalog filter to use if filter not specified in exposure; if blank then filter must be specified in exposure
config.astrometry.solver.defaultFilter=''

# Padding to add to 4 all edges of the bounding box (pixels)
# 	Valid Range = [0,inf)
config.astrometry.solver.pixelMargin=50

# Mapping of camera filter name: reference catalog filter name; each reference filter must exist
config.astrometry.solver.filterMap={}

# Use the parity (flip / handedness) of the image from the input exposure's WCS headers?
config.astrometry.solver.useWcsParity=True

# Maximum CPU time to spend solving, in seconds
# 	Valid Range = [0.0,inf)
config.astrometry.solver.maxCpuTime=0.0

# Polynomial order of SIP distortion terms
# 	Valid Range = [2,inf)
config.astrometry.solver.sipOrder=4

# Use the RA,Dec center information from the input exposure's WCS headers?
config.astrometry.solver.useWcsRaDecCenter=True

# Matching threshold for Astrometry.net solver (log-odds)
# 	Valid Range = [13.815510558,inf)
config.astrometry.solver.matchThreshold=27.631021115928547

# Use the pixel scale from the input exposure's WCS headers?
config.astrometry.solver.useWcsPixelScale=True

# Maximum number of stars to use in Astrometry.net solving
# 	Valid Range = [10,inf)
config.astrometry.solver.maxStars=50

# List of flags which cause a source to be rejected as bad
config.astrometry.solver.badFlags=['slot_Centroid_flag', 'base_PixelFlags_flag_edge', 'base_PixelFlags_flag_saturated', 'base_PixelFlags_flag_crCenter']

# maximum number of iterations of match sources and fit WCSignored if not fitting a WCS
# 	Valid Range = [1,inf)
config.astrometry.solver.maxIter=5

# Compute polynomial SIP distortion terms?
config.astrometry.solver.calculateSip=True

# When useWcsRaDecCenter=True, this is the radius, in degrees, around the RA,Dec center specified in the input exposure's WCS to search for a solution.
# 	Valid Range = [0.0,inf)
config.astrometry.solver.raDecSearchRadius=1.0

# Retrieve all available fluxes (and errors) from catalog?
config.astrometry.solver.allFluxes=True

# The match and fit loop stops when maxMatchDist minimized:  maxMatchDist = meanMatchDist + matchDistanceSigma*stdDevMatchDistance  (where the mean and std dev are computed using outlier rejection); ignored if not fitting a WCS
# 	Valid Range = [0,inf)
config.astrometry.solver.matchDistanceSigma=2.0

# Deblend sources?
config.doDeblend=True

# Name of coadd
config.coaddName='deep'

# Match sources to reference catalog?
config.doMatchSources=False

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
config.measurement.plugins['base_PixelFlags'].masksFpAnywhere=['CLIPPED']

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
config.measurement.plugins['base_ClassificationExtendedness'].fluxRatio=0.925

config.measurement.plugins.names=['base_CircularApertureFlux', 'base_InputCount', 'base_NaiveCentroid', 'base_PixelFlags', 'base_SkyCoord', 'base_PsfFlux', 'base_Variance', 'base_GaussianCentroid', 'base_SdssCentroid', 'base_GaussianFlux', 'base_SdssShape', 'base_ClassificationExtendedness']
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
config.deblend.maskLimits={}

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
config.deblend.maxFootprintArea=1000000

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
config.deblend.propagateAllPeaks=True

# Mask name for footprints not deblended, or None
config.deblend.notDeblendedMask='NOT_DEBLENDED'

# Only deblend the brightest maxNumberOfPeaks peaks in the parent (<= 0: unlimited)
config.deblend.maxNumberOfPeaks=0

# When splitting stray flux, clip fractions below this value to zero.
config.deblend.clipStrayFluxFraction=0.01

# Name of field in schema with number of deblended children
config.setPrimaryFlags.nChildKeyName='deblend_nChild'

