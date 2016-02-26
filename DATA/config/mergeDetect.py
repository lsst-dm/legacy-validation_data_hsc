import lsst.pipe.tasks.multiBand
assert type(config)==lsst.pipe.tasks.multiBand.MergeDetectionsConfig, 'config is of type %s.%s instead of lsst.pipe.tasks.multiBand.MergeDetectionsConfig' % (type(config).__module__, type(config).__name__)
# Minimum distance from closest peak to create a new one (in arcsec).
config.minNewPeak=1.0

# When adding new catalogs to the merge, all peaks less than this distance  (in arcsec) to an existing peak will be flagged as detected in that catalog.
config.maxSamePeak=0.3

# Priority-ordered list of bands for the merge.
config.priorityList=['HSC-I', 'HSC-R', 'HSC-Z', 'HSC-Y', 'HSC-G']

# Always keep this many peaks in each family
# 	Valid Range = [1,inf)
config.cullPeaks.rankSufficient=20

# Keep peaks with less than this normalized rank that also match the rankConsidered condition.
# 	Valid Range = [0.0,inf)
config.cullPeaks.rankNormalizedConsidered=0.7

# Keep peaks with less than this rank that also match the rankNormalizedConsidered condition.
# 	Valid Range = [1,inf)
config.cullPeaks.rankConsidered=30

# Always keep peaks detected in this many bands
# 	Valid Range = [1,inf)
config.cullPeaks.nBandsSufficient=2

# Name of coadd
config.coaddName='deep'

