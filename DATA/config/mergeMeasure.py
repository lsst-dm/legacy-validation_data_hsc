import lsst.pipe.tasks.multiBand
assert type(config)==lsst.pipe.tasks.multiBand.MergeSourcesConfig, 'config is of type %s.%s instead of lsst.pipe.tasks.multiBand.MergeSourcesConfig' % (type(config).__module__, type(config).__name__)
# Priority-ordered list of bands for the merge.
config.priorityList=['HSC-I', 'HSC-R', 'HSC-Z', 'HSC-Y', 'HSC-G']

# Name of coadd
config.coaddName='deep'

