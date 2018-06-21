"""
Configure to force using astrometry_net instead of directly
loading Pan-STARRS1 reference catalogs
"""

from lsst.pipe.tasks.setConfigFromEups import setPhotocalConfigFromEups, setAstrometryConfigFromEups

from lsst.meas.extensions.astrometryNet import LoadAstrometryNetObjectsTask
config.measureCoaddSources.match.refObjLoader.retarget(LoadAstrometryNetObjectsTask)

menu = { "ps1*": {}, # No changes
         "sdss*": { "filterMap": {"y": "z"} }, # No y-band, use z instead
         "2mass*": { "filterMap": {ff:"J" for ff in 'grizy'} }, # No optical; use J
       }

# Both configs need to have the filterMap set
setAstrometryConfigFromEups(config.measureCoaddSources.match.refObjLoader, menu)

config.measureCoaddSources.match.refObjLoader.filterMap["y"] = "z"
