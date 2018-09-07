# We do not have transmission curves attached to our validation repos yet
config.processCcd.isr.doAttachTransmissionCurve = True
# these commissioning data do not have the correct header info to apply the stray light correction
config.processCcd.isr.doStrayLight = False
