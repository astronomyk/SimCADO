###############################################################################
# SpectralGrating
#
# DESCRIPTION
#
# "The universe is a pretty big place. If it's just us, seems like an awful
#   waste of space." - Carl Sagan, Contact
#
# This is a place holder for the SpectralGrating object. The main point of this
# object to imitate the use of a grating. The LightObject should be turned on
# its side, so that the spectral dimension is splayed out over the x,y plane.
#
# We can achieve this in two ways:
# - rotate the input cube and then use a PlaneEffect to translate the photons
#   from each spectral slice to the correct pixels, or
# - generate an ADC-like delta PSF for all the 1nm wide slices and apply each
#   PSF to the LightObject data.
#
# The SpectralGrating object should include information for: where photons of
# each wavelength will land, what percentage of photons get there, the scatter
#
# Talk to Wolfgang about what else needs to go where
#
#
#
#
#
#
