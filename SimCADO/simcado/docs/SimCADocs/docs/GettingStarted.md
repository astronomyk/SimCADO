# Getting Started with SimCADO
SimCADO can be super easy to use, or super complicated. The level of complexity is completely up to the user. Regardless of your intended use for SimCADO, it's probably a good idea to at least have a vague understanding of what is going on during a simulation.

## Behind the scenes of SimCADO
SimCADO uses 4 main classes during a simulation: 
* `Source` holds spatial and spectral information about the astronomical source of photons, e.g. galaxy, star cluster, etc.
* `OpticalTrain` contains information on the various elements along the optical path, e.g. mirrors reflectivity curves, PSFs, instrumental distortion, etc.
* `Detector` represents the focal plane detector array and contains information on the electronic characteristics of the detector chips and their physical positions.
* `UserCommands` is a dictionary of all the important keywords needed by SimCADO to run the simultationm, e.g. `OBS_EXPTIME` (exposure time) or `INST_FILTER_TC` (filter curve)

## Source
Create a Source object by calling >>> my_src = sim.Source(...)

The Source object holds 6 arrays

`.x, .y` - [each a numpy.ndarray]. Coordinates for each point source in the image in units of [arcsec] from the focal plane centre
lam - [numpy.ndarray]. An array with the centre of the wavelength bins in [um] for each unique spectrum
spectra - [numpy.ndarray]. An (n, m) array holding n spectra, each with m values. Default units are [ph/s]
Note - lam and spectra should use use equal bin width. Variable bin widths leads to unpredictable results

ref - [numpy.ndarray]. An array to connect the point source at x[i], y[i] to a unique spectrum at spectra[j], i.e. ref[i] = j
weight - [numpy.ndarray]. If two sources share the same spectrum, but are at different distances or have different luminosities a scaling factor can be specified to the spectrum when applied to each specific point source
The Source object can be created by reading in a previously generated Source FITS file, or by passing the arrays needed to construct the Source object from scratch:

sim.Source(filename=None, lam=None, spectra=None, x=None, y=None, ref=None, weight=None, **kwargs)

If filename is specified, the other arguments are ignored. If filename == None, it is expected that at least lam, spectra, x, y, ref are specified. weight is optional. Default is weight[:] = 1

Optional keywords can be specified:

units [default ph/s] is the units for the spectra, i.e. x phontons per second per spectral bin.
Use the built in Source functions from the module optic_utils
Here we create a Source object from a 1E4 Msun open cluster in the LMC (50 kpc) with the function

>>> sim.source.source_1E4_Msun_cluster()


## OpticalTrain
## Detector
## UserCommands

## Simulating with SimCADO
## The quick, the dirty and the ugly
### Creating `Source` objects

## Building your own simulation run
### Creating `UserCommand` objects
### Creating `OpticalTrain` objects
### Creating `Detector` array objects
### Running the simulation

## Things to watch out for





