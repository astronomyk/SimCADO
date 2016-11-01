# Welcome to SimCADocs
The (slowly expanding) documentation base for SimCADO

## SimCADO in a nutshell
SimCADO is a python package designed to simulate the effects of the Atmosphere, E-ELT, and MICADO instrument on incoming light. The current version (v0.2) can simulate the MICADO imaging modi (4mas and 1.5mas per pixel in the wavelength range 0.7µm to 2.5µm).


### Downloading and Installing
The quick way:

    $ pip install http://homepage.univie.ac.at/kieran.leschinski/SimCADO/SimCADO-0.2dev.zip

For slightly more options, see the [Downloads](Download.md) section


### Reference Material
* The inner workings of SimCADO are described in detail in [Leschinski et al. (2016)](https://arxiv.org/pdf/1609.01480v1.pdf)

* The current status of MICADO is described in [Davies et al. (2016)](https://arxiv.org/pdf/1607.01954.pdf)

## Running a simulation in 3 lines

The easiest way to run a simulation is to create, or load, a Source object and then call the `.run()` command. If you specify a filename, the resulting image will be output to a FITS file under that name. If you do not specify a filename, the output will be returned to the console/notebook as an `astropy.io.fits.HDUList` object.

To begin, we will import the simcado module (assuming it is already installed).

    >>> import simcado as sim

At the very least, we need to create a `Source` object which contains both spatial and spectral information on our object of interest. Here we use the built-in command `.source.source_1E4_Msun_cluster()` to create a `Source`-object for a 10000-Msun stellar cluster. (See [Creating Sources](examples/Source.md) for more information).

    >>> src = sim.source.source_1E4_Msun_cluster()

We now pass the `source` object through SimCADO. This is as easy as calling `.run()`. If we specify a `filename`, SimCADO will write the output to disk in the form of a FITS file. If no `filename` is given, then SimCADO returns an `astropy.io.fits` object to the console/notebook.

    >>> sim.run(src, filename="my_first_sim.fits")

That's it. Of course SimCADO can also go in the other direction, providing many more levels of complexity, but for that the reader is directed to the examples pages and/or the [API](API/_build/index.html) documentation

## SimCADO building blocks
For a brief explanation of how SimCADO works and which classes are relevant, please see either the [Getting Started](GettingStarted.md) or [SimCADO in depth](deep_stuff/SimCADO.md) section.

## Contact

For questions and complaints alike, please contact the authors:

* [kieran.leschinski@univie.ac.at]()
* [oliver.czoske@univie.ac.at]()

**Developers (Vienna):** Kieran Leschinski, Oliver Czoske

**Data Flow Team Leader (Gronigen):** Gijs Verdoes Kleijn

**MICADO home office (MPE):** http://www.mpe.mpg.de/ir/micado
