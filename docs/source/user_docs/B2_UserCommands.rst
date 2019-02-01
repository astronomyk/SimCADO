.. todo:: Write this doc page

Things to be added / updated

----

UserCommands can still be called without referencing any file is particular
This returns an object with all default parameters, i.e. lots of Nones

----

Load a package (assuming already downloaded) by passing::

    >>> sim.UserCommands(instrument="MICADO")

Load the observing mode with ``mode="MODE_<name>"``. See further down the page

+ Reference downloading packages
+ Reference setting instruments and modes

----

Single parameters can be updated by using a dictionary call:
UserCommands[KEYWORD] = my_value
This won't work for setting instruments or modes, or filters

----

Multiple Parameters can be updated with .update(), specifying
a filename,
a string, or
a dictionary

----

Prepackaged instruments can be loaded by using .set_instrument(). Check to see
which packages are downloaded by calling simcado.server.get_local_package()

+ Reference downloading instrument packages

List which observation modes are available with UserCommands.list("modes").
Load a mode with UserCommands.set_mode().

----

Making your own mode is as simple as creating an ASCII file with the relevant
KEYWORD VALUE pairs in SExtractor format ("KEYWORD  value  # Comment") and
reading this file in with UserCommands.update(filename)

+ Reference making your own instrument package
+ Reference dumping .default.config
+ .. todo:: reinstate dump_config(filename)
+ .. todo:: write() method for UserCommands

----

Selecting a filter can be done with .select_filter(). Both the simple name or
a file name/path are accepted. Check which filters are available with
.list("filters")

----

PSFs can be set by passing the name of a PSF from the local PSF database to
UserCommands["SCOPE_PSF_FILE"].
The PSFs names can be listed with simcado.server.list_psfs()

+ Reference Downloading packages # PSF files

Alternatively, user PSF files can be used by passing the file path to
UserCommands["SCOPE_PSF_FILE"].

