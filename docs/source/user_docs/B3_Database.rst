.. todo:: write this doc page

SimCADO can access a series of instrument, and telescope packages. These
packages are put together by Kieran. These are accessed by the simcado.server
module.
Additionally there are various PSF packages available, as well as basic data for
building the source module objects

----

To make use of these packages with simcado, the user needs to set up a
directory for saving and tracking these packages. This is called with
simcado.server.set_up_local_package_directory(path_name)
When starting an ipython session, simcado should be told about this path by
setting the rc value::

    >>> import simcado as sim
    >>> sim.rc.__rc__["FILE_LOCAL_DOWNLOADS_PATH"] = "path/to/simcado/downloads"

By default this is set to "D:/simcado_downloads/"

----

Avoid this by setting the value in the .simcadorc file
The file is found in the package directory, which can be found with::

    >>> simcado.__pkg_dir__

----

Find which packages are available on the server with .get_server_package()
Display remote and local packages with .list_instruments(), .list_psfs(), or
.list_source_packages().
If you like scrollling, use list_all()

----

Download a specific package with .download_package(package_name)
It doesn't matter which part a package belongs to, .download_package() searches
all database files on the server and returns the first matching entry it finds

----

List which packages you already have with .get_local_packages()

----

Load a package into simcado by calling a UserCommands object with the package
name specified::

    >>> sim.UserCommands(instrument="MICADO")

As long as the MICADO package has been downloaded, and simcado knows where to
look, it should find the package data.

----

The same goes for downloading PSFs: ``.download_package(psf_name)``

----

Making your own packages

+ Reference needed data formats
+ Reference package structure
+ Reference dumping
