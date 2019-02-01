Changes between SimCADO v0.4 and SimCADO 1.0
============================================

There have been several major changes in the recent release of SimCADO. We have
tried to maintain backwards compatibility, however in some cases this just
wasn't feasible. These changes are listed below.

If you find a change which isn't listed here, please send a brief email to
kieran.leschinski@univie.ac.at, and we will update this list.

UserCommands
------------

* ``UserCommands`` **can still be called** directly from the top level
  (``simcado.UserCommands``), however the code for the class has been moved to
  ````simcado.commands.user_commands.UserCommands``
* ``UserCommands`` **no longer** contains subcategory dictionaries,
  e.g. ``<UserCommands>.inst``. Instead, subcategory dictionaries can be
  accessed by using the category stub as a dictionary item call::

      cmd = simcado.UserCommands()
      atmo_dict = cmd["ATMO"]
      inst_dict = cmd["INST"]

* An empty ``UserCommand`` initialisation will return a ``UserCommand`` object
  with **NO** instrument specific values.::

      >>> cmd = simcado.UserCommands()
      >>> print(cmd["INST_FILTER_TC"])
      None

  In other words, it will be a skeleton
  object. To load instrument specific parameters, a downloaded instrument
  package is needed. See the new `UserCommands <B2_UserCommands.rst>`_ and
  `Database <B3_Database.rst>`_ sections

Database
--------
The instrument package server was not included in the original version of
SimCADO. As such all of this is new.
See the new `UserCommands <B2_UserCommands.rst>`_ and
`Database <B3_Database.rst>`_ sections for how to use instrument packages.




