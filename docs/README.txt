Notes on compiling these docs
=============================
They are an absolute pain in the posterior

For whatever reason, numpydoc is broken for versions > 0.6.0, and this only
works together with sphinx <= 1.6.7, hence we need to re-install these specific
versions of the packages::

    $ pip install sphinx ==1.6.7 --force-reinstall
    $ pip install numpydoc ==0.6.0 --force-reinstall

Also, with the current directory tree the package path is three levels higher::

    sys.path.insert(0, os.path.abspath('../../..'))

Make sure the version.py file isn't in the .gitignore

To exclude modules that we don't want in the apidocs, just import a package
that isn't in the ``requirements.readthedocs.txt`` file. E.g.

    import pytest




