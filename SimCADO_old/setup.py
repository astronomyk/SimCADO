#!/usr/bin/env python
"""SimCADO: A python package to simulate MICADO

"""

## BEFORE importing distutils, remove MANIFEST. distutils doesn't properly
## update it when the contents of directories change.

#import os
#if os.path.exists('MANIFEST'):
#    os.remove('MANIFEST')

# Is this the version number scheme that we want?
MAJOR = 0
MINOR = 1
MICRO = 0
VERSION = '%d.%d.%d' % (MAJOR, MINOR, MICRO)

## Is this needed?
#def write_version_py(filename='version.py'):
#    cnt = """
## THIS FILE GENERATED FROM SIMCADO SETUP.PY
#__version__ = '%(version)s'
#"""
#    a = open(filename, 'w')
#    try:
#        a.write(cnt % {'version': VERSION})
#    finally:
#        a.close()

from distutils.core import setup

def setup_package():
    # Rewrite the version file everytime
    #write_version_py()

    setup(name = 'SimCADO',
          version = VERSION,
          description = "MICADO Instrument simulator",
          author = "Kieran Leschinski, Oliver Czoske",
          author_email = "kleschinski@univie.ac.at, oczoske@univie.ac.at",
          url = "none",
          package_dir={'SimCADO': ''}, # SimCADO is in the top-level dir
          packages=['SimCADO'],
          scripts = [],
          package_data = {'SimCADO': ['data/*']},
          )
    
    
if __name__ == '__main__':
    setup_package()

