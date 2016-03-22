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
ATTR = 'dev'
VERSION = '%d.%d%s' % (MAJOR, MINOR, ATTR)

## Is this needed?
def write_version_py(filename='SimCADO/version.py'):
    cnt = """
# THIS FILE GENERATED FROM SIMCADO SETUP.PY
version = '{}'
"""
    a = open(filename, 'w')
    try:
        a.write(cnt.format(VERSION))
    finally:
        a.close()

from distutils.core import setup

def setup_package():
    # Rewrite the version file everytime
    write_version_py()

    setup(name = 'SimCADO',
          version = VERSION,
          description = "MICADO Instrument simulator",
          author = "Kieran Leschinski, Oliver Czoske",
          author_email = "karl.jansky@univie.ac.at,oliver.czoske@univie.ac.at",
          url = "none",
          package_dir={'SimCADO': 'SimCADO'}, 
          packages=['SimCADO'],
          scripts = ['scripts/simcado.py'],
          package_data = {'SimCADO': ['data/*']},
          )
    
    
if __name__ == '__main__':
    setup_package()

