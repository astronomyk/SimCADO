#!/usr/bin/env python
"""
SimCADO: A python package to simulate MICADO
"""

## BEFORE importing distutils, remove MANIFEST. distutils doesn't properly
## update it when the contents of directories change.

#import os
#if os.path.exists('MANIFEST'):
#    os.remove('MANIFEST')

def get_old_version(filename='simcado/version.py'):
    f = open(filename, "r")
    for i in range(3): vers = f.readline()
    print(vers)
    vers = vers.replace("dev", "").replace("'","").split(".")[-1]
    return int(float(vers))

    
    
# Is this the version number scheme that we want?
MAJOR = 0
MINOR = 2
TINY = get_old_version() + 1
ATTR = 'dev'
VERSION = '%d.%d.%d%s' % (MAJOR, MINOR, TINY, ATTR)



## Is this needed?
def write_version_py(filename='simcado/version.py'):
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
          url = "http://homepage.univie.ac.at/kieran.leschinski/",
          package_dir={'simcado': 'simcado'}, 
          packages=['simcado', 'simcado.tests'],
          scripts = ['scripts/simcado.py'],
          include_package_data=True,
          package_data = {'simcado': ['data/*', 'docs/*']},
          )
    
    
if __name__ == '__main__':
    setup_package()

