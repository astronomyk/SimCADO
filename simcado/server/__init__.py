## Author : Kieran
## Date-created :  2018-11-09
## Date-modified : 2018-11-09
## Changes :
##    2018-11-09 Added keywords INST_PKG_LOCAL_PATH and INST_PKG_SERVER_PATH
##


from os.path import dirname
from inspect import getfile, currentframe

# Package directory
PKG_DIR = dirname(getfile(currentframe()))
__PKG_DIR__ = PKG_DIR       # for backwards compatibility


