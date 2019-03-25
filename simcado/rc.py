import inspect
import os
from os.path import join

# Search path for finding files
__pkg_dir__ = os.path.dirname(inspect.getfile(inspect.currentframe()))
__data_dir__ = join(__pkg_dir__, "data")
__search_path__ = ['./', __pkg_dir__, __data_dir__]
