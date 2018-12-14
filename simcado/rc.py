import os
import inspect
import yaml

__pkg_dir__ = os.path.dirname(inspect.getfile(inspect.currentframe()))
__data_dir__ = os.path.join(__pkg_dir__, "data")

# For utils.find_file()
__search_path__ = ['./', __pkg_dir__, __data_dir__]

# load in settings from rc file ".simcadorc"
with open(os.path.join(__pkg_dir__, ".simcadorc"), "r") as rc_file:
    __rc__ = yaml.load(rc_file)

