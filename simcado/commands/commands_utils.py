import os
from collections import OrderedDict

# import shutil
# import logging
# from ..utils import find_file


def read_config(config_str):
    """
    Read in a SimCADO configuration file

    The configuration file is in SExtractor format:
       'PARAMETER    Value    # Comment'

    Parameters
    ----------
    config_str : str
        the filename of the .config file

    Returns
    -------
    config_dict : dict (collections.OrderedDict)
        A dictionary with keys 'PARAMETER' and values 'Value'.

    Notes
    -----
    The values of the dictionary are strings and will have to be converted to
    the appropriate data type as they are needed.
    """

    if isinstance(config_str, str):
        if os.path.exists(config_str):
            with open(config_str, 'r') as fp1:
                lines = fp1.readlines()
        else:
            lines = [line for line in config_str.split("\n") if len(line) > 0]
    else:
        raise ValueError("config_str must be a filename or multi-line string")

    config_dict = lines_to_dict(lines)

    return config_dict


def lines_to_dict(lines):
    """
    Create a OrderedDict from a list of lines

    Parameters
    ----------
    lines : list
        A list of strings read in from a multi-line string or a file-like object
    Returns
    -------
    config_dict : OrderedDict

    """

    import re

    if not isinstance(lines, (tuple, list)):
        raise ValueError("input must be a list of strings")

    # remove lines that are all spaces or spaces + '#'
    # these are the regular expressions
    isempty = re.compile(r'^\s*$')
    iscomment = re.compile(r'^\s*#')

    # Read the file into a list of strings
    config_dict = OrderedDict()

    for line in lines:
        if isempty.match(line):
            continue
        if iscomment.match(line):
            continue

        line = line.rstrip()             # remove trailing \n
        content = line.split('#', 1)[0]  # remove comment
        try:
            param, value = content.split(None, 1)
        except:
            raise ValueError("Line is missing value: \n {}".format(line))

        value = str_to_python_type(value)
        config_dict[param] = value

    return config_dict


def update_config(new_config, old_config):
    """
    Update a SimCADO configuration dictionary

    A configuration file in the SExtractor format::

        'PARAMETER    Value    # Comment'

    Parameters
    ----------
    new_config : str
        the filename of the .config file

    Returns
    -------
    config_dict : dict
        A dictionary with keys 'PARAMETER' and values 'Value'.

    Notes
    -----
    the values of the dictionary are strings and will have
    to be converted to the appropriate data type as they are needed.

    """

    if isinstance(new_config, str):
        old_config.update(read_config(new_config))
    elif isinstance(new_config, dict):
        old_config.update(new_config)
    else:
        raise ValueError("new_config must be `str` or `dict`")

    return old_config


def str_to_python_type(input_str):
    if isinstance(input_str, str):
        # Convert to number if possible
        try:
            input_str = float(input_str.strip())
        except ValueError:
            input_str = input_str.strip()

            # Convert string "none" to python None
            if input_str.strip().lower() == "none":
                input_str = None
            # Convert string booleans to python booleans
            elif input_str.strip().lower() == "true":
                input_str = True
            elif input_str.strip().lower() == "false":
                input_str = False

    return input_str


def convert_dict_strings_to_python_types(dic):
    """Converts all str python types to corresponding python type """

    for key in dic:
        dic[key] = str_to_python_type(dic[key])

    return dic


def is_item_subcategory(item, dic):
    return item.upper() in set([key.split("_")[0] for key in dic])


def get_subcategory(item, dic):
    return {key : dic[key] for key in dic if item.upper() in key.split("_")[0]}


def extract_filter_name(path_name):
    filt_name = os.path.basename(path_name).split(".")[0]
    return filt_name.replace("TC_filter_", "")


################################################################################
#                          Dump files to disk                                  #
################################################################################
#
# def dump_mirror_config(path=None, what="scope"):
#     """
#     Dump the EC_mirrors_scope.tbl or the EC_mirrors_ao.tbl to disk
#
#     Parameters
#     ----------
#     path : str, optional
#         path where the mirror configuration file is to be saved
#     what : str, optional
#         ["scope", "ao"] dump the mirror configuration for either the telescope
#         or the AO module
#     """
#
#     if what.lower() == "scope":
#         print("Dumping telescope mirror configuration.")
#         fname = find_file("EC_mirrors_scope.tbl")
#     elif what.lower() == "ao":
#         print("Dumping AO mirror configuration.")
#         fname = find_file("EC_mirrors_ao.tbl")
#
#     if path is None:
#         with open(fname, "r") as fd1:
#             print(fd1.read())
#     else:
#         path = os.path.dirname(path)
#         shutil.copy(fname, path)
#
#
# def dump_chip_layout(path=None):
#     ## OC, 2016-08-11: changed parameter from 'dir' (redefines built-in)
#     """
#     Dump the FPA_chip_layout.dat file to a path specified by the user
#
#     Parameters
#     ----------
#     path : str, optional
#         path where the chip layout file is to be saved
#     """
#     fname = find_file("FPA_chip_layout.dat")
#
#     if path is None:
#         with open(fname, "r") as fd1:
#             print(fd1.read())
#     else:
#         path = os.path.dirname(path)
#         shutil.copy(fname, path)
#         logging.debug("Printed chip layout to file: "+path+"/"+fname)
#
#
# def dump_defaults(filename=None, selection="freq"):
#     ## OC, 2016-08-11: changed parameter from 'type' to 'selection' as
#     ##    'type' redefines built-in
#     """
#     Dump the frequent.config file to a path specified by the user
#
#     Parameters
#     ----------
#     filename : str, optional
#         path or filename where the .config file is to be saved
#     selection : str, optional
#         ["freq", "all"] amount of keywords to save. "freq" only prints the most
#         frequently used keywords. "all" prints all of them
#     """
#
#     from .. import rc
#
#     if "freq" in selection.lower():
#         fname = "frequent.config"
#     elif "all" in selection.lower():
#         fname = "default.config"
#
#     if filename is None:
#         gname = os.path.join(rc.__data_dir__, fname)
#         with open(gname, "r") as fd1:
#             print(fd1.read())
#         return None
#     else:
#         path, gname = os.path.split(filename)
#         if path == "":
#             path = "."
#
#         if gname == "":
#             gname = fname
#         shutil.copy(os.path.join(rc.__data_dir__, fname),
#                     os.path.join(path, gname))