import os
import yaml

YAMLS_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                          "../yamls/"))


def _tiny_yaml_dict():
    text = """
# ATMOSPHERE
object : atmosphere
name : armazones

properties :
    temperature : 0         # [-270..270] deg C
    pressure : 0.6          # [0..1] bar

effects :
-   name : super_psf
    class : GaussianDiffractionPSF
    z_order : 0
    kwargs :
        diameter : 39
"""
    return yaml.load(text)


def _inst_yaml_dict():
    text = """
# INSTRUMENT OPTICS
object : instrument
name : micado_wide_field
z_order : 3
inst_pkg_name : micado

properties :
    temperature : -190
    plate_scale : 0.004

effects :
-   name : micado_surface_list
    class : SurfaceList
    kwargs :
        file_name : micado_mirror_list.tbl

-   name : micado_adc
    class : AtmosphericDispersion
    kwargs :
        zenith_distance : 30
        reverse_shifts : True

-   name : pupil_mask
    class : ApertureList
    kwargs :
        file_name : aperture_list.tbl        
        
    """
    return yaml.load(text)
