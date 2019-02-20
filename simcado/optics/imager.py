import os
import warnings
from copy import deepcopy

from astropy.io import ascii as ioascii
from astropy.io import fits
from astropy.table import Table, vstack
from astropy import units as u

import simcado.optics.effects.effects
from ..utils import find_file
from . import effects as eff_mod
from .image_plane import ImagePlane

from ..source.source2 import Source


class OpticalTrain:
    def __init__(self, cmds=None):
        self.observation_dict = None
        self.optics_manager = None
        self.radiometry_table = None
        self.fov_manager = None
        self.image_plane = None
        self.yaml_docs = None

        if cmds is not None:
            self.load(cmds)

    def load(self, user_commands):
        # UserCommands contains the filenames of which yaml files to use as
        # well what the observational parameters will be

        self.observation_dict = user_commands.cmds
        self.yaml_docs = user_commands.yaml_docs
        self.optics_manager = OpticsManager(user_commands.yaml_docs)

    def observe(self, source):
        # prepare the observation
        self.optics_manager.update(self.observation_dict)
        self.fov_manager = FOVManager(self.optics_manager,
                                      self.observation_dict)
        self.image_plane = ImagePlane(self.optics_manager.image_plane_header)

        source.rotate(self.observation_dict["PUPIL_ANGLE_OFFSET"])
        source = source * self.optics_manager.throughput
        source.append(self.optics_manager.background_source)

        # run the observation
        for fov in self.fov_manager.fovs_list:
            fov.extract_from(source)
            for effect in self.optics_manager.sorted_effects:
                fov = effect.apply_to(fov)
            self.image_plane.add(fov)


class OpticsManager:
    def __init__(self, yaml_docs=None):

        self.optical_elements = None
        self.radiometry_table = None

        if yaml_docs is not None:
            self.load_effects(yaml_docs)

    def load_effects(self, yaml_docs):
        self.optical_elements = [OpticalElement(dic) for dic in yaml_docs]

    def make_radiometry_table(self, filename):
        self.radiometry_table = Table

    @property
    def sorted_effects(self):
        sorted_effects_list = []
        return sorted_effects_list

    @property
    def image_plane_header(self):
        header = fits.Header()
        return header

    @property
    def background_source(self):
        spec = self.radiometry_table.emission
        bg_src = Source
        return Source


class FOVManager:
    def __init__(self, optics_manager, observation_dict):
        self.fovs_list = []

    def generate_fovs_list(self, optics_manager, observation_dict):
        fovs_list = [None]
        return fovs_list


class OpticalElement:
    def __init__(self, yaml_dict=None):
        self.meta = None
        self.properties = None
        self.effects = []

        if isinstance(yaml_dict, dict):
            self.meta = {key : yaml_dict[key] for key in yaml_dict
                         if key not in ["properties", "effects"]}
            if "properties" in yaml_dict:
                self.properties = yaml_dict["properties"]
            if "effects" in yaml_dict:
                self.effects_dicts = yaml_dict["effects"]
                self.make_effects(yaml_dict["effects"])

    def make_effects(self, effects_dicts):
        for effdic in effects_dicts:
            self.effects += [make_effect(effdic, **self.properties)]

    @property
    def surface_list(self):
        surf_list = [effect for effect in self.effects
                     if isinstance(effect,
                                   simcado.optics.effects.effects.SurfaceList)]
        return surf_list

    @property
    def mask_list(self):
        mask_list = [effect for effect in self.effects
                     if isinstance(effect, eff_mod.MaskList)]
        return mask_list


def make_effect(effect_dict, **super_kwargs):
    effect_meta_dict = {key : effect_dict[key] for key in effect_dict
                        if key not in ["class", "kwargs"]}
    effect_class_name = effect_dict["class"]
    effect_cls = getattr(eff_mod, effect_class_name)

    effect_kwargs = effect_dict["kwargs"]
    effect_kwargs.update(super_kwargs)

    effect = effect_cls(**effect_kwargs)
    effect.meta.update(effect_meta_dict)

    return effect


class RadiometryTable:
    surfaces_list = []


class DetectorArray:
    detectors_list = []

