import os
import warnings
from copy import deepcopy

import numpy as np

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
        self.optics_manager = OpticsManager(user_commands.yaml_docs,
                                            **self.observation_dict)

    def observe(self, orig_source):
        # prepare the observation
        self.optics_manager.update(self.observation_dict)
        self.fov_manager = FOVManager(self.optics_manager)
        self.image_plane = ImagePlane(self.optics_manager.image_plane_header)

        # Make a FOV list - z_order = 0..99
        # Make a image plane - z_order = 100..199
        # Apply Source altering effects - z_order = 200..299
        # Apply FOV specific (3D) effects - z_order = 300..399
        # Apply FOV-independent (2D) effects - z_order = 400..499

        source = deepcopy(orig_source)
        for effect in self.optics_manager.source_effects:
            source = effect.apply_to(source)

        for fov in self.fov_manager.fovs:
            fov.extract_from(source)
            for effect in self.optics_manager.fov_effects:
                fov = effect.apply_to(fov)
            self.image_plane.add(fov)

        for effect in self.optics_manager.image_plane_effects:
            self.image_plane = effect.apply_to(self.image_plane)


class OpticsManager:
    def __init__(self, yaml_docs=[], **kwargs):
        self.optical_elements = [OpticalElement({"name": "misc"})]
        self.meta = {}
        self.meta.update(kwargs)

        if yaml_docs is not None:
            self.load_effects(yaml_docs)

    def load_effects(self, yaml_docs):
        if isinstance(yaml_docs, dict):
            yaml_docs = [yaml_docs]
        self.optical_elements += [OpticalElement(dic) for dic in yaml_docs]

    def add_effect(self, effect, ext=0):
        if isinstance(effect, eff_mod.Effect):
            self.optical_elements[ext].add_effect(effect)

    def update(self, obs_dict):
        self.meta.update(obs_dict)

    def get_all(self, class_type):
        effects = []
        for opt_el in self.optical_elements:
            effects += opt_el.get_all(class_type)

        return effects

    def get_z_order_effects(self, z_level):
        effects = []
        for opt_el in self.optical_elements:
            effects += opt_el.get_z_order_effects(z_level)

        return effects

    @property
    def image_plane_header(self):
        detector_lists = self.get_all(eff_mod.DetectorList)

        if len(detector_lists) != 1:
            warnings.warn("None or more than one DetectorList found. Using the"
                          " first instance.{}".format(detector_lists))

        pixel_scale = self.meta["SIM_DETECTOR_PIX_SCALE"] * u.arcsec
        header = detector_lists[0].image_plane_header(pixel_scale)

        return header

    @property
    def image_plane_effects(self):
        imp_effects = []
        for opt_el in self.optical_elements:
            imp_effects += opt_el.get_z_order_effects([400, 499])

        return imp_effects

    @property
    def fov_effects(self):
        fov_effects = []
        for opt_el in self.optical_elements:
            fov_effects += opt_el.get_z_order_effects([300, 399])

        return fov_effects

    @property
    def source_effects(self):
        src_effects = []
        for opt_el in self.optical_elements:
            src_effects += opt_el.get_z_order_effects([200, 299])

        return src_effects

    @property
    def background_source(self):
        return bg_src

    @property
    def radiometry_table(self):
        return rad_table

    def __add__(self, other):
        self.add_effect(other)

    def __getitem__(self, item):
        if isinstance(item, eff_mod.Effect):
            effects = []
            for opt_el in self.optical_elements:
                effects += opt_el.get_all(item)
            return effects

    def __repr__(self):
        msg = "\nOpticsManager contains {} OpticalElements \n" \
              "".format(len(self.optical_elements))
        for ii, opt_el in enumerate(self.optical_elements):
            msg += '[{}] "{}" contains {} effects \n' \
                   ''.format(ii, opt_el.meta["name"], len(opt_el.effects))

        return msg


class OpticalElement:
    def __init__(self, yaml_dict=None):
        self.meta = {"name": "<empty>"}
        self.properties = {}
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

    def add_effect(self, effect):
        if isinstance(effect, eff_mod.Effect):
            self.effects += [effect]

    def get_all(self, effect_class):
        return [eff for eff in self.effects if isinstance(eff, effect_class)]

    def get_z_order_effects(self, z_level):
        if isinstance(z_level, int):
            zmin = z_level
            zmax = zmin + 99
        elif isinstance(z_level, (tuple, list)):
            zmin, zmax = z_level[:2]
        else:
            zmin, zmax = 0, 500

        effects = []
        for eff in self.effects:
            z = eff.meta["z_order"]
            if isinstance(z, (list, tuple)):
                if any([zmin <= zi <= zmax for zi in z]):
                    effects += [eff]
            else:
                if zmin <= z <= zmax:
                    effects += [eff]

        return effects

    @property
    def z_orders(self):
        z_orders = []
        for eff in self.effects:
            z = eff.meta["z_order"]
            if isinstance(z, (list, tuple)):
                z_orders += z
            else:
                z_orders += [z]

        return z_orders

    @property
    def surface_list(self):
        surf_list = [effect for effect in self.effects
                     if isinstance(effect, (eff_mod.SurfaceList,
                                            eff_mod.TERCurve))]
        return surf_list

    @property
    def mask_list(self):
        mask_list = [effect for effect in self.effects
                     if isinstance(effect, eff_mod.ApertureList)]
        return mask_list

    def __add__(self, other):
        self.add_effect(other)

    def __getitem__(self, item):
        return self.get_all(item)

    def __repr__(self):
        msg = '\nOpticalElement : "{}" contains {} Effects: \n' \
              ''.format(self.meta["name"], len(self.effects))
        eff_str = "\n".join(["[{}] {}".format(i, eff.__repr__())
                             for i, eff in enumerate(self.effects)])
        return msg + eff_str


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


class FOVManager:
    def __init__(self, effects=[], **kwargs):
        self.meta = {}
        self.meta.update(kwargs)
        self.effects = effects
        self._fovs_list = []

    def generate_fovs_list(self):

        # 1. Find the Wavelength range
            # Build from edges of throughput curve,

        # 2. Find the wavelength bins
            # If TraceList and Aperture list, then Spectroscopy
            # TraceList
                # for each trace dlam along the trace centre in increments
                #   of SIM_SUB_PIXEL_FRACTION
                # Must be accompanied by an ApertureList

            # If not, then imaging
            # PSF core increase (atmo, ncpas)
                # If from a files, what is the bin size?
                # If analytic, dlam between a FWHM or SIM_SUB_PIXEL_FRACTION
            # ADC + AD shifts
                # dlam between shift of SIM_SUB_PIXEL_FRACTION

        # 3. Find the spatial range
            # If Spectroscopy
            # ApertureList
                # For each Trace set the sky header to the aperture footprint
                #   plus any shifts from AtmosphericDispersion
                # Set the Image plane footprint centred on the image plane
                #   position

            # If Imaging
            # DetectorList, or ApertureMask, plus any shift from
            #   AtmosphericDispersion

        fovs_list = []
        return fovs_list

    @property
    def fovs(self):
        self.generate_fovs_list()
        return self._fovs_list
