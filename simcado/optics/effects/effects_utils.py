from copy import deepcopy

import numpy as np
from astropy import units as u

import simcado.optics.effects.spectroscopy_effects
import simcado.optics.effects.ter_curves
from ... import utils
from ..radiometry_utils import empty_surface_list
from .. import effects as efs


def combine_surface_effects(surface_effects):
    surflist_list = [eff for eff in surface_effects
                     if isinstance(eff, efs.SurfaceList)]
    surf_list = [eff for eff in surface_effects
                 if isinstance(eff, simcado.optics.effects.ter_curves.TERCurve)]

    if len(surflist_list) == 0:
        tbl = empty_surface_list()
        tbl.meta["name"] = "Radiometry Table"
        surflist_list += [tbl]

    new_surflist = deepcopy(surflist_list[0])
    for surflist in surflist_list[1:]:
        new_surflist.add_surface_list(surflist)

    for surf in surf_list:
        new_surflist.add_surface(surf, surf.meta["name"])

    return new_surflist


def get_all_effects(effects, effect_class):
    return [eff for eff in effects if isinstance(eff, effect_class)]


def make_effect(effect_dict, **super_kwargs):
    effect_meta_dict = {key : effect_dict[key] for key in effect_dict
                        if key not in ["class", "kwargs"]}
    effect_class_name = effect_dict["class"]
    effect_cls = getattr(efs, effect_class_name)

    effect_kwargs = {}
    if "kwargs" in effect_dict:
        effect_kwargs = effect_dict["kwargs"]
    effect_kwargs.update(effect_meta_dict)
    effect_kwargs.update(super_kwargs)

    effect = effect_cls(**effect_kwargs)
    effect.meta.update(effect_meta_dict)

    return effect


def is_spectroscope(effects):
    has_trace_lists = sum([isinstance(eff,
                                      simcado.optics.effects.spectroscopy_effects.SpectralTraceList)
                           for eff in effects])
    has_apertures = sum([isinstance(eff, (
    simcado.optics.effects.spectroscopy_effects.ApertureList,
    simcado.optics.effects.spectroscopy_effects.ApertureMask))
                         for eff in effects])

    return bool(has_apertures and has_trace_lists)
