from copy import deepcopy

from simcado.optics import effects as efs
from simcado.optics.radiometry_utils import empty_surface_list


def combine_radiometry_effects(surfaces):
    surflist_list = [eff for eff in surfaces if isinstance(eff,
                                                           efs.SurfaceList)]
    surf_list = [eff for eff in surfaces if isinstance(eff, efs.TERCurve)]

    if len(surflist_list) == 0:
        tbl = empty_surface_list()
        tbl.meta["name"] = "Radiometry Table"
        surflist_list += [tbl]

    rad_table = deepcopy(surflist_list[0])
    for surflist in surflist_list[1:]:
        rad_table.add_surface_list(surflist)

    for surf in surf_list:
        rad_table.add_surface(surf, surf.meta["name"])

    return rad_table


def _get_all(effects, effect_class):
    return [eff for eff in effects if isinstance(eff, effect_class)]


def make_effect(effect_dict, **super_kwargs):
    effect_meta_dict = {key : effect_dict[key] for key in effect_dict
                        if key not in ["class", "kwargs"]}
    effect_class_name = effect_dict["class"]
    effect_cls = getattr(efs, effect_class_name)

    effect_kwargs = effect_dict["kwargs"]
    effect_kwargs.update(super_kwargs)

    effect = effect_cls(**effect_kwargs)
    effect.meta.update(effect_meta_dict)

    return effect


def is_spectroscope(effects):
    has_trace_lists = sum([isinstance(eff, efs.TraceList) for eff in effects])
    has_apertures = sum([isinstance(eff, (efs.ApertureList, efs.ApertureMask))
                         for eff in effects])

    return bool(has_apertures and has_trace_lists)
