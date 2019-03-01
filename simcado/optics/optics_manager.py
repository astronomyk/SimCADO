import warnings

from simcado.optics import effects as efs
from simcado.optics.effects.effects_utils import combine_radiometry_effects
from simcado.optics.optical_element import OpticalElement


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
        if isinstance(effect, efs.Effect):
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
        detector_lists = self.get_all(efs.DetectorList)
        header = detector_lists[0].image_plane_header

        if len(detector_lists) != 1:
            warnings.warn("None or more than one DetectorList found. Using the"
                          " first instance.{}".format(detector_lists))

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
        src_effects = [self.radiometry_table]
        for opt_el in self.optical_elements:
            src_effects += opt_el.get_z_order_effects([200, 299])

        return src_effects

    @property
    def image_plane_setup_effects(self):
        implane_setup_effects = []
        for opt_el in self.optical_elements:
            implane_setup_effects += opt_el.get_z_order_effects([100, 199])

        return implane_setup_effects

    @property
    def fov_setup_effects(self):
        fovmanager_effects = [self.radiometry_table]
        for opt_el in self.optical_elements:
            fovmanager_effects += opt_el.get_z_order_effects([0, 99])

        return fovmanager_effects

    @property
    def background_source(self):
        bg_src = None
        return bg_src

    @property
    def radiometry_table(self):
        surface_like_effects = []
        for opt_el in self.optical_elements:
            surface_like_effects += opt_el.ter_list

        rad_table = combine_radiometry_effects(surface_like_effects)

        return rad_table

    def __add__(self, other):
        self.add_effect(other)

    def __getitem__(self, item):
        if isinstance(item, efs.Effect):
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