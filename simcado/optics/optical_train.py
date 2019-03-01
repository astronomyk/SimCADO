from copy import deepcopy

from .optics_manager import OpticsManager
from .fov_manager import FOVManager
from .image_plane import ImagePlane


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

    def observe(self, orig_source, **kwargs):
        # prepare the observation
        self.optics_manager.update(**self.observation_dict, **kwargs)
        self.fov_manager = FOVManager(self.optics_manager.fov_setup_effects,
                                      **self.optics_manager.meta)
        self.image_plane = ImagePlane(self.optics_manager.image_plane_header,
                                      **self.optics_manager.meta)

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


