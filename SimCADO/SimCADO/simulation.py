"""
simualtion.py
"""

from astropy.io import ascii
import simcado as sim

def run(src, cmds=None, opt_train=None, fpa=None, 
        small_fov=True, filename=None, return_internals=False):
    """
    Run a MICADO simulation with default parameters
    """
    
    if cmds is None: 
        cmds = commands.UserCommands()
        cmds["INST_FILTER_TC"] = "J"
    
    if small_fov:
        chip_layout = """#  id    x_cen    y_cen   x_len   y_len
                         #       arcsec   arcsec   pixel   pixel
                             0        0        0    1024    1024"""
        cmds["FPA_CHIP_LAYOUT"] = chip_layout
        
    if opt_train is None:
        opt_train = optics.OpticalTrain(cmds)
    if fpa is None:
        fpa = detector.Detector(cmds)
     
    print(fpa.layout)
    src.apply_optical_train(opt_train, fpa)
    
    if filename is not None:
        fpa.read_out(filename=filename, to_disk=True)
    else: 
        hdu = fpa.read_out()
        if return_internals:
            return hdu, (cmds, opt_train, fpa)
        else:
            return hdu