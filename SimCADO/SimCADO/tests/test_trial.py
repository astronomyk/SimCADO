# content of test_sample.py
def func(x):
    return x + 1

def test_right_answer():
    assert func(3) == 4

def run_simcado():
    import simcado as sim
    import os
    cmd = sim.UserCommands()
    cmd["OBS_EXPTIME"] = 3600
    cmd["OBS_NDIT"] = 1
    cmd["INST_FILTER_TC"] = "J"

    src = sim.optics_utils.source_1E4_Msun_cluster()
    opt = sim.OpticalTrain(cmd)
    fpa = sim.Detector(cmd)

    src.apply_optical_train(opt, fpa)
    fpa.read_out("my_output.fits")
    
    is_there = os.path.exists("my_output.fits")
    os.remove("my_output.fits")
    
    return is_there

def test_run_simcado():
    assert run_simcado() == True
