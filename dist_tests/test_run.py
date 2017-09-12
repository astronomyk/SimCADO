import numpy as np
import simcado as sim

def test_full_sim_with_single_star_all_defaults():
    cmd = sim.UserCommands()
    opt = sim.OpticalTrain(cmd)
    fpa = sim.Detector(cmd)
    
    src = sim.source.star(mag=20, filter_name="Ks")
    src.apply_optical_train(opt, fpa)
    
    hdu = fpa.read_out()
    
    # "It’s over nine thousand!" - Vegeta Breigh
    assert np.max(hdu[0].data) > 9000

    
def test_quick_sim_with_single_star_all_defaults():
    src = sim.source.star(mag=20, filter_name="Ks")
    hdu = sim.run(src)
    
    assert np.max(hdu[0].data) > 9000
    
    
def test_centre_chip_with_single_star_all_defaults():
    src = sim.source.star(mag=20, filter_name="Ks")
    hdu = sim.run(src, detector_layout="centre")
    
    assert np.max(hdu[0].data) > 9000
    

def test_full_sim_with_single_star_double_corr_readout():
    cmd = sim.UserCommands()
    opt = sim.OpticalTrain(cmd)
    fpa = sim.Detector(cmd)
    
    src = sim.source.star(mag=20, filter_name="Ks")
    src.apply_optical_train(opt, fpa)
    hdu = fpa.read_out(read_out_type="non_destructive")
    
    # "It’s over nine thousand!" - Vegeta Breigh
    assert np.max(hdu[0].data) > 9000

    
def test_flux_from_non_destructive_is_same_as_superfast():

    cmd = sim.UserCommands()
    cmd["FPA_READ_OUT_SCHEME"] = "double_corr"
    opt = sim.OpticalTrain(cmd)
    fpa = sim.Detector(cmd)
    
    src = sim.source.star(mag=20, filter_name="Ks")
    src.apply_optical_train(opt, fpa)
    hdu_nd = fpa.read_out(read_out_type="non_destructive")
    hdu_sf = fpa.read_out(read_out_type="superfast")
    
    sig_nd = np.max(hdu_nd[0].data) - np.min(hdu_nd[0].data)
    sig_sf = np.max(hdu_sf[0].data) - np.min(hdu_sf[0].data)

    # is the difference between the signal part fo the two read-outs < 10%?   
    assert sig_nd - sig_sf < 0.1 * sig_sf 



    
if __name__ == "__main__":
    test_full_sim_with_single_star_and_all_defaults()
    test_quick_sim_with_single_star_all_defaults()
    test_centre_chip_with_single_star_all_defaults()