import simcado as sim
cmd = sim.UserCommands()
opt = sim.OpticalTrain(cmd)

def test_gen_telescope_shake():
    gauss = opt._gen_telescope_shake()
    assert type(gauss == sim.psf.PSF)
    assert gauss.params["pix_res"] == cmd["SCOPE_JITTER_FWHM"]
    
    
def test_gen_master_psf(self, psf_type="Airy"):
    psf = opt._gen_master_psf(psf_type="Airy")