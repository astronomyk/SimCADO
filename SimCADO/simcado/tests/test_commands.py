def test_load_UserCommands():
    import simcado as sim
    cmd = sim.UserCommands()
    assert type(cmd) == sim.commands.UserCommands
    
def test_all_defaults_paths_exist():
    import simcado as sim