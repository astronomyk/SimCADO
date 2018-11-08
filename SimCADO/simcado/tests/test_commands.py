import pytest
from simcado.commands import UserCommands

def test_load_UserCommands():
    cmd = UserCommands()
    assert type(cmd) == UserCommands

def test_update():
    cmd = UserCommands()
    cmd['OBS_EXPTIME'] = 30
    assert cmd.cmds['OBS_EXPTIME'] == 30
    
def test_wrong_keyword():
    cmd = UserCommands()
    with pytest.raises(KeyError):
        cmd.update({'NO_EXISTE' : 30})
