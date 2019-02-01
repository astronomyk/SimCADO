import pytest

from simcado.commands import user_commands as usr_cmds

# .. todo:: finish the tests for checking instrument and filter


@pytest.fixture(scope="class")
def empty_cmds():
    cmd = usr_cmds.UserCommands()
    return cmd


@pytest.mark.usefixtures("empty_cmds")
class TestUserCommandsInit:
    def test_empty_arguments_produces_object(self, empty_cmds):
        assert isinstance(empty_cmds, usr_cmds.UserCommands)

    def test_str_nones_are_converted_to_nones(self, empty_cmds):
        assert empty_cmds.cmds["OBS_RA"] is None

    def test_str_booleans_are_converted_to_booleans(self, empty_cmds):
        assert empty_cmds.cmds["SIM_LOGGING"] is False

    def test_str_floats_are_converted_to_floats(self, empty_cmds):
        assert empty_cmds.cmds["SIM_SIM_MESSAGE_LEVEL"] == 3.

    def test_multiline_str_accepted_as_argument_for_filename(self):
        file_str = "SIM_LOGGING True \n SIM_SIM_MESSAGE_LEVEL None # comment"
        cmds = usr_cmds.UserCommands(file_str)
        assert isinstance(cmds, usr_cmds.UserCommands)

    def test_multiline_str_foreign_keyword_ignored_but_throw_warning(self):
        file_str = "BOGUS my_bogus.txt"
        cmds = usr_cmds.UserCommands(file_str)
        assert isinstance(cmds, usr_cmds.UserCommands)
        # Can only be fixed when filename is read in by self.update,
        # NOT self.cmds.update
        # assert 0


@pytest.mark.usefixtures("empty_cmds")
class TestUserCommandsGettersAndSetters:
    def test_individual_keywords_can_be_updated_via_call(self, empty_cmds):
        empty_cmds["OBS_RA"] = 266.4167
        assert empty_cmds["OBS_RA"] == 266.4167

    def test_individual_values_are_converted_to_none_bool_float(self, empty_cmds):
        empty_cmds["OBS_RA"] = "none"
        empty_cmds["OBS_DEC"] = "False"
        empty_cmds["OBS_SEEING"] = "2.0"
        assert empty_cmds["OBS_RA"] is None
        assert empty_cmds["OBS_DEC"] is False
        assert empty_cmds["OBS_SEEING"] == 2.

    def test_subcategory_dicts_are_updated_when_called(self, empty_cmds):
        empty_cmds["OBS_SEEING"] = "1.2"
        assert empty_cmds["OBS_SEEING"] == 1.2

    def test_foreign_keywords_returns_warning_and_are_ignored(self, empty_cmds):
        empty_cmds["BOGUS"] = "bogus.txt"
        assert "BOGUS" not in empty_cmds.keys


class TestUserCommandsUpdate:
    def test_dict_updated_with_filename(self):
        pass

    def test_dict_updated_with_string(self):
        pass

    def test_dict_updated_with_new_dict(self):
        pass

    def test_foreign_keywords_returns_warning_but_not_exception(self):
        pass


@pytest.mark.usefixtures("empty_cmds")
class TestDeepcopyObject:
    def test_object_returned_is_not_old_object(self, empty_cmds):
        from copy import deepcopy
        new_cmd = deepcopy(empty_cmds)
        assert new_cmd is not empty_cmds


class TestSetInstrument:
    pass


class TestSetMode:
    pass


class TestSetFilterName:
    pass




