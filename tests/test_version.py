import pymerra2


class TestPyMerra2:

    def test_version_definition(self):
        version = pymerra2.__version__
        assert version >= 0.2

    def test_variables(self):
        assert isinstance(pymerra2.variables.var_list, dict)

    def test_variable_keys(self):
        for i in pymerra2.variables.var_list:
            assert {"esdt_dir", "collection", "merra_name"}.issubset(list(pymerra2.variables.var_list[i].keys()))
