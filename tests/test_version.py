import pymerra2 as pm2


def test_version_definition():
    version = pm2.__version__
    assert version >= 0.2

def test_modules():
    pass


def test_variables():
    assert type(pm2.variables.var_list) == dict


if __name__ == '__main__':
    pass
