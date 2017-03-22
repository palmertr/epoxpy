import pytest


class BaseTest:

    @pytest.fixture(autouse=True)
    def initdir(self, tmpdir):
        tmpdir.chdir()

    @pytest.fixture
    def a(self):
        from epoxpy.lib import A
        return A()

    @pytest.fixture
    def b(self):
        from epoxpy.lib import B
        return B()

    @pytest.fixture
    def c(self):
        from epoxpy.lib import C
        return C()
