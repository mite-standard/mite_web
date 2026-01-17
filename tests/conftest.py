import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--slow", action="store_true", default=False, help="run slow tests"
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: mark test as slow-running")


def pytest_runtest_setup(item):
    if "slow" in item.keywords and not item.config.getoption("--slow"):
        pytest.skip("test requires '--slow' option to run")
