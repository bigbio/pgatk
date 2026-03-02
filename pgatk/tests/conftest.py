"""Shared pytest fixtures for pgatk tests."""

from pathlib import Path

import pytest

# Absolute paths to key directories, resolved at import time.
_PACKAGE_DIR = Path(__file__).resolve().parent.parent          # pgatk/pgatk/
TESTDATA_DIR = _PACKAGE_DIR / "testdata"
CONFIG_DIR = _PACKAGE_DIR / "config"


@pytest.fixture
def testdata_dir():
    """Return the absolute path to the testdata directory."""
    return TESTDATA_DIR


@pytest.fixture
def config_dir():
    """Return the absolute path to the config directory."""
    return CONFIG_DIR


@pytest.fixture(autouse=True)
def _chdir_to_package_root(monkeypatch):
    """Set CWD to the pgatk package directory for backward compatibility.

    Existing tests reference files with relative paths like
    ``'testdata/test.vcf'`` and ``'config/ensembl_config.yaml'``.
    This fixture ensures those paths resolve correctly regardless of
    where pytest is invoked from.

    Uses monkeypatch.chdir for guaranteed cleanup even if teardown is interrupted.
    """
    monkeypatch.chdir(_PACKAGE_DIR)
