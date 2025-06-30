import os
from pathlib import Path
from typing import Any

import pytest
from pytest import FixtureRequest

from tests.utils import load_yaml


@pytest.fixture
def solv_test_path(request: FixtureRequest) -> Path:
    return Path(os.path.dirname(request.module.__file__))


@pytest.fixture
def solv_fixtures_path(solv_test_path: Path) -> Path:
    return solv_test_path / "fixtures"


@pytest.fixture
def create_solv_config_filename(solv_fixtures_path: Path) -> Path:
    return solv_fixtures_path / "solv.yml"


@pytest.fixture
def create_solv_config_data(create_solv_config_filename: Path) -> dict[str, Any]:
    return load_yaml(create_solv_config_filename)
