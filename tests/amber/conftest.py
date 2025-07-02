import os
from pathlib import Path

import pytest
from pytest import FixtureRequest

from mdmix.core.models import Residue
from tests.utils import load_json


@pytest.fixture
def amber_test_path(request: FixtureRequest) -> Path:
    return Path(os.path.dirname(request.module.__file__))


@pytest.fixture
def amber_fixtures_path(amber_test_path: Path) -> Path:
    return amber_test_path / "fixtures"


@pytest.fixture
def off_file_path(amber_fixtures_path: Path) -> Path:
    return amber_fixtures_path / "ANTWAT20.off"


def get_residue_from_fixture(filename: Path) -> Residue:
    current_dir = Path(__file__).parent / "fixtures"
    data = load_json(current_dir / filename)
    return Residue(**data)


ANT_RESIDUE = get_residue_from_fixture(Path("ant_residue.json"))
ANT_RESIDUE_NOH = get_residue_from_fixture(Path("ant_residue_noh.json"))
