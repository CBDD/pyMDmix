from pathlib import Path
from typing import Any, Callable, Iterable

import pytest

from mdmix.plugins.solvent import CreateSolventRequest, Solvent, SolventInMemoryRepository, SolventRepository
from tests.utils import load_yaml

RepositoryFactory = Callable[[], SolventRepository]


FIXTURES_PATH = Path(__file__).parent / "fixtures"

CREATE_SOLVENT_REQUEST_FILENAME = FIXTURES_PATH / "solv.yml"

CREATE_SOLVENT_REQUEST_DATA = load_yaml(CREATE_SOLVENT_REQUEST_FILENAME)

CREATE_SOLVENT_REQUEST = CreateSolventRequest(**CREATE_SOLVENT_REQUEST_DATA)

BAD_CREATE_SOLVENT_REQUEST_WRONG_FILE = CreateSolventRequest(**CREATE_SOLVENT_REQUEST_DATA)
BAD_CREATE_SOLVENT_REQUEST_WRONG_FILE.solvents[0].off_file = Path("fakefile.off")


@pytest.fixture
def create_solv_config_data() -> dict[str, Any]:
    return CREATE_SOLVENT_REQUEST_DATA


def get_in_memory_repository(data: Iterable[Solvent] = []) -> SolventInMemoryRepository:
    return SolventInMemoryRepository(data)


@pytest.fixture
def in_memory_repository() -> SolventRepository:
    return get_in_memory_repository()


@pytest.fixture(params=[get_in_memory_repository])
def repository(request: pytest.FixtureRequest) -> SolventRepository:
    factory: RepositoryFactory = request.param
    return factory()
