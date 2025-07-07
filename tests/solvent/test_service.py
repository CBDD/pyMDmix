import pytest

from mdmix.plugins.solvent import CreateSolventRequest, SolventException, SolventRepository, SolventService

from .conftest import (
    BAD_CREATE_SOLVENT_REQUEST_WRONG_FILE,
    CREATE_SOLVENT_REQUEST,
)


@pytest.mark.parametrize(
    "create_request",
    [pytest.param(CREATE_SOLVENT_REQUEST)],
    # more to be added, so we can test different situations
)
def test_create_solvent(create_request: CreateSolventRequest, repository: SolventRepository) -> None:
    service = SolventService(repository)
    response = service.create(create_request)
    assert response.ids


@pytest.mark.parametrize(
    "ignore_errors",
    [
        pytest.param(True, id="ignore errors"),
        pytest.param(False, id="raise errors"),
    ],
)
@pytest.mark.parametrize(
    "create_request, exception",
    [
        pytest.param(BAD_CREATE_SOLVENT_REQUEST_WRONG_FILE, SolventException),
        # TODO: add more failure modes
    ],
)
def test_create_solvent_failure(
    ignore_errors: bool,
    create_request: CreateSolventRequest,
    exception: type[Exception],
    repository: SolventRepository,
) -> None:
    service = SolventService(repository)
    create_request.ignore_errors = ignore_errors
    if ignore_errors:
        response = service.create(create_request)
        assert len(response.ids) == 0
    else:
        with pytest.raises(exception):
            service.create(create_request)
