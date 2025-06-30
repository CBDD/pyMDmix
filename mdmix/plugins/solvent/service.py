import logging
from typing import Iterable

from .models import CreateSolventRequest, CreateSolventResponse, DeleteSolventRequest, DeleteSolventResponse, Solvent
from .repositories import DEFAULT_REPOSITORY, SolventRepository

logger = logging.getLogger("Solvent")


class SolventService:
    def __init__(self, repository: SolventRepository = DEFAULT_REPOSITORY):
        self.repository = repository

    def list(self) -> Iterable[Solvent]:
        return self.repository.get_all()

    def info(self, name: str) -> str:
        solvent = self.repository.get_by_id(name)
        return solvent.toRecord()

    def create(self, request: CreateSolventRequest) -> CreateSolventResponse:
        results = self.repository.bulk_create(
            request.solvents,
            ignore_errors=request.ignore_errors,
            update_existing=request.update_existing,
        )
        return CreateSolventResponse(ids=results)

    def delete(self, request: DeleteSolventRequest) -> DeleteSolventResponse:
        results = self.repository.bulk_delete(
            request.ids,
            ignore_errors=request.ignore_errors,
            ignore_missing=request.ignore_missing,
        )
        return DeleteSolventResponse(ids=results)


def get_service() -> SolventService:
    return SolventService()
