import logging

from .exceptions import SolventNotFound
from .models import CreateSolventRequest, DeleteSolventRequest
from .repositories import DEFAULT_REPOSITORY, SolventRepository

logger = logging.getLogger("Solvent")


class SolventService:
    def __init__(self, repository: SolventRepository = DEFAULT_REPOSITORY):
        self.repository = repository

    def list(self):
        for solvent in self.repository.get_all():
            print(solvent.toLine())

    def info(self, name):
        try:
            solvent = self.repository.get_by_id(name)
            print(solvent.toRecord())
        except SolventNotFound as e:
            logger.error(e)

    def create(self, request: CreateSolventRequest):
        results = self.repository.bulk_create(
            request.solvents,
            ignore_errors=request.ignore_errors,
            update_existing=request.update_existing,
        )
        print(f"Created {len(results)} solvents:")
        print("\n\t -".join(results))

    def delete(self, request: DeleteSolventRequest):
        results = self.repository.bulk_delete(
            request.ids,
            ignore_errors=request.ignore_errors,
            ignore_missing=request.ignore_missing,
        )
        print(f"Deleted {len(results)} solvents:")
        print("\n\t -".join(results))


def get_service():
    return SolventService()
