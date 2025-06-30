from typing import Iterable, Protocol

import yaml

from .exceptions import SolventException, SolventExists, SolventNotFound
from .models import CreateSolventRequest, Solvent, SolventDefinition


class SolventRepository(Protocol):
    def get_all(self) -> Iterable[Solvent]: ...

    def get_by_id(self, id: str) -> Solvent: ...

    def create(self, request: SolventDefinition, update_existing: bool = False) -> str: ...

    def delete(self, id: str, ignore_missing: bool = True) -> str | None: ...

    def bulk_create(
        self,
        requests: list[SolventDefinition],
        ignore_errors: bool = False,
        update_existing: bool = False,
    ) -> list[str]: ...

    def bulk_delete(self, ids: list[str], ignore_errors: bool = False, ignore_missing: bool = True) -> list[str]: ...


class SolventInMemoryRepository:
    def __init__(self, initial_solvents: Iterable[Solvent] = []) -> None:
        self.data = {solvent.name: solvent for solvent in initial_solvents}

    def get_all(self) -> Iterable[Solvent]:
        return self.data.values()

    def get_by_id(self, id: str) -> Solvent:
        solvent = self.data.get(id)
        if solvent is None:
            raise SolventNotFound(id)
        return solvent

    def create(self, request: SolventDefinition, update_existing: bool = False) -> str:
        if not update_existing and request.name in self.data:
            raise SolventExists(request.name)
        self.data[request.name] = Solvent.from_solvent_definition(request)
        return request.name

    def delete(self, id: str, ignore_missing: bool = True) -> str | None:
        try:
            self.data.pop(id)
            return id
        except KeyError as e:
            if not ignore_missing:
                raise SolventNotFound from e
        return None

    def bulk_create(
        self,
        requests: list[SolventDefinition],
        ignore_errors: bool = False,
        update_existing: bool = False,
    ) -> list[str]:
        results = []
        for request in requests:
            try:
                results.append(self.create(request, update_existing=update_existing))
            except Exception as e:
                if not ignore_errors:
                    raise SolventException from e
        return results

    def bulk_delete(self, ids: list[str], ignore_errors: bool = False, ignore_missing: bool = True) -> list[str]:
        results = []
        for id in ids:
            try:
                results.append(self.delete(id, ignore_missing=ignore_missing))
            except Exception as e:
                if not ignore_errors:
                    raise SolventException from e
        return [r for r in results if r is not None]


def get_repository() -> SolventRepository:
    repository = SolventInMemoryRepository()
    with open("data/solvents/config.yml", "r") as f:
        data = yaml.full_load(f)
    request = CreateSolventRequest(**data)
    repository.bulk_create(request.solvents)
    return repository


DEFAULT_REPOSITORY = get_repository()
