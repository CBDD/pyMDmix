import typer

from .cli import app as solvent_app
from .exceptions import SolventException, SolventExists, SolventNotFound
from .models import CreateSolventRequest, Solvent, SolventDefinition
from .repositories import DEFAULT_REPOSITORY as SOLVENT_DEFAULT_REPOSITORY
from .repositories import SolventInMemoryRepository, SolventRepository, get_repository


def register(core_app: typer.Typer):
    core_app.add_typer(solvent_app, name="solvent")


__all__ = [
    "register",
    "SolventException",
    "SolventExists",
    "SolventNotFound",
    "CreateSolventRequest",
    "Solvent",
    "SolventDefinition",
    "SOLVENT_DEFAULT_REPOSITORY",
    "SolventInMemoryRepository",
    "SolventRepository",
    "get_repository",
]
