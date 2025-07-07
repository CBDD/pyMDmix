import logging
from pathlib import Path

import typer
import yaml

from .models import CreateSolventRequest, DeleteSolventRequest
from .service import get_service

logger = logging.getLogger("Solvent")

app = typer.Typer(name="solvent")


@app.command("list")
def list_solvents() -> None:
    service = get_service()
    solvents = service.list()
    print("Available solvents:")
    for solvent in solvents:
        print(solvent.toLine())


@app.command("info")
def info_solvent(name: str) -> None:
    service = get_service()
    info = service.info(name)
    print(info)


@app.command("create")
def create_solvent(config_file: Path, ignore_errors: bool = False, update_existing: bool = False) -> None:
    service = get_service()
    with open(config_file, "r") as f:
        data = yaml.load(f, yaml.FullLoader)
    request = CreateSolventRequest(**data, ignore_errors=ignore_errors, update_existing=update_existing)
    response = service.create(request)
    print(f"Created {len(response.ids)} solvents:")
    print("\n\t -".join(response.ids))


@app.command("delete")
def delete_solvent(names: list[str], ignore_errors: bool = False, ignore_missing: bool = True) -> None:
    request = DeleteSolventRequest(ids=names, ignore_errors=ignore_errors, ignore_missing=ignore_missing)
    service = get_service()
    response = service.delete(request)
    print(f"Deleted {len(response.ids)} solvents:")
    print("\n\t -".join(response.ids))
