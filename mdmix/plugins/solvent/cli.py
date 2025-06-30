import logging
from pathlib import Path

import typer
import yaml

from .models import CreateSolventRequest, DeleteSolventRequest
from .service import get_service

logger = logging.getLogger("Solvent")

app = typer.Typer(name="solvent")


@app.command("list")
def list_solvents():
    service = get_service()
    service.list()


@app.command("info")
def info_solvent(name: str):
    service = get_service()
    service.info(name)


@app.command("create")
def create_solvent(config_file: Path, ignore_errors: bool = False, update_existing: bool = False):
    service = get_service()
    with open(config_file, "r") as f:
        data = yaml.load(f, yaml.FullLoader)
    request = CreateSolventRequest(**data, ignore_errors=ignore_errors, update_existing=update_existing)
    service.create(request)


@app.command("delete")
def delete_solvent(names: list[str], ignore_errors: bool = False, ignore_missing: bool = True):
    request = DeleteSolventRequest(ids=names, ignore_errors=ignore_errors, ignore_missing=ignore_missing)
    service = get_service()
    service.delete(request)
