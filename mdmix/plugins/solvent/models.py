from collections import defaultdict
from functools import cached_property
from io import StringIO
from pathlib import Path
from typing import Self

import yaml
from pydantic import BaseModel, model_validator

from mdmix.amber.offparser import OFFFile
from mdmix.core.models import Residue

DEFAULT_WATER_MODEL = "TIP3P"


class Probe(BaseModel):
    name: str
    mask: str
    types: list[str]

    @cached_property
    def res_atom_mapping(self) -> dict[str, list[str]]:
        "Split string of format RESNAME@ATOMNAME,ATOMNAME;RESNAME@ATOMNAMES and return a dictionary {RES:[ATOMNAME, ATOMNAME..], RES:[ATOMNAMES]}"
        raw = self.mask[1:] if self.mask.startswith(":") else self.mask
        out: dict[str, list[str]] = defaultdict(list)
        for entry in raw.split(";"):
            parts = entry.split("@")
            out[parts[0]] += [atom.strip() for atom in parts[1].split(",")] if len(parts) == 2 else ["all"]
        return out


class SolventDefinition(BaseModel):
    name: str
    off_file: Path
    description: str
    box_unit: str
    water_model: str
    probes: list[Probe]
    frcmods: list[str]


class CreateSolventRequest(BaseModel):
    solvents: list[SolventDefinition]
    ignore_errors: bool = False
    update_existing: bool = False


class DeleteSolventRequest(BaseModel):
    ids: list[str]
    ignore_errors: bool = False
    ignore_missing: bool = True


class CreateSolventResponse(BaseModel):
    ids: list[str]


class DeleteSolventResponse(BaseModel):
    ids: list[str]


class Solvent(BaseModel):
    name: str
    data: str
    probes: dict[str, Probe]
    box_unit: str
    # TODO: check how (or if) this is being used in the original mdmix implementation
    frcmods: dict[str, str] = {}
    water_model: str = DEFAULT_WATER_MODEL
    description: str = ""

    def __str__(self) -> str:
        return f"{self.name}"

    def toLine(self) -> str:
        return f"{self.name}({self.water_model}): {self.description}"

    def toRecord(self) -> str:
        data = {self.name: self.model_dump(exclude={"name", "data"})}
        return yaml.dump(data)

    @cached_property
    def off_file(self) -> OFFFile:
        return OFFFile(StringIO(self.data))

    @cached_property
    def residues(self) -> dict[str, Residue]:
        return {res_name: self.off_file.get_residue(res_name) for res_name in self.off_file.get_residues(self.box_unit)}

    # TODO
    # @cached_property
    # def com_probes(self) -> dict[str, ???]:
    #     pass

    @cached_property
    def types(self) -> set[str]:
        return {probe_type for probe in self.probes.values() for probe_type in probe.types}

    @cached_property
    def volume(self) -> float:
        return self.off_file.get_box_size(self.box_unit)

    @model_validator(mode="after")
    def validate_solvent(self) -> Self:
        assert self.box_unit in self.off_file.units, f"Box unit {self.box_unit} not present in the off file"
        missing_residues = set(self.residues) - set(self.off_file.units)
        assert len(missing_residues) == 0, f"Residue definitions not present in the off file for {missing_residues}"
        assert self.volume > 0
        return self


def SolventFactory(definition: SolventDefinition) -> Solvent:
    with open(definition.off_file, "r") as f:
        data = f.read()

    frcmods: dict[str, str] = {}
    for filename in definition.frcmods:
        with open(filename, "r") as f:
            frcmods[filename] = f.read()

    return Solvent(
        name=definition.name,
        data=data,
        probes={probe.name: probe for probe in definition.probes},
        box_unit=definition.box_unit,
        frcmods=frcmods,
        water_model=definition.water_model,
        description=definition.description,
    )

