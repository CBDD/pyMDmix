from pathlib import Path

from pydantic import BaseModel

DEFAULT_WATER_MODEL = "TIP3P"


class Probe(BaseModel):
    name: str
    mask: str
    types: list[str]


class SolventDefinition(BaseModel):
    name: str
    off_file: Path
    description: str
    box_unit: str
    water_model: str
    probes: list[Probe]
    frcmods: list[str]

    def read_off(self) -> bytes:
        with open(self.off_file, "br") as f:
            return f.read()


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
    data: bytes
    probes: list[Probe]
    box_unit: str
    frcmods: list[str] = []
    water_model: str = DEFAULT_WATER_MODEL
    description: str = ""

    def __str__(self) -> str:
        return f"{self.name}"

    def toLine(self) -> str:
        return f"{self.name}({self.water_model}): {self.description}"

    def toRecord(self) -> str:
        return "\n".join(
            [
                f"{self.name}:",
                f"\t",
            ]
        )

    @staticmethod
    def from_solvent_definition(solvent_definition: SolventDefinition) -> "Solvent":
        return Solvent(
            name=solvent_definition.name,
            data=solvent_definition.read_off(),
            probes=solvent_definition.probes,
            box_unit=solvent_definition.box_unit,
            frcmods=solvent_definition.frcmods,
            water_model=solvent_definition.water_model,
            description=solvent_definition.description,
        )


"""
previously, a config file for a solvent would be an ini file like this:
[GENERAL]
# solvation name (ex: ION)
name = ANT 
info = Acetonitrile 20%% mixture
# path to off file
objectfile = ANTWAT20.off 
# Name of the solvation box  unit in object file(ex: IONWAT20)
boxunit = ANTWAT20
watermodel = TIP3P

[PROBES]
# map probe names with residue@atoms (ie. NEG=COO@O1,O2)
# probe names must be unique
WAT=WAT@O
N=ANT@N1
C=ANT@C3

[TYPES]
WAT=Wat
N=Acc
C=Hyd


now it's a yaml like that:

solvents:
  - name: ANT
    description: Acetonitrile 20% mixture
    objectfile: ANTWAT20.off
    boxunit: ANTWAT20
    watermodel: TIP3P
    probes:
      - name: WAT
        mask: WAT@O
        types: [Wat]
      - name: N
        mask: ANT@N1
        types: [Acc]
      - name: C
        mask: ANT@C3
        types: [Hyd]
    frcmod_paths: []
"""
