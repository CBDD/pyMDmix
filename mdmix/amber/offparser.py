import tempfile
from collections import defaultdict
from contextlib import contextmanager
from enum import StrEnum
from functools import cached_property
from pathlib import Path
from typing import IO, Generator

from mdmix.core.models import Atom, Residue


class OFFSectionType(StrEnum):
    ATOMS = "atoms"
    ATOMSPERTINFO = "atomspertinfo"
    BOUNDBOX = "boundbox"
    CHILDSEQUENCE = "childsequence"
    CONNECT = "connect"
    CONNECTIVITY = "connectivity"
    HIERARCHY = "hierarchy"
    NAME = "name"
    POSITIONS = "positions"
    RESIDUECONNECT = "residueconnect"
    RESIDUES = "residues"
    RESIDUESPDBSEQUENCENUMBER = "residuesPdbSequenceNumber"
    SOLVENTCAP = "solventcap"
    VELOCITIES = "velocities"


class OFFFile:
    def __init__(self, file: Path | str | IO[str]):
        self.data: str
        if isinstance(file, (str, Path)):
            with open(file, "r") as f:
                self.data = f.read()
        else:
            self.data = file.read()

    @cached_property
    def sections(self) -> dict[tuple[str, OFFSectionType], list[str]]:
        lines = iter(self.data.splitlines())
        matching_section: tuple[str, OFFSectionType] | None = None
        skip_lines = False
        all_sections: dict[tuple[str, OFFSectionType], list[str]] = defaultdict(list)
        expected_sections = []
        while (line := next(lines, None)) is not None:
            line = line.strip()
            if line == "":
                continue
            if not line.startswith("!"):
                if skip_lines:
                    continue
                if matching_section is not None:
                    all_sections[matching_section].append(line)
                else:
                    expected_sections.append(line.replace('"', ""))
                continue
            if line.startswith("!"):
                if line.startswith("!!index"):
                    skip_lines = False
                    matching_section = None
                else:
                    path = line.split(" ")[0].split(".")
                    if len(path) < 4:
                        raise ValueError
                    entry, unit_name, unit_type, section_type, *_ = path
                    if entry != "!entry":
                        raise ValueError
                    if unit_name not in expected_sections:
                        raise ValueError
                    if unit_type != "unit":
                        skip_lines = True
                        matching_section = None
                        continue
                    if section_type not in OFFSectionType:
                        raise ValueError
                    skip_lines = False
                    matching_section = (unit_name, OFFSectionType(section_type))

        return all_sections

    @cached_property
    def units(self) -> dict[str, dict[OFFSectionType, list[str]]]:
        all_units: dict[str, dict[OFFSectionType, list[str]]] = defaultdict(dict)
        for (unit_name, unit_type), data in self.sections.items():
            all_units[unit_name][unit_type] = data
        return all_units

    def get_coords(self, unit_name: str) -> list[tuple[float, float, float]]:
        lines = self.sections[(unit_name, OFFSectionType.POSITIONS)]
        splits = (line.split(None, 3) for line in lines)
        coords = [(float(x), float(y), float(z)) for x, y, z in splits]
        return coords

    def get_connectivity(self, unit_name: str) -> list[tuple[int, int]]:
        lines = self.sections.get((unit_name, OFFSectionType.CONNECTIVITY))
        if lines is None:
            return []
        splits = (line.split() for line in lines)
        connectivity = [(int(atom1), int(atom2)) for atom1, atom2, *_ in splits]
        reversed_pairs = [(atom2, atom1) for atom1, atom2 in connectivity]
        connectivity.extend(reversed_pairs)
        return connectivity

    def get_atoms(self, unit_name: str, skip_hydrogens: bool = False) -> list[Atom]:
        lines = self.sections[(unit_name, OFFSectionType.ATOMS)]
        splits = (line.split() for line in lines)
        return [
            Atom(
                name=name.replace('"', ""),
                id=int(id),
                type=type.replace('"', ""),
                element=int(element),
                charge=float(charge),
            )
            for name, type, _, _, _, id, element, charge in splits
            if not skip_hydrogens or int(element) != 1
        ]

    def get_residue(self, res_name: str, skip_hydrogens: bool = False) -> Residue:
        return Residue(
            name=res_name,
            atoms=self.get_atoms(res_name, skip_hydrogens=skip_hydrogens),
            connectivity=self.get_connectivity(res_name),
            coordinates=self.get_coords(res_name),
        )

    def get_residues(self, unit_name: str, unique: bool = True) -> list[str]:
        lines = self.sections[(unit_name, OFFSectionType.RESIDUES)]
        splits = (line.split() for line in lines)
        residues = (residue.replace('"', "") for residue, *_ in splits)
        return list(set(residues)) if unique else list(residues)

    def get_residue_count(self, unit_name: str, res_name: str) -> int:
        return self.get_residues(unit_name, unique=False).count(res_name)

    def get_atom_count(self, unit_name: str, res_name: str, atom_name: str) -> int:
        # First the number of residues
        residues = self.get_residue_count(unit_name, res_name)
        atoms = self.get_atoms(res_name, skip_hydrogens=False)
        atom_names = [atom.name for atom in atoms if atom.name == atom_name]
        return len(atom_names) * residues

    def get_box_dimensions(self, unit_name: str) -> tuple[float, float, float]:
        lines = self.sections[(unit_name, OFFSectionType.BOUNDBOX)]
        _, _, x, y, z = lines
        return float(x), float(y), float(z)

    def get_box_size(self, unit_name: str) -> float:
        x, y, z = self.get_box_dimensions(unit_name)
        return x * y * z

    def write(self, file: Path | str | IO[str]) -> None:
        if isinstance(file, (Path, str)):
            with open(file, "w") as f:
                f.write(self.data)
        else:
            file.write(self.data)

    @contextmanager
    def tempfile(self) -> Generator[str, None, None]:
        with tempfile.NamedTemporaryFile("w") as f:
            self.write(f.file)
            yield f.name
