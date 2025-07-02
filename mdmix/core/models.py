from functools import cached_property
from typing import Any

import numpy
from numpy.typing import NDArray
from pydantic import BaseModel


class Atom(BaseModel):
    id: int
    name: str
    type: str
    element: int
    charge: float

    def __str__(self) -> str:
        parts = [f"{key}: {value}" for key, value in self.model_dump(exclude={"charge"}).items()]
        parts.append(f"charge: {self.charge:.4f}")
        return "\n".join(parts)

    def __repr__(self) -> str:
        return self.name

    def __eq__(self, other: Any) -> bool:
        return other.name == self.name if isinstance(other, Atom) else other == self.name


class Residue(BaseModel):
    name: str
    atoms: list[Atom]
    connectivity: list[tuple[int, int]]
    coordinates: list[tuple[float, float, float]]

    @cached_property
    def xyz(self) -> NDArray[numpy.float64]:
        return numpy.array(self.coordinates, dtype=numpy.float64)

    @cached_property
    def charge(self) -> float:
        return round(sum(atom.charge for atom in self.atoms), 4)

    @cached_property
    def atoms_by_id(self) -> dict[int, Atom]:
        return {atom.id: atom for atom in self.atoms}

    @cached_property
    def atoms_by_name(self) -> dict[str, Atom]:
        return {atom.name: atom for atom in self.atoms}

    @cached_property
    def center(self) -> NDArray[numpy.float64]:
        return self.xyz.mean(axis=0)  # type: ignore

    def __str__(self) -> str:
        parts = [f"RESIDUE NAME: {self.name}", "ATOMS:"]
        parts += [f"\t{atom}" for atom in self.atoms]
        return "\n".join(parts)

    def __repr__(self) -> str:
        return self.name

    # originally the __eq__ operator
    def eq(self, other: Any) -> bool:
        return other.name == self.name if isinstance(other, Residue) else other == self.name


class Probe(BaseModel):
    name: str
    residue: Residue
    atoms: list[Atom]
    types: list[str]
    probability: float

    @cached_property
    def mask(self) -> NDArray[numpy.bool]:
        return numpy.array([at.name in self.atoms for at in self.residue.atoms])

    def istype(self, type: str) -> bool:
        return type in self.types

    def __repr__(self) -> str:
        return self.name

    # originally the __eq__ operator
    def eq(self, other: Any) -> bool:
        return other.name == self.name if isinstance(other, Probe) else other == self.name
