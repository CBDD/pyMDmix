import itertools
from pathlib import Path

import pytest

from mdmix.amber.offparser import OFFFile, OFFSectionType
from mdmix.core.models import Residue

from .conftest import ANT_RESIDUE, ANT_RESIDUE_NOH


@pytest.fixture
def off_file(off_file_path: Path) -> OFFFile:
    off_file = OFFFile(off_file_path)
    return off_file


def test_parse_off_file(off_file: OFFFile) -> None:
    EXPECTED_UNITS = ["ANT", "ANTWAT20", "WAT"]
    EXPECTED_SECTIONS = itertools.product(EXPECTED_UNITS, OFFSectionType)

    assert off_file.sections
    for section in EXPECTED_SECTIONS:
        assert section in off_file.sections

    assert all(unit in off_file.units for unit in EXPECTED_UNITS)


@pytest.mark.parametrize(
    "res_name, expected_residue, skip_hydrogens",
    [
        pytest.param("ANT", ANT_RESIDUE, False, id="ANT"),
        pytest.param("ANT", ANT_RESIDUE_NOH, True, id="ANT(no hydrogens)"),
    ],
)
def test_get_residue(res_name: str, expected_residue: Residue, skip_hydrogens: bool, off_file: OFFFile) -> None:
    # indirectly testing also get_coords, get_connectivity and get_atoms
    # TODO: implement dedicated tests later
    residue = off_file.get_residue(res_name, skip_hydrogens=skip_hydrogens)
    assert residue == expected_residue


@pytest.mark.parametrize(
    "unit_name, expected_res_names, unique",
    [
        pytest.param("ANT", ["ANT"], True, id="ANT"),
        pytest.param("ANTWAT20", ["ANT", "WAT"], True, id="ANTWAT20"),
        pytest.param("WAT", ["WAT"], True, id="WAT"),
        pytest.param("ANTWAT20", ["ANT"] * 17 + ["WAT"] * 209, False, id="ANTWAT20 - no unique"),
    ],
)
def test_get_residues(unit_name: str, expected_res_names: list[str], unique: bool, off_file: OFFFile) -> None:
    residues = off_file.get_residues(unit_name, unique)
    assert sorted(residues) == sorted(expected_res_names)


@pytest.mark.parametrize(
    "unit_name, res_name, expected_value",
    [
        pytest.param("ANT", "ANT", 1, id="single present"),
        pytest.param("ANT", "WAT", 0, id="not present"),
        pytest.param("ANTWAT20", "ANT", 17, id="multiple ANT"),
        pytest.param("ANTWAT20", "WAT", 209, id="multiple WAT"),
    ],
)
def test_get_residue_count(unit_name: str, res_name: str, expected_value: int, off_file: OFFFile) -> None:
    count = off_file.get_residue_count(unit_name, res_name)
    assert count == expected_value


@pytest.mark.parametrize(
    "unit_name, res_name, atom_name, expected_value",
    [
        pytest.param("ANT", "ANT", "N1", 1, id="single present"),
        pytest.param("ANT", "WAT", "H1", 0, id="residue not present"),
        pytest.param("ANT", "ANT", "F1", 0, id="atom not present"),
        pytest.param("ANTWAT20", "ANT", "N1", 17, id="multiple ANT"),
        pytest.param("ANTWAT20", "WAT", "H2", 209, id="multiple WAT"),
    ],
)
def test_get_atom_count(unit_name: str, res_name: str, atom_name: str, expected_value: int, off_file: OFFFile) -> None:
    count = off_file.get_atom_count(unit_name, res_name, atom_name)
    assert count == expected_value


ANTWAT20_DIMENSIONS = (19.861128, 19.861128, 19.861128)


def test_get_box_dimensions(off_file: OFFFile) -> None:
    # not the greatest way of doing this, but it works
    assert off_file.get_box_dimensions("ANTWAT20") == ANTWAT20_DIMENSIONS


def test_get_box_size(off_file: OFFFile) -> None:
    x, y, z = ANTWAT20_DIMENSIONS
    assert off_file.get_box_size("ANTWAT20") == x * y * z


def test_write(off_file: OFFFile) -> None:
    # indirectly testing write()
    with off_file.tempfile() as filename:
        with open(filename, "r") as f:
            assert f.read() == off_file.data
    assert not Path(filename).exists()
