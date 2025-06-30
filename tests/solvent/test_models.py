from typing import Any

from mdmix.plugins.solvent.models import CreateSolventRequest


def test_parse_create_solvent_request(create_solv_config_data: dict[str, Any]):
    parsed = CreateSolventRequest(**create_solv_config_data)
    assert len(parsed.solvents) > 0
