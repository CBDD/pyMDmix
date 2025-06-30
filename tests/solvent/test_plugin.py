from mdmix.core.plugin import Plugin
from mdmix.plugins import solvent


def test_solvent_plugin_is_a_plugin() -> None:
    assert isinstance(solvent, Plugin)
