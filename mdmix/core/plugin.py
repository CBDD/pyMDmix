import importlib
from typing import Protocol, runtime_checkable

import typer


@runtime_checkable
class Plugin(Protocol):
    def register(self, core_app: typer.Typer) -> None: ...


def load_plugin(module: str) -> Plugin:
    plugin = importlib.import_module(module)
    if not isinstance(plugin, Plugin):
        raise TypeError(f"The module {module} does not implement the proper Plugin interface")
    return plugin
