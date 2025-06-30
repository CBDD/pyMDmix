from pathlib import Path
from typing import Any, Collection

import yaml
from pydantic_settings import BaseSettings, PydanticBaseSettingsSource, YamlConfigSettingsSource


class AppConfig(BaseSettings):
    version: str = "1.0.0"
    debug: bool = False
    # List of plugin modules to load
    plugins: list[str] = [
        "mdmix.plugins.solvent"
    ]
    logging_level: str = "WARNING"

    @classmethod
    def settings_customise_sources(
        cls,
        settings_cls: type[BaseSettings],
        init_settings: PydanticBaseSettingsSource,
        env_settings: PydanticBaseSettingsSource,
        dotenv_settings: PydanticBaseSettingsSource,
        file_secret_settings: PydanticBaseSettingsSource,
    ) -> tuple[PydanticBaseSettingsSource, ...]:
        return (YamlConfigSettingsSource(settings_cls),)

    class Config:
        env_prefix = "MDMIX_"  # Optional: can load from environment variables with this prefix
        case_sensitive = False
        use_enum_values = True


def dict_merge(destination: dict[str, Any], source: dict[str, Any]) -> dict[str, Any]:
    for key, value in source.items():
        if isinstance(value, dict):
            destination[key] = dict_merge(destination.get(key, {}), value)
        elif isinstance(value, (list, tuple)):
            destination[key] += value
        elif isinstance(value, set):
            destination[key] |= value
        elif isinstance(value, Collection):
            raise TypeError(f"unhandled collection type: {value.__class__.__name__}")
        else:
            destination[key] = value
    return destination


def get_config_dict(overrides: list[Path]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for filename in overrides:
        with open(filename, "r") as f:
            incoming_data = yaml.load(f, yaml.FullLoader) or {}
        dict_merge(result, incoming_data)
    return result
