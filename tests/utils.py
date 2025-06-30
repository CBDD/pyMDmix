from pathlib import Path
from typing import Any

import yaml


def load_yaml(path: Path) -> dict[str, Any]:
    with open(path, "r") as f:
        return yaml.load(f, yaml.FullLoader)
