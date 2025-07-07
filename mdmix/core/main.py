import typer

from .config import AppConfig
from .plugin import load_plugin


class MDMixApp:
    def display_header(self) -> None:
        print("mdmix 0.3.0 - to be implemented")

    def __init__(self, config: AppConfig = AppConfig()):
        self.config = config
        self.typer_app = typer.Typer(name="mdmix")
        for module in self.config.plugins:
            plugin = load_plugin(module)
            plugin.register(self.typer_app)

    def run(self) -> None:
        self.display_header()
        self.typer_app()


def main() -> None:
    config = AppConfig()
    MDMixApp(config).run()


if __name__ == "__main__":
    main()
