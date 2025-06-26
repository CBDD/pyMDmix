
import typer
from typing import Optional, List
from pathlib import Path
from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich import print as rprint
import logging
from abc import ABC, abstractmethod

# Initialize console for rich output
console = Console()

# Base CLI app
app = typer.Typer(
    name="pymdmix3",
    help="PyMDMix3 - Advanced Molecular Dynamics Mixed Solvent Simulations",
    rich_markup_mode="rich"
)


class CLIPlugin(ABC):
    """Base class for CLI plugins"""
    
    def __init__(self, name: str, help_text: str):
        self.name = name
        self.help_text = help_text
        self.app = typer.Typer(help=help_text)
        
    @abstractmethod
    def register_commands(self):
        """Register plugin commands"""
        pass
        
    def attach_to_main(self, main_app: typer.Typer):
        """Attach plugin to main app"""
        self.register_commands()
        main_app.add_typer(self.app, name=self.name)


class DisplayUtils:
    """Utilities for rich console output"""
    
    @staticmethod
    def show_header(title: str, subtitle: str = ""):
        """Display a styled header"""
        if subtitle:
            content = f"[bold cyan]{title}[/bold cyan]\n[dim]{subtitle}[/dim]"
        else:
            content = f"[bold cyan]{title}[/bold cyan]"
        console.print(Panel(content, expand=False))
        
    @staticmethod
    def show_error(message: str):
        """Display error message"""
        console.print(f"[bold red]Error:[/bold red] {message}")
        
    @staticmethod
    def show_success(message: str):
        """Display success message"""
        console.print(f"[bold green]Success:[/bold green] {message}")
        
    @staticmethod
    def show_warning(message: str):
        """Display warning message"""
        console.print(f"[bold yellow]Warning:[/bold yellow] {message}")
        
    @staticmethod
    def show_info(message: str):
        """Display info message"""
        console.print(f"[bold blue]Info:[/bold blue] {message}")
        
    @staticmethod
    def create_table(title: str, columns: List[str]) -> Table:
        """Create a styled table"""
        table = Table(title=title, show_header=True, header_style="bold magenta")
        for col in columns:
            table.add_column(col)
        return table


# Logging configuration
def setup_logging(verbose: bool = False):
    """Setup logging configuration"""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


# Main app callbacks
@app.callback()
def main_callback(
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Enable verbose output"),
    version: bool = typer.Option(False, "--version", help="Show version")
):
    """
    PyMDMix3 - Advanced Molecular Dynamics Mixed Solvent Simulations
    
    A modern Python implementation for mixed-solvent MD simulations.
    """
    if version:
        from .. import __version__
        console.print(f"PyMDMix3 version {__version__}")
        raise typer.Exit()
        
    setup_logging(verbose)


# Main command
@app.command()
def info():
    """Show PyMDMix3 information and status"""
    from .. import __version__, __description__
    
    DisplayUtils.show_header("PyMDMix3", __description__)
    
    info_table = Table(show_header=False, box=None)
    info_table.add_column("Property", style="cyan")
    info_table.add_column("Value", style="white")
    
    info_table.add_row("Version", __version__)
    info_table.add_row("Python", "3.8+")
    info_table.add_row("License", "GPL-3.0")
    
    console.print(info_table)
    console.print("\nUse [bold]pymdmix3 --help[/bold] to see available commands.")
