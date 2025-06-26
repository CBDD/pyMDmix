import typer
from pathlib import Path
from typing import Optional, List
import yaml

from cli.base import CLIPlugin, DisplayUtils, console
from core.models import (
    ProjectRequest, SystemRequest, ReplicaRequest, 
    SimulationRequest, RestraintMode, EnsembleType,
    SimulationEngine, QueueSystem
)
from core.services import ProjectSetupService


class ProjectPlugin(CLIPlugin):
    """Project management commands"""
    
    def __init__(self):
        super().__init__("project", "Project setup and management commands")
        self.service = ProjectSetupService()
        
    def register_commands(self):
        """Register project commands"""
        
        @self.app.command("create")
        def create_project(
            config_file: Path = typer.Argument(..., 
                help="YAML configuration file for project",
                exists=True,
                file_okay=True,
                dir_okay=False,
                readable=True
            ),
            output_dir: Optional[Path] = typer.Option(None, "--output", "-o",
                help="Output directory (default: current directory)"
            ),
            dry_run: bool = typer.Option(False, "--dry-run",
                help="Show what would be created without actually creating files"
            )
        ):
            """Create a new PyMDMix project from YAML configuration"""
            
            DisplayUtils.show_header("PyMDMix3 Project Creation", 
                                   f"Configuration: {config_file.name}")
            
            try:
                # Load configuration
                with open(config_file, 'r') as f:
                    config_data = yaml.safe_load(f)
                    
                # Create request from config
                request = self._parse_project_config(config_data)
                
                if output_dir:
                    request.output_dir = output_dir
                    
                # Show project summary
                self._show_project_summary(request)
                
                if dry_run:
                    DisplayUtils.show_info("Dry run mode - no files created")
                    return
                    
                # Confirm creation
                if not typer.confirm("Create project with these settings?"):
                    raise typer.Abort()
                    
                # Create project
                with console.status("[bold green]Creating project..."):
                    response = self.service.setup_project(request)
                    
                # Show results
                self._show_creation_results(response)
                
            except Exception as e:
                DisplayUtils.show_error(f"Project creation failed: {e}")
                raise typer.Exit(1)
                
        @self.app.command("setup")
        def setup_project(
            name: str = typer.Argument(..., help="Project name"),
            structure_file: Path = typer.Argument(..., 
                help="Structure file (OFF or PDB)",
                exists=True
            ),
            solvents: List[str] = typer.Option(..., "--solvent", "-s",
                help="Solvents to use (can specify multiple)"
            ),
            replicas: int = typer.Option(3, "--replicas", "-r",
                help="Number of replicas per solvent"
            ),
            time: float = typer.Option(20.0, "--time", "-t",
                help="Simulation time in nanoseconds"
            ),
            temperature: float = typer.Option(300.0, "--temp",
                help="Temperature in Kelvin"
            ),
            restraints: RestraintMode = typer.Option(RestraintMode.FREE, "--restraints",
                help="Restraint mode"
            ),
            output_dir: Optional[Path] = typer.Option(None, "--output", "-o",
                help="Output directory"
            ),
            force_fields: List[str] = typer.Option(
                ["leaprc.protein.ff14SB", "leaprc.water.tip3p"],
                "--forcefield", "-ff",
                help="Force field files"
            )
        ):
            """Quick project setup with command-line parameters"""
            
            DisplayUtils.show_header("PyMDMix3 Quick Project Setup", name)
            
            # Determine structure type
            structure_type = "off" if structure_file.suffix.lower() == ".off" else "pdb"
            
            # Build request
            request = ProjectRequest(
                name=name,
                description=f"PyMDMix3 project for {structure_file.stem}",
                system=SystemRequest(
                    name=structure_file.stem,
                    structure_file=structure_file,
                    structure_type=structure_type,
                    force_fields=force_fields
                ),
                replicas=ReplicaRequest(
                    count=len(solvents) * replicas,
                    solvents=solvents * replicas  # Replicate solvent list
                ),
                simulation=SimulationRequest(
                    total_time=time,
                    temperature=temperature,
                    restraint_mode=restraints
                ),
                output_dir=output_dir
            )
            
            # Show summary
            self._show_project_summary(request)
            
            # Confirm
            if not typer.confirm("Create project?"):
                raise typer.Abort()
                
            # Create project
            with console.status("[bold green]Creating project..."):
                response = self.service.setup_project(request)
                
            self._show_creation_results(response)
            
        @self.app.command("template")
        def generate_template(
            output_file: Path = typer.Option("project_template.yaml", "--output", "-o",
                help="Output file for template"
            )
        ):
            """Generate a template YAML configuration file"""
            
            template = {
                "project": {
                    "name": "my_project",
                    "description": "Mixed-solvent MD simulation project"
                },
                "system": {
                    "name": "protein_system",
                    "structure_file": "path/to/structure.off",
                    "structure_type": "off",  # or "pdb"
                    "force_fields": [
                        "leaprc.protein.ff14SB",
                        "leaprc.water.tip3p"
                    ],
                    "extra_files": [],  # Additional parameter files
                    "box_buffer": 12.0
                },
                "replicas": {
                    "count": 3,
                    "solvents": ["WAT", "ETA", "WAT"],  # One per replica
                    "naming_scheme": "{solvent}_{index}"
                },
                "simulation": {
                    "total_time": 20.0,  # nanoseconds
                    "timestep": 2.0,  # femtoseconds
                    "temperature": 300.0,  # Kelvin
                    "pressure": 1.0,  # bar
                    "ensemble": "NPT",  # or "NVT"
                    "restraint_mode": "FREE",  # or "HA", "BB"
                    "restraint_force": 5.0,
                    "output_frequency": 5000,
                    "log_frequency": 1000,
                    "engine": "AMBER"  # or "NAMD", "OPENMM"
                },
                "queue_system": "SLURM"  # or "SGE", "LOCAL", etc.
            }
            
            with open(output_file, 'w') as f:
                yaml.dump(template, f, default_flow_style=False, sort_keys=False)
                
            DisplayUtils.show_success(f"Template saved to: {output_file}")
            DisplayUtils.show_info("Edit the template and use 'pymdmix3 project create' to create your project")
            
    def _parse_project_config(self, config: dict) -> ProjectRequest:
        """Parse YAML configuration into ProjectRequest"""
        
        # Extract sections
        project_cfg = config.get('project', {})
        system_cfg = config.get('system', {})
        replica_cfg = config.get('replicas', {})
        sim_cfg = config.get('simulation', {})
        
        # Build request
        return ProjectRequest(
            name=project_cfg.get('name', 'pymdmix_project'),
            description=project_cfg.get('description', ''),
            system=SystemRequest(**system_cfg),
            replicas=ReplicaRequest(**replica_cfg),
            simulation=SimulationRequest(**sim_cfg),
            queue_system=QueueSystem(config.get('queue_system', 'LOCAL'))
        )
        
    def _show_project_summary(self, request: ProjectRequest):
        """Display project configuration summary"""
        
        summary_table = DisplayUtils.create_table("Project Configuration", 
                                                ["Setting", "Value"])
        
        summary_table.add_row("Project Name", request.name)
        summary_table.add_row("System", request.system.name)
        summary_table.add_row("Structure", str(request.system.structure_file))
        summary_table.add_row("Total Replicas", str(request.replicas.count))
        summary_table.add_row("Solvents", ", ".join(set(request.replicas.solvents)))
        summary_table.add_row("Simulation Time", f"{request.simulation.total_time} ns")
        summary_table.add_row("Temperature", f"{request.simulation.temperature} K")
        summary_table.add_row("Ensemble", request.simulation.ensemble.value)
        summary_table.add_row("Restraints", request.simulation.restraint_mode.value)
        
        console.print(summary_table)
        
    def _show_creation_results(self, response):
        """Display project creation results"""
        
        if response.status == "success":
            DisplayUtils.show_success(f"Project created successfully!")
        else:
            DisplayUtils.show_warning(f"Project created with issues")
            
        # Summary
        console.print(f"\nProject Directory: [bold]{response.project_dir}[/bold]")
        console.print(f"Total Replicas: {response.total_replicas}")
        console.print(f"Successful: [green]{response.successful_replicas}[/green]")
        
        if response.successful_replicas < response.total_replicas:
            failed = response.total_replicas - response.successful_replicas
            console.print(f"Failed: [red]{failed}[/red]")
            
        # Show any errors
        if response.system.errors:
            DisplayUtils.show_error("System setup errors:")
            for error in response.system.errors:
                console.print(f"  - {error}")
                
        # Show failed replicas
        failed_replicas = [r for r in response.replicas if r.status == "failed"]
        if failed_replicas:
            DisplayUtils.show_error("Failed replicas:")
            for replica in failed_replicas:
                console.print(f"  - {replica.name}: {replica.error}")
                
        # Next steps
        console.print("\n[bold]Next Steps:[/bold]")
        console.print("1. Change to project directory:")
        console.print(f"   cd {response.project_dir}")
        console.print("2. Run system preparation:")
        console.print("   ./scripts/setup_project.sh")
        console.print("3. Submit simulations to queue or run locally")