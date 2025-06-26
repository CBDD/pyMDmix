from pathlib import Path
from typing import List, Dict, Optional, Any
import logging
from datetime import datetime
import shutil
import yaml

from core.models import (
    ProjectRequest, SystemRequest, ReplicaRequest, SimulationRequest,
    ProjectResponse, SystemResponse, ReplicaResponse,
    SolventInfo, AnalysisRequest
)

from core.template_engine import ScriptTemplateEngine
from solvents import SolventDatabase


class SystemSetupService:
    """Service for setting up molecular systems"""
    
    def __init__(self, template_engine: ScriptTemplateEngine):
        self.template_engine = template_engine
        self.logger = logging.getLogger("SystemSetupService")
        
    def setup_system(self, request: SystemRequest, output_dir: Path) -> SystemResponse:
        """Setup molecular system and create necessary scripts"""
        
        self.logger.info(f"Setting up system: {request.name}")
        
        # Create system directory
        system_dir = output_dir / "system"
        system_dir.mkdir(parents=True, exist_ok=True)
        
        files_created = {}
        warnings = []
        errors = []
        
        try:
            # Copy structure file
            dest_file = system_dir / f"system.{request.structure_file.suffix[1:]}"
            shutil.copy2(request.structure_file, dest_file)
            files_created['structure'] = dest_file
            
            # Create tLEaP script
            leap_script = self._create_leap_script(request, system_dir)
            files_created['leap_script'] = leap_script
            
            # Create execution script
            exec_script = self._create_execution_script(request, system_dir)
            files_created['execution_script'] = exec_script
            
            # Copy additional files
            for extra_file in request.extra_files:
                dest = system_dir / extra_file.name
                shutil.copy2(extra_file, dest)
                files_created[f'extra_{extra_file.stem}'] = dest
                
            status = "success"
            
        except Exception as e:
            self.logger.error(f"System setup failed: {e}")
            errors.append(str(e))
            status = "failed"
            
        return SystemResponse(
            name=request.name,
            status=status,
            files_created=files_created,
            warnings=warnings,
            errors=errors
        )
        
    def _create_leap_script(self, request: SystemRequest, output_dir: Path) -> Path:
        """Create tLEaP script for system preparation"""
        
        script_lines = [
            f"# PyMDMix3 System Preparation Script",
            f"# System: {request.name}",
            f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            ""
        ]
        
        # Load force fields
        for ff in request.force_fields:
            if ff.startswith('leaprc'):
                script_lines.append(f"source {ff}")
            else:
                script_lines.append(f"loadAmberParams {ff}")
                
        # Load additional parameter files
        for extra_file in request.extra_files:
            if extra_file.suffix == '.frcmod':
                script_lines.append(f"loadAmberParams {extra_file.name}")
            elif extra_file.suffix == '.off':
                script_lines.append(f"loadOff {extra_file.name}")
                
        script_lines.append("")
        
        # Load structure
        if request.structure_type == "off":
            script_lines.extend([
                f"loadOff system.off",
                "mol = first",  # Get first unit
                "check mol"
            ])
        else:  # PDB
            script_lines.extend([
                f"mol = loadPdb system.pdb",
                "check mol"
            ])
            
        script_lines.extend([
            "",
            "# Save dry system",
            f"saveAmberParm mol {request.name}_dry.prmtop {request.name}_dry.inpcrd",
            f"savePdb mol {request.name}_dry.pdb",
            f"saveOff mol {request.name}_dry.off",
            "",
            "quit"
        ])
        
        # Write script
        script_path = output_dir / "prepare_system.leap"
        with open(script_path, 'w') as f:
            f.write('\n'.join(script_lines))
            
        return script_path
        
    def _create_execution_script(self, request: SystemRequest, output_dir: Path) -> Path:
        """Create bash execution script"""
        
        script_content = f"""#!/bin/bash
# PyMDMix3 System Preparation
# System: {request.name}

echo "=== PyMDMix3 System Preparation ==="
echo "System: {request.name}"
echo "Date: $(date)"

cd "$(dirname "$0")"

if ! command -v tleap &> /dev/null; then
    echo "ERROR: tLEaP not found. Please load AMBER."
    exit 1
fi

echo "Running tLEaP..."
tleap -f prepare_system.leap

if [ -f "{request.name}_dry.prmtop" ]; then
    echo "SUCCESS: System prepared"
    ls -la *.prmtop *.inpcrd *.pdb *.off 2>/dev/null
else
    echo "ERROR: System preparation failed. Check leap.log"
    exit 1
fi
"""
        
        script_path = output_dir / "run_preparation.sh"
        with open(script_path, 'w') as f:
            f.write(script_content)
        script_path.chmod(0o755)
        
        return script_path


class ReplicaSetupService:
    """Service for setting up simulation replicas"""
    
    def __init__(self, template_engine: ScriptTemplateEngine, 
                 solvent_db: SolventDatabase):
        self.template_engine = template_engine
        self.solvent_db = solvent_db
        self.logger = logging.getLogger("ReplicaSetupService")
        
    def setup_replicas(self, request: ReplicaRequest, system_name: str,
                      simulation_request: SimulationRequest,
                      output_dir: Path) -> List[ReplicaResponse]:
        """Setup all simulation replicas"""
        
        replicas_dir = output_dir / "replicas"
        replicas_dir.mkdir(parents=True, exist_ok=True)
        
        responses = []
        
        for i in range(request.count):
            solvent_name = request.solvents[i]
            replica_name = request.naming_scheme.format(
                solvent=solvent_name, 
                index=i+1
            )
            
            response = self._setup_single_replica(
                replica_name, solvent_name, system_name,
                simulation_request, replicas_dir
            )
            responses.append(response)
            
        return responses
        
    def _setup_single_replica(self, name: str, solvent_name: str, 
                            system_name: str, sim_request: SimulationRequest,
                            replicas_dir: Path) -> ReplicaResponse:
        """Setup a single replica"""
        
        self.logger.info(f"Setting up replica: {name}")
        
        replica_dir = replicas_dir / name
        files = {}
        
        try:
            # Create directory structure
            for subdir in ["system", "minimize", "equilibrate", "production", "logs"]:
                (replica_dir / subdir).mkdir(parents=True, exist_ok=True)
                
            # Get solvent info
            solvent = self.solvent_db.get_solvent(solvent_name)
            if not solvent:
                raise ValueError(f"Solvent not found: {solvent_name}")
                
            # Create solvation script
            solv_script = self._create_solvation_script(
                solvent, system_name, replica_dir / "system"
            )
            files['solvation_script'] = solv_script
            
            # Create MD input files
            md_files = self._create_md_inputs(
                sim_request, replica_dir
            )
            files.update(md_files)
            
            # Create run scripts
            run_scripts = self._create_run_scripts(
                name, replica_dir
            )
            files.update(run_scripts)
            
            status = "success"
            error = None
            
        except Exception as e:
            self.logger.error(f"Replica setup failed: {e}")
            status = "failed"
            error = str(e)
            
        return ReplicaResponse(
            name=name,
            solvent=solvent_name,
            status=status,
            directory=replica_dir,
            files=files,
            error=error
        )
        
    def _create_solvation_script(self, solvent: SolventInfo, 
                                system_name: str, output_dir: Path) -> Path:
        """Create solvation tLEaP script"""
        
        script_lines = [
            f"# PyMDMix3 Solvation Script",
            f"# Solvent: {solvent.display_name}",
            f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            "",
            "# Load force fields",
            "source leaprc.protein.ff14SB",
            "source leaprc.water.tip3p",
            ""
        ]
        
        # Load water model if different
        if solvent.water_model and solvent.water_model != "TIP3P":
            script_lines.append(f"source leaprc.water.{solvent.water_model.lower()}")
            
        # Load solvent parameters (simplified - in real implementation would handle OFF files)
        script_lines.extend([
            "",
            "# Load system",
            f"mol = loadOff ../../system/{system_name}_dry.off",
            "",
            "# Solvate",
            f"solvatebox mol {solvent.water_model}BOX 12.0",
            "",
            "# Add ions to neutralize",
            "addions mol Na+ 0",
            "addions mol Cl- 0",
            "",
            "# Save solvated system",
            "saveAmberParm mol solvated.prmtop solvated.inpcrd",
            "savePdb mol solvated.pdb",
            "",
            "quit"
        ])
        
        script_path = output_dir / "solvate.leap"
        with open(script_path, 'w') as f:
            f.write('\n'.join(script_lines))
            
        return script_path
        
    def _create_md_inputs(self, sim_request: SimulationRequest, 
                         replica_dir: Path) -> Dict[str, Path]:
        """Create MD input files"""
        
        files = {}
        
        # Minimization
        min_vars = {
            'minsteps': 5000,
            'ntpr': 100,
            'restraintmask': '@CA,C,N' if sim_request.restraint_mode != 'FREE' else '',
            'restraint_wt': 10.0
        }
        
        min_content = self.template_engine.create_amber_minimization(min_vars)
        if min_content:
            min_file = replica_dir / "minimize" / "min.in"
            min_file.write_text(min_content)
            files['minimize_input'] = min_file
            
        # Equilibration (5 stages)
        for stage in range(1, 6):
            eq_vars = {
                'final_temp': sim_request.temperature,
                'restraintmask': '@CA,C,N' if sim_request.restraint_mode != 'FREE' else '',
                'restraint_wt': 10.0 / stage,  # Reduce restraints
                'nsteps': 25000 if stage < 5 else 50000
            }
            
            eq_content = self.template_engine.create_amber_equilibration(stage, eq_vars)
            if eq_content:
                eq_file = replica_dir / "equilibrate" / f"eq{stage}.in"
                eq_file.write_text(eq_content)
                files[f'equilibrate_{stage}_input'] = eq_file
                
        # Production
        prod_vars = {
            'nsteps': int(sim_request.total_time * 500000),  # 500k steps = 1ns
            'temp0': sim_request.temperature,
            'pres0': sim_request.pressure,
            'timestep': sim_request.timestep / 1000,  # fs to ps
            'freq': sim_request.output_frequency,
            'ntpr': sim_request.log_frequency,
            'restraintmask': self._get_restraint_mask(sim_request.restraint_mode),
            'restraint_wt': sim_request.restraint_force if sim_request.restraint_mode != 'FREE' else 0
        }
        
        prod_content = self.template_engine.create_amber_production(
            sim_request.ensemble.value, prod_vars
        )
        if prod_content:
            prod_file = replica_dir / "production" / "prod.in"
            prod_file.write_text(prod_content)
            files['production_input'] = prod_file
            
        return files
        
    def _get_restraint_mask(self, mode: str) -> str:
        """Get restraint mask based on mode"""
        if mode == 'HA':
            return '!@H='
        elif mode == 'BB':
            return '@CA,C,N'
        else:
            return ''
            
    def _create_run_scripts(self, replica_name: str, replica_dir: Path) -> Dict[str, Path]:
        """Create execution scripts for the replica"""
        
        files = {}
        
        # Master run script
        master_script = f"""#!/bin/bash
# PyMDMix3 Replica Execution Script
# Replica: {replica_name}

echo "=== PyMDMix3 Replica: {replica_name} ==="
echo "Date: $(date)"

cd "$(dirname "$0")"

# Step 1: Solvation
echo "Step 1: Solvating system..."
cd system
tleap -f solvate.leap
if [ ! -f "solvated.prmtop" ]; then
    echo "ERROR: Solvation failed"
    exit 1
fi
cd ..

# Step 2: Minimization
echo "Step 2: Running minimization..."
cd minimize
pmemd.cuda -O -i min.in -p ../system/solvated.prmtop -c ../system/solvated.inpcrd \\
    -o min.out -r min.rst7 -inf min.mdinfo
cd ..

# Step 3: Equilibration
echo "Step 3: Running equilibration..."
cd equilibrate
for i in {{1..5}}; do
    echo "  Stage $i..."
    if [ $i -eq 1 ]; then
        pmemd.cuda -O -i eq${{i}}.in -p ../system/solvated.prmtop \\
            -c ../minimize/min.rst7 -o eq${{i}}.out -r eq${{i}}.rst7 \\
            -x eq${{i}}.nc -inf eq${{i}}.mdinfo
    else
        prev=$((i-1))
        pmemd.cuda -O -i eq${{i}}.in -p ../system/solvated.prmtop \\
            -c eq${{prev}}.rst7 -o eq${{i}}.out -r eq${{i}}.rst7 \\
            -x eq${{i}}.nc -inf eq${{i}}.mdinfo
    fi
done
cd ..

# Step 4: Production
echo "Step 4: Running production MD..."
cd production
pmemd.cuda -O -i prod.in -p ../system/solvated.prmtop \\
    -c ../equilibrate/eq5.rst7 -o prod.out -r prod.rst7 \\
    -x prod.nc -inf prod.mdinfo
cd ..

echo "Replica {replica_name} completed at: $(date)"
"""
        
        script_path = replica_dir / "run_replica.sh"
        with open(script_path, 'w') as f:
            f.write(master_script)
        script_path.chmod(0o755)
        files['run_script'] = script_path
        
        return files


class ProjectSetupService:
    """Main service for project setup orchestration"""
    
    def __init__(self):
        self.template_engine = ScriptTemplateEngine()
        self.solvent_db = SolventDatabase()
        self.system_service = SystemSetupService(self.template_engine)
        self.replica_service = ReplicaSetupService(self.template_engine, self.solvent_db)
        self.logger = logging.getLogger("ProjectSetupService")
        
    def setup_project(self, request: ProjectRequest) -> ProjectResponse:
        """Setup complete PyMDMix project"""
        
        self.logger.info(f"Setting up project: {request.name}")
        
        # Determine output directory
        if request.output_dir:
            project_dir = request.output_dir / request.name
        else:
            project_dir = Path.cwd() / request.name
            
        project_dir.mkdir(parents=True, exist_ok=True)
        
        # Setup system
        system_response = self.system_service.setup_system(
            request.system, project_dir
        )
        
        # Setup replicas
        replica_responses = self.replica_service.setup_replicas(
            request.replicas, request.system.name,
            request.simulation, project_dir
        )
        
        # Create project-level scripts
        scripts = self._create_project_scripts(request, project_dir)
        
        # Save project configuration
        self._save_project_config(request, project_dir)
        
        # Determine overall status
        successful_replicas = sum(1 for r in replica_responses if r.status == "success")
        status = "success" if successful_replicas == len(replica_responses) else "partial"
        
        return ProjectResponse(
            project_name=request.name,
            status=status,
            project_dir=project_dir,
            system=system_response,
            replicas=replica_responses,
            total_replicas=len(replica_responses),
            successful_replicas=successful_replicas,
            scripts_created=scripts
        )
        
    def _create_project_scripts(self, request: ProjectRequest, 
                               project_dir: Path) -> Dict[str, Path]:
        """Create project-level scripts"""
        
        scripts_dir = project_dir / "scripts"
        scripts_dir.mkdir(exist_ok=True)
        
        scripts = {}
        
        # Master setup script
        setup_script = f"""#!/bin/bash
# PyMDMix3 Project Setup
# Project: {request.name}

echo "=== PyMDMix3 Project Setup ==="
echo "Project: {request.name}"
echo "Replicas: {request.replicas.count}"
echo ""

cd "$(dirname "$0")/.."

# System preparation
echo "Preparing system..."
cd system
./run_preparation.sh
cd ..

# Setup all replicas
echo "Setting up replicas..."
for replica in replicas/*/; do
    if [ -d "$replica" ]; then
        echo "Processing $(basename "$replica")..."
        cd "$replica/system"
        tleap -f solvate.leap
        cd ../../..
    fi
done

echo "Setup completed!"
"""
        
        setup_path = scripts_dir / "setup_project.sh"
        with open(setup_path, 'w') as f:
            f.write(setup_script)
        setup_path.chmod(0o755)
        scripts['setup'] = setup_path
        
        return scripts
        
    def _save_project_config(self, request: ProjectRequest, project_dir: Path):
        """Save project configuration to YAML"""
        
        config_file = project_dir / "project_config.yaml"
        
        config_data = {
            'project': request.dict(),
            'created_at': datetime.now().isoformat(),
            'pymdmix_version': '3.0.0'
        }
        
        with open(config_file, 'w') as f:
            yaml.dump(config_data, f, default_flow_style=False)


class SolventManagementService:
    """Service for solvent database management"""
    
    def __init__(self):
        self.solvent_db = SolventDatabase()
        self.logger = logging.getLogger("SolventManagementService")
        
    def list_solvents(self) -> List[SolventInfo]:
        """List all available solvents"""
        solvent_names = self.solvent_db.list_solvents()
        solvents = []
        
        for name in solvent_names:
            solvent = self.solvent_db.get_solvent(name)
            if solvent:
                solvents.append(SolventInfo(
                    name=solvent.name,
                    display_name=solvent.display_name,
                    description=solvent.description,
                    density=solvent.density,
                    molecular_weight=solvent.molecular_weight,
                    water_model=solvent.water_model,
                    ionic=solvent.ionic,
                    residue_counts=solvent.residue_counts,
                    probes=solvent.probes
                ))
                
        return solvents
        
    def get_solvent_info(self, name: str) -> Optional[SolventInfo]:
        """Get detailed information about a solvent"""
        solvent = self.solvent_db.get_solvent(name)
        
        if solvent:
            return SolventInfo(
                name=solvent.name,
                display_name=solvent.display_name,
                description=solvent.description,
                density=solvent.density,
                molecular_weight=solvent.molecular_weight,
                water_model=solvent.water_model,
                ionic=solvent.ionic,
                residue_counts=solvent.residue_counts,
                probes=solvent.probes
            )
            
        return None
        
    def add_solvent_from_config(self, config_file: Path) -> bool:
        """Add solvent from configuration file"""
        try:
            return self.solvent_db.load_from_config_file(config_file)
        except Exception as e:
            self.logger.error(f"Failed to add solvent: {e}")
            return False