from pathlib import Path
from typing import List, Dict, Optional, Any, Union
from pydantic import BaseModel, Field, field_validator
from enum import Enum
from datetime import datetime


# Enums for validation
class RestraintMode(str, Enum):
    FREE = "FREE"
    HA = "HA"  # Heavy atoms
    BB = "BB"  # Backbone atoms


class EnsembleType(str, Enum):
    NVT = "NVT"
    NPT = "NPT"


class SimulationEngine(str, Enum):
    AMBER = "AMBER"
    NAMD = "NAMD"
    OPENMM = "OPENMM"


class QueueSystem(str, Enum):
    SLURM = "SLURM"
    SGE = "SGE"
    LOCAL = "LOCAL"
    CTE = "CTE"
    DAINT = "DAINT"


# Request Models (for CLI input)
class SystemRequest(BaseModel):
    """Request model for system setup"""
    name: str = Field(..., description="System name")
    structure_file: Path = Field(..., description="Path to OFF or PDB file")
    structure_type: str = Field("off", description="Structure file type (off/pdb)")
    force_fields: List[str] = Field(
        default_factory=lambda: ["leaprc.protein.ff14SB", "leaprc.water.tip3p"],
        description="Force field files"
    )
    extra_files: List[Path] = Field(default_factory=list, description="Additional parameter files")
    box_buffer: float = Field(12.0, description="Solvation buffer in Angstroms")
    
    @field_validator('structure_file')
    @classmethod
    def validate_structure_file(cls, v):
        if not v.exists():
            raise ValueError(f"Structure file not found: {v}")
        return v
    
    @field_validator('structure_type')
    @classmethod
    def validate_structure_type(cls, v):
        if v.lower() not in ['off', 'pdb']:
            raise ValueError("Structure type must be 'off' or 'pdb'")
        return v.lower()


class ReplicaRequest(BaseModel):
    """Request model for replica configuration"""
    count: int = Field(3, ge=1, description="Number of replicas")
    solvents: List[str] = Field(..., description="List of solvents for replicas")
    naming_scheme: str = Field("{solvent}_{index}", description="Replica naming pattern")
    
    @field_validator('solvents')
    @classmethod
    def validate_solvents(cls, v, info):  # Change parameter name
        # Access other field values through info.data
        if info.data and 'count' in info.data:
            count = info.data['count']
            if len(v) != count:
                # If single solvent provided, replicate for all replicas
                if len(v) == 1:
                    return v * count
                raise ValueError(f"Number of solvents ({len(v)}) must match replica count ({count})")
        return v



class SimulationRequest(BaseModel):
    """Request model for simulation parameters"""
    total_time: float = Field(20.0, gt=0, description="Total simulation time in ns")
    timestep: float = Field(2.0, gt=0, description="Integration timestep in fs")
    temperature: float = Field(300.0, gt=0, description="Temperature in K")
    pressure: float = Field(1.0, description="Pressure in bar")
    ensemble: EnsembleType = Field(EnsembleType.NPT, description="Thermodynamic ensemble")
    restraint_mode: RestraintMode = Field(RestraintMode.FREE, description="Restraint mode")
    restraint_force: float = Field(5.0, ge=0, description="Restraint force constant")
    output_frequency: int = Field(5000, gt=0, description="Trajectory output frequency")
    log_frequency: int = Field(1000, gt=0, description="Log output frequency")
    engine: SimulationEngine = Field(SimulationEngine.AMBER, description="MD engine")


class ProjectRequest(BaseModel):
    """Complete project request model"""
    name: str = Field(..., description="Project name")
    description: str = Field("", description="Project description")
    system: SystemRequest
    replicas: ReplicaRequest
    simulation: SimulationRequest
    output_dir: Optional[Path] = Field(None, description="Output directory")
    queue_system: QueueSystem = Field(QueueSystem.LOCAL, description="Queue system")
    

# Response Models (for service output)
class SystemResponse(BaseModel):
    """Response model for system setup results"""
    name: str
    status: str
    files_created: Dict[str, Path]
    warnings: List[str] = Field(default_factory=list)
    errors: List[str] = Field(default_factory=list)


class ReplicaResponse(BaseModel):
    """Response model for replica status"""
    name: str
    solvent: str
    status: str
    directory: Optional[Path] = None
    files: Dict[str, Path] = Field(default_factory=dict)
    error: Optional[str] = None


class ProjectResponse(BaseModel):
    """Response model for complete project setup"""
    project_name: str
    status: str
    created_at: datetime = Field(default_factory=datetime.now)
    project_dir: Path
    system: SystemResponse
    replicas: List[ReplicaResponse]
    total_replicas: int
    successful_replicas: int
    scripts_created: Dict[str, Path] = Field(default_factory=dict)
    
    def get_summary(self) -> Dict[str, Any]:
        """Get project summary"""
        return {
            "project": self.project_name,
            "status": self.status,
            "total_replicas": self.total_replicas,
            "successful": self.successful_replicas,
            "failed": self.total_replicas - self.successful_replicas,
            "directory": str(self.project_dir)
        }


class SolventInfo(BaseModel):
    """Solvent information response"""
    name: str
    display_name: str
    description: str = ""
    density: float
    molecular_weight: float
    water_model: str = "TIP3P"
    ionic: bool = False
    residue_counts: Dict[str, int] = Field(default_factory=dict)
    probes: Dict[str, Dict[str, Any]] = Field(default_factory=dict)


class AnalysisRequest(BaseModel):
    """Request model for analysis parameters"""
    project_dir: Path
    replicas: List[str] = Field(default_factory=list, description="Specific replicas to analyze")
    probes: List[str] = Field(default_factory=list, description="Specific probes to analyze")
    cutoff: float = Field(3.5, description="Probe cutoff distance")
    grid_spacing: float = Field(0.5, description="Grid spacing for density calculation")
    output_format: str = Field("dx", description="Output format for grids")


# Template Variables Models
class TemplateVariables(BaseModel):
    """Base model for template variables"""
    pass


class MinimizationTemplateVars(TemplateVariables):
    """Variables for minimization templates"""
    maxcyc: int = 5000
    ncyc: int = 2500
    ntpr: int = 100
    restraint_wt: float = 10.0
    restraintmask: str = ""


class EquilibrationTemplateVars(TemplateVariables):
    """Variables for equilibration templates"""
    nsteps: int = 25000
    timestep: float = 2.0
    temp0: float = 300.0
    tempi: float = 100.0
    ntpr: int = 1000
    ntwx: int = 5000
    restraint_wt: float = 5.0
    restraintmask: str = ""
    stage: int = 1


class ProductionTemplateVars(TemplateVariables):
    """Variables for production templates"""
    nsteps: int = 500000
    timestep: float = 2.0
    temp0: float = 300.0
    pres0: float = 1.0
    ntpr: int = 1000
    ntwx: int = 5000
    ioutfm: int = 1  # NetCDF format
    iwrap: int = 1  # Wrap coordinates
    restraint_wt: float = 0.0
    restraintmask: str = ""


# Configuration Models (for YAML)
class YAMLConfig(BaseModel):
    """Base configuration model for YAML files"""
    
    class Config:
        extra = "allow"
        use_enum_values = True
    
    def to_yaml(self, path: Path):
        """Save configuration to YAML file"""
        import yaml
        with open(path, 'w') as f:
            yaml.dump(self.dict(), f, default_flow_style=False)
    
    @classmethod
    def from_yaml(cls, path: Path):
        """Load configuration from YAML file"""
        import yaml
        with open(path, 'r') as f:
            data = yaml.safe_load(f)
        return cls(**data)


class ProjectConfig(YAMLConfig):
    """Project configuration for YAML storage"""
    project: ProjectRequest
    created_at: datetime = Field(default_factory=datetime.now)
    pymdmix_version: str = "3.0.0"
