Introduction

MD settings refer to all parameters the program needs to set up a Replica (being a Replica an independent MD set up containing input files - topology and coordinates of the solvated system - and md configuration files - for minimization, equilibration and production steps -. pyMDMix needs a MD Settings Configuration File to create Replicas. This file will contain all parameters related to the simulation set up. Including:

Solvent: Add a truncated octahedron box with the solvent mixture selected to the system to simulate.
Simulation parameters: temperature, length in nanoseconds, trajectory writing frequency, MD simulation program (currently supported simulation programs are AMBER or NAMD.), steps per job, etc...
Simulation restraints: It is possible to automatically set positional restraints over the protein atoms. The residues to restrain are selected automatically or by using masks. Two options are predefined to select atoms in each residue that will be restrained: heavy atoms (HA option - all non-hydrogen atoms) or backbone atoms (BB option). A force value will be needed when restraining positions.
Residue masks: Specify with residues will be selected for applying positional restraints to (RESTRMASK) or residues over which the trajectory will be aligned at analysis stage (ALIGNMASK). If none of these options is specified at replica creation time, the program will try to detect the protein residues automatically.
The full list of options that can be customized is quite large and very often many parameters are rarely modified. For this reason, most of the options are pre-defined with a default value in pyMDMix. The user will not have to worry about many settings. In fact, the simplest and minimmum MD Settings Configuration File (MSC file) that can be written includes only the solvent to use:

[MDSETTINGS]
SOLVENT = ETA
With this input, a Replica of the selected project system will be solvated with an ethanol/water mixture called ETA (name should exist in Solvents Database). All MD parameters will take default values (e.g. 20 ns, 300K, without restraints, etc...). Each extra option given in this configuration file will override the defaults and create a unique set of settings for this replica (see examples below).

Default values can be found and modified inside two different files which are read in a hierarchical manner at replica creation time. Each parameter is accompanied with a short description. Find a copy of the file at the bottom of this page.

Package default settings ($PYTHON_INSTALL_DIR/pyMDMix/data/defaults/md-settings.cfg): These default parameters are valid for all users in the platform. Modifications at this level will alter the defaults for all the program installation.
User default settings ($HOME/.mdmix/md-settings.cfg): A copy of the package default settings file will be saved inside the user's home directory the first time pyMDMix is executed. Uncomment and modify any line you wish to alter. The options here defined will override the package defaults and will be valid only for current user projects.
Any option found in these two files can be included in the MSC file at replica creation time. Any option defined in MSC file has the highest priority and will override any user or package defaults. For instance, md-settings.cfg includes trajfrequency option (more precisely, found with INI syntax: int-trajfrequency). We could increase this value for the replica with ethanol we want to create with this MSC file:

[MDSETTINGS]
solvent = ETA
trajfrequency = 1000
MD Settings Configuration File (MSC file) for single replica definition

At this point, we can think of this file as the settings definition for creating a single Replica (see next how this same file can be used to define multiple replicas at once). The file which should be written by the user and adapted to the simulation requirements should is formed of  a single section named MDSETTINGS followed by a list of option=value pairs.

Example 1: Simple replica setup

We want to simulate the system using ethanol mixture saved  in solvent database with name ETA for 50 nanoseconds at 298K temperature.

[MDSETTINGS]
solvent = ETA
nanos = 50
temp = 298
Example 2: Replica setup with restraints

Same setup as before but now we want to apply positional restraints over all heavy atoms of the protein with a force of 0.5 kcal/mol·A2 to limit protein motion.

[MDSETTINGS]
solvent = ETA
nanos = 50
temp = 298
restr = HA
force = 0.5
Notice that we did not define what residues to restrain. The program will identify all protein residues automatically (at this point, it is important that non standard residues that you wish to be restrained are identified in the System Configuration file).

[MDSETTINGS]
solvent = ETA
nanos = 50
temp = 298
restr = HA
force = 0.5
restrmask = 3-150,155-290
In this last example we want to apply restraints to particular residues.

Example 3: Replica setup with longer simulation jobs

By default, the production stage is planned to include 1ns for each job run. That is, each nanosecond of simulation will be run independently and will output one trajectory file. If we want to include more nanoseconds per job run, we should modify nvt_prod_steps (nstlim in AMBER configuration files or numsteps in NAMD configuration). The default timestep is 2fs and default nvt_prod_steps is 500.000 (500.000 x 2 = 1.000.000 fs = 1 ns). For running 5ns per job run, nvt_prod_steps should be set to 2.500.000. At the same time, we want to ouput trajectory frames every 4ps (= every 2000 steps).

[MDSETTINGS]
solvent = ETA
nanos = 50
temp = 298
restr = HA
force = 0.5
restrmask = 3-150,155-290
nvt_prod_steps = 2500000
trajfrequency = 2000
MD Settings Configuration File (MSC file) for multiple replica definition

It is possible to define multiple settings inside the same MSC file in two different ways:

Denining multiple MDSETTINGS entries
Using multivalue special options
1. Definining multiple MDESTTINGS sections

MSC File template

# pyMDMix MD Settings Configuration File (MSC File)
# All non used options should be commented/removed
[MDSETTINGS]
###########
# GENERAL #
###########
# Comma separated list of solvent box names to be used in the project (e.g. MAM,ETA,WAT). NO DEFAULT. MANDATORY
SOLVENTS =
# Number of replicas for each solvent (read documentation for advanced options). E.g: 3, WAT:1. DEFAULT:3
#NREPL =
# Number of nanoseconds to run (advanced configuration in the documentation). E.g: 20, WAT//10. DEFAULT:20
#NANOS =
# Temperature of each replica (again advanced configuration assigning independent temperatures is possible). E.g.: 300. DEFAULT:300K
#TEMP =
# Restraining scheme (OPTIONS: HA, FREE, BB for heavyatoms, no-restraints, back-bone atoms). DEFAULT: FREE. Also accepts advanced configuration.
#RESTR=
# Restraining force when RESTR!=FREE. In kcal/mol.A^2. Default=0.0. Advanced configuration is possible.
#FORCE=
# Residue mask of residues where restraints should be applied (default: auto detect). E.g.: 1-100
#RESTRMASK = auto
# Residue mask of residues to which backbone atoms we should align the trajectories in analysis process (default: same as restrmask: all protein). E.g.: 1-100
#ALIGNMASK = auto
#####################
# SETTINGS OVERRIDE #
#####################
# All parameters present in replica-settings.cfg can be here overriden.
# For instance, if we want to modify the trajectory writing frequency:
# TRAJFREQUENCY =  
# Modify number of steps per run file (step).
# We might be interested in increase the file size to include several nanosecods
# instead of only 1 (the default; 500000 steps with md_timestep of 2fs)
# NVT_PROD_STEPS = 1000000
# NPT_PROD_STEPS = 1000000
# ETC...
Package/User Default MD settings

pyMDMix is distributed with a file which defines all default MD parameters all projects and users will be using if no other instructions are given. This file is located at the python module installation directory ($PYTHON_INSTALL_DIR/pyMDMix/data/defaults/md-settings.cfg) and will be copied to the user's home directory at first execution as well ($HOME/.mdmix/md-settings.cfg). The file will be similar to this one:

[GENERAL]
## The following options are only ckecked for their type
## use type-name if you want to enforce type conversion on a parameter
## example:
## int-param1 = 10 creates a variable param1 of type int with value 10
## By contrast, param1 = 10 gives a variable param1 of type str with value '10'
## Possible types: int, float, bool, list.
## list type will chop the string by commas. Eg. list-ff=a,b will become a list ff=['a','b']
# Set simulation options
 mdnetcdf = 1 # 1 (write trajectory in nc format) or 0 (write in ascii format)
 restrMode = FREE # Restraining scheme: FREE, HA (heavy atoms) or BB (backbone only)
 float-restrForce = 0.0 # Restraining force if applicable. Default 0 kcal/mol.A^2
 int-nanos = 20 # Production length in nanoseconds. Default: 20ns
 float-temp = 300 # Simulation temperature. Default = 300K
 mdProgram = AMBER # Default simulation program. Options: AMBER or NAMD currently
 int-trajfrequency = 500 # Trajectory writing frequency = 1000 snapshots per nanosecond = int-production_nsteps / int-trajfrequency
 int-minsteps = 5000 # Number of minimization steps to run
 int-heating_steps = 100000 # Heating steps for each file. 100.000 steps = 200ps
 float-parm_heating_tempi = 100 # Start heating at 100 K
 int-npt_eq_steps = 500000 # 1ns equilibration at NPT
 int-nvt_prod_steps = 500000 # 1ns production files = nvt_prod_steps*(md_timestep/10e6)
 float-md_timestep = 2 # 2 fs timestep
 int-namd_heating_steps = 500000 # 1ns equilibration total time to increase temperature from float-heating_tempi to float-temp in NAMD
 list-FF = leaprc.ff99SB, leaprc.gaff # Default forcefield files to load when opening tLeap

# DEFINE REFERENCE STRUCTURE FOR RESTRAINED SIMULATIONS
# If 1: Use output of minimization as reference structure to set positional restraints.
# If 0: (default) use initial input structure (CRD file) as reference state.
 int-minimizationAsRef = 0
As you may probably have realized, some of the options are prefixed with float- int- or list-. These prefix will not be part of the parameter and internally identify the value type. Be careful to keep these prefix if the parameter should be modified.