A Replica is considered a single independent MD setup, with the input files for running the simulations, topologies and coordinates for an already solvated system (system for which the solvent mixture has been added). The user can move this replica folder to the computation facility and return it to the project folder whenever the simulation is finished. The replica will identify if the simulation is completed and will orchestrate analysis. (Visit Getting started pictogram for a visual summary of what a Replica is)

Replica Settings configuration files

Once a project is started and a system has been added, to create replicas (independent MDs runs), a Settings configuration file should be written specifying the solvents you wish to simulate and the number of independent replicas to create.  As described in Getting Started, replicas of a system are created from a particular set of MD settings, the most important being the solvents and number of replicas to run for each solvent. This file will contain definitions for the replicas to be create as well as parameters for the MD (see MD Settings).

It consists of a unique [MDSETTINGS] section with a single mandatory entry: SOLVENTS. The rest of options, if not given, are taken from the defaults.

Settings configuration file template:

# pyMDMix MD Settings Configuration File
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
For instance, if we want to create 3 Replicas of a particular system using ETA mixture and 3 replicas with pure water as solvent (check Solvents page for details), this would be the configuration file we should use:

[MDSETTINGS]
SOLVENTS = ETA, WAT
NREPL = 3
Notice we didn't give any length of simulation, temperature nor any md parameter. They will be taken from the defaults (check MD settings section for learning about default parameters). The default parameters include 20ns as standard simulation lenght. If we wish to make longer simulations we have two options:

Modify the defaults
Only modify length for this particular replicas
For modifying the length of simulation for these particular replicas (i.e. to 40ns), we should use the following configuration file:

[MDSETTINGS]
SOLVENTS = ETA, WAT
NREPL = 3
NANOS = 40
Most common options that are usually changed during replica creation (nanoseconds, temperature, restraining scheme and force), have their special treatment in the configuration file. A complete input where all these options are modified, could be the following:

[MDSETTINGS]
SOLVENTS = ETA, WAT
NREPL = 3
NANOS = 40
TEMP = 298
RESTRMODE = HA
RESTRFORCE = 0.05
With this input, we will generate 3 replicas for ETA mixture, 3 replicas for WAT solvent box each of them with restraints in all heavy atoms with a force of 0.05 kcal/molA^2 that will run for 40ns at 298K temperature.

Create replicas in a project

Once we have the Settings file ready (e.g. the file above saved as mysettings.cfg) with the replicas we want to add, it's time to create the folders and MD input files. To do so, move to an existing Project folder which should contain the system to solvate (e.g. here named SystemName) and type:

mdmix add replicas -f mysettings.cfg -sys SystemName
This will create a folder for each of the replicas (the Replica folder).

Replica Folder

The folder created will contain a file ended with mrepl which is internally used by pyMDMix to identify the folder and contains all replica parameters. This file should not be removed!
