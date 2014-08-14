(visit Getting started pictogram for a summary on the Project creation workflow)

The pyMDMix project consists of different systems and individual simulation set ups for the systems (the Replicas). Once a project is started, all program commands should be issued inside the project folder.

Creating an empty project

mdmix create project -n MyProject
This will write a new folder under the current working directory with name MyProject. This new folder is considered the

PROJECT FOLDER
from now on and across all documentation. To work with the project, we should move to the recently created directory.

The project folder contains a file named MyProject.mproj. This file is internally used by the project and should never be removed.

Adding systems

The project is prepared for dealing with more than one system. This option has been implemented thinking in different states of a same macromolecule that could be simulated independently applying restraints. In the majority of situations, we will be interested in adding just one system. Check System preparation section for details and examples on how to create the system and prepare the System configuration file. This system is a macromolecule prepared for a simulation in vacuo. Solvent will be added when creating replicas.

In this example, we have created a system identified with name SystemA inside configuration file  mysystem.cfg. To add it to current project:

mdmix add system -f mysystem.cfg
This command will save a file SystemA.msys in the project folder for internal use. Don't remove this file. Remember this command should be executed inside the project folder.

Creating replicas

Once a system is added to the project, we can start creating Replicas of the system. To do so, we have to prepare a Settings configuration file. This file contains the set of parameters needed to create the solvated system and the MD input files. NREPL option allows you to create several replicas sharing the same settings at once. Check The Replicas and MD Settings for details on default parameters and examples.

We have decided the replicas to create and saved all the settings in replicas_settings.cfg file. To finally create the MD inputs needed to run the simulations, execute:

mdmix add replicas -f replicas_settings.cfg -sys SystemA
New folders will be created inside MD folder. Each folder corresponds to one of the replicas. Take now these folders to run the simulation!

Creating projects in one step

A more comfortable way of creating a project from scratch when you already decided the system and what replicas to simulate is to join Settings configuration file and System configuration file sections into one unique file (e.g. joining previous two example cfg files into SystemA_replicas.cfg). The project would be created with:

mdmix create project -n MyProject -f SystemA_replicas.cfg
The separate method of creation is interesting when you decide to add more replicas to an already created project.

Adding groups

During analysis process you might want to do joint analysis of a group of replicas that share same simulation conditions. For instance, if we have simulated a system with ethanol using restraints and without restraints, we will be interested in analyzing these two groups separately. To do so, we should create a group:

mdmix add group -gn ethanol_restrained -s ETA_1 ETA_2 ETA_3
mdmix add group -gn ethanol_free -s ETA_4 ETA_5 ETA_6
We have created two groups: ethanol_restrained and ethanol_free. The first one is composed by ETA_1 to ETA_3 and the second from ETA_4 to ETA_6 replica names. Check Analysis Guide section for details on how to use this groups for joint analysis.

Get project description

By calling the following command inside a project folder, a summary of the project status will be printed. Information on the systems linked, the replicas created or groups defined will be printed. More over, a summary of the different replica parameters and the simulation status will be checked and printed (that is, the program will check if all simulation steps are complete and the trajectory has been aligned).

> mdmix info project
------------------------------
PROJECT pep_project INFO
------------------------------
SYSTEMS:
========
SYSTEM:pep
REPLICAS
========
REPLICA:ETA_1 system:pep_ETA nanos:1 restrMode:HA restrForce:0.010 Min:True Eq:True Prod:TrueAlign:True
REPLICA:WAT_1 system:pep_WAT nanos:1 restrMode:HA restrForce:0.010 Min:True Eq:True Prod:TrueAlign:True
GROUPS
======
project_group: ETA_1, WAT_1
------------------------------
