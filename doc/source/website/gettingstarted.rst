x allows an easy set up of several simulations for the same system under different conditions: solvent, temperature, restraining schemes, etc. Moreover, after simulations are done, many analysis tools will help you to quality check and extract useful information from these simulations in a aqueous-organic environments.

Project setup up workflow

Manual system preparation. Generate an AMBER object file with Leap (or PDB file) which is the main input to pyMDMix. Check out System preparation.
Create an empty pyMDMix project and add the system to study.
Create replicas of the system: A replica is considered as a single independent MD system and will contain all input files needed to run the simulation (including commands to be executed and queue input files if requested). Chose solvent and simulation conditions (check MD settings for details).
Run the simulations. This is done by the user in his/her computation facility. The execution commands may be in this point adapted to your cluster specifications (i.e. modify queue input files and so on). This step is out of control of the program as every cluster has specific operation rules.
Bring back results to the project and start analysis.
Analyze results! Common analysis involves the generation of density maps for each solvent and probe. These maps are transformed to free energies. Check out Analysis Guide for more information and details.
pyMDMix diagram

To better understand the program operation, this picture illustrates some of the main objects the program internally manages and how the user interacts with them through a pyMDMix project. This graphic will help you understand many parts of the documentation:

allobjectsoverview

Default MDMix simulation protocol

When building the project, input files to run the simulations will be automatically generated. The default simulation protocol implemented in the package is the following:

Minimization: 5000 steps.
Equilibration: Heating from 100K to the final temperature in 800ps in NVT ensemble + 1ns of NPT equilibration at the final temperature. Langevin thermostat. 2fs timestep.
Production: NVT ensemble with Langevin Thermostat, PME for electrostatic calculations. Timestep of 2fs. Each file will correspond to 1ns of simulation: writing frequency of 500 and number of steps of 500.000.
This default parameters can be modified or adapted at each program setup (see MD setting section).

