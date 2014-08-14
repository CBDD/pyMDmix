When simulation is finished, we can start processing results. The process commonly involves the following steps:

Align the trajectory of different replicas to a common reference structure
Compute the density maps
Convert densities of independent replicas to a single energy map for each probe
Before entering into the analysis process, let me introduce you the selection algebra for selecting what replicas are going to be analyzed and other common options to the majority of analysis commands.

Replica selection for analysis commands

Most of pyMDMix analysis commands require the user to specify which replicas are going to be analyzed. The construction of such command has a shape similar to this one:

mdmix analyze [command] [selection] [command-options]
Where command depends on the analysis stage (alignment, density calculation, energy conversion, etc.) and command-options are command specific options that the user can always check by getting help (mdmix analyze [command] -h). Common command-options include the number of processors to use for the analysis or the nanoseconds to analyze (if we are interested in a sub-selection of the trajectory).

Selection methods:

all: All replicas in the project will be processed.
bysolvent -s [solvent list]: Giving a list of solvent box names, all replicas belonging to those solvents will be analyzed.
byname -s [replica name list]: Give a list of replica names and only those replicas will be analyzed
group -s [group name]: Give a name of a group and replicas belonging to that group will be analyzed. The group must be defined beforehand.
Examples:

We want to align trajectory of all replicas:

mdmix analyze align all
We want to calculate density maps for replicas with MAM solvent mixture:

mdmix analyze density bysolvent -s MAM
We want to convert to energies all density maps in replicas ETA_1 and ETA_2 specifically:

mdmix analyze energy byname -s ETA_1 ETA_2
We want to calculate density maps for a group we created before:

mdmix analyze density group -s ethanol_free
Step selection and number of processors to use

The simulation input files will be prepared to write one file for each nanosecond of simulation by default. It will be possible therefor to select a subset of nanoseconds to analyze by using -N option. If you change the default settings and produce files with larger number of nanoseconds, instead of nanosecond selection, will be step selection (as the program is only able to work with whole files and not portions of it). In any case, the command will be the same. To check what analysis commands include this possibility, call the command with -h flag to get help.

For instance, we might be interested in aligning the trajectory only on the first four nanoseconds (or steps) for all replicas belonging to ETA solvent (because the rest is still running):

mdmix analyze align bysolvent -s ETA -N1:4
Being the colon a range maker: all numbers from 1 to 4 in this case.

Besides the step selection, some analysis commands allow parallel execution using multiple processors (e.g. density command). In this case, we should tell the program how many processors to use with -C option.

mdmix analyze density bysolvent -s ETA -C8
This command will analyze all replicas belonging to ETA solvent, all known nanoseconds, using 8 processors.

Trajectory alignment

TODO

Density calculation

TODO

Energy conversion

TODO

Hot spots calculation

TODO

Plotting commands

TODO
