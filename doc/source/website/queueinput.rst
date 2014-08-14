pyMDMix is prepared for writing queue input files for each job to be executed during the simulation run (minimization, equilibration and production steps input files). As every queueing system has a particular configuration, it is not possible to provide a universal automatic input writing system.

The solution in pyMDMix comes with the queue template files. The user can generate a queue template which should contain a series of magic keywords that the program will use to fill in the appropriate commands for each simulation step.

The queue template file

To write scripts for grid distributed systems (queueing systems), the user should first save a template file inside the pyMDMix user directory (located at $HOME/.mdmix/) with suffix '_queue_temp.txt'. As an example, this is the template file for the SGE queueing system already distributed with the software:

#!/bin/csh
######################################
# pyMDMix Queue input template #
# http://mdmix.sourceforge.net #
# Usage: qsub FILE_NAME.csh #
######################################
###########################################################
# To use with SGE: Sun Grid Engine Queue Manager #
# -q SGE queue name. #
# -cwd Working in the current working directory. #
# -N SGE job name. #
# -o Job standard output file. #
# -e Job standard error file. #
###########################################################
#$ -q QUEUENAME
#$ -cwd
#$ -N {jobname}
#$ -o {jobname}.out
#$ -e {jobname}.err
###############################################################
{precommands}
{commands}
{postcommands}
qsub {next}
Adapt it to your needs, all keywords in curly brackets should be present in all templates you generate. These keywords are used by pyMDMix to fill in with corresponding simulation job information. All queue files will automatically chained (first file will submit to queue the second, and so on...).

Creating queue inputs

The prefix name of the saved file will be used by the program to identify the queue template. E.g. if we save a file with name MYQUEUE_queue_temp.txt, the queue name will be identified as MYQUEUE by pyMDMix.

Move to a project folder to create the queue input for all replicas:

mdmix queue list  # Will show all installed queues (template files)
mdmix queue create -n MYQUEUE # Will use MYQUEUE_queue_temp.txt to generate inputs
