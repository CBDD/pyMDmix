##  ------------------------ pyMDMix -----------------------------------
##                  http://mdmix.sourceforge.net
##  --------------------------------------------------------------------
## 
##  Software for preparation, analysis and quality control
##  of solvent mixtures molecular dynamics.
## 
##  Copyright (C) 2014 dalvarez
## 
##  This program is free software: you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
## 
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
## 
##  You should have received a copy of the GNU General Public License
##  along with this program.  If not, see <http://www.gnu.org/licenses/>.
## 
##  Please cite your use of pyMDMix in published work:
## 
##              TOBEPUBLISHED
## 
##  --------------------------------------------------------------------
"""
Analysis module introduces Analysis Class to manage replicas, trajectories and Actions.
"""
__author__="dalvarez"
__date__ ="$19-mar-2014 0:57:29$"

#import tools as T
import logging
import os

#from Solvents import SolventManager
from .Replicas import Replica
from .Actions import *
import multiprocessing as multi
#from GridsManager import Grid, GridSpace

class ActionsManagerError(Exception):
    pass

class ActionsManager(object):
    def __init__(self, ncpus=1, parallel=False):
        self.log = logging.getLogger("AnalysisManager")
        self.ncpus = ncpus
        if ncpus > 1: parallel = True
        self.parallel = parallel
        self.replicas = []
        self.actions = []
        self.jobs = {}
        self.results = {}

    def addReplicas(self, replicas):
        "Add replicas to current analysis"
        if not isinstance(replicas, list): replicas = [replicas]

        # Check al memebers are replicas
        for i in replicas:
            if not isinstance(i, Replica):
                raise ActionsManagerError("Expected type Replica but got %s instead"%(type(i)))

        # Add to current
        self.replicas.extend(replicas)

    def addActions(self, analysis):
        "Add analysis actions to run on each replica"
        if not isinstance(analysis, list): analysis = [analysis]

        # Check all memebers are Actions
        actions = []
        for act in analysis:
            if isinstance(act, str):
                a = globals().get(act)
            elif hasattr(DensityGrids,'action_name'):
                # Assume its an action class
                a = act
            else:
                raise ActionsManagerError("Unrecognized Action type: %s, %s"%(act,type(act)))
            if a: actions.append(a)
            else: raise ActionsManagerError("Action %s not found"%act)
        self.actions.extend(actions)

    def prepareRun(self, **kwargs):
        "Join replicas and actions"
        self.log.info("Preparing analysis...")
        for replica in self.replicas:
            self.jobs[replica] = [action(replica, **kwargs) for action in self.actions]

    def __submitReplica(self, replica, **kwargs):
        "Submit trajectory to diferent actions for same replia *replica*"
        traj = replica.getTrajectory(**kwargs)
        actions = self.jobs[replica]

        # Start workers and queues (one queue for each action
        # Start nworkers accordint to self.ncpus and trying to be
        # equitative between number of actions
        nworkers = [self.ncpus/len(actions) for _ in actions]
        for i in range(self.ncpus%len(actions)): nworkers[i] += 1 # add one more for each remaning cpu
        queues = [multi.Queue(n*5) for n in nworkers] # Maximum queue size is the number of workers times 5, to not overload and still keep them busy
        processes = []
        [processes.extend(act.prepareParallelWorkers(queues[i], nworkers[i])) for i,act in enumerate(actions)]

        # Start processes
        [p.start() for p in processes]

        i=0
        for step in traj:
            self.log.info("Working on file %s"%step.fname)
            for snaps in step:
                i+=1
                [q.put((snaps,i)) for q in queues]

        # Wait for jobs to complete
        for i,nw in enumerate(nworkers):
#            print "Stopping workers"
            [queues[i].put(None) for _ in range(nw)]
#        [p.terminate() for p in processes]
        [p.join() for p in processes]

        # Finally process results
        self.results[replica] = dict([(a.action_name, a.calcResults()) for a in actions])

    def run(self, **kwargs):
        """
        Execute each job as a multiprocessing.Process as each one should submit
        command by command to the global executor class
        """
        for repl in self.replicas:
            self.log.info('-'*50)
            self.log.info("Running %s analysis..."%repl.name)
            self.__submitReplica(repl, **kwargs)
            self.log.info('-'*50)

    def processResults(self, **kwargs):
        for repl, results in self.results.items():
            self.log.info("Processing results for replica %s"%repl.name)
            for actname, actres in results.items():
                self.log.info("Action %s..."%actname)
                process = globals().get(actname+'_postprocess')
                process(results=actres, replica=repl, **kwargs)




if __name__ == "__main__":
    print("""
import pyMDMix
import pyMDMix.Analysis as A
p=pyMDMix.loadProject()
anal=A.ActionsManager(3)
anal.addActions('Residence')
anal.addReplicas(p.replicas.values())
anal.prepareRun(spherecenter=[10,6,8])
anal.run()
anal.processResults()
""")
