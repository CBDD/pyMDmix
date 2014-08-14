#  ------------------------ pyMDMix ----------------------------------- 
#                  http://mdmix.sourceforge.net
#  -------------------------------------------------------------------- 
# 
#  Software for preparation, analysis and quality control
#  of solvent mixtures molecular dynamics.
# 
#  Copyright (C) 2014 daniel
# 
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
# 
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
# 
#  Please cite your use of pyMDMix in published work:
# 
#              TOBEPUBLISHED
# 
#  -------------------------------------------------------------------- 
"""
Action to calculate residence plots
Report residue ids occupying certain hotspot or coordinate+tolerance regions.
"""

__author__="dalvarez"
__date__ ="$17-Jun-2014 15:20$"

import os
import os.path as osp
import numpy as npy
import multiprocessing
import logging
import scipy
if scipy.__version__ > "0.11.0":
    from scipy.spatial import KDTree as KDTree
else:
    from scipy.spatial import KDTree

import pyMDMix
import pyMDMix.tools as T

def simplifyNestedList(inlist, l=[]):
    "Return a plain list from nested list"
    for el in inlist:
        if isinstance(el, list):
            simplifyNestedList(el, l)
#            [l.append(i) for i in el]
        else:
            l.append(el)
    return l

class ResidenceError(Exception):
    pass

class Residence_Worker(multiprocessing.Process):
    def __init__(self, snapQueue, hotspot_coords, tolerance, maskNoH, mapIndexToResID,
                        mapResIDToResName, trackResNames, results, lock):
        multiprocessing.Process.__init__(self)
        self.snapQueue = snapQueue
        self.maskNoH = maskNoH
        self.hotspot_coords = hotspot_coords
        self.tolerance = tolerance
        self.mapIndexToResID = mapIndexToResID
        self.mapResIDToResName = mapResIDToResName
        self.trackResNames = trackResNames
        self.results = results
#        print 'thread dict id %d'%(id(self.results))
        self.lock = lock

    def run(self):
        "Same as Bfactor run but adapted to multiprocessing threading"
        while True:
            # Get next task or exit if None
            snapshot = self.snapQueue.get()
            if snapshot is None: break
            frame, framenum = snapshot
            framecoords = frame[self.maskNoH,:]

            # Fetch ids near hotspot_coords using a KDTree
            tree = KDTree(framecoords)    
            ids = simplifyNestedList(tree.query_ball_point(self.hotspot_coords, self.tolerance),[])
            del(tree)

            if npy.any(ids):
                # Save unique IDs for tracked residues
                resids = npy.unique(map(self.mapIndexToResID.get, ids))
                resnames = map(self.mapResIDToResName.get, resids)
                if self.trackResNames:
                    m = npy.array([r in self.trackResNames for r in resnames])
                    if npy.any(m):
                        if len(m) > 1:
                            resids = resids[m]
                        else:
                            # m is true but only one element. It will make npy crash if taking the mask
                            resids = resids
                    else:
                        # Noone of the ids is a tracked residue. Add a zero!
                        resids = npy.array([0])
                self.results.update({framenum:resids.tolist()})
#                    print id(self.results)
            else:
                # Not occupied
#                with self.lock: 
                self.results.update({framenum:[0]})
#                print "Added frame %d to results. empty."%framenum
        return self.results

class Residence(object):
    action_name = "Residence"
    def __init__(self, replica, hotspot=None, spherecenter=None, tolerance=0.5,
                        trackResNames=False, stepselection=[], parallel=True,
                        outfilename='', *args, **kwargs):
        """
        Parse a trajectory and check if a region is occupied or not and by what residue ID and name.
        
        The region must be defined using a sphere: with a tolerance around the center coordinate in Angstrom.

        Args:
            replica         (ReplicaInfo)       What replica are we considering
            hotspot         (HotSpot instance)  Hotspot to study.
            spherecenter    (npy.array)         3D coordinates of a cartesian point to consider the center of the sphere.
            tolerance       (float)             Radius around the hotspot coordinates or sphere center to track if residues are inside or not.
            trackResNames   (list)              Residue names to track if they are inside the hotspot.
                                                ATTENTION: If no list is given. Then all residue names can be occupying the hotspot (even protein residues).
            parallel        (bool)              Use parallel calculation.
            outfilename     (str)               File name where to store the resutls in ASCII format. 
                                                If not given, a default 'occupancy_results.txt' will be used and saved in folder of execution.
        Returs:

            Once all trajectory is parsed, execute calcResults() and getResults() to obtain a dictionary
            with the resulting occupied frames and name mappings.

            Read calcResults() documentation for returning dict format.
        """
        self.log = logging.getLogger("Residence")
        if not isinstance(replica, pyMDMix.Replica): raise DensityError, "replica argument of wrong type."
        self.replica= replica
        self.pdb = replica.getPDB()
        self.parallel = parallel
        self.hotspot = hotspot
        self.spherecenter = spherecenter
        self.tolerance = tolerance
        self.trackResNames = trackResNames
        self.outfilename = outfilename
	self.hotspot_coords = []
        self.workerList = []
        self.results = {}
        self.mapIndexToResID = {}

        # Algned trajetory is needed
        if not replica.isAligned(stepselection):
            raise ResidenceError, "Cannot calculate residence plots over non-aligned trajectory"

        self.setup()

    def setup(self):
        "Get a mask for the atoms to track in bfactor calculation"
        # Set coordinates
        if isinstance(self.spherecenter, npy.ndarray) and self.spherecenter.size == 3:
            #User gave already a mask. Use that one!
            self.hotspot_coords = [self.spherecenter]
            self.log.info("Residence: Tracking occupancy with sphere method. Center: %s"%self.hotspot_coords)
        elif self.hotspot:
            import pyMDMix.HotSpotsManager as HM
            if isinstance(self.hotspot, HM.HotSpot):
                self.hotspot_coords = self.hotspot.coordList
                self.log.info("Residence: Tracking occupancy of hotspot %s"%self.hotspot)
            else:
                raise ResidenceError, "Wrong 'hotspot' argument type"
        else:
            raise ResidenceError, "Hotspot+tolerance or Spherecenter+tolerance must be given."

        # Print info
        self.log.info("Residence: Resnames to track: %s"%self.trackResNames)
        self.log.info("Tolerance: %.2f"%self.tolerance)

        
        # Set maps to residuenames and ids
        self.maskNoH = ~self.pdb.maskH()
        self.mapIndexToResID = dict(enumerate(npy.array(self.pdb['residue_number'], dtype=int)[self.maskNoH]))
        self.mapResIDToResName = dict(zip(npy.array(self.pdb['residue_number'], dtype=int)[self.maskNoH], 
                                        npy.array(self.pdb['residue_name'])[self.maskNoH]))
        # Add a dummy residue with id 0 to identify non-occupied frames
        self.mapResIDToResName[0] = 'NO_RESIDENCE'
        
        # If trackResName not given, issue a warning! All residues will be included in the tracking! Even protein ones.
        if not self.trackResNames:
            self.log.warn("Residence: No trackResNames attribute given. Will track any residue falling into the hotspot (even protein ones!).")
        
    def run(self,frameandnum):
        "For each snapshot just save the coordinates of the atoms of interest in a file"
        frame, framenum = frameandnum
        framecoords = frame[self.maskNoH,:]

        # Fetch ids near hotspot_coords using a KDTree
        tree = KDTree(framecoords)  
        ids = simplifyNestedList(tree.query_ball_point(self.hotspot_coords, self.tolerance),[])
        del(tree)

        if npy.any(ids):
            # Save unique IDs for tracked residues
            resids = map(self.mapIndexToResID.get, ids)
            resnames = map(self.mapResIDToResName.get, resids)
            if self.trackResNames:
                m = [r in self.trackResNames for r in resnames]
                if npy.any(m):
                    if len(m) > 1:
                        resids = npy.array(resids)[m].tolist()
                    else:
                        # m is true but only one element. It will make npy crash if taking the mask
                        resids = npy.array(resids).tolist()
                else:
                    # Noone of the ids is a tracked residue. Add a zero!
                    resids = [0]
            self.results[framenum] = npy.unique(resids).tolist()
            return (1, frameandnum)
        else:
            # Not occupied
            self.results[framenum] = [0]
            return (0, frameandnum)

    def calcResults(self):
        """Process results"""
        self.log.info("Residence: processing results...")
        # Save results to a file
        if self.outfilename: self.outfilename = os.path.basename(self.outfilename)
        else: self.outfilename = 'occupancy_results.txt'
        self.outfilename=osp.join(self.replica.path, self.outfilename)
        self.log.info("Residence: Writting ASCII results under replica %s folder with filename %s ..."%(self.replica.name, self.outfilename))
        out = open(self.outfilename, 'w')

        # First lines for ID - Residue Name identification
        # Build a map: {RESNAME:[ID,ID,ID], ...}
        allids = set()
        [[allids.add(i) for i in framids] for framids in self.results.values()]
        resnames = map(self.mapResIDToResName.get, allids)
        finalmap = {}
        for i,idx in enumerate(allids):
            name = resnames[i]
            if not finalmap.has_key(name): finalmap[name] = []
            finalmap[name].append(idx)
        
        # Write map to first lines in file preceded with a #
        out.write("# Residence results study for replica %s"%self.replica.name)
        if self.trackResNames: out.write(". Tracked residues: %s\n"%self.trackResNames)
        else: out.write('\n')
        out.write("# RESNAME-RESID MAP\n")
        for k, v in finalmap.iteritems(): out.write("# %s = %s\n"%(k, ','.join(map(str,v))))
        out.write("# DATA\n")
        
        # Finally write frames and ids
        for frameid in sorted(self.results.keys()):
            out.write("%d\t%s\n"%(frameid, '\t'.join(map(str, self.results[frameid]))))
        
        out.close()        
        self.log.info("Results writen to %s"%os.path.abspath(self.outfilename))
        
        # Add name to resid map to results
        self.results['map'] = finalmap
        return dict(self.results)

    def getResults(self):
        return dict(self.results)

    def prepareParallelWorkers(self, dataQueue, nworkers):
        "Will return as many Processes as probes to study. Each will need a snapshot to run. So when sending snapshots to the task queue, repeat each snap as needed"
        # Instantiate workers
        self.lock = multiprocessing.Lock()
#        global results
#        results = self.results
        # Replace self.results dict with a Manager.dict for synchronous work
        self.results = multiprocessing.Manager().dict()
        for _ in range(nworkers):
            self.workerList.append(Residence_Worker(dataQueue, self.hotspot_coords, self.tolerance, self.maskNoH, 
                                                self.mapIndexToResID, self.mapResIDToResName, 
                                                self.trackResNames, self.results, self.lock))
        return self.workerList

def Residence_postprocess(results, replica, **kwargs):
    """
    Plot results of a residence calculation for replica *replica*.

    :arg dict ResidenceResults: Result from a :class:`Action.Residence` call.
    :arg replica: Replica were to save results.
    :type replica: :class:`Replicas.Replica`
    """
    from pyMDMix.Plotter import Plot
    plot = Plot()
    replica.go()
    plot.plotResidenceResults(results, outfilename='residence.png', **kwargs)
#    replica.go()
#    if not osp.exists(replica.densityfolder): os.mkdir(replica.densityfolder)
#    for probe, grid in results.iteritems():
#        grid.writeDX(osp.join(replica.folder, probe+'.dx'))
    pyMDMix.browser.goback()