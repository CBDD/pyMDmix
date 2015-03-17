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
Action to calculate density grids
"""
import os
import os.path as osp
import logging
import tempfile
import numpy as npy
import Biskit as bi

import multiprocessing

import pyMDMix
import pyMDMix.tools as T
from pyMDMix.GridsManager import NewGrid

class DensityError(Exception):
    pass

class DensityGrids(object):
    action_name = "DensityGrids"
    def __init__(self, replica, probeselection=False, includeCOM=False, 
                        onlyCOM=False, subregion=None, stepselection=[], 
                        reference=False, outprefix='', *args, **kwargs):
        """Calculate density grids from a previously aligned trajectory for all probes in replica and its center of mass

        :arg list probeselection: List of probe names to calculate for replica solvent. If not given, all probes will be calculated.
        :arg bool includeCOM: Calculate center of mass coordinates as an extra probe.
        :arg bool onlyCOM: Consider only COM probes for calculation. Any probeselection list will be ignored as well as the probelist in the solvent.
        :arg tuple subregion: coordinates of minimum and maximum point to consider. 
        :arg str reference: Path to a reference PDB file over which to construct the initial 
                grid where counts will be added. If not given, the default replica reference PDB file will be used.
                Important to give this argument if trajectory was aligned using a reference PDB other than the replica's one.
        """
        self.log = logging.getLogger("DensityGrids")
        self.replica = replica

        if not isinstance(replica, pyMDMix.Replica): raise DensityError, "replica argument of wrong type."
                
        self.log.info("Setting up density grids calculation for replica %s"%replica.name)
        if not replica.isAligned(stepselection):
            raise DensityError, "Cannot calculate density over non-aligned trajectory"
        
        self.solvent = replica.getSolvent()
        if not self.solvent:
            raise DensityError, "Cannot fetch solvent %s from the database! Make sure there are no conflicting files"%(replica.solvent)
        self.pdb = replica.getPDB()
        
        self.probeselection= probeselection
        self.includeCOM= includeCOM
        self.onlyCOM = onlyCOM
        if outprefix: self.outprefix = outprefix
        else: self.outprefix = ''

        # Wrok only on subregion?
        if subregion:
            if len(subregion) != 2:
                raise DensityError, "subregion argument should have two xyz coordinates: min and max. E.G.: ((x0,y0,z0),(x1,y1,z1))"
            if len(subregion[0]) != 3 or len(subregion[1] != 3):
                raise DensityError, "subregion argument should have two xyz coordinates: min and max. E.G.: ((x0,y0,z0),(x1,y1,z1))"
        self.subregion = subregion
        
        # To be set un setup()
        self.probes = None
        self.container = None
        self.countGrids = {}
        
        # Reference given or use replica's one?
        self.ref = osp.join(self.replica.path,self.replica.ref)
        if reference:
            if osp.exists(reference): 
                self.log.info("Using PDB file %s as reference to construct density grid"%reference)
                self.ref=osp.abspath(reference)
            else:
                raise DensityError, "Reference PDB file %s not found."%reference
            
        self.setup()        
        
    
    def setup(self):
        self.__setProbes()
        self.prepareGrids()
        f = self.container.getIndexFunction()
        def fx(x):
            r = f(x)
            if not r: return [-1,-1,-1]
            else: return r
        self.indexFunction = fx

        # Set subregion function?
        self.__setSubregion()
        self.log.info("Ready to calculate")

    def __setSubregion(self):
        if self.subregion:
            def subfx(xyz):
                "x is an Nx3 numpy array"
                x,y,z = xyz.T
                validx = (x < self.subregion[0][0]) * (x > self.subregion[1][0])
                validy = (y < self.subregion[0][1]) * (y > self.subregion[1][1])
                validz = (z < self.subregion[0][2]) * (z > self.subregion[1][2])
                validcoords = validx*validy*validz
                return xyz[validcoords]
            self.subregionfx = subfx
        else:
            self.subregionfx = None
            
    def __setProbes(self):
        if self.onlyCOM:
            self.log.info("Density calculation over COM probes ONLY")
            self.probes = self.solvent.comprobes
        elif self.probeselection:
            self.log.info("User list of probes for density calculation: %s"%self.probeselection)
            l = []
            for p in self.probeselection:
                if p in self.solvent.probelist: l.append(p)
                elif p in self.solvent.comprobes: 
                    l.append(p)
                    self.includeCOM = True
                else: self.log.debug("Probe %s not in current solvent. skipping..."%p)
            if not l:
                raise DensityError, "No probes selected"
            self.probes = l
        else:
            l = self.solvent.probelist
            if self.includeCOM: l+=self.solvent.comprobes
            self.probes = l
            
        self.log.info("Calculating density for probes: %s"%self.probes)
        
    def __calcGridDimensionsAndOrigin(self):
        "From reference PDB calculate dimensions of the grid to be calculated and the origin"
        refpdb = bi.PDBModel(self.ref)
        maxdist = npy.sqrt(((refpdb.xyz.max(axis=0) - refpdb.xyz.min(axis=0))**2).sum())
        dimensions = npy.round(maxdist)
        origin = refpdb.xyz.mean(axis=0) - (dimensions/2.)
        
        # Move origin to nearest half integer
        mask = npy.abs(origin - npy.trunc(origin)) < 0.5 # True when first decimal is bigger than 0.5
        origin = npy.trunc(origin)
        origin[mask == 0]+=0.5
        
        return dimensions, origin

    def __calcSubRegionGrid(self):
        "Calc dimensions and origin of subregion grid"
        min, max = npy.array(self.subregion)
        d = max-min # dmension + 2.5Angstroms buffer
        origin = min
        
        # Move origin to nearest half integer
        mask = npy.abs(origin - npy.trunc(origin)) < 0.5 # True when first decimal is bigger than 0.5
        origin = npy.trunc(origin)
        origin[mask == 0]+=0.5
                
        return d, origin

    def prepareGrids(self):
        "Build empty grid knowing size and spacing centered at origin. Make npy.memmaps for each probe."
        if self.subregion:
            dimensions, origin = self.__calcSubRegionGrid()
        else:
            dimensions, origin = self.__calcGridDimensionsAndOrigin()
        self.log.debug("Creating grid template with dimensions: %s and origin: %s"%(dimensions,origin))
        self.container = NewGrid(dimensions, origin=origin, dtype=npy.uint32)

        for probe in self.probes:
            tmp = tempfile.mktemp(prefix='mdmix_mmap_')
            cmmap = npy.memmap(tmp, mode='w+', dtype='uint32',shape=self.container.data.shape)
            self.countGrids[probe] = cmmap

    def prepareParallelWorkers(self, dataQueue, nworkers):
        "Will return as many Processes as probes to study. Each will need a snapshot to run. So when sending snapshots to the task queue, repeat each snap as needed"
        # Substitute grids' data by memmap array/files
        # Instantiate workers with memmap arguments and probe info
        workerList = []
        for _ in range(nworkers):
            workerList.append(CountGridsWorker(dataQueue, self.probes, self.pdb,
                                self.countGrids, self.indexFunction, self.subregionfx))
        return workerList

    def calcResults(self):
        "Create Grid for each resulting counts grid"
        self.results = {}
        self.log.info("Processing density results for replica %s"%self.replica.name)
        self.replica.go()
        if not osp.exists(self.replica.densityfolder): os.mkdir(self.replica.densityfolder)
        for probe, data in self.countGrids.iteritems():
            g = self.container.copy()
            g.update(data)
            g.setProbe(probe)
            g.setType('MDMIX_DENS')
            out = osp.join(self.replica.densityfolder, '{}{}_{}.dx'.format(self.outprefix, self.replica.name, probe))
            self.log.info("Writing density grid %s"%out)
            g.writeDX(out)
            # Remove tempfile
            os.remove(data.filename)
            self.results[probe] = g
        
        pyMDMix.browser.goback()
        return self.results

#    def run(self, frameandnum):
#        T.EXECUTOR.submitCmd(fx=self.__run, args=frameandnum, kwargs={})

    def run(self, frameandnum):
        snapshot, num = frameandnum
        self.pdb.setXyz(snapshot)
        self.log.debug("+ Frame %i"%num)
        # Add counts to corresponding grid
        for probe in self.probes:
            coords = self.pdb.getProbeCoords(probe)
            if self.subregionfx: coords = self.subregionfx(coords)
            idx = npy.apply_along_axis(self.indexFunction, 1, coords)
            for idx in idx:
                if npy.all(idx == [-1,-1,-1]): continue
                self.countGrids[probe][tuple(idx)] += 1

class CountGridsWorker(multiprocessing.Process):
    def __init__(self, snapQueue, probelist, solvatedPdb, outCountGrids, indexFunction, subregionfx):
        multiprocessing.Process.__init__(self)
        self.snapQueue = snapQueue
        self.probelist = probelist
        self.pdb = solvatedPdb.clone()
        self.outCountGrids = outCountGrids
        self.toIndex = indexFunction
        self.subregionfx = subregionfx

    def run(self):
        while True:
            # Get next task or exit if None
            snapshot = self.snapQueue.get()
            if snapshot is None: break
            frame = snapshot[1]
#            print "++ Frame %i"%frame
            snapshot = snapshot[0]
            self.pdb.setXyz(snapshot)
            # Add counts to corresponding grid
            for probe in self.probelist:
                coords = self.pdb.getProbeCoords(probe)
                if self.subregionfx: coords = self.subregionfx(coords)
                idx = npy.apply_along_axis(self.toIndex, 1, coords)
                for idx in idx:
                    if npy.all(idx == [-1,-1,-1]): continue
                    try:
                        self.outCountGrids[probe][tuple(idx)] += 1
                    except IndexError:
                        raise IndexError, "{} {} {} grid shape: {}".format(idx, tuple(idx), probe, self.outCountGrids[probe].shape)

class CountProteinWorker(multiprocessing.Process):
    def __init__(self, snapQueue, lock, mask, solvatedPdb, outCountGrids, indexFunction, subregionfx):
        multiprocessing.Process.__init__(self)
        self.snapQueue = snapQueue
        self.lock = lock
        self.mask = mask
        self.pdb = solvatedPdb.clone()
        self.outCountGrids = outCountGrids
        self.toIndex = indexFunction
        self.subregionfx = subregionfx

    def run(self):
        while True:
            # Get next task or exit if None
            snapshot = self.snapQueue.get()
            if snapshot is None: break
            frame = snapshot[1]
#            print "++ Frame %i"%frame
            snapshot = snapshot[0]
            self.pdb.setXyz(snapshot)
            # Add counts to corresponding grid
            # Add counts to corresponding grid
            coords = self.pdb.xyz[self.mask]
            if self.subregionfx: coords = self.subregionfx(coords)
            idx = npy.apply_along_axis(self.toIndex, 1, coords)
            with self.lock:
                for idx in idx:
                    if npy.all(idx == [-1,-1,-1]): continue
                    self.outCountGrids[tuple(idx)] += 1

class DensityProtein(object):
    action_name = "DensityProtein"
    def __init__(self, replica, subregion=None, stepselection=[], outprefix='', *args, **kwargs):
        """Calculate density grids from a previously aligned trajectory for all atoms in protein.

        :arg replica: Replica to work on
        :type replica: :class:`Replicas.Replica`
        :arg tuple subregion: coordinates of minimum and maximum point to consider. 
        """
        self.log = logging.getLogger("DensityProtein")
        if not isinstance(replica, pyMDMix.Replica): raise DensityError, "replica argument of wrong type."
        self.replica = replica
        self.log.info("Setting up density grids calculation BY MASK for replica %s"%replica.name)
        if not replica.isAligned(stepselection):
            raise DensityError, "Cannot calculate density over non-aligned trajectory"
        
        self.pdb = replica.getPDB()
        
        if outprefix: self.outprefix = outprefix
        else: self.outprefix = ''

        # Wrok only on subregion?
        if subregion:
            if len(subregion) != 2:
                raise DensityError, "subregion argument should have two xyz coordinates: min and max. E.G.: ((x0,y0,z0),(x1,y1,z1))"
            if len(subregion[0]) != 3 or len(subregion[1] != 3):
                raise DensityError, "subregion argument should have two xyz coordinates: min and max. E.G.: ((x0,y0,z0),(x1,y1,z1))"
        self.subregion = subregion
        
        # To be set in setup()
        self.container = None
        self.countGrid = None
        
        self.setup()
    
    def setup(self):
        self.prepareGrid()
        f = self.container.getIndexFunction()
        def fx(x):
            r = f(x)
            if not r: return [-1,-1,-1]
            else: return r
        self.indexFunction = fx

        # Set protein mask
        self.mask = self.pdb.soluteMask

        # Set subregion function?
        self.__setSubregion()
        self.log.info("Ready to calculate")

    def __setSubregion(self):
        if self.subregion:
            def subfx(xyz):
                "x is an Nx3 numpy array"
                x,y,z = xyz.T
                validx = (x < self.subregion[0][0]) * (x > self.subregion[1][0])
                validy = (y < self.subregion[0][1]) * (y > self.subregion[1][1])
                validz = (z < self.subregion[0][2]) * (z > self.subregion[1][2])
                validcoords = validx*validy*validz
                return xyz[validcoords]
            self.subregionfx = subfx
        else:
            self.subregionfx = None
                    
    def __calcGridDimensionsAndOrigin(self):
        "From reference PDB calculate dimensions of the grid to be calculated and the origin"
        from pyMDMix.PDB import SolvatedPDB
        refpdb = SolvatedPDB(osp.join(self.replica.path, self.replica.ref)).getSolute()
        maxdist = npy.sqrt(((refpdb.xyz.max(axis=0) - refpdb.xyz.min(axis=0))**2).sum())
        dimensions = npy.round(maxdist) + 5 #Add 5 Ansgroms
        origin = refpdb.xyz.mean(axis=0) - (dimensions/2.)
        
        # Move origin to nearest half integer
        mask = npy.abs(origin - npy.trunc(origin)) < 0.5 # True when first decimal is bigger than 0.5
        origin = npy.trunc(origin)
        origin[mask == 0]+=0.5
                
        return dimensions, origin

    def __calcSubRegionGrid(self):
        "Calc dimensions and origin of subregion grid"
        min, max = npy.array(self.subregion)
        d = max-min # dmension + 2.5Angstroms buffer
        origin = min
        
        # Move origin to nearest half integer
        mask = npy.abs(origin - npy.trunc(origin)) < 0.5 # True when first decimal is bigger than 0.5
        origin = npy.trunc(origin)
        origin[mask == 0]+=0.5
                
        return d, origin

    def prepareGrid(self):
        "Build empty grid knowing size and spacing centered at origin. Make npy.memmaps for each probe."
        if self.subregion:
            dimensions, origin = self.__calcSubRegionGrid()
        else:
            dimensions, origin = self.__calcGridDimensionsAndOrigin()
        self.log.debug("Creating grid template with dimensions: %s and origin: %s"%(dimensions,origin))
        self.container = NewGrid(dimensions, origin=origin, dtype=npy.uint32)
        tmp = tempfile.mktemp(prefix='mdmix_mmap_')
        cmmap = npy.memmap(tmp, mode='w+', dtype='uint32',shape=self.container.data.shape)
        self.countGrid = cmmap

    def prepareParallelWorkers(self, dataQueue, nworkers):
        "Will return as many Processes as probes to study. Each will need a snapshot to run. So when sending snapshots to the task queue, repeat each snap as needed"
        # Substitute grids' data by memmap array/files
        # Instantiate workers with memmap arguments and probe info
        workerList = []
        from multiprocessing import Lock
        lock = Lock()
        for _ in range(nworkers):
            workerList.append(CountProteinWorker(dataQueue, lock, self.mask, self.pdb,
                                self.countGrid, self.indexFunction, self.subregionfx))
        return workerList

    def calcResults(self):
        "Create Grid for each resulting counts grid"
        self.results = {}
        self.log.info("Processing density results for replica %s"%self.replica.name)
        self.replica.go()
        if not osp.exists(self.replica.densityfolder): os.mkdir(self.replica.densityfolder)
        g = self.container.copy()
        g.update(self.countGrid)
        g.setProbe('PROT')
        g.setType('MDMIX_DENS')
        out = osp.join(self.replica.densityfolder, '{}{}_{}.dx'.format(self.outprefix, self.replica.name, 'protein'))
        self.log.info("Writing density grid %s"%out)
        g.writeDX(out)
        # Remove tempfile
        os.remove(self.countGrid.filename)
        self.results['prot'] = g
        
        pyMDMix.browser.goback()
        return self.results

    def run(self, frameandnum):
        snapshot, num = frameandnum
        self.pdb.setXyz(snapshot)
        self.log.debug("+ Frame %i"%num)
        # Add counts to corresponding grid
        coords = self.pdb.xyz[self.mask]
        if self.subregionfx: coords = self.subregionfx(coords)
        idx = npy.apply_along_axis(self.indexFunction, 1, coords)
        for idx in idx:
            if npy.all(idx == [-1,-1,-1]): continue
            self.countGrids[tuple(idx)] += 1


class DensityGridsAllHA(object):
    action_name = "DensityGridsAllHA"
    def __init__(self, replica, subregion=None, stepselection=[], outprefix='', *args, **kwargs):
        """Calculate density grids from a previously aligned trajectory for all heavy atoms and center of mass

        :arg list probeselection: List of probe names to calculate for replica solvent. If not given, all probes will be calculated.
        :arg bool includeCOM: Calculate center of mass coordinates as an extra probe.
        :arg bool onlyCOM: Consider only COM probes for calculation. Any probeselection list will be ignored as well as the probelist in the solvent.
        :arg tuple subregion: coordinates of minimum and maximum point to consider. 
        """
        self.log = logging.getLogger("DensityGridsAllHA")
        self.replica = replica
        self.log.info("Setting up density grid calculation for all heavy atoms in replica %s"%replica.name)
        if not replica.isAligned(stepselection):
            raise DensityError, "Cannot calculate density over non-aligned trajectory"
        
        self.solvent = replica.getSolvent()
        if not self.solvent:
            raise DensityError, "Cannot fetch solvent %s from the database! Make sure there are no conflicting files"%(replica.solvent)
        self.pdb = replica.getPDB()

        if outprefix: self.outprefix = outprefix
        else: self.outprefix = ''

        # Wrok only on subregion?
        if subregion:
            if len(subregion) != 2:
                raise DensityError, "subregion argument should have two xyz coordinates: min and max. E.G.: ((x0,y0,z0),(x1,y1,z1))"
            if len(subregion[0]) != 3 or len(subregion[1] != 3):
                raise DensityError, "subregion argument should have two xyz coordinates: min and max. E.G.: ((x0,y0,z0),(x1,y1,z1))"
        self.subregion = subregion
        
        # To be set un setup()
        self.probes = None
        self.container = None
        self.countGrids = {}
        
        if not isinstance(replica, pyMDMix.Replica): raise DensityError, "replica argument of wrong type."
        self.setup()
    
    def setup(self):
        self.setHAinfo()
        self.prepareGrids()
        f = self.container.getIndexFunction()
        def fx(x):
            r = f(x)
            if not r: return [-1,-1,-1]
            else: return r
        self.indexFunction = fx

        # Set subregion function?
        self.__setSubregion()
        self.log.info("Ready to calculate")

    def setHAinfo(self):
        "Prepare a list with all Heavy atoms to track and center of mass for each residue in the solvent"
        self.hainfo = {}
        self.probes = []
        for res in self.solvent.residues:
            if res == 'WAT' or res == 'HOH': continue # Don't track water
            self.hainfo[res.name] = dict((a.name,a.id) for a in res.atoms if a.element != 1) # take non-hydrogen atoms only (id and name)
            self.probes.extend(['%s_%s'%(res.name, a.name) for a in res.atoms if a.element != 1])
            self.probes.extend(['%s_%s'%(res.name, 'COM')])

    def __setSubregion(self):
        if self.subregion:
            def subfx(xyz):
                "x is an Nx3 numpy array"
                x,y,z = xyz.T
                validx = (x < self.subregion[0][0]) * (x > self.subregion[1][0])
                validy = (y < self.subregion[0][1]) * (y > self.subregion[1][1])
                validz = (z < self.subregion[0][2]) * (z > self.subregion[1][2])
                validcoords = validx*validy*validz
                return xyz[validcoords]
            self.subregionfx = subfx
        else:
            self.subregionfx = None
                    
    def __calcGridDimensionsAndOrigin(self):
        "From reference PDB calculate dimensions of the grid to be calculated and the origin"
        refpdb = bi.PDBModel(osp.join(self.replica.path,self.replica.ref))
        maxdist = npy.sqrt(((refpdb.xyz.max(axis=0) - refpdb.xyz.min(axis=0))**2).sum())
        dimensions = npy.round(maxdist)
        origin = refpdb.xyz.mean(axis=0) - (dimensions/2.)
        
        # Move origin to nearest half integer
        mask = npy.abs(origin - npy.trunc(origin)) < 0.5 # True when first decimal is bigger than 0.5
        origin = npy.trunc(origin)
        origin[mask == 0]+=0.5
        
        return dimensions, origin

    def __calcSubRegionGrid(self):
        "Calc dimensions and origin of subregion grid"
        min, max = npy.array(self.subregion)
        d = max-min # dmension + 2.5Angstroms buffer
        origin = min
        
        # Move origin to nearest half integer
        mask = npy.abs(origin - npy.trunc(origin)) < 0.5 # True when first decimal is bigger than 0.5
        origin = npy.trunc(origin)
        origin[mask == 0]+=0.5
                
        return d, origin

    def prepareGrids(self):
        "Build empty grid knowing size and spacing centered at origin. Make npy.memmaps for each probe."
        if self.subregion:
            dimensions, origin = self.__calcSubRegionGrid()
        else:
            dimensions, origin = self.__calcGridDimensionsAndOrigin()
        self.log.debug("Creating grid template with dimensions: %s and origin: %s"%(dimensions,origin))
        self.container = NewGrid(dimensions, origin=origin, dtype=npy.uint32)

        # Set one grid for each HA and COM
        for probe in self.probes:
            tmp = tempfile.mktemp(prefix='mdmix_mmap_')
            cmmap = npy.memmap(tmp, mode='w+', dtype='uint32',shape=self.container.data.shape)
            self.countGrids[probe] = cmmap

    def prepareParallelWorkers(self, dataQueue, nworkers):
        "Will return as many Processes as probes to study. Each will need a snapshot to run. So when sending snapshots to the task queue, repeat each snap as needed"
        # Substitute grids' data by memmap array/files
        # Instantiate workers with memmap arguments and probe info
        workerList = []
        for _ in range(nworkers):
            workerList.append(CountGridsWorkerHA(dataQueue, self.probes, self.pdb, self.hainfo,
                                self.countGrids, self.indexFunction, self.subregionfx))
        return workerList

    def calcResults(self):
        "Create Grid for each resulting counts grid"
        self.results = {}
        self.log.info("Processing density results for replica %s"%self.replica.name)
        self.replica.go()
        if not osp.exists(self.replica.densityfolder): os.mkdir(self.replica.densityfolder)
        for probe, data in self.countGrids.iteritems():
            g = self.container.copy()
            g.update(data)
            g.setProbe(probe)
            g.setType('MDMIX_DENS')
            out = osp.join(self.replica.densityfolder, '{}{}_{}.dx'.format(self.outprefix, self.replica.name, probe))
            self.log.info("Writing density grid %s"%out)
            g.writeDX(out)
            # Remove tempfile
            os.remove(data.filename)
            self.results[probe] = g
        
        pyMDMix.browser.goback()
        return self.results

#    def run(self, frameandnum):
#        T.EXECUTOR.submitCmd(fx=self.__run, args=frameandnum, kwargs={})

    def run(self, frameandnum):
        snapshot, num = frameandnum
        self.pdb.setXyz(snapshot)
        self.log.debug("+ Frame %i"%num)
        # Add counts to corresponding grid
        for probe in self.probes:
            if 'COM' in probe: coords = self.pdb.getProbeCoords(probe)
            else:
                # get all coords for HA name
                res, atname = probe.split('_')
                ati = self.hainfo[res][atname] -1
                coords = npy.array([rxyz[ati] for rxyz in self.pdb.iterResidues(res)])
            if self.subregionfx: coords = self.subregionfx(coords)
            idx = npy.apply_along_axis(self.indexFunction, 1, coords)
            for ix in idx:
                if npy.all(ix == [-1,-1,-1]): continue
                self.countGrids[probe][tuple(ix)] += 1

class CountGridsWorkerHA(multiprocessing.Process):
    def __init__(self, snapQueue, probelist, solvatedPdb, hainfo, outCountGrids, indexFunction, subregionfx):
        multiprocessing.Process.__init__(self)
        self.snapQueue = snapQueue
        self.probelist = probelist
        self.hainfo = hainfo
        self.pdb = solvatedPdb.clone()
        self.outCountGrids = outCountGrids
        self.toIndex = indexFunction
        self.subregionfx = subregionfx

    def run(self):
        while True:
            # Get next task or exit if None
            snapshot = self.snapQueue.get()
            if snapshot is None: break
            frame = snapshot[1]
#            print "++ Frame %i"%frame
            snapshot = snapshot[0]
            self.pdb.setXyz(snapshot)
            # Add counts to corresponding grid
            for probe in self.probelist:
                if 'COM' in probe: coords = self.pdb.getProbeCoords(probe)
                else:
                    # get all coords for HA name
                    res, atname = probe.split('_')
                    ati = self.hainfo[res][atname] -1
                    coords = npy.array([rxyz[ati] for rxyz in self.pdb.iterResidues(res)])
                if self.subregionfx: coords = self.subregionfx(coords)
                idx = npy.apply_along_axis(self.toIndex, 1, coords)
                for ix in idx:
                    if npy.all(ix == [-1,-1,-1]): continue
                    try:
                        self.outCountGrids[probe][tuple(ix)] += 1
                    except IndexError:
                        raise IndexError, "{} {} {} grid shape: {}".format(ix, tuple(idx), probe, self.outCountGrids[probe].shape)


def DensityGrids_postprocess(results, replica, **kwargs):
    """
    Save results of a density calculation for replica *replica*.

    :arg dict densityResults: Result from a :class:`Action.Density` call.
    :arg replica: Replica were to save results.
    :type replica: :class:`Replicas.Replica`
    """
    replica.go()
    if not osp.exists(replica.densityfolder): os.mkdir(replica.densityfolder)
    for probe, grid in results.iteritems():
        grid.writeDX(osp.join(replica.folder, probe+'.dx'))
    pyMDMix.browser.goback()
    
def DensityGridsAllHA_postprocess(results, replica, **kwargs):
    """
    Save results of a density calculation for replica *replica*.

    :arg dict densityResults: Result from a :class:`Action.Density` call.
    :arg replica: Replica were to save results.
    :type replica: :class:`Replicas.Replica`
    """
    replica.go()
    if not osp.exists(replica.densityfolder): os.mkdir(replica.densityfolder)
    for probe, grid in results.iteritems():
        grid.writeDX(osp.join(replica.folder, probe+'.dx'))
    pyMDMix.browser.goback()