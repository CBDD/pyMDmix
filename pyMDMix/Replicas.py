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
##              TOBEPUBLISHED^
## 
##  --------------------------------------------------------------------

"""
This module contains the main class :class:`Replica` for storage and manipulation
of a single simulation run.

This object contains all the MD parameters and system information needed to generate
the input files to run molecular dynamics.

Replica initialization and MD input creation
--------------------------------------------
The :class:`Replica` has to be configured using two sorts of information::
    - System related information, through a :class:`Systems.SolvatedSystem` object
    - MD related parameters, through a :class:`MDSettings.MDSettings` object

Once a replica has been created, a pickled version will be automatically saved inside the replica folder (with extension ``.mrepl``) and continuously updated. At any time,
you can re-load any replica using ``fromfile`` parameter:

    >>> previousreplica = R.Replica(fromfile='savedpickle.mrepl')
    
"""

__author__="dalvarez"
__date__ ="$20-ene-2014 22:56:58$"

import os
import os.path as osp
import tempfile
import logging
import tools as T

from Structures import FileLock

class ReplicaError(Exception):
    pass

class BadAttribute(ReplicaError):
    pass

class BadFile(ReplicaError):
    pass

class Replica(object):
    """
    Class to contain an independent simulation run (a Replica).
    Create folders, input files and control completeness of the different steps.

        :Arguments:
            -name (str): Replica name. Name should be given now or later with setName(name).
            -solvent (str):   Solvent name. It must exist in Solvent database.
            -fromfile (str):  Load existing Replica instance from pickled file. This instantiation way has priority.
            -path (str):  Path to replica folder structure. If not given now, can be assigned with :meth:`setPath`.
            -restrMask (str): Amber format mask to select residues and atoms to be restrained (if needed). If emtpy and restrains are requested, an automatic mask will be calculated from the pdb
            -alignMask (str): Amber format mask to select atoms and residues over which trajectory should be aligned. If empty, an automatic mask will be calcualted.
            -refPdb (str): Path to a PDB file used as reference for trajectory alignment.

        :Keywords:
            This provides a very flexible attribute assignment system.
            - Every pair key=value will be assigned as an attribute to current replica.
            - Pair values not given will take default values from Global Settings or User Settings specifications.

    """
    def __init__(self, system=None, mdsettings=None, name=None, fromfile=None, top=None, crd=None, **kwargs):
        """
        Constructor method for Replica objects from SolvatedSystem and MDSettings objects

        :arg system: SolvatedSystem or Empty System with macromolecule to be solvated. Can be given later.
        :type: :class:`Systems.System` or :class:`Systems.SolvatedSystem`
        :arg mdsettings: MD settings: solvent, restraints, length, etc.
        :type: :class:`MDSettings`

        :arg str name:   Replica name. Name should be given now or later with setName(name).
        :arg str fromfile:  Load existing Replica instance from pickled file. This instantiation way has priority.
        :arg str top: Amber Topology file path to be loaded.
        :arg str crd: Amber Coordinate file to be loaded. A solvated system will be created from top and crd.
        
        """
        if fromfile:
            # Load existing replica file
            self.load(fromfile)
            # Update paths if changed
            path = osp.split(osp.abspath(fromfile))[0]
            if path != self.path: self.setPath(path)
            # For back compatibility
            if not hasattr(self,'iwrap'): self.iwrap=1    
        else:
            # Folder structure and control
            self.__folderscreated = False
            self.__mdinput = False
            self.__backdir = None

            # Set filenames/formats/replica info
            self.name = name                    #: Name to identify the replica
            self.setName(name)

            # Init system
            import Systems
            if system:
                # Solvated system as input to create the replica
                if not isinstance(system, Systems.SolvatedSystem):  raise ReplicaError, "system should be SolvatedSystem instance. Got %s instead."%(type(system))
                self.system = system                #: SolvatedSystem. System with solvent box already created.
            elif top and crd:
                # Input as PRMTOP and PRMCRD, create solvated system from them
                self.system = Systems.SolvatedSystem(top=top, crd=crd)
            self.top = None                      #: Amber PRMTOP file with solvated system. To be set from input system.
            self.pdb = None                      #: PDB file generated from PRMTOP and PRMCRD. To be set from input system.
            self.crd = None                      #: Amber PRMCRD file with solvated system
            self.ref = None
            self.path = None                    #: Path to replica folder
            self.replFilePath = None            #: path to replica pickled file, should be placed inside replica directory

            # ADOPT ATTRIBUTES FROM MDSETTINGS
            # If not defined in kwargs
            if not mdsettings:
                from MDSettings import MDSettings
                mdsettings = MDSettings(**kwargs)
                
            for k,v in mdsettings.__dict__.iteritems():
                if k == 'name': continue
                if k == 'FF': continue
                setattr(self, k, v)

            # Add some operation attrs
            # Set names and paths if replicapath specified or name given,
            # then autocreate the paths with this root name
            self.extension = {} # Store extension format for different folders
            self.restrPdb = None    # For NAMD runs
            self.attached = {}      # Attached data
            self.create = False     
            self.grids = None

    def __str__(self):
        """Print some info about the replica"""
        info = """
REPLICA INFORMATION:
-------------
Replica name: {name}
Solvent: {solvent}

------- System: {system}
        -------
        top: {top}
        crd: {crd}
        pdb: {pdb}
        ref: {ref}

------- Files, folders and formats:
        ---------------------------
        Replica path: {path}
        MD production folder: {mdfolder}
        Minimization folder: {minfolder}
        Equilibration folder: {eqfolder}
        Aligned Trajectory folder: {alignfolder}
        Density grids folder: {densityfolder}
        Energy grids folder: {energyfolder}

        File extensions already detected: {extension}
        Write center trajectory in NetCDF?: {mdnetcdf}
        Can we read netcdf files from python?: {analyzenetcdf}

        MD Output filename template: '{mdoutfiletemplate}'
        Equilibration Output filename template: '{eqoutfiletemplate}'

------- SIMULATION:
        ------------
        MD program: {mdProgram}
        Temperature: {temp}
        Nanoseconds: {nanos}
        Use restraints?: {hasRestraints}

        """
        if self.hasRestraints:
            info+="""
------- RESTRAINTS:
        -----------
        Schema: {restrMode}
        Force: {restrForce}

        """
        return info.format(**self.__dict__)

    def __repr__(self):
        return "%s Replica"%self.name

    def __getstate__(self):
        d = self.__dict__.copy()
        del d['log']
        return d

    def __setstate__(self, d):
        d['log'] = logging.getLogger("Replica (%s)"%d['name'])
        self.__dict__.update(d)

    def __calcNumberSnapshots(self):
        self.nsnaps = self.ntrajfiles*(self.prod_steps/float(self.trajfrequency))

    def __calcExpectedNTrajFiles(self):
        "Calculate from self.nanos and depending on steps per file and timestep, the total number of trajetory files expected"
        # convert nanosecond to femtosecond (1e6), divide by num of fs per step to obtain number of steps
        # Finally divide num of steps needed by num of steps per file
        self.ntrajfiles=int((self.nanos*1e6/self.md_timestep)/self.prod_steps)

    def __fixTopology(self, prmtop):
        """
        Fix prmtop by removing SCEE and SCNB sections, they will take default values.
        File will be copied with same name and suffix ``_back`` to prevent file loss.
        Fixed topology will overwrite old filename.
        
        :arg str prmtop: File path to amber parm7 topology file.
        """
        import shutil
        shutil.copy(prmtop, prmtop+'_back')
        top = open(prmtop+'_back','r')
        out = open(prmtop,'w')
        line=top.readline()
        while line:
                if not '%FLAG SCEE_SCALE_FACTOR' in line:
                        out.write(line)
                else:
                        while not '%FLAG SOLTY ' in line:
                                line = top.readline()
                        out.write(line)
                line = top.readline()
        out.close()
        top.close()
        return True

    def desc(self):
        "Return a summary description string"
        s='REPLICA:%s\tsystem:%s\tnanos:%i'%(self.name, self.system.name, self.nanos)
        if self.hasRestraints: s+='\trestrMode:%s\trestrForce:%.3f'%(self.restrMode, self.restrForce)
        else: s+='\trestrMode:FREE'
        s+='\tMin:%s\tEq:%s\tProd:%s\tAlign:%s'%(self.isMinimizationFinished(), self.isEquilibrationFinished(),
                                                self.isProductionFinished(), self.isAligned())
        return s+'\n'

    def asMDSettings(self):
        "Return MDSettings object with settings for current replica"
        from MDSettings import MDSettings
        d = self.__dict__.copy()
        [d.pop(i) for i in ('log', 'top', 'crd', 'pdb','ref','path','attached',
                        'extension','replFilePath','system','restrPdb','create',
                        'hasRestraints', '_Replica__folderscreated','ntrajfiles',
                        'analyzenetcdf','_Replica__mdinput')]
        return MDSettings(**d)

    def fetchGrids(self, prefix=None, suffix=None, **kwargs):
        """
        Explore replica folder looking for known grid formats and fecthing type/probes.
        Additional filters are allowed: if *prefix* or *suffix* are given, only grids
        matching this starting or ending of filename will be returned (file extension is not included in suffix).
        
        :arg str prefix: String with filename prefix.
        :arg str suffix: String with filename suffix.
        """
        if not prefix: prefix = ''
        if not suffix: suffix = ''
        self.log.debug("Fetching replica %s grids (prefix: %s suffix: %s)"%(self.name, prefix, suffix))
        import GridsManager as GM
        grids = []
        for root, dir, files in os.walk(self.path):
            for f in files:
                fname, ext = osp.splitext(f)
                ext = ext.lstrip('.')
                if ext in ('dx','xplor','cns') and fname.startswith(prefix) and fname.endswith(suffix):
                    self.log.debug(" -> %s"%(f))
                    grids.append(GM.Grid(osp.join(root, f)))
        self.grids = grids
        return grids
        
    def getChecker(self, **kwargs):
        "Get MD checker according to the simulation program used"
        # Set checker according to mdProgram
        if self.mdProgram == 'AMBER':
            from Amber import AmberCheck
            return AmberCheck(self, **kwargs)
        if self.mdProgram == 'OPENMM':
            from OpenMM import OpenMMCheck
            return OpenMMCheck(self, **kwargs)
        elif self.mdProgram == 'NAMD':
            from NAMD import NAMDCheck
            return NAMDCheck(self, **kwargs)

    def getSolvent(self):
        "Return solvent instance assigned to current replica"
        import Solvents
        return Solvents.getSolvent(self.solvent)

    def getProbes(self):
        "Return a list with all possible probes to be calculated from solvent"
        solv = self.getSolvent()
        probes = solv.probelist
        probes.extend(solv.comprobes.keys())
        return probes

    def getTrajectory(self, stepselection=[], usealigned=True, framestep=1, frameselection=[]):
        """
        Get Trajectory object for current replica. If stepselection is given, only 
        files corresponding to those steps will be returned. If usealigned is True,
        will try to use aligned trajectory if it exists.

        :arg list stepselection: List of ints identifying steps to return as trajectory. Will return all if empty.
        :arg bool usealigned: Use aligned trajectory if possible.

        :return: :class:`Trajectory.Trajectory` object
        """
        from Trajectory import Trajectory
        
        if not stepselection: stepselection = range(1, self.ntrajfiles+1)
        else: self.log.info("Trajectory selected steps: %s"%stepselection)
        if usealigned:
            if not self.isAligned(stepselection): usealigned = False
        
        if usealigned:
            path = self.alignpath
            checkext = self.checkAlignExtension
            self.log.debug("Using aligned trajectory")
        else:
            if not self.isProductionFinished(stepselection):
                raise ReplicaError, "Cannot retrieve trajectory for non-finished steps: %s"%stepselection
            path = self.mdpath
            checkext = self.checkProductionExtension
            self.log.debug("Using not aligned trajectory")
        
        # Build File list to parse
        flist = []
        stepselection.sort()
        self.log.debug("Selected steps for trajectory: %s"%stepselection)
        for step in stepselection:
            extension = checkext(step)[step]
            if not extension:
                self.log.warn("File for step %i not found"%step)
                continue
            f = self.mdoutfiletemplate.format(step=step, extension=extension)
            flist.append(osp.join(path, f))

        self.log.debug("Filelist for trajectory: %s"%flist)
        #Build Trajectory object
        self.log.debug("Replica %s Trajectory: stepselection - %s ; framestep - %i "%(self.name, stepselection, framestep))
        return Trajectory(flist, self.getPDB(), step=framestep, frameselection=frameselection)

    def getPDB(self):
        "Return a SolvatedPDB with replica pdb"
        from PDB import SolvatedPDB
        p = SolvatedPDB(pdb=osp.join(self.path, self.pdb), solvent=self.solvent, extraResidues=self.system.extraResList)
#        p.setSolvent(self.solvent)
        p.fixNumbering()
        return p
    
    def getGridsByType(self, grid_type=None, **kwargs):
        """Return grids found in current replica.
        If *type* is given, will only return grids matching the type selected.
        If *type* is **None**, will return any grid found.
        :kwargs: Other params to pass to :meth:`fetchGrids`
        :return: A nested dictionary will be returned with this shape: ``{type: {probe: grid}}``.
        """
        if not self.grids: self.fetchGrids(**kwargs)
        d = {}
        grids_by_type = {g.type: [] for g in self.grids}
        for g in self.grids:
            grids_by_type[g.type].append(g)
        d = {
            gtype: {g.probe: g for g in glist} if grid_type is None or grid_type == gtype else {}
            for gtype, glist in grids_by_type.items()
        }
        return d

    def getGridsByProbe(self, probelist, **kwargs):
        """
        Return grids corresponding to selected probes. Inside each probe entry,
        multiple grid types can be returned if found.

        :arg list probelist: list of probe names to fetch grids
        :kwargs: Other params to pass to :meth:`fetchGrids`
        :return: dictionary with shape ``{probe:[grid, grid,..],}``
        """
        if not isinstance(probelist, list): probelist = [probelist]
        if not self.grids: self.fetchGrids(**kwargs)
        d = {}
        for g in self.grids:
            if set([g.probe]) & set(probelist):
                if not d.has_key(g.probe): d[g.probe] = []
                d[g.probe].append(g)
        return d
        
    def checkProductionExtension(self, steps=[]):
        """
        Check file extension for production trajectories in 'mdfolder'.
        If *steps* list is given, check only those steps files.

        :return: dictionary with {step:extension} format. If any format is not recognized or the file is missing **None** will be output.
        """
        import re
        self.go()
        self.log.debug("Checking production extension %s"%self.name)
        steps = steps or range(1, self.ntrajfiles+1)
        if not isinstance(steps, list): steps = [steps]
        result = {}

        files = os.listdir(self.mdfolder)
        for i in steps:
            n = self.mdoutfiletemplate.format(step=i, extension='')
            ext = re.compile(n+'(\w+)')
            tmp = []
            for f in files:
                m = ext.search(f)
                if m: tmp.append(m.groups()[0])
            match = set(tmp) & set(self.avail_trajext)
            if match:
                match = list(match)
                if len(match) > 1: self.log.warn("More than one valid extension found for production trajectory step %i: %s. Using first one: %s"%(i, match, match[0]))
                match = match[0]
            else: match = None
            result[i] = match

        T.BROWSER.goback()
        return result

    def checkAlignExtension(self, steps=[]):
        """
        Check file extension for aligned trajectories in 'alignfolder'.
        If *steps* list is given, check only those steps files.

        :return: dictionary with {step:extension} format. If any format is not recognized or file does not exists **None** will be output.
        """
        import re
        self.go()
        self.log.debug("Checking align extension %s"%self.name)
        steps = steps or range(1, self.ntrajfiles+1)
        if not isinstance(steps, list): steps = [steps]
        result = {}

        if not osp.exists(self.alignfolder):
            result = dict([(i, None) for i in steps])
            return result
        
        files = os.listdir(self.alignfolder)
        for i in steps:
            n = self.mdoutfiletemplate.format(step=i, extension='')
            ext = re.compile(n+'(\w+)')
            tmp = []
            for f in files:
                m = ext.search(f)
                if m: tmp.append(m.groups()[0])
            match = set(tmp) & set(self.avail_trajext)
            if match:
                match = list(match)
                if len(match) > 1: self.log.warn("More than one valid extension found for aligned trajectory step %i: %s. Using first one: %s"%(i, match, match[0]))
                match = match[0]
            else: match = None
            result[i] = match
            
        T.BROWSER.goback()
        return result

    def isAligned(self, stepselection=[]):
        """
        Return True if the trajectory has been already aligned to the reference structure.
        """
        if not stepselection: stepselection=range(1,self.ntrajfiles+1)
        if not isinstance(stepselection, list): stepselection = [stepselection]
        
        exts = self.checkAlignExtension(stepselection) # Check all file extensions in align folder
        if sum([el != None for el in exts.values()]) != len(stepselection):
            return False
        return True

    def isProductionFinished(self, stepselection=[], warn=True):
        """
        Return True if MD production stage has been completed
        
        :arg list stepselection: Selection of steps to be checked. Default is all.
        :arg bool warn: Print a warning message when some file is missing or uncomplete.
        
        :returns: Bool indicating if steps are correctly finished.
        """
        if not stepselection: stepselection=range(1,self.ntrajfiles+1)
        if not isinstance(stepselection, list): stepselection = [stepselection]
        check = self.getChecker(warn=warn)
        return check.checkProduction(stepselection=stepselection)

    def lastCompletedProductionStep(self, startstep=1):
        """
        Check each productio step and return the number of the last incompleted step.
        Useful to track progress or in energy conversion to take volume from last completed step when analyzing
        incompleted runs. Will return zero if no production step is complete.
        
        :arg int startstep: Check production steps starting with this number to reduce function timing when we already know info.
        
        :returns: int
        """
        last = 0
        check = self.getChecker(warn=False)
        for step in range(startstep,self.ntrajfiles+1):
            if check.checkProduction(stepselection=[step]): last+=1
            else: break
        if last == 0: self.log.warn("No output files found in production folder OR no production step is still complete.")
        return last

    def isEquilibrationFinished(self):
        """
        Return True if MD equilibration stage has been completed
        """
        check = self.getChecker(warn=False)
        return check.checkEquilibration()
    
    def isMinimizationFinished(self):
        """
        Return True if MD minimization stage has been completed
        """
        check = self.getChecker(warn=False)
        return check.checkMinimization()
        
    def setNanos(self, nanos):
        "Change number of nanoseconds for current replica"
        self.nanos = int(nanos)
        self.__calcExpectedNTrajFiles()
        self.__calcNumberSnapshots()
        self.log.debug("Changed nanos to %i, expectedtrajfiles: %i"%(self.nanos, self.ntrajfiles))
        if self.__folderscreated: self.write()

    def setProductionSteps(self, steps):
        """
        Change number of production steps per job 
        E.g. at 2fs per step, 500000 steps is 1ns, half nanosecond will be 250000 steps.
        Remember to rewrite MDInput and Queue files if needed.
        """
        self.prod_steps = steps
        self.__calcExpectedNTrajFiles()
        self.write()

    def setOutFileTemplate(self, outfiletemplate):
        """
        Set/Modify output filename template for current replica. All filename templates must include {nano} and {extension}.
        E.g.: md{nano}.{extension}
        """
        self.outfiletemplate = outfiletemplate
        self.log.debug("Changed tempalte to %s"%outfiletemplate)
        if self.__folderscreated: self.write()

    def go(self):
        "Move to replica folder if created"
        if self.__folderscreated or osp.exists(self.path):
            self.__backdir = T.BROWSER.cwd
            T.BROWSER.gotoReplica(self)
            return True
        return False

    def goback(self):
        if self.__backdir: 
            T.BROWSER.chdir(self.__backdir)
            return True
        return False

    def setName(self, name):
        "Set replica name. Adapt logger."
        self.name = name
        self.log = logging.getLogger("Replica (%s)"%name)
        if self.__folderscreated: self.write()

    def createMDInput(self, **kwargs):
        """
        Create MD input config files for the program selected (AMBER or NAMD).
        """
        # Check folder was created and select appropriate writer
        if not self.__folderscreated: self.createFolder()
        if self.mdProgram == 'AMBER':
            from Amber import AmberWriter as writer
        elif self.mdProgram == 'NAMD':
            from NAMD import NAMDWriter as writer
        elif self.mdProgram == 'OPENMM':
            from OpenMM import OpenMMWriter as writer

        else:
            raise ReplicaError, "MD Program not recognized: %s"%self.mdprog

        # Write commands file and replica config input files
        self.go()
        w = writer(self)
        w.writeCommands()
        w.writeReplicaInput()

        self.__mdinput = True
        if self.__folderscreated: self.write()
        self.goback()
        return True

    def createQueueInput(self, queue, **kwargs):
        """
        Write queue input files using templates.

        :arg str queue: Name of the template file defining a queue system. File ``queue_queue_temp.txt`` should exist in users mdmix home directory or package templates.
        """
        if not self.folderscreated(): return False
        import QueueWriting as Q
        self.log.info("Writing Queue %s input files for replica %s"%(queue,self.name))
        self.go()
        queue = Q.QueueInputWriter(queue, **kwargs)
        queue.write(self, **kwargs)
        self.goback()
        return True
    
    def folderscreated(self):
        "Return **True** if replica directory structure is created."
        return self.__folderscreated

    def mdinputwritten(self):
        "Return **True** if replica MD input files are writen."
        return self.__mdinput

    def iscreated(self):
        "Return **True** if replica folder and MD inputs have been written"
        return self.folderscreated() and self.mdinputwritten()

    def createFolder(self, where=False, fixtop=True, **kwargs):
        """
        Create directory tree for current replica. :attr:`path` should have been set
        with :meth:`setPath`. Copy inside the :attr:`top`, :attr:`crd` files if given or generate them from
        a object file :attr:`off`. Create also replica pickle file.

        :arg str where: Path where folder structure should be created. If **None**, use current folder. Path must exist.
        :arg bool fixtop: When saving amber parm7 topology file, remove SCEE and SCNB sections from it. When loading some solvent boxes that include tailored parameters in new Amber programs (> 9.0), they may make the program crash because no SCNB and SCEE scaling factors were specifically given. If the sections are removed, they all take default values.

        Tree structure
        ::
            replica.name/
                replica.top # will be copied if existent or created from off if given
                replica.crd # will be copied if existent
                replica.pdb # will be copied if existent
                replica.name.mrepl
                replica.minfolder/
                replica.eqfolder/
                replica.mdfolder/

        """
        if not self.name:
            raise ReplicaError, "Unnamed replica folder can not be created."
        
        if self.system and self.eqfolder and self.mdfolder:
            import distutils.dir_util as du
            import shutil

            self.log.info("Creating folder structure for replica %s"%self.name)


            pwd = T.BROWSER.cwd
            if where and osp.exists(where): where=osp.abspath(where)
            else: where = T.BROWSER.getcwd()
            T.BROWSER.chdir(where)

            du.create_tree(self.name, [self.minfolder+os.sep, self.eqfolder+os.sep,
                                        self.mdfolder+os.sep],verbose=True)


            # Save top, crd and pdb files for system in current created folder
            basenames = '{0}_{1}'.format(self.system.name, self.name)
            topcrdok = self.system.saveTopCrd(osp.join(self.name, basenames))
            pdbok = self.system.savePDB(osp.join(self.name, basenames+'.pdb'))

            if topcrdok and pdbok:
                self.top = basenames+'.prmtop'
                self.crd = basenames+'.prmcrd'
                self.pdb = basenames+'.pdb'
            else:
                raise ReplicaError, "Error saving system top, crd or pdb files"

            # update replica path and save replica file
            T.BROWSER.chdir(self.name)
            self.setPath(T.BROWSER.getcwd())
            self.__folderscreated = True

            # If fixtop, will remove SCEE and SCNB entries from topology file
            # Needed for some boxes
            if fixtop: self.__fixTopology(self.top)

            # Save reference pdb from pdb
            refpdb = self.system.ref
            self.ref = basenames+'_ref.pdb'
            refpdb.writePdb(self.ref)
            self.log.debug("Created replica: %s"%str(self))
            self.write() # write project file
            T.BROWSER.chdir(pwd)
#            self.log.info("Created folder structure for replica %s"%self.name)
        else:
            raise ReplicaError, "Folder names or replica name not set. Cannot create folders."

    def createAll(self, queue=False, **kwargs):
        """
        Create Folders and MDinput, if queue is given, create also queue input for specified queue name.
        """
        cwd = T.BROWSER.cwd
        self.createFolder(**kwargs)
        self.createMDInput(**kwargs)
        if queue: self.createQueueInput(queue,**kwargs)
        T.BROWSER.chdir(cwd)

    def importData(self, **kwargs):
        """
        Import existing data into current replica. Useful when analyzing data from external simulations
        not run under pyMDMix.

        :keywords:
            Give key=value pairs to be imported where:
            - **value**: absolute path of existing folder containing data to link or a existing file.
            - **key**: repica folder or file where to link data to:

                FILES:
                    - *pdb*:       System PDB
                    - *top*:       System Amber Topology
                    - *crd*:       System Amber Coordinates
                    - *solvent*:   Simulated solvent

                FOLDERS
                    - *mdfolder*:      Production trajectory and output files
                    - *eqfolder*:      Equilibration folder
                    - *alignfolder*:   Aligned trajectory folder
                    - *densityfolder*: Containing density grids
                    - *energyfolder*:  Contraining energy converted grids

            All keys are optional. Only keys assigned will be imported.

            Example::~
                >>> replica = Replica('ETA', name='test')
                >>> replica.importData(pdb='/oldfolder/system.pdb', crd='/oldpath/system.crd', top='/oldpath/system.top',
                >>>                    mdfolder='/oldpath/production', eqfolder='/oldpath/equilibration')

            E.g. /oldfolder/production content will be linked into :attr:`mdfolder` folder.
        """
        import os
        if not self.__folderscreated: self.createFolder()
        self.go()
        self.log.debug("Importing data to %s"%T.BROWSER.getcwd())
        for key, value in kwargs.iteritems():
            # Check attribute exists
            if not hasattr(self, key):
                raise BadAttribute, "%s attribute not in Replica object"%key 
            # File pair
            if osp.isfile(value):
                fname = osp.basename(value)
                if osp.exists(fname): os.remove(fname)
                os.symlink(value, fname)
                setattr(self, key, fname)
                self.log.debug("Linked file %s to replica file %s"%(value, fname))
            # Dir pair
            elif osp.isdir(value):
                # Get replica folder that will contain linked files/subdirs
                destpath = getattr(self, key)
                # Get first level of files/dir structure at source folder
                # and link each file and subfolder
                firstlevelpath = os.walk(value).next()
                main = firstlevelpath[0]
                for f in T.simplifyNestedList(firstlevelpath[1:], []):
                    dest = osp.join(destpath, f)
                    if osp.exists(dest): os.remove(dest)
                    if not osp.exists(destpath): os.mkdir(destpath)
                    ori = osp.join(main, f)
                    self.log.debug("Linking %s to %s"%(ori, dest))
                    os.symlink(ori, dest)
                self.log.debug("Linked folder %s content to replica folder %s"%(value, destpath))
            else:
                raise BadFile, "Not a valid path: %s"%value
        self.write()

    def setPath(self, path, update=True):
        """
        Set replica path. If update is True, update subfolder paths.
        """
        self.path = path
        if update: self.updatePaths()
        elif self.__folderscreated: self.write()

    def updatePaths(self):
        "Update replicaPaths using :attr:`path` as base path"
        self.minpath = self.path+os.sep+self.minfolder
        self.eqpath = self.path+os.sep+self.eqfolder
        self.mdpath = self.path+os.sep+self.mdfolder
        self.alignpath = self.path+os.sep+self.alignfolder
        self.densitypath = self.path+os.sep+self.densityfolder
        self.energypath = self.path+os.sep+self.energyfolder
        name = self.name or 'replica'
        self.replFilePath = self.path+os.sep+name+'.mrepl'
        self.log.debug("Updating replica paths")
        if self.__folderscreated: self.write()
#        self.energypathavg = self.energypath+os.sep+'avg'
#        self.correctedpath = self.replPath+os.sep+self.correctfolder

    def attach(self, object, attachname, desc=''):
        """
        Attach an object to this replica. This will create a pickle of the object
        with a temporary filename and store the pickle file name with a *name* and *description*
        in current replica attribute :attr:`Replica.attached` dictionary.

        :arg any object: Object to pickle and link to the replica
        :arg str attachname: Name to assign to the attachment
        :arg str desc: Description. Optional.
        """
        # Go to replica folder if not there
        T.BROWSER.gotoReplica(self)
        # If attachname present, remove it first
        if self.attached.has_key(attachname): self.dettach(attachname)

        # Obtain temp name, save pickle and store info in self.attached
        attachfile = fname = os.path.basename(tempfile.mktemp(prefix=self.name+'_%s_'%attachname))
        T.dump(object, attachfile)
        self.attached[attachname] = {'file':attachfile, 'desc':desc}

        if self.__folderscreated: self.write()

        # Goback to previous folder
        T.BROWSER.goback()
        
    def dettach(self, attachname):
        """
        Remove attachement with name **attachname**
        """
        attachment = self.attached.get(attachname)
        if attachment:
            T.BROWSER.gotoReplica(self)
            os.remove(attachment['file'])
            self.attached.pop(attachname)
            T.BROWSER.goback()
            if self.__folderscreated: self.write()
            return True
        return False

    def getAttached(self, attachname):
        """
        Load attached object. See ::meth::`attach` for more info.

        :returns: Unpickled object.
        """
        val = self.attached.get(attachname)
        if val:
            fname = val['file']
            T.BROWSER.gotoReplica(self)
            if not osp.exists(fname):
                T.BROWSER.goback()
                raise BadFile, "%s file not found in replica folder."%fname
            obj = T.load(fname)
            T.BROWSER.goback()
            return obj
        else: return False

    def runAlignment(self, ncpus=1, steps=[], waitend=True, **kwargs):
        if kwargs.get('run') and not self.isProductionFinished(steps): raise ReplicaError, "Cannot align replica because production stage is not completed."
        from Align import Align
        Align(self, steps=steps, nthreads=ncpus, waitend=waitend, **kwargs)     
    def runcppDensity(self, ncpus, waitend=True, **kwargs):
        if kwargs.get('run') and not self.isAligned(): raise ReplicaError, "Cannot calculate density because alignment is not completed."
        from Actions.Density import DensityGrids, cppDensity
        samplegrid = DensityGrids(self, probeselection=kwargs['probelist'], outprefix=kwargs['outprefix'], includeCOM=kwargs['includeCOM'],
                            onlyCOM=kwargs['onlyCOM'], stepselection=kwargs['nanosel'], reference=kwargs['ref'])
        samplegrid.prepareGrids()
        origin = samplegrid.container.origin
        dimensions = samplegrid.container.shape
        cppDensity(self, nthreads=ncpus, waitend=waitend,  griddimensions = dimensions, gridorigin=origin,**kwargs)
    def calcEnergy(self, **kwargs):
        "Convert density to energies. Give in ``\*\*kwargs`` all parameters to :meth:`Energy.EnergyConversion.convert`."
        from Energy import EnergyConversion
        econv = EnergyConversion()
        econv.convert(self, **kwargs)

    def write(self):
        "Save object __dict__ to pickled file."
        self.log.debug("Writing pickle")
        with FileLock(self.replFilePath) as lock:
            d = self.__dict__.copy()
            del d['log']
            del d['grids']
            T.dump(d, self.replFilePath)

    def load(self, replfile=None):
        "Load existing project from pickled file"
        f = replfile or self.replFilePath
        if not osp.exists(f): raise BadFile, "File %s not found."%f
        with FileLock(f) as lock:
            d = T.load(f)
            d['log'] = logging.getLogger("Replica (%s)"%d['name'])
            d['grids'] = None
            if not d.has_key('prod_steps'): d['prod_steps'] = d['nvt_prod_steps']
            self.__dict__.update(d)


def renameReplicaList(repls):
    """
    Set Name accoding to solvent and number of replicas. E.g. if the list contains 3 replicas with ethanol and 2 with water,
    the names will be: ``['ETA_1','ETA_2','ETA_3','WAT_1','WAT_2']``. This function does not return anything. Names are changed in 
    the replicas objects.
    
    :arg list repls: List of :class:`Replica` objects to be renamed
    """
    import collections
    solvs = collections.Counter()
    for r in repls:
        solvs[r.solvent] += 1
        n = '%s_%i'%(r.solvent, solvs[r.solvent])
        r.setName(n)
        
def loadReplica(replicafile=None):
    """
    Load existing replica. 
    If *replicafile* path is not given, will try to load any ``\*.mrepl`` file present in current folder.

    :returns: :class:`~Replicas.Replica` object
    """
    if not replicafile:
        import glob
        files = glob.glob('*.mrepl')
        if len(files) > 1:
            raise ReplicaError,"More than one project file in current folder. Please remove the invald one."
        elif not files:
            raise ReplicaError,"No file found with extension *.mrepl in current folder and no path was given."
        replicafile = files[0]
    return Replica(fromfile=replicafile)


import Biskit.test as BT

class Test(BT.BiskitTest):
    """Test"""

    def test_NewReplica(self):
        """Create new replica test"""
        self.r1 = Replica()
        self.r2 = Replica(name='testReplica', )

    def test_PepReplica(self):
        rfile = T.testRoot('pep','pep_mdmix','MD','')
#
#    def cleanUp(self):
#        T.tryRemove( self.f_out )

if __name__ == '__main__':
    BT.localTest()
