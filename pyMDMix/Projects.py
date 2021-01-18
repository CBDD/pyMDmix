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

__author__="dalvarez"
__date__ ="$20-ene-2014 22:52:30$"

import os
import os.path as osp
import logging

import Systems
from MDSettings import MDSettings
import settings as S
import tools as T
from Replicas import Replica
from Structures import FileLock


class ProjectError(Exception):
    pass

class BadAttribute(ProjectError):
    pass

class BadFile(ProjectError):
    pass

class Project(object):
    def __init__(self, name='mdmix_project', fromfile=None, replPaths='MD', **kwargs):
        """
        Constructor.
        """
        if fromfile:
            if osp.exists(fromfile):
                self.load(fromfile)
                self.projFilePath = osp.abspath(fromfile)
                # Update paths if projectPath and projFilePath do not match
                dir = osp.split(self.projFilePath)[0]
                if self.projectPath != dir:
                    _updated=True
                    self.updatePath(dir)
                else:
                    _updated=False
                T.BROWSER.setHome(self.projectPath)
                T.BROWSER.goHome()
                self.fetchSystems()
                self.fetchReplicas()
                if _updated: self.updateReplicaPaths()
            else:
                raise BadFile, "Project File %s not found"%fromfile
        else:
            import collections
            
            self.version = S.VERSION
            self.projName = name                 #: Name of the project
            self.projFileName = name+'.mproj'    #: Name of the file that will be generated to save all project information
            self.projectPath = osp.join(os.getcwd(), name)
            self.projFilePath = osp.join(self.projectPath, self.projFileName)

            self.log = logging.getLogger("Project %s"%name)
            self.__folderscreated = False
            
            # Path to store all replicas / and look for them there
            self.replPaths = replPaths

            self.replicas = {}
            self.systems = {}
            self.replicagroups = {}
            self.solventReplicas = {}
            self.solventCounter = collections.Counter()


    def __getitem__(self, replicaname):
        "Return replica by name"
        return self.replicas.get(replicaname)
    
    def __str__(self):
        "Do some pretty printing of most important attributes"
        sysinfo = ''
        replinfo= ''
        s ="\n"+"-"*30
        s+= "\nPROJECT %s INFO\n"%self.projName
        s+="-"*30
        s+='\n'
        s+= 'SYSTEMS:\n'
        s+= '========\n'
        for sys in self.systems.values(): sysinfo+=str(sys)
        if not sysinfo: sysinfo='NO SYSTEM FOUND\n'
        s+= sysinfo
        s+= '\n'
        s+= 'REPLICAS\n'
        s+= '========\n'
        for repl in self.replicas.values(): replinfo+=repl.name+'\n'
        if not replinfo: replinfo='NO SYSTEM FOUND\n'
        s+= replinfo
        s+='\n'
        if self.replicagroups:
            s+='GROUPS\n'
            s+='======\n'
            for n, g in self.replicagroups.iteritems():
                s+="%s: %s\n"%(n,', '.join(g))
        s+='\n'
        s+="-"*30
        s+='\n'
        return s

    def __repr__(self):
        return "pyMDMix Project: %s"%(self.projName)

    def __createSystem(self, system):
        "Save a reference pdb for the system and write system file in current folder"
        system.write()
        self.log.info("Writing system reference file: %s"%system.name+'_ref.pdb')
        system.ref.writePdb(system.name+'_ref.pdb')

    def longdesc(self):
        "Do some pretty printing of most important attributes"
        sysinfo = ''
        replinfo= ''
        s ="\n"+"-"*30
        s+= "\nPROJECT %s INFO\n"%self.projName
        s+="-"*30
        s+='\n'
        s+= 'SYSTEMS:\n'
        s+= '========\n'
        for sys in self.systems.values(): sysinfo+=str(sys)
        if not sysinfo: sysinfo='NO SYSTEM FOUND\n'
        s+= sysinfo
        s+= '\n'
        s+= 'REPLICAS\n'
        s+= '========\n'
        for repl in self.replicas.values(): replinfo+=repl.desc()
        if not replinfo: replinfo='NO SYSTEM FOUND\n'
        s+= replinfo
        s+='\n'
        if self.replicagroups:
            s+='GROUPS\n'
            s+='======\n'
            for n, g in self.replicagroups.iteritems():
                s+="%s: %s\n"%(n,', '.join(g))
        s+='\n'
        s+="-"*30
        s+='\n'
        return s

    def addNewSystems(self, systems):
        """
        Add systems to current project.

        :arg list systems: List of :class:`Systems.System` instances. Make sure System names are unique.
        """
        if not self.__folderscreated: raise ProjectError, "Can not add Systems if Project folder is not created"
        if not isinstance(systems, list): systems = [systems]
        T.BROWSER.goHome()
        for s in systems:
            if not isinstance(s, Systems.System):
                self.log.warn("Invalid type. system should be Systems.System instance.")
            if s.name not in self.systems.keys():
                self.__createSystem(s)
            else:
                self.log.warn("System %s already exists. Skipping..."%s.name)
                continue
        T.BROWSER.goback()
        self.fetchSystems()

    def applyReplicas(self, action, replicalist=[], **kwargs):
        "General method to call method *action* with arguments in **kwargs** to all replicas in replicalist. If list is empty, apply to all known replicas"
        if not replicalist: replicalist = self.fetchReplicas().values()
        process = []
        
        # Set executor number of threads here and
        # don't let replicas to alter it by setting ncpus to false
        ncpus = kwargs.get('ncpus')
        kwargs['ncpus'] = False
        ncpus = ncpus or 1
        T.EXECUTOR.changeNthreads(ncpus)
        
        for r in replicalist:
            if not isinstance(r, Replica): 
                self.log.warn("Unexpected type %s. Expected Replica type."%type(r))
                continue
                
            try: meth = getattr(r, action)
            except: raise ProjectError, "Replica do not have method with name %s"%action
            
            import multiprocessing as multi
            self.log.debug("Submitting replica %s action %s"%(r.name, action))
            process.append(multi.Process(target=meth, kwargs=kwargs))
        
        # Start jobs
        [p.start() for p in process]
        [p.join() for p in process]

    def getSystem(self, sysname):
        "Return System object with name *sysname*. Must exists in current project folder."
        return self.fetchSystems().get(sysname)
    
    def getReplica(self, replicaname):
        "Return System object with name *sysname*. Must exists in current project folder."
        return self.fetchReplicas().get(replicaname)

    def createReplicas(self, systemname, settings, **kwargs):
        """
        Create Replicas of system with name *systemname* using MDSettings in settings list.
        Replicas will be automatically named according to the solvent and number of replicas for same solvent (e.g. ETA_1, ETA_2, WAT_1, etc.)

        :arg str systemname: Name of system to use. Must exist in current project (add it before calling this method with :meth:`addNewSystems`)
        :arg list setttings: List of :class:`MDSettings` objects.
        """
        if not self.__folderscreated: raise ProjectError, "Can not create Replicas if Project folder is not created"
        if not systemname in self.systems.keys(): raise BadAttribute, "System name %s not found in current project"%systemname
        if not isinstance(settings, list): settings = [settings]

        # Check types
        for sets in settings:
            if not isinstance(sets, MDSettings): raise BadAttribute, "Expected MDSettings type, but %s given"%(type(sets))

        # Build replicas and add names
        T.BROWSER.goHome()
        self.log.info("Creating replicas for system %s..."%systemname)
        repls = self.getSystem(systemname) + settings
        if not isinstance(repls, list): repls=[repls]
        self.__nameNewReplicas(repls)

        # Create replicas inside self.replPaths
        T.BROWSER.chdir(self.replPaths)
        [r.createAll(**kwargs) for r in repls]
        T.BROWSER.goback()

        self.fetchReplicas()
        self.log.info("DONE")

    def createQueueInputs(self, queue, **kwargs):
        "Create Queue input for all replicas in the project"
        self.applyReplicas('createQueueInput',queue=queue,**kwargs)

    def __nameNewReplicas(self, repls):
        "Put automatic names to replicas according to the known solvents"
        self.fetchReplicas() # Set to zero and count existing solvents
        c = self.solventCounter.copy()
        if not isinstance(repls, list): repls = [repls]
        for r in repls:
            c[r.solvent] += 1
            r.setName("{0}_{1}".format(r.solvent, c[r.solvent]))

    def fetchReplicas(self):
        "Fetch all replicas pickled files inside `ProjectFolder/replPaths`"
        repls = []
        cwd = T.BROWSER.cwd
        if cwd != T.BROWSER.home: T.BROWSER.goHome()
        T.BROWSER.chdir(self.replPaths)
        for root, dir, files in os.walk(os.curdir):
            for f in files:
                if f.endswith('mrepl'): 
                    r = Replica(fromfile=osp.join(root, f))
                    repls.append(r)
        T.BROWSER.chdir(cwd)
        
        # Check no duplicate names exist
        names = []
        duplicate = []
        for i,r in enumerate(repls): 
            if not r.name in names: names.append(r.name)
            else:
                # Duplicate name
                self.log.warning("Detected a duplicated replica with name %s in path %s"%(r.name,r.path))
                self.log.warning("Skipping this replica. Make sure the replica file (*.mrepl) are located in the right directory and no duplicate names exist")
                duplicate.append(i)
                
        self.replicas = dict([(r.name, r) for i,r in enumerate(repls) if not i in duplicate])
        self.solventCounter.clear()
        for r in repls: self.solventCounter[r.solvent] += 1
        return self.replicas

    def fetchSystems(self):
        "Fetch pickled systems at project home folder"
        cwd = T.BROWSER.cwd
        change = cwd != T.BROWSER.home
        if change: T.BROWSER.goHome()
        systems=[]
        for root, dir, files in os.walk(os.curdir):
            for f in files:
                if f.endswith('msys'): systems.append(Systems.System(fromfile=osp.join(root, f)))
        if change: T.BROWSER.chdir(cwd)
        self.systems = dict([(r.name, r) for r in systems])
        return self.systems

    def changeFilePath(self, projfile):
        self.projFilePath = osp.abspath(projfile)

    def updatePath(self, path=None):
        """
        Update main project folder path with *path* or current working directory if *path* not given.
        All project replica's paths will be also updated.

        This method to work requires that the expected project file (:attr:`Project.projFileName`) is placed inside the current working directory.
        """
        path = path or T.BROWSER.getcwd()
            
        # Check expected project file exists inside *path*
        if not osp.exists(osp.join(path, self.projFileName)):
            raise ProjectError, "Trying to set project path to folder %s that does not contain expected project file %s"%(path, self.projFileName)

        self.projectPath = path
        self.projFilePath = osp.join(self.projectPath, self.projFileName)
        self.updateReplicaPaths()

        # Finally update BROWSER home
        T.BROWSER.setHome()
        self.write()

    def updateReplicaPaths(self):
        "Update replica paths with current project main path information"
        # Update replica paths
        for r in self.replicas.values():
            r.setPath(osp.join(self.projectPath,self.replPaths,r.name))
            r.write()

    def createProjectFolder(self):
        "Create folder structure for current project name and update paths"
        import distutils.dir_util as du
        if not self.__folderscreated:
            self.log.info("Creating project folder")
            self.__folderscreated= True
            du.create_tree(self.projName,[self.replPaths+os.sep], verbose=True)
            T.BROWSER.chdir(self.projName)
            self.projFilePath = self.projFileName
            self.write()
            self.updatePath()

    def listGroups(self):
        return self.replicagroups.keys()

    def createGroup(self, groupname, replicanames):
        """
        Create a group of replicas for joint analysis

        :arg str groupname: Name to identify the group
        :arg list replicanames: Replica names to add to group
        """
        g = []
        for r in replicanames:
            if self.replicas.get(r): g.append(r)
            else: self.log.warn("Replica name %s not in current project"%r)
        self.replicagroups[groupname] = g
        self.write()

    def getGroup(self, groupname):
        """
        Get a list with replicas belonging to the group *groupname*

        :arg str groupname: Name of the group to retrieve

        :returns: List with :class:`~Replicas.Replica` instances or False if group does not exists
        """
        groupnames = self.replicagroups.get(groupname)
        if not groupname: return False
        
        return [self.replicas.get(r) for r in groupnames]

    def removeGroup(self, groupname):
        """
        Remove group.
        :arg str groupname: group name to be removed from project
        """
        if self.replicagroups.has_key(groupname):
            self.replicagroups.pop(groupname)
            self.write()
            return True
        return False

    def extendSimulations(self, replicalist, nanos):
        """
        Extend simulation for existing replicas.
        
        :arg list replicalist: List of :class:`Replicas.Replica` instances or strings to extend.
        :arg int nanos: Number of nanoseconds to add to current replica nanoseconds. 
        """
        if not isinstance(replicalist, list): replicalist = [replicalist]
        for r in replicalist:
            if isinstance(r, str): r = self.getReplica(r)
            if not isinstance(r, Replica): raise BadAttribute, "Expected argument of type Replica, not %s"%(type(r))
            self.log.info("Extending replica %s simulation %i nanoseconds to %s"%(r.name, nanos, r.nanos+nanos))
            r.setNanos(r.nanos+nanos)
            r.createMDInput()
        return
    
    def alignReplicas(self, replicalist, ncpus=1, steps=[], waitend=True, **kwargs):
        """
        Run alignment process on replicas in replicalist
        
        :arg list replicalist: List of replica instances
        :arg int ncpus: Threads to start
        :arg list steps: Steps to align. If emtpy: analyze all.
        :arg bool waitend: Wait until all alignment process is done before exiting the method        
        """
        self.log.info("Running alignment for replicas %s"%replicalist)
        self.applyReplicas('runAlignment',replicalist, steps=steps, ncpus=ncpus, **kwargs)
        if waitend: T.EXECUTOR.waitJobCompletion()
        
    def calc_cppdensityReplicas(self, replicalist, ncpus=1, waitend=True, **kwargs):
        """
        Run density calculation process on replicas in replicalist
        
        :arg list replicalist: List of replica instances
        :arg int ncpus: Threads to start (normal will be 1)
        :arg list steps: Steps to use for the density calculation. If emtpy: analyze all.
        :arg bool waitend: Wait until all alignment process is done before exiting the method        
        """
        self.log.info("Running cpp density calculation for replicas %s"%replicalist)
        self.applyReplicas('runcppDensity',replicalist, ncpus=ncpus, **kwargs)
        if waitend: T.EXECUTOR.waitJobCompletion()
    
    def write(self):
        "Save object __dict__ to pickled file."
        with FileLock(self.projFilePath) as lock:
            d = self.__dict__.copy()
            del d['log']
            d['systems'] = {}    # Empty systems and replicas
            d['replicas'] = {}
            T.dump(d, self.projFilePath)

    def load(self, projfile=None):
        "Load existing project from pickled file"
        f = projfile or self.projFilePath
        with FileLock(f) as lock:
            d = T.load(f)
            d['log'] = logging.getLogger("Project %s"%d['projName'])
            self.__dict__.update(d)

#    def importReplica(self, replica):
#        """
#        Import existing replica object. Will copy all information in :attr:`replPaths` directory.
#
#        :arg str replica: Name of the replica to import
#        :arg str solvent: Solvent used.
#
#        :replicaargs:
#            - top, crd, pdb: existing file paths to be imported
#            - mdfolder, eqfolder, minfolder: existing paths
#            - other configuration options desrcibing imported replica (restrMode, mdProgram, etc. see :class:`~Replicas.Replica`)
#
#        """
#        from Replicas import Replica
#        r = Replica(name=name, solvent=solvent)
#        self.addNewReplica(r)
#        T.BROWSER.goMD()
#        r.createFolder()
#        r.importData(**replicaargs)
#        T.BROWSER.goHome()
#        self.updateReplica(r)
#        self.write()

#    def updateReplica(self, replica):
#        """Update replica. When a replica is modified, it is recomended to run this method to update project file information.
#        :arg replica: replica instance to be updated in current project
#        :type replica: :class:`~Replicas.Replica`
#        """
#        if self.replicas.has_key(replica.name): self.replicas[replica.name] = replica
#        else: self.log.warn("Replica with name %s not in current project. Call :meth:`Project.addNewReplica` method")
#        self.write()
#


### Functions to create project from configfile
### and load existing one
def createProject(projectConfigFile, name):
    """
    Auxiliary function to build a project from a project configuration file (PCF)

    :arg str projectConfigFile: path to config file containing SYSTEM and MDSETTINGS sections
    :arg str name: Name to project to be created
    """
    from Systems import parseSystemConfigFile
    from MDSettings import parseSettingsConfigFile
    if not osp.exists(projectConfigFile): raise BadFile, "File %s not found."%projectConfigFile
    print "Parsing System information..."
    sys = parseSystemConfigFile(projectConfigFile)
    print "Parsing md settings for replica creation..."
    sets = parseSettingsConfigFile(projectConfigFile)
    print "Creating project %s"%name
    project = Project(name=name)
    project.createProjectFolder()
    project.addNewSystems(sys)
    project.createReplicas(sys.name, sets)
    return project

def loadProject(projectfile=None):
    """
    Load existing project.
    If projecfile path is not given, will try to load any \*.mproj file present in current folder.
    """
    if not projectfile:
        import glob
        files = glob.glob('*.mproj')
        if not files: raise ProjectError,"No project file found in current folder. Make sure you are in a pyMDMix project folder."
        if len(files) > 1:
            raise ProjectError,"More than one project file in current folder. Please remove the invald one."
        projectfile = files[0]
    return Project(fromfile=projectfile)

def returnMDMixProject(parserargs):
    if parserargs.debug: level='DEBUG'
    else: level='INFO'
    import pyMDMix
    pyMDMix.setLogger(level=level)
    try:
        p = pyMDMix.loadProject()
        return p
    except:
        return False

def returnMDMixProjectOrFail(parserargs):
    #When command is different to CREATE PROJECT or INFO, this program should be executed in project folder
    #Let's try to load a project or exit
    p = returnMDMixProject(parserargs)
    from pyMDMix import MDMixError
    if not p: raise MDMixError, 'No project file found in current folder. Make sure you are in a pyMDMix project folder.'
    return p


###TESTING
import Biskit.test as BT

class Test(BT.BiskitTest):
    """Test"""
    def test_CREATEpepProject(self):
        """Create complete test project"""
        import shutil
        print T.testRoot()
        cfg = T.testRoot('pep','pep_amber_mdmix.cfg')
        off = T.testRoot('pep','pep.off')
        self.f_out = T.tempDir()
        T.BROWSER.chdir(self.f_out)
        shutil.copy(cfg, 'pep_amber_mdmix.cfg')
        shutil.copy(off, 'pep.off')
        proj = createProject('pep_amber_mdmix.cfg', 'pep_test_project')
        self.f_out += 'pep_test_project'

    def cleanUp(self):
        T.tryRemove( self.f_out, tree=1 )

if __name__ == '__main__':
    BT.localTest()
