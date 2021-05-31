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
Defines two classes to work and store structure and amber parameters for a particular biomolecule and other related information.

At the end, to define the system in the MD, we need an amber parameter file (PRMTOP) and a coordinates file (PRMCRD)
idependently of the simulation program (AMBER or NAMD) chosen.

Two different classes have been implemented:
    - :class:`Systems.System`
    - :class:`Systems.SolvatedSystem`
    
"""

__author__="dalvarez"
__date__ ="$14-mar-2014 20:13:59$"

import os
import os.path as osp
import logging

import Solvents
import tools as T
import settings as S
from Structures import FileLock

class SystemError(Exception):
    pass

class BadFile(SystemError):
    pass

class System(object):
    """
    """
    def __init__(self, name=None, fromfile=None, amberPDB=None, amberOFF=None, unitName=None,
                    extraResList=[], FF=S.DEF_AMBER_FF, **kwargs):
        """
        Initialize a system object.

        :arg str name: Identified name for the system. If None, a random name will be generated ``mdmix_system_XXX``.
        :arg str fromfile: Path to already pickled System object to be loaded.
        :arg str amberPDB: Path to PDB file already prepared to amber standards and ready for topology generation.
        :arg str amberOFF: Path to Amber Object File (OFF) containing the system already parameterized and stored as a Leap Unit.
        :arg str uniName: Name of the unit containing the system in the OFF file. If not given, the first name found in the file will be used.
        :arg list extraResList: Non standard residue names to consider as part of the macromolecule.
        :arg list FF: Forcefields used to parameterize the system in the OFF file or in the preparation of the PDB (atom names should coincide). If non-standard residues are present, any frcmod file used should also be given for correctly parameterizing them.

        """
        if fromfile:
            self.load(fromfile)
        else:
            if not name:
                # Generate random string
                import random, string
                randstr = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(3))
                name = 'mdmix_system_'+randstr
            
            self.name = name
            self.extraResList = extraResList
            self.sysFilePath = self.name+'.msys'

            # Try to find path of extra FF files if given
            self.FF = []
            for f in FF:
                if osp.exists(f): self.FF.append(osp.abspath(f))
                else: self.FF.append(f) # it might be a parm file located in amber standard folders (to be checked later)

            # Files to be generated
            self.amberOFF = None #: :class:`~OFFManager.OFFManager` object to store the amber object file with the system. To be set with :meth:`setOFF`
            self.unitName = None #: To be set with :meth:`setOFF`
            self.ref = None #: PDB of the macromolecule, generated from :attr:`amberOFF` and :attr:`extraResList`

            # Internal control attrs
            self.create = False
            self.__filecreated = False

            # OFF constructor has priority
            if amberOFF:
                self.setOFF(amberOFF, unitName)
            elif amberPDB:
                self.setAmberPDB(amberPDB, **kwargs)
            else:
                pass # Wait for incorporation of data later

    def __repr__(self):
        return "%s System"%self.name

    def __str__(self):
        s="SYSTEM:%s\n"%self.name
        if self.FF: s+="\tFF:%s\n"%self.FF
        if self.extraResList: s+="\tExtraRes:%s\n"%self.extraResList
        return s

    def __add__(self, other):
        from MDSettings import MDSettings
        from Replicas import Replica

        if isinstance(other, MDSettings): other = [other]
        if isinstance(other, list):
            # Check elements are MDSettings objects
            out = []
            for el in other:
                if isinstance(el, MDSettings):
                    sys = self.solvate(el.solvent)
                    out.append(Replica(sys, el))

            if len(out) == 1: return out[0]
            return out

#    def __getstate__(self):
#        d = self.__dict__.copy()
#        del d['log']
#        return d
#
#    def __setstate__(self, d):
#        d['log'] = logging.getLogger("System (%s)"%d['name'])
#        self.__dict__.update(d)

    def __initCreate(self):
        if not self.create:
            from Amber import AmberCreateSystem
            self.create = AmberCreateSystem(FFlist=self.FF, informative=False)
            if self.amberOFF:
                tmpoff = self.amberOFF.writeTmp()
                self.create.loadOff(tmpoff) # will write to temporary file location
                self.amberOFF.cleanTmp()

    def __cleanCreate(self):
        if self.create:
            self.create.leap.close()
            self.create = None
            #if self.amberOFF:
            #    self.amberOFF.cleanTmp()

    def setAmberPDB(self, amberPDB, FF=[], **kwargs):
        """
        Set up an amber object file departing from a PDB file already prepared to amber standards and ready for topology generation. Despite a good PDB file is expected, some
        cleaning process will take place:
            - Cap N and C terminus
            - Rename HIS residues according to protonation
            - Rename CYS to CYX if needed, and disulfide bridges built in tLeap
            - Removal of all hydrogen atoms to let tLeap add them back.
            
        :arg str amberPDB: Path to PDB file complying with amber specifications.
        :arg list FF: Forcefield names and modification files that might be needed to correctly interpret the PDB by tLeap.

        :kwargs: Check :meth:`~AutoPrepare.AmberPDBCleaner.cleanPDB` method for extra parameters that are accepted to control automatic cleaning of the PDB.

        :return: Nothing will be returned, the object file generated will be stored in :attr:`amberOFF` as an :class:`~OFFManager.OFFManager` object.
        """
        if not osp.exists(amberPDB):
            raise BadFile, "PDB file %s does not exists"

        FF = FF or self.FF

        self.__initCreate()
        self.create.createOFF(self.name, amberPDB, extraff=FF, **kwargs)
        outoff = self.name+'.lib'
        self.setOFF(outoff ,'sys')


    def setOFF(self, amberOFF, unitname=None):
        """
        Save Amber Object file as :class:`OFFManager.OFFManager` object in :attr:`amberOFF`. :attr:`unitName` will be set.
        
        :arg str amberOFF: Path to amber object file.
        :arg str unitname: Name of the unit we should use in the future. If not given, automatically take the first unit found in the file.
        """
        from OFFManager import OFFManager
        
        # Special case: test
        # Grab off from testing directory
        if amberOFF.lower() == 'test': 
            print "WORKING WITH A TEST SYSTEM!"
            amberOFF = T.testRoot('pep','pep.off')
            self.name = 'testsystem'
            self.sysFilePath = self.name+'.msys'
        
        self.amberOFF = OFFManager(amberOFF)
        if not unitname:
            units = self.amberOFF.getUnits()
            unitname =  units[0]
        self.unitName = unitname
        self.__createRef()

    def solvate(self, solvent, suffix=None, tmp=True):
        """
        Solvate system in OFF using solvent *solvent*. A PRMTOP and PRMCRD files will be saved and paths stored to :attr:`top` and :attr:`crd`.

        :arg solvent: Solvent Instance or Name to fetch the Solvent Instance in the database.
        :type: str or :class:`~Solvents.Solvent`
        :arg str suffix: Suffix for output files. Prefix will be the system name. E.g.: ``systemname_suffix.prmtop`` and ``systemname_suffix.prmcrd``.
        :arg bool tmp: Work in temporary folder.
        """
        if not self.amberOFF:
            raise SystemError, "Can not solvate if no Amber Object File is assigned."
        
        if isinstance(solvent, str): solvent = Solvents.getSolvent(solvent)
        if not solvent:
            raise SystemError, 'Invalid solvent instance or name.'

        # Build a loger just for this method
        log = logging.getLogger("SystemLogger")
        log.info("Solvating %s with solvent mixture %s"%(self.name, solvent.name))

        suffix = suffix or solvent.name
        name = '{0}_{1}'.format(self.name, suffix)
        prmtop = name+'.prmtop'
        prmcrd = name+'.prmcrd'

        if tmp:
            if isinstance(tmp, str): path = tmp
            else: path = T.tempDir()
            prmtop = osp.join(path, prmtop)
            prmcrd = osp.join(path, prmcrd)

        # Initiate AmberCreateSystem with loaded AmberOFF
        # Solvate, neutralize
        self.__initCreate()
        self.create.solvateOrganic(self.unitName, solvent=solvent)  # It will work even if its water
        self.create.saveAmberParm(self.unitName, prmtop, prmcrd)
        self.__cleanCreate()

        s = SolvatedSystem(name, prmtop, prmcrd, solvent=solvent.name, ref=self.ref)
        return s

    def writeOFF(self, fname):
        "Write Object File to disk file with name *fname*"
        if self.amberOFF:
            return self.amberOFF.write(fname)

    def __createRef(self):
        """
        Create reference pdb and store it in self.ref
        """
        import Biskit as bi
        # Initiate AmberCreateSystem with loaded AmberOFF
        self.__initCreate()
        self.create.saveAmberParm(self.unitName, 'tmp.top', 'tmp.crd')
        self.create.ambpdb('tmp.top','tmp.crd','tmp.pdb')
        self.ref = bi.PDBModel('tmp.pdb')
        self.ref.source.set_string("System %s reference"%self.name)
        # Remove tmp files
        [os.remove('tmp.'+ext) for ext in ('top','crd','pdb')]
        self.__cleanCreate()

    def copy(self):
        """
        Return a copy of the system object so it can be altered.
        sysFilePath will be reset and __creted
        """
        import copy
        return copy.deepcopy(self)

    def write(self, outfile=None):
        "Save object __dict__ to pickled file."
        outfile = outfile or self.sysFilePath
        with FileLock(outfile) as lock:
            d = self.__dict__.copy()
            d['create'] = None #Dont pickle Amber class if loaded
#            del d['log']
            T.dump(d, outfile)
            self.sysFilePath = osp.abspath(outfile)
            self.__sysfile = True

    def load(self, sysfile=None):
        "Load existing project from pickled file"
        f = sysfile or self.sysFilePath
        if not osp.exists(f): raise BadFile, "File %s not found."%f
        with FileLock(f) as lock:
            d = T.load(f)
#            d['log'] = logging.getLogger("System (%s)"%d['name'])
            self.__dict__.update(d)


class SolvatedSystem(System):
    """
    Subclass of System. Contains an already solvated system and thus, should be instantiated using a PARMTOP and PARMCRD file.
    """
    def __init__(self, name=None, top=None, crd=None, pdb=None, solvent=None, ref=None, **kwargs):
        """
        Solvated System constructor

        :arg str name: Name of solvated system.
        :arg str top: File path to amber prmtop
        :arg str crd: File path to amber prmcrd
        :arg str pdb: File path to pdb file created from prmtop and prmcrd files. If none, will be automatically created
        :arg str solvent: Solvent used in solvating the System
        :arg ref: PDBModel that will be used as reference for trajectory alignment, common to all solvated systems coming from same System.
        :type ref: :class:`Biskit.PDBModel`
        """
        System.__init__(self, name=name, **kwargs)
        self.top = None #: Amber PRMTOP file. To be set from :meth:`setTOPCRD` or from :attr:`setOFF
        self.crd = None #: Amber PRMCRD file`. Same as PRMTOP
        self.pdb = None #: PDB created from PRMTOP and PRMCRD
        self.ref = ref  #: Reference PDB, PDB before solvation that will be used as reference for al solvated systems

        # Store temprorary files paths
        self.tmp_top = None
        self.tmp_crd = None
        self.tmp_pdb = None
        self.solvent = solvent #: If System contains solvated system, identify the solvent
        
        if top and crd: self.setTopCrd(top, crd)       
        
    def __repr__(self):
        return "%s SolvatedSystem"%self.name

    def __add__(self, other):
        from MDSettings import MDSettings
        from Replicas import Replica
        
        if isinstance(other, MDSettings): other = [other]
        if isinstance(other, list):
            # Check elements are MDSettings objects
            out = []
            for el in other:
                if isinstance(el, MDSettings):
                   out.append(Replica(self, el))

            if len(out) == 1: return out[0]
            return out

    def setTopCrd(self, top, crd):
        """
        Set PRMTOP and PRMCRD for the system. In this case, the system might be already solvated.

        :arg str top: path to PRMTOP file
        :arg str crd: path to PRMCRD file
        """
        if not osp.exists(top): raise BadFile, "File %s not found."%top
        if not osp.exists(crd): raise BadFile, "File %s not found."%crd
        self.top = open(top,'r').read()
        self.crd = open(crd,'r').read()
        self.setPDBfromTOPCRD()

        if self.ref:
            from PDB import SolvatedPDB
            self.ref = SolvatedPDB(self.ref)
        else:
            # Reference not given, create one from self.pdb
            self.ref = self.getSolvatedPDB().getSolute()

        if not self.solvent:
            self.solvent = self.getSolvatedPDB().solvent
            
    def saveTopCrd(self, prefix):
        """
        Save top disk file system topology and coordinates (amber PRMTOP and PRMCRD).
        :arg str prefix: Prefix name to output files. Extensions will be ``.prmtop`` and ``.prmcrd``.
        """
        if self.top and self.crd:
            open(prefix+'.prmtop','w').write(self.top)
            open(prefix+'.prmcrd','w').write(self.crd)
            return osp.exists(prefix+'.prmtop') and osp.exists(prefix+'.prmcrd')
        return False

    def savePDB(self, fname):
        """
        Save PDB created from topology and coordinates to disk file.
        """
        if self.pdb:
            open(fname, 'w').write(self.pdb)
            return osp.exists(fname)
        return False

    def getSolvatedPDB(self):
        """
        Return a :class:`~PDB.SolvatedPDB` object from the PDB file generated from TOP and CRD
        """
        from PDB import SolvatedPDB
        pdb = self.getTmpPdbFile()
        o = SolvatedPDB(pdb, self.extraResList)
        self.cleanTmp()
        return o

    def getTmpTopCrdFiles(self):
        """
        Save Top and Crd as temporary files

        :return: tuple (prmtop tempfile path, prmcrd tempfile path)
        """
        tmp = T.tempfile.mktemp()
        self.tmp_top = tmp+'.prmtop'
        self.tmp_crd = tmp+'.prmcrd'
        if self.saveTopCrd(tmp):
            return self.tmp_top, self.tmp_crd

    def getTmpPdbFile(self):
        """
        Save pdb in temporary file and return file path
        """
        if self.pdb:
            tmp = T.tempfile.mktemp()+'.pdb'
            self.tmp_pdb = tmp
            open(self.tmp_pdb, 'w').write(self.pdb)
            return self.tmp_pdb

    def cleanTmp(self):
        "Remove temporary files"
        if self.tmp_top: os.remove(self.tmp_top)
        if self.tmp_crd: os.remove(self.tmp_crd)
        if self.tmp_pdb: os.remove(self.tmp_pdb)
        self.tmp_top = None
        self.tmp_crd = None
        self.tmp_pdb = None
    
    def setPDBfromTOPCRD(self):
        "Save a PDB file from the TOP and CRD files in attributes."
        import time
        from Amber import AmberCreateSystem

        self.getTmpTopCrdFiles()
        self.tmp_pdb = tmp = T.tempfile.mktemp()+'.pdb'

        create = AmberCreateSystem(informative=False)
        create.ambpdb(self.tmp_top, self.tmp_crd, self.tmp_pdb)
        time.sleep(1) # delay to allow ampdb work
        
        self.pdb = open(self.tmp_pdb,'r').read()
        self.cleanTmp()

    def solvate(self, **kwargs):
        pass

        
def loadSystem(systemfile=None):
    """
    Load existing system. 
    If *systemfile* path is not given, will try to load any ``\*.mrepl`` file present in current folder.

    :returns: :class:`~Systems.System` object
    """
    if not systemfile:
        import glob
        files = glob.glob('*.msys')
        if len(files) > 1:
            raise SystemError,"More than one system file in current folder. Please give as argument the file to load."
        elif not files:
            raise SystemError,"No file found with extension *.msys in current folder and no path was given."
        systemfile = files[0]
    return System(fromfile=systemfile)


def parseSystemConfigFile(projectConfigFile):
    """
    Auxiliary function to build a System from a System configuration file (SCF)
    """
    from Parsers import SystemConfigFileParser
    if not osp.exists(projectConfigFile): raise BadFile, "File %s not found."%projectConfigFile
    sys = SystemConfigFileParser().parse(projectConfigFile)
    return sys


import Biskit.test as BT

class Test(BT.BiskitTest):
    """Test"""

    def test_SystemOFF(self):
        """Create System from OFF"""
        off = T.testRoot('pep','pep.off')
        self.f_out =  T.tempDir()+os.sep+'SystemOFF'
        if not os.path.exists(self.f_out): os.mkdir(self.f_out)
        T.BROWSER.chdir(self.f_out)
        self.r1 = System(amberOFF=off)
        solvatedsys = self.r1.solvate('WAT', tmp=self.f_out)

    def test_SystemPDB(self):
        pdb = T.testRoot('pep','pep.pdb')
        self.f_out =  T.tempDir()+os.sep+'SystemPDB'
        if not os.path.exists(self.f_out): os.mkdir(self.f_out)
        T.BROWSER.chdir(self.f_out)
        self.r1 = System(amberPDB=pdb)
        solvatedsys = self.r1.solvate('WAT', tmp=self.f_out)
        
    def test_SolvatedSystem(self):
        top = T.testRoot('pep','pep.prmtop')
        crd = T.testRoot('pep','pep.prmcrd')
        self.f_out =  T.tempDir()+os.sep+'SolvatedSystem'
        if not os.path.exists(self.f_out): os.mkdir(self.f_out)
        T.BROWSER.chdir(self.f_out)
        solvatedsys = SolvatedSystem(name='pep', top=top, crd=crd)
        solvatedsys.setPDBfromTOPCRD()
        pdb = solvatedsys.getSolvatedPDB()

    def cleanUp(self):
        if self.f_out: T.tryRemove( self.f_out, tree=1 )


if __name__ == '__main__':
    BT.localTest()
