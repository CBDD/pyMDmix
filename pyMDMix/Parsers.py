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

import logging
import sys
import os.path as osp
import ConfigParser

import settings as S
from Solvents import SolventManager
import Systems
import MDSettings

class ParserError(Exception):
    pass
class SystemParserError(ParserError):
    pass
class BadFile(ParserError):
    pass
class MDSettingsParserError(ParserError):
    pass
class BadSolvent(MDSettingsParserError):
    pass

class SystemConfigFileParser(object):
    "Read a config file and generate a System"
    def __init__(self):
        self.log = logging.getLogger('SystemConfigFileParser')

    def parse(self, configFile):
        "Parse all options in configFile"
        if not osp.exists(configFile): raise BadFile, "File %s not found."%configFile
        self.__configHandle = ConfigParser.ConfigParser()
        self.__configHandle.read(configFile)

        ########################################################################
        ######## READ PROJECT SECTION AND CREATE PROJECT INSTANCE #########
        ########################################################################
        
        try:
            fileSection = dict(self.__configHandle.items('SYSTEM'))
        except:
            raise SystemParserError, "Section SYSTEM missing."

        parms = {} # Store parameters to be passed to System

#        projName = fileSection.get('projname', None)
#        if not projName:
#            projName = "mdmixProjName%0.4i"%(random.randint(1,9999))
#            self.log.warn("No SystemName given! Will generate a random name: %s"%projName)
#        parms.update({'name':projName})

        # Build can be from an object file or a PDB
        name = fileSection.get('name') or 'system'
        off = fileSection.get('off')
        pdb = fileSection.get('pdb')
        top = fileSection.get('top')
        crd = fileSection.get('crd')
        if off:
            # If OFF given, this has priority
            if pdb: self.log.warn("Ignoring PDB entry in config file. OFF entry has priority and will be used.")
            
            # CONSIDER SPECIAL SYSTEM CASE FOR TESTING
            
            if off == 'test': parms.update({'amberOFF':off})
            else:
                if not osp.exists(off): 
                    # Try to find off in folder where config file is found
                    offpath = osp.join(osp.split(configFile)[0], off)
                    if not osp.exists(offpath): raise BadFile, "Entry OFF pointing to file %s not found."%off
                    else: off = offpath
                parms.update({'amberOFF':osp.abspath(off)})

            parms.update({'unitName':fileSection.get('uname')})
        elif pdb:
            if not osp.exists(pdb): raise BadFile, "Entry PDB pointing to file %s not found."%pdb
            parms.update({'amberPDB':osp.abspath(pdb)})
        elif top and crd:
            # Directly create Solvated system
            if not osp.exists(top) or not osp.exists(crd): raise BadFile, "PRMTOP file %s or PRMCRD file %s not found."%(top,crd)
            parms.update({'top':top, 'crd':crd})

        # Optional: specify extra residue names to keep in reference structure (useful for modified residues)
        extrares = fileSection.get('extrares', [])
        if extrares:
            extrares = [l.strip() for l in extrares.split(',')]
            self.log.info("Given extra residue list: %s"%extrares)
        parms.update({'extraResList':extrares})

        # Reatrin mask is common
        restrmask = fileSection.get('restrmask') or 'auto'
        alignmask = fileSection.get('alignmask') or 'auto'
        # If extra forcefield parameters are to be loaded
        extraFF = fileSection.get('extraff', [])
        if extraFF: extraFF = [f.strip() for f in extraFF.split(',')]
        parms.update({'restrainMask':restrmask, 'alignMask':alignmask, 'FF':extraFF})

        # Return System
        if not off and (pdb and crd):
            return Systems.SolvatedSystem(name=name, **parms)
        else:
            return Systems.System(name=name, **parms)

class MDSettingsConfigFileParser(object):
    "Read a replica configuration file which should contain a [REPLICA] section and create and add ReplicaInfo to the project."
    def __init__(self):
        self.log = logging.getLogger('MDSettingsConfigFileParser')

    def __checkSolventList(self, l):
        "Check any of the solvent numbers given is diferent from de available ones"""
        available = SolventManager().listSolvents()
        missing = set(l) - set(available)
        if missing: raise BadSolvent, "Solvent list invalid: %s"%missing
        return True

    def __splitPerReplicaSlash(self, string, valControl=None, valFormat=str):
        "ValControl, accepted values"
        if not string:
            # Empty field or not assigned Return None for each replica
            replInfo = {}
            for solv, nrepl in self.solv_nrepl.iteritems():
                replInfo[solv] = {}
                for i in range(1, nrepl+1): replInfo[solv][i] = None
            return replInfo

        # Else treat the string
        string = string.strip()
        tempD = {}
        if ',' in string:
            split = string.split(',')
            for piece in split:
                if '/' in piece:
                    p = piece.split('/') # Here we have 3 parts: solvent, replica/s, value
                    solv = p[0].strip()
                    repl = p[1].strip()
                    value = valFormat(p[2].strip())
                    if valControl and value not in valControl:
                        self.log.error("(__splitPerReplicaSlash) value '%s' not accepted. Possible values are: %s"%(value, valControl))
                        sys.exit(1)

                    expectedNrepl = self.solv_nrepl[solv]
                    if not tempD.has_key(solv): tempD[solv] = {}

                    # repl can be empty if all replicas should be modified
                    # it can contain a single integer or
                    # it can contain a range (1-5)
                    if not repl:
                        r = range(1, expectedNrepl+1)
                    elif '-' in repl:
                        s = repl.split('-')
                        r = range(int(s[0]), int(s[1])+1)
                    else:
                        r = [int(repl)]

                    # Obtain combinations
                    for i in range(1, expectedNrepl+1):
                        if i in r:
                            tempD[solv][i] = value

                else:
                    value = valFormat(piece)
                    if valControl and value not in valControl:
                        self.log.error("(__splitPerReplicaSlash) value '%s' not accepted. Possible values are: %s"%(value, valControl))
                        sys.exit(1)
                    tempD['COMMON'] = value
        else:
            value = valFormat(string)
            if valControl and value not in valControl:
                self.log.error("(__splitPerReplicaSlash) value '%s' not accepted. Possible values are: %s"%(value, valControl))
                sys.exit(1)
            tempD['COMMON'] = value
#            print 'Common',value

        # Assign COMMON value to replicas not specified
        replInfo = {}
        for solv, nrepl in self.solv_nrepl.iteritems():
            replInfo[solv] = {}
            for i in range(1, nrepl+1):
                if tempD.has_key(solv) and tempD[solv].has_key(i):
                    replInfo[solv][i] = tempD[solv][i]
                else:
                    cval = tempD.get('COMMON')
#                    if cval == None:
#                        self.log.error("(__splitPerReplicaSlash) No common value given. Could not assign %s_%i replica."%(solv,i))
#                        sys.exit(1)
#                        continue
                    replInfo[solv][i] = cval

        return replInfo

    def __splitPerReplicaColon(self, string, valFormat=int):
        # Split in parts and divide name:value
        # format: NREPL = 3, WAT:1, ETA:6
        temp_replInfo = {}
        string = string.strip()
        if ',' in string:
            split = string.split(',')
            for piece in split:
                if ':' in piece:
                    p = piece.split(':')
                    temp_replInfo[p[0].strip()] = valFormat(p[1])
                else:
                    temp_replInfo['COMMON'] = valFormat(piece)
        else:
            temp_replInfo['COMMON'] = valFormat(string)
        # Assign a value to each solvent
        replInfo = {}
        common = temp_replInfo.get('COMMON')
        for solv in self.solvents:
            if solv in temp_replInfo.keys():
                replInfo[solv] = temp_replInfo[solv]
            else:
                if not common: raise MDSettingsParserError, "Common value not given in differential assignment and solvent %s cannot be assigned"%solv
                replInfo[solv] = common
        return replInfo

    def parse(self, configfile):
        """
        Parse replica configuration file and return a list of Replica objects.

        :arg str configfile: Path to replica configuration file.

        :return: Replica objects as a list
        """
        import SettingsParser as P
        import difflib

        settingsInstances = []

        if not osp.exists(configfile): raise BadFile, "Config file does not exist: %s"%configfile
        self.__configHandle = ConfigParser.ConfigParser()
        self.__configHandle.read(configfile)

        #################################################################
        ########                READ MDSETTINGS INFO                #####
        #################################################################

        # Make sure there are MDSETTINGS sections in file
        sections = self.__configHandle.sections()
        mdsections = [s for s in sections if s.startswith('MDSETTINGS')]
        nummdsections = len(mdsections)
        if not nummdsections: raise MDSettingsParserError, "No sections found starting with MDSETTINGS name."

        # Visit all mdsettings sections
        for section in mdsections:
            fileSection = dict(self.__configHandle.items(section))

            solvents = fileSection.get('solvents') or fileSection.get('solvent')
            if not solvents: raise MDSettingsParserError, "SOLVENT(S) option missing in replica config file."
            solvents = [el.strip() for el in solvents.split(',')]
            self.__checkSolventList(solvents)
            self.solvents = solvents

            # Number of replicas for the solvents
            # In this version it is possible to specify different number of replicas for different solvents
            # If 1 int is given, all replicas share the same number
            # Else, with a comma, separate SOLVENT:NUMREPL, if only an int is given, replicas not specified
            # are considered as common.
            #
            # Examples: NREPL = 3                 --  All solvents will have 3 replicas
            #           NREPL = 3, WAT:1          --  All solvents 3 replicas except water, with 1 only
            #           NREPL = 3, ETA:6, WAT:2   --  Ethanol 6 replicas, Water 2 and the rest 3 replicas
            #

            nrepl = fileSection.get('nrepl','%i'%S.DEF_NREPLICAS)
            self.solv_nrepl = self.__splitPerReplicaColon(nrepl)

            # Restraining mode // Force // Nanos // Temperature parsing per replica
            #
            # Here per-replica specifications can be given.
            # Format:
            #   RESTR = BB              -- All replicas restrain in the backbone
            #   RESTR = HA, WAT//FREE   -- All WAT replicas without restraints, the rest restrains in the heavy atoms
            #   RESTR = HA, ETA/1/BB, ETA/1/HA, ETA/2/FREE  -- Ethanol replica 1 with BB restraint, 1 with HA and 2 FREE, rest is HA
            #   RESTR = HA, MAM/1-3/HA, MAM/4-5/FREE        -- MAM replicas 1 to 3 HA restrained. 4 to 5 FREE. Rest HA.
            #
            # Numbering of the replicas starts at 0. The numbers must coincide within the total NREPL number given.
            # SAME FORMAT FOR NANOS, TEMPERATURE AND FORCE
            restr = fileSection.get('restr') or None
            if restr: restr = restr.upper().strip()
            splitedRestr = self.__splitPerReplicaSlash(restr, valControl=['BB','HA','FREE'])
            force = fileSection.get('force') or None
            splitedForce = self.__splitPerReplicaSlash(force, valFormat=float)
            nanos = fileSection.get('nanos') or None
            splitednanos = self.__splitPerReplicaSlash(nanos, valFormat=int)
            temp = fileSection.get('temp') or None
            splitedTemp = self.__splitPerReplicaSlash(temp, valFormat=float)

            # Optional TOP and CRD files to create the replica directly without project
            #        top = fileSection.get('top') or None
            #        crd = fileSection.get('crd') or None
            #        off = fileSection.get('off') or None
            #        extraResidues= fileSection.get('extrares')
            #        if extraResidues: [el.strip() for el in extraResidues.split(',')]
            #        else: extraRediues = []
            #        refPdb= fileSection.get('ref') or None
            restrMask= fileSection.get('restrmask') or ''
            alignMask=fileSection.get('alignmask') or ''

            # Options in fileSection to be ignored as they where already treaten
            mainopts = ['restr','solvents','solvent','nanos','temp', 'force', 'nrepl', 'restrmask','alignmask']

            ### CHECK FOR CONFIGURATION PARAMETERS PRESENT IN CONFIG FILES
            ### THAT CAN BE MODIFIED
            m = P.SettingsManager(S.CFG_MD_DEFAULT, S.CFG_MD_USER, createmissing=True  )
            m.collectSettings()
            settingKeys = m.settings2dict().keys()
            extracfg = {}
            for k,v in fileSection.iteritems():
                if v is None: continue
                if k in mainopts: continue
                # Will try to do fuzzy comparison to identify what config parameter should be modified
                # to give some flexibility to the user (specially for case matching
                bestmatch = difflib.get_close_matches(k, settingKeys, 1, 0.8)
                if bestmatch:
                    setting = m.settings[bestmatch[0]]
                    extracfg.update({setting.name:setting.vtype(v)})
                else: raise MDSettingsParserError, "Attribute %s not present in md-settings. Make sure the spelling is correct"%k

            for solv, nrepl in self.solv_nrepl.iteritems():
                for i in range(1, nrepl+1):
                    replicaRestrMode = splitedRestr[solv][i]
                    replicaRestrForce = splitedForce[solv][i]
                    replicaNanos = splitednanos[solv][i]
                    replicaTemp = splitedTemp[solv][i]
                    settingsInstances.append(MDSettings.MDSettings(solvent=solv, nanos=replicaNanos, restrMode=replicaRestrMode,
                                    restrForce=replicaRestrForce, temp=replicaTemp, restrMask=restrMask, alignMask=alignMask, **extracfg))

        return settingsInstances


import Biskit.test as BT
import tools as T
class Test(BT.BiskitTest):
    """Test"""

    def test_PepSystemParser(self):
        """Parse system config file"""
        cfg = T.testRoot('pep','pep_amber_mdmix.cfg')
        system = SystemConfigFileParser().parse(cfg)

    def test_PepSettingsParser(self):
        """Parse settings config file"""
        cfg = T.testRoot('pep','pep_amber_mdmix.cfg')
        settings = MDSettingsConfigFileParser().parse(cfg)
        assert len(settings) == 2
        assert settings[0].trajfrequency == 2000

    def test_PepSettingsParserMulti(self):
        """Parse settings config file with multiple MDSETTINGS entries"""
        cfg = T.testRoot('pep','multisettings.cfg')
        settings = MDSettingsConfigFileParser().parse(cfg)
        assert len(settings) == 5
        assert settings[-1].trajfrequency == 500

if __name__ == "__main__":
    BT.localTest()

#
#import Biskit.test as BT
#import tools as T
#class Test1(BT.BiskitTest):
#    """Test"""
#    def test_PepSystemParser(self):
#        """Parse system config file"""
#        cfg = T.testRoot('pep','pep_amber_mdmix.cfg')
#        system = SystemConfigFileParser().parse(cfg)
#
#    def test_PepSettingsParser(self):
#        """Parse settings config file"""
#        cfg = T.testRoot('pep','pep_amber_mdmix.cfg')
#        settings = MDSettingsConfigFileParser().parse(cfg)
#        assert len(settings) == 2
#        assert settings[0].trajfrequency == 2000
#
#    def test_PepSettingsParserMulti(self):
#        """Parse settings config file with multiple MDSETTINGS entries"""
#        cfg = T.testRoot('pep','multisettings.cfg')
#        settings = MDSettingsConfigFileParser().parse(cfg)
#        assert len(settings) == 5
#        assert settings[-1].trajfrequency == 500
##
##if __name__ == '__main__':
##    BT.localTest()
#
