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
This module provides one important class: :class:`MDSettings`. When instantiated, this object
will contain all parameters needed to create input files for running the Molecular Dynamics simualtions. These parameters include
the length (in nanoseconds), restraining schemes (either flexible, heavy atoms restrained or backbone restrained), trajectory wirting frequency,
simulation steps, etc...

A default parameter configuration has been proposed in the MD Settings File (MSF) :ref:`mdsettings`.
All parameters there found are automatically assigned to each MDSettings object attribute. If MDSettings is instantiated with
any of the options given as parameter, the attibutes from the settings file are overriden:

    >>> import pyMDMix
    >>>
    >>> # Use all default values given in the md-settings file
    >>> defaultsettings = pyMDMix.MDSettings('WAT')
    >>> defaultreplica.nanos
    20
    >>> defaultreplica.md_timestep
    2
    >>>
    >>> customreplica = R.Replica(nanos=40, md_timestep=1) # 40ns at 1fs timestep
    >>> customreplica.md_timestep
    1
    >>>
    >>> # Applying restraints over heavy atoms (HA) over all protein residues (automatic mask) with a force of 0.1 kcal/mol.A^2
    >>> restrainedReplica = R.Replica('MAM', name='restrainTest', restrMode='HA', restrForce=0.1)
    >>> print restrainedReplica
    REPLICA INFORMATION:
    -------------
    Replica name: restrainTest
    Solvent: MAM
    ------- Amber files:
            -----------
            Amber PRMTOP: None
            Amber PRMCRD: None
            PDB file: None
    (...)
    ------- SIMULATION:
            ------------
            MD program: AMBER
            Temperature: 300.0
            Nanoseconds: 20
            Use restraints?: True
    ------- RESTRAINTS:
            -----------
            Schema: HA
            Force: 0.1
    >>>

To see all settings that a replica has finally adopted, use :meth:`~Replica.printSettings` method::

    >>> restrainedReplica.printSettings()
    Replica restrainTest settings:
    ------------------------
    md_timestep: 2.0
    parm_heating_tempi: 100.0
    restrMode: HA
    minsteps: 5000
    nanos: 20
    amber_solvate_buffer: 14.0
    (...)
    >>>

As mentioned here :ref:`settingsmodule`, the user should find an editable copy of the ``md-settings.cfg`` file in his home directory (``~/.mdmix/md-settings.cfg``).
Edit this file to include changes that will affect all future program runs. Otherwise, if you are interested in modifying the parameters for concrete replicas,
do it in replica creation time.
"""

__author__="dalvarez"
__date__ ="$16-mar-2014 21:03:48$"

import sys
import difflib
import SettingsParser as P
import settings as S

class MDSettingsError(Exception):
    pass

class MDSettings(object):
    def __init__(self, solvent=None, name=None, restrMask='', alignMask='', **kwargs):
        """
        Constructor method for MDSettings objects.

        :arg str name:   Settings name.
        :arg str solvent:   Solvent name. It must exist in Solvent database.
        :arg str restrMask: Amber format mask to select residues and atoms to be restrained (if needed). If emtpy and restrains are requested, an automatic mask will be calculated from the pdb
        :arg str alignMask: Amber format mask to select atoms and residues over which trajectory should be aligned. If empty, an automatic mask will be calcualted.

        :keywords: This provides a very flexible attribute assignment system.
            - Every pair key=value will be assigned as an attribute to current MDSettings.
            - Pair values not given will take default values from Global Settings or User Settings specifications.
        """
        # Set filenames/formats/replica info
        self.solvent = solvent              #: Solvent mixture name (e.g. ETA)
        self.analyzenetcdf = S.NETCDF_AVAIL #: Is NetCDF installed? If so, use it for trajectory storage
        self.restrMask = restrMask          #: Mask of residues ids to apply restrains to (it doesn't need atomic information, just residue numbers)
        self.alignMask = alignMask          #: Mask of residue ids to align to when imaging, centering and aligning trajectory in analysis process

        # ADOPT ATTRIBUTES FROM KWARGS
        # ONLY IF THEY ARE PRESENT IN SETTINGS
        # FINALY GET DEFAULT ATTRIBUTES FROM SETTINGS
        # ONLY ADOPT NOT DEFINED IN KWARGS
        m = self.__getSettings()
        settingKeys = m.settings2dict().keys()
        for k,v in kwargs.iteritems():
            if v is None: continue
            if k in settingKeys: setattr(self, k, v)
            else:
                # Will try to do fuzzy comparison
                # to give some flexibility to the user (specially for case matching
                bestmatch = difflib.get_close_matches(k, settingKeys, 1, 0.8)
                if bestmatch: setattr(self, bestmatch[0], v)
                else: print >> sys.stderr, "Attribute %s not in md-settings. skipping..."%k
        m.updateNamespace(self.__dict__, keepdefined=True)

        # Set name (automatic name if not given)
        if name: self.name = name
        else: self.name = 'Settings_{solvent}_{restrMode}'.format(**self.__dict__)

        # Add some operation attrs
        self.hasRestraints = self.restrMode != 'FREE'
        # Calculate expected number of trajectory files
        self.__calcExpectedNTrajFiles()
        self.__calcNumberSnapshots()

    def __str__(self):
        """Print settings"""
        return self.getSettingsStr()

    def __repr__(self):
        return "%s MD Settings"%self.name

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __getSettings(self):
        "Get default settings"
        m = P.SettingsManager(S.CFG_MD_DEFAULT, S.CFG_MD_USER, createmissing=True)
        m.collectSettings()
        return m

    def getSettingsStr(self):
        "Print all settings"
        s = "MD settings {0}:\n-------------------------------\n".format(self.name)
        m = self.__getSettings()
        sdict = m.settings2dict().keys()
        for k in sdict:
            s+="{0}: {1}\n".format(k, getattr(self, k))
        return s

    def updateFromSettings(self):
        "Auxiliary function to update object with attributes from configuration files."
        m = P.SettingsManager(S.CFG_MD_DEFAULT, S.CFG_MD_USER, createmissing=True)
        m.updateNamespace(self.__dict__, keepdefined=True)
        if self.__folderscreated: self.write()

    def __calcNumberSnapshots(self):
        self.nsnaps = self.nanos*(self.prod_steps/float(self.trajfrequency))

    def __calcExpectedNTrajFiles(self):
        "Calculate from self.nanos and depending on steps per file and timestep, the total number of trajetory files expected"
        # convert nanosecond to femtosecond (1e6), divide by num of fs per step to obtain number of steps
        # Finally divide num of steps needed by num of steps per file
        self.ntrajfiles=int((self.nanos*1e6/self.md_timestep)/self.prod_steps)


### Functions to create project from configfile
### and load existing one
def parseSettingsConfigFile(settingsConfigFile, noSolvent=False):
    """
    Auxiliary function to build replicas from a replica configuration file (RCF)

    :arg str replicaConfigFile: Path to config file containing REPLICA section.
    :arg bool makefolders: Create folders for each replica in current path? Automatic names will be given to each replica.

    :return: List of :class:`~Replicas.Replica` objects constructed with all parameters in the RCF
    """
    import os.path as osp
    from Parsers import MDSettingsConfigFileParser
    
    if not osp.exists(settingsConfigFile): raise BadFile, "File %s not found."%settingsConfigFile
    if noSolvent: sets = MDSettingsConfigFileParser().parseNoSolvent(settingsConfigFile)
    else: sets = MDSettingsConfigFileParser().parse(settingsConfigFile)
    return sets


import pyMDMix.test as BT
import pyMDMix.tools as T

class Test(BT.BiskitTest):
    """Test"""
    def test_MDSettings(self):
        """Test MDSettings object"""
        defaults = MDSettings()
        watdefaults = MDSettings(solvent='WAT')
        self.assertEqual(defaults.trajfrequency, watdefaults.trajfrequency)
        self.assertEqual(defaults.restrMode, watdefaults.restrMode)
        
    def test_MDSettings_CFG(self):
        etaha = MDSettings(solvent='ETA', restrMode='HA', restrForce=0.1)
        etacfg = parseSettingsConfigFile(T.testRoot('mdsettings','sets1.cfg'))
        self.assertEqual(etaha, etacfg[0])
    
    def test_MDSettings_snapshots(self):
        sets = MDSettings(prod_steps=1000000, trajfrequency=1000) #1Milion = 2ns at 1000 steps = 1000 snaps per nano x 20ns default 
        self.assertEqual(sets.nsnaps, 20000)
        
    def test_MDSettings_ntrajfiles(self):
        sets = MDSettings(prod_steps=1000000, nanos=50) # 50ns and 2ns per file: 25 files?
        self.assertEqual(sets.ntrajfiles, 25)
#    def cleanUp(self):
#        T.tryRemove( self.f_out )

if __name__ == '__main__':
    BT.localTest()
