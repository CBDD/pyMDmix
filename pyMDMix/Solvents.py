##
## pyMDMix --- http://mdmix.sourceforge.net
## Software for preparation, analysis and quality control
## of solvent mixtures molecular dynamics
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 3 of the
## License, or any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
## General Public License for more details.
##
## You find a copy of the GNU General Public License in the file
## license.txt along with this program; if not, write to the Free
## Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
##
## Please cite your use of pyMDMix in published work:
##
##    TOBEPUBLISHED
##

"""

This module provides main :class:`Solvent` object class and
a :class:`SolventManager` class to create/remove solvents
from the solvent database.

The easiest way to create a new solvent is through a configuration file
where all parameters can be assigned.

Second way is to instantiate the solvent directly giving all required
parameters to the constructor method.

Instantiating and working with Solvent objects
----------------------------------------------
Solvent objects might be instantiated giving all required parameters to the constructor.
Here I exemplify it's usage and basic attributes::

    >>> import pyMDMix.tools as T
    >>> import pyMDMix.Solvents as S
    >>>
    >>> off_file = T.testRoot('ETAWAT20.off')   # Amber object file with box unit
    >>> name = 'ETA'
    >>> probesmap = {'OH':'ETA@O1', 'CT':'ETA@C2'}   # Will link probe name OH with residue name ETA atom name O1
    >>> typesmap = {'OH':'Don,Acc', 'CT':'Hyd'}     # Link probe OH with types Donor and Acceptor. Link CT with hydrophobic. Names are arbitrary.
    >>>
    >>> boxunit = 'ETAWAT20' # As named inside the off file
    >>> info = 'Test direct instantiation of a Solvent object'
    >>>
    >>> # Create instance
    >>> solv = S.Solvent(name=name, info=info, off=off_file, probesmap=probesmap, typesmap=typesmap, boxunit=boxunit)
    >>> print solv
    SOLVENT: ETA
    INFO: Test direct instantiation of a Solvent object
    BOXUNIT: ETAWAT20
    >>> print solv.probes # print configured probes
    [CT, OH]
    >>> print solv.types  # set object
    set(['Acc', 'Don', 'Hyd'])
    >>>
    >>> # Calculate the probability of finding atom O1 of residue ETA in a grid voxel of 0.5x0.5x0.5 Angstroms
    >>> print solv.getProbability('ETA','O1',voxel=0.5**3)
    0.0005320194076762995


Adding new solvents to databases
--------------------------------
The Solvent object must be configured using :meth:`SolventManager.createSolvent`
method giving a configuration file as argument.
Here is an exaple of valid configuraiton file with all available options commented:

    .. literalinclude:: solvent_template.cfg
        :language: bash

Configuring and adding a new solvent into the default database::

    >>> import tools as T
    >>> from Solvents import SolventManager
    >>> SM = SolventManager()
    >>>
    >>> # Read configuration and create object.
    >>> # Object file path in the configuration file must be correct or
    >>> # errors will arise
    >>> configfile = T.testRoot('solvent_template.cfg')
    >>> print configfile
    /Users/dalvarez/Dropbox/WORK/pyMDMix/pyMDMix/data/test/solvent_template.cfg
    >>> newsolvent = SM.createSolvent(configfile)
    >>> print newsolvent
    SOLVENT: ETA
    INFO: Ethanol 20% mixture
    BOXUNIT: ETAWAT20
    >>>
    >>> # Add to the database
    >>> SM.saveSolvent(newsolvent)
    ETA saved to database /Users/dalvarez/Dropbox/WORK/pyMDMix/pyMDMix/data/solventlib/SOLVENTS.db
    >>> SM.listSolvents()
    ['PYZ', 'ISX', 'TFE', 'CLE', 'MSU', 'IMZ', 'ANT', 'ION', 'ETA', 'MOH', 'ISO5', 'ETAA', 'ISO', 'WAT', 'COM', 'PYR', 'MAM']

Configure and save object in a specific, empty database::

    >>> # This action will copy all solvents in the default db to the new db and add the new solvent
    >>> customdb = 'mycustomdb.db'
    >>> SM.saveSolvent(newsolvent, db=customdb)
    >>> SM.listSolvents(customdb)
    ['PYZ', 'ISX', 'TFE', 'CLE', 'MSU', 'IMZ', 'ANT', 'ION', 'ETA', 'MOH', 'ISO5', 'ETAA', 'ISO', 'WAT', 'COM', 'PYR', 'MAM']
    >>> # Optionally, an empty database can be created and only add the new solvent
    >>> SM.saveSolvent(newsolvent, db='otherdb.db', createEmpty=True)
    >>> SM.printSolvents('otherdb.db')
    ['ETA']



"""

__author__="dalvarez"
__date__ ="$16-ene-2014 17:09:33$"

import os.path as osp
import logging
import settings as S
import tools as T
import OFFManager as O
from containers import Probe


class SolventParserError(Exception):
    pass

class SolventManagerError(Exception):
    pass

class BadFile(SolventManagerError):
    pass

class MissingSection(SolventParserError):
    pass

class MissingOption(SolventParserError):
    pass

class MappingError(SolventParserError):
    pass

class Solvent(object):
    "Solvent class for storing information on solvent mixtures for simulations"
    def __init__(self, name, info, offpath, boxunit, probesmap, typesmap, frcmodpaths=[], watermodel=S.DEF_AMBER_WATBOX, *args, **kwargs):
        """
        Constructor of a Solvent object.
        It expects some mandatory fields and some optional extra fields that can be assigned through kwargs.

        :arg str name: Name identifying the current solvent mixture.
        :arg str info: Some string describing the solvent mixture. Will be used for printing the solvent.
        :arg str off:  Filename that must exist. Will be used to fetch all information about the mixture and
                        also to later solvate the systems for simulation. Make sure this file is correct and the units
                        correctly working when setting up systems with it.
        :arg str boxunit: String specifying the name of the Leap unit containing the mixture inside the objectfile. Mandatory, specially important when more than 1 units are saved inside same object file.
        :arg probesmap: Map probe names to residues and atom names in the solvent box. Check documentation at the website or the documents for more information.
        :type probesmap: dict
        :arg dict typesmap: Dictionary mapping chemical types to the probes in *probesmap*.

        :arg str watermodel: If the solvent box contains waters, specify the water model used. Example; TIP3P, TIP4P... By default, it will be assigned to :attr:`S.DEF_WATER_BOX`

        :keywords: key=value pairs that will be set as Solvent attributes. Useful for adding extra information in the solvent instance for easy access from custom functions.

        *probemap* and *typesmap* Examples::

            >>> probemap = {'OH':'ETA@O1'} # Identify atom O1 of residue ETA as probe OH.
            >>> typesmap = {'OH':'Don,Acc'} # Assign probe OH the chemical types Don and Acc (for Donor and Acceptor).

        *kwargs* assignment example::

            >>> solvent = Solvent(name='mysolvent',info='custom solvent',off='path/to/objectfile.off',..., myspecialattr = 300)
            >>> print solvent.myspecialattr
            300
                    
        """
        self.name = name
        self.info = info
        self.boxunit = boxunit
        self.watermodel = watermodel
        self.offpath = offpath
        self.frcmodpaths = frcmodpaths
        self.probesmap = probesmap
        self.typesmap = typesmap

        # To be set in self.__setup()
        self.comprobes = []
        self.probes = None
        self.probelist = None
        self.types = None
        self.frcmod = {}
        self.off = None
        self.volume = None
        self.residues = None
        
        # ADD other attributes from options file
        # can be used in specific actions later
        for key, val in kwargs.iteritems():
            if val: print "Adding attribute,vals: %s, %s"%(key,val)
            setattr(self, key, val)

        # Parse probesmap and typesmap
        # and process OFF file to fecth all information from the solventbox
        self.__setup(probesmap, typesmap)

    def __str__(self):
        s = "SOLVENT: {name}\nINFO: {info}\nBOXUNIT: {boxunit}"
        s = s.format(**self.__dict__)
        s += '\nPROBES: %s'%(', '.join(self.probelist))
        return s

    def __repr__(self):
        return self.name+' Solvent Instance'

    def __buildMaps(self, probesmap, typesmap):
        """
        Build a :class:`Probe` list from the :attr:`probesmap` and :attr:`typesmap given in the
        configuration file.
        """
        ###################################################
        #     Set up probes, types and interrelations     #
        ###################################################
        probelist = []
        for probe, mask in probesmap.iteritems():
            d = T.amberMaskToDict(mask)
            
            # Check all residues in the mapping are present in the solvent box
            resnames = set(d.keys())
            if not resnames <= set(map(lambda x: x.name, self.residues)): # all resnames should be present in the off
                raise MappingError, "In PROBES section, mapping to unkown residues: %s"%(resnames - self.residues)

            # Check all atoms are present in each residue
            # TODO keep loop structure for multiple residues although
            # one probe can currently map to only one residue
            for res in resnames:
                resatoms = set([at.name for at in self.getResidue(res).atoms])
                atoms = set(d[res])
                if not atoms <= resatoms:
                    raise MappingError, "In PROBES section, probe %s is mapping \
                                        to unkown atom names for residue %: %s"%(probe, res, atoms - resatoms)

            res = self.getResidue(d.keys()[0])
            atoms = d[res.name]
            # Create Probe instance
            # Add type information and expected probability
            probability = self.getProbability(res.name, atoms)
            probename = '{0}_{1}'.format(self.name, probe)
            P = Probe(name=probename, residue=res, atoms=atoms, type=None, probability=probability)
            P.type = typesmap[probe].split(',')
            probelist.append(P)
        return probelist

    def __setup(self, probesmap, typesmap):
        self.__parseOffandFrcmod()
        self.probes = self.__buildMaps(probesmap, typesmap)
        self.probelist = [p.name for p in self.probes]
        self.types = set(T.simplifyNestedList([p.type for p in self.probes]))
        self.isIonic()
        self.__setCOMs()

    def __setCOMs(self):
        "Set center of mass probe for each residue"
        self.comprobes = dict([('%s_COM'%res.name,res) for res in self.residues])

    def isIonic(self):
        """Check wether the solvent box contains charged residues"""
        # Check if residues are charged
        charged_res = sum([r.charge != 0.0 for r in self.residues])
        self.totalcharge = sum([r.charge for r in self.residues])
        self.isionic = charged_res > 0
        return self.isionic

    def getProbeByName(self, name):
        """
        Returns :class:`~pyMDMix.containers.Probe` instance with name *name*.
        If name is a COM probe, will return the linked residue insetad of a probe instance.
        Return **False** if *name* not found.
        """
        if name in self.probelist:
            return self.probes[self.probelist.index(name)]
        elif name in self.comprobes.keys():
            #It is a COM probe
            return self.comprobes[name]
        return False

    def getProbesByType(self, type):
        """
        Returns a list of :class:`~pyMDMix.containers.Probe` instances that match type *type*.
        """
        return [p for p in self.probes if type in p.type]

    def getResidue(self, resname):
        """
        Return :class:`~pyMDMix.containers.Residue` object with name *resname*
        """
        for r in self.residues:
            if r.name == resname: return r
        return False

    def getProbeProbability(self, probename):
        """
        Return probe probability
        """
        if probename in self.probelist: return self.getProbeByName(probename).p
        elif probename in self.comprobes.keys(): return self.getProbability(self.comprobes[probename].name)
        return False

    def getNumRes(self, res):
        "Count number of residues with name *res* in current box"
        offparser = O.OFFManager(offString=self.off)
        return offparser.getNumRes(self.boxunit, res)

    def getProbability(self, res, atoms=[], voxel=None):
        """
        Obtain expected probability for residue atoms to fall into the voxel volume.
        If atoms list is empty, the probability of the residue to be in a voxel volume will
        be returned.

        :arg str res: Residue name.
        :arg list atoms:    Atom name list. If not given, will consider 1 (= probability of the residue)
        :arg float voxel:   Volume of the voxel.
                            If not given, it is automatically calculated from defaults.

        :returns: Probability of finding any of the atoms in a voxel.
        :type: float
        """
        # Calculate number of residues
        # And number of atoms inside each residue
        if not isinstance(atoms, list): atoms=[atoms]
        nres = self.getNumRes(res)
        numatoms = len(atoms) or 1  # Any residue will not contain duplicated names so just the length will make it
        if not voxel:
            # Calculate voxel volume from default grid spacing
            voxel = S.GRID_SPACING**3
        return (nres*numatoms*voxel)/self.volume

    def writeOff(self, name):
        """Write object file of current solvent to filname *name*"""
        return open(name, 'w').write(self.off)

    def __parseOffandFrcmod(self):
        "Parse off file to fetch info about the solvent box"
        ###################################################
        #     Read frcmod and off files and fecth info    #
        ###################################################
        
        # If frcmod files given, read in content
        # Path should be correct. Will be stored in self.frcmod dictionary with the filename as key.
        if self.frcmodpaths:
            [self.frcmod.update({f:open(f,'r').read()}) for f in self.frcmodpaths]

        # Read the object file and fetch all relevant information
        self.off = open(self.offpath, 'r').read()
        offparser = O.OFFManager(offString=self.off)

        # Check self.boxunit exists
        allunits = offparser.getUnits()
        if not self.boxunit in allunits:
            raise SolventParserError, "Main box unit %s not present in object file %s!"%(self.boxunit, self.offpath)

        self.volume = offparser.getVolume(self.boxunit)
        if not self.volume:
            raise SolventParserError, "Main box unit %s does not have box information! \
                            Are you sure this is the pre-equilibrated solvent box?"%self.boxunit

        # Now check all residues in the boxunit are also inside the object file as
        # separate units. Will save the list as a set for easy comparison.
        unitresidues = set(offparser.getResidueList(self.boxunit, unique=True))
        if not unitresidues <= set(allunits):  # Check residues set is a subset of allunits or equal
            missingres = unitresidues - set(allunits)
            raise SolventParserError, "Objectfile %s does not contain units for the residues %s present in \
                                   main solvent box %s!"%(self.offpath, ','.join(missingres), self.boxunit)
        else:
            # Correct. Fecth residue instances from off
            self.residues = [offparser.getResidue(r) for r in unitresidues]
        
        # Store atomic information for each residue
#        self.atoms = {}
#        [self.atoms.update({res:offparser.getAtoms(res, skipH=False)}) for res in self.residues]


class SolventManager(object):
    "Class to manage solvent creation/removal and database manipulation"
    def __init__(self):
        "Manage solvents and solvent databases"
        self.log = logging.getLogger("SolventManager")

    def __str__(self):
        strout='-'*50
        strout+='\n'
        strout+="SOLVENT DB in USE: %s\nPackage SOLVENT DB: %s\nSOLVENTS Directory: %s\nSolvents:"%(self.__getDatabase(), S.SOLVENTDB, S.SOLVENTS_DIR)
        # Fetch solvents
        lib = self.getDatabase()
        solv=[]
        for s, info in lib.iteritems():
            solv.append(info.__str__().replace('\n','\t'))
        strout+='\n\t'+'\n\t'.join(solv)+'\n'
        strout+='-'*50
        strout+='\n'
        return strout
    
    def __repr__(self):
        return "SolventManager Instance"

    def __getSection(self, dict, name):
        """
        Mandatory get an key named *name* from de dictionary *dict* or raise exception.
        
        :arg dict dict: Dictionary which should contain section named :attr:`name`
        :arg str name: Mandatory section to get
        
        :returns: Dictionary with the section content.
        :rtype: dict
        
        :raise MissingSection: Section is not present in :attr:`dict`
        """
        sect = dict.get(name)
        if not sect: raise MissingSection, "%s section missing in Solvent Config file."%name
        return sect
    
    def __getOption(self, dict, name):
        """
        Mandatory get a key named *name* from de dictionary *dict* or raise exception.
        
        :arg dict dict: Dictionary which should contain` section named :attr:`name`
        :arg str name: Mandatory option to get

        :returns: :class:`SettingsParser.Settings` object with the options content.

        :raises MissingOption: Option is not present in :attr:`dict
        """
        setting = dict.get(name)
        if not setting: raise MissingOption, "%s option missing in Solvent Config file."%name
        return setting

    def __todict(self, settingoption):
        """
        Convert a setting section into dictionary fromat"
        """
        m = dict([[k, v.value] for k,v in settingoption.iteritems()])
        return m

    def __parseConfig(self, configFile):
        """
        Read in a solvent configuration file and process the data to
        create a new solvent in the main or local solvent database.

        Templates for the configuration file can be obtained through
        the user interface.

        :arg str configFile: Solvent configuration filename to read and parse.
                It should contain all mandatory fields.
        """
        import SettingsParser as P
        file = T.absfile(configFile)
        config = P.SettingsParser(file)
        result = config.parse(keepsections=True)

        d = {}

        # GENERAL SECTION
        generalSection = self.__getSection(result,'GENERAL')
        name = self.__getOption(generalSection, 'name').value
        info = self.__getOption(generalSection, 'info').value
        off = self.__getOption(generalSection, 'objectfile').value    
        boxunit =  self.__getOption(generalSection, 'boxunit').value
        if generalSection.get('watermodel'):
            watermodel = generalSection.get('watermodel').value
        else:
            watermodel = S.DEF_AMBER_WATBOX

        d.update({'name':name, 'info':info, 'boxunit':boxunit, 'watermodel':watermodel})

        # Check off file exists and can find it
        if not osp.exists(off):
            # Trye to search it in configfile directory
            off = osp.join(osp.dirname(file), off)
            if not osp.exists(off): raise P.InvalidPath

        # This one is optional. If present, will be loaded along the off file in tleap
        # Usually parameters will be saved inside the off file.
        # If frcmod was given, convert to list and check all files exist
        frcmod = generalSection.get('frcmod')
        frcmodlist = []
        if frcmod:
            frcmod.typeCast(list)
            frcmod = frcmod.value
            frcmodlist = []
            for f in frcmod:
                if not osp.exists(f):
                    # Trye to search it in configfile directory
                    f = osp.join(osp.dirname(file), f)
                    if not osp.exists(f): raise P.InvalidPath, "%f frcmod file not found"%f
                    frcmodlist.append(f)
                else:
                    frcmodlist.append(f)

        # PROBES SECTION
        probesmap = self.__todict(self.__getSection(result,'PROBES'))

        # TYPES section -- Will contain chemical type mapping to probenames
        typesmap = self.__todict(self.__getSection(result,'TYPES'))

        d.update({'offpath':off, 'probesmap':probesmap, 'typesmap':typesmap, 'frcmodpaths':frcmodlist})

        # Check typesMapping and probesMapping have same keys
        if not set(probesmap.keys()) == set(typesmap.keys()):
            typesSet = set(typesmap.keys())
            probesSet = set(probesmap.keys())
            miss = typesSet - probesSet
            self.log.error("Missmatch between probes in [PROBES] section and chemical types in [TYPES] section: %s"%miss)
            raise MappingError

        # OPTIONAL SECTIONS BEGIN HERE
        # CORRECTION SECTION IF MODULES ARE IMPLEMENTED THESE SECTIONS ARE NECESSARY
        # TODO fix this mess of correction parsing
        correctionSect = result.get('CORRECTION')
        corrections = {}
        if correctionSect:
            corr = self.__todict(correctionSect)
            for type, val in corr.iteritems():
                type = type.upper()
                corrections[type] = dict([[el.strip().upper() for el in probe.split(':')] for probe in val.strip().split(',')])
            # transform to float the numbers
            for k,val in corrections.iteritems():
                corrections[k] = dict(zip(val.keys(),map(float, val.values())))
            del correctionSect, type, k, val

        d.update({'corrections':corrections})

        # EXTRA SECTION. OPTIONAL.
        # Specific pairs attribute to add to the solvent instances.
        # Implemented for further development.
        extra = result.get('EXTRA') or {}
        if extra:
            self.log.info("Solvent has EXTRA field. Will add each flag to attributes.")
            extra = self.__todict(extra)
        d.update(extra)
        
        return d

    def createSolvent(self, configfile):
        """
        Create a Solvent isntance from the information in *configfile* configuration file.
        Examples of this file are available at template folder :func:`T.templatesRoot`

        :arg str configfile: Filename with solvent configuration file

        :return solvent: Return a Solvent object configured.
        :rtype: :class:`Solvent`
        """
        if osp.exists(configfile):
            options = self.__parseConfig(configfile)
            return Solvent(**options)
        else:
            raise BadFile, "File %s does not exist"%configfile

    def __getDatabase(self, db=None, createEmpty=False):
        """
        Select and return database (which is actually a dictionary type) for which the user has writing permisions.
        By default, the package database will be used but if the user cannot write there,
        a new database will be created under the users home directory: $HOME/.mdmix/SOLVENTS.db
        which will also contain all solvents from the package db.

        :arg None|str db: If *db* != None:
                                - Check if this DB filename is writtable and return full path. Create it otherwise.
                            If None:
                                - Check if user has writing permisions in package DB and return the path.
                                - If not, create a user-spacific database and return the path.
        :arg bool createEmpty:   If the database has to be created, ignore solvents from the package DB.
                            Will create an empty DB.

        :return: Path to a db for which the user should have writting permisions.
        :type: string
        """
        import shutil
        # Try to use specified db
        if db:
            if osp.exists(db) and T.filePermission(db)['WRITE']: return T.absfile(db)
            else:
                if createEmpty:  T.dump({}, db)
                else:  shutil.copy(S.SOLVENTDB, db)
                return T.absfile(db)
        else:
            # USE USER DATABASE IF IT EXISTS
            if osp.exists(S.USER_SOLVENT_DB): return T.absfile(S.USER_SOLVENT_DB)
            else: db = S.SOLVENTDB

        # Create USER DB from package if user cannot write to package directory
        if not T.filePermission(db)['WRITE']:
            # User cannot write to package DB so create a custom one
            if createEmpty: T.dump({}, S.USER_SOLVENT_DB)
            else:           shutil.copy(db, S.USER_SOLVENT_DB)  # Make a copy to user's home mdmix folder
            db = S.USER_SOLVENT_DB

        return T.absfile(db)

    def getDatabase(self, db=None):
        "Open the database / unpickle it."
        return T.load(self.__getDatabase(db))

    def getDefaultDatabasePath(self):
        return self.__getDatabase()

    def saveSolvent(self, solvent, db=None, createEmpty=False):
        """
        Save a Solvent isntance in the database :attr:`db` or default DB locations.
        Selection of database is done in :meth:`self.__getDatabase`.
        
        :arg solvent: Solvent object to save.
        :type solvent: :class:`Solvent`
        :arg str db:        Database where to save the solvent.
                        If None, default ones will be used (package DB if
                        the user can write there, or user DB otherwise).
        :arg bool createEmpty:   If new database, create empty.
                            If False, copy data from package DB to the new DB.
        """
        db = self.__getDatabase(db, createEmpty)
        solvLib = T.load(db)
        solvLib[solvent.name] = solvent
        T.dump(solvLib, db)
        self.log.info("Solvent %s saved to database %s"%(solvent.name, db))

    def removeSolvent(self, solvName, db=None):
        """
        Remove solvent from database
        Same as saveSolvent, :attr:`db` will be chosen automatically if None.

        :arg str solvName: Solvent name. 
        :arg str db: Database path. If None, automatically detect.

        :raises SolventManagerError: if *db* does not contain *solvName*.
        """
        solvents = self.listSolvents(db)
        if solvName in solvents:
            solvLib.pop(solvName)
            T.dump(solvLib, db)
            self.log.info("Removed solvent %s from database %s"%(solvName, db))
            return True
        else:
            raise SolventManagerError, "DB %s does not contain solvent name %s"%(db, solvName)

    def getSolvent(self, name, db=None):
        """
        Fetch solvent from the database by name.

        :arg str name: Solvent name.
        :arg str db: Database path. If None, automatically detect.

        :Return: Solvent object from the database
        :rtype: :class:`Solvent` instance or False if not found.
        """
        lib = self.getDatabase(db)
        o = lib.get(name)
#        print o.__dict__.keys()
#        S = Solvent(**o.__dict__) # Create a new instance from old one just in case code has changed 
        return o

    def fetchSolventByProbe(self, probename):
        "Giving a probe name, return the corresponding solvent"
        possible = []
        for solv in self.getDatabase().values():
            if probename in solv.probes or probename in solv.comprobes: possible.append(solv)
        
        if not possible:
            self.log.warn("No solvent found for probe name %s"%probename)
            return False
        elif len(possible) > 1:
            self.log.warn("More than one solvent found for probe named %s: %s. Returning %s."%(probename, possible, possible[0].name))
        return possible[0]
    
    def listSolvents(self, db=None):
        """
        Fetch solvent names from the database.

        :arg str db: Database path. If None, automatically detect.

        :Return: Name list
        :rtype: list
        """
        lib = self.getDatabase(db)
        return lib.keys()

    def printSolvents(self, db=None):
        """
        Like list solvents but will print to screen information about the solvents.
        """
        lib = self.getDatabase(db)
        for s, info in lib.iteritems():
            print s,info.info, info.boxunit


# AUXILIARY FUNCTION TO GET SOLVENT INSTANCES FROM DEFAULT DATABASES
# WITHOUT NEED TO INSTANTATE THE MANAGER
def getSolvent(name):
    """
    Axuxiliary module function to get :class:`Solvent` instance *name* from the
    default databases.

    Useful to get solvents without instantiating :class:`SolventManager` from code.
    """
    M = SolventManager()
    return M.getSolvent(name)

###TESTING
import Biskit.test as BT

class Test(BT.BiskitTest):
    """Test"""

    def test_SolventCreation(self):
        """SolventManager Solvent creation test"""

        f_in = osp.join(T.testRoot(), 'solvents','solvent_template.cfg')
        self.f_out =  osp.join(T.tempDir(), 'solvent.db')

        self.m = SolventManager()
        solv = self.m.createSolvent(f_in)
#        print solv.getProbeByName('CT').mask

        self.m.saveSolvent(solv, db=self.f_out)

        self.assertEqual(str(solv) ,str(self.m.getSolvent('ETA', self.f_out)))
#        self.assertEqual(solv, getSolvent('ETA'))
        #        self.assertEqual( r, 42) ## from 'int-testparam = 42' in settings.cfg


    def cleanUp(self):
        T.tryRemove( self.f_out )

if __name__ == '__main__':
    BT.localTest()
