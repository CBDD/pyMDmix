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
This module provides a reader for Amber OFF file format.

Example::
    >>> import os.path as osp
    >>> import pyMDMix.tools as T
    >>> import pyMDMix.OFFManager as O
    >>>
    >>> f_in = osp.join(T.testRoot(), 'ETAWAT20.off')
    >>> m = O.OFFManager(offFile=f_in)
    >>>
    >>> print m.getUnits() # Get unit names present in the OFF file
    ['ETA','ETAWAT20','WAT']
    >>> print m.getResidueList('ETAWAT20', unique=True) # Get residues inside unit 'ETAWAT20'
    ['WAT','ETA']
    >>> print m.getVolume('ETAWAT20') # Volume of the box
    7988.43038
    >>> print m.getNumRes('ETAWAT20', 'ETA') # Number of 'ETA' residues inside 'ETAWAT20' unit.
    17


"""

__author__="dalvarez"
__date__ ="$17-ene-2014 1:07:04$"

import os.path as osp
from containers import Residue, Atom

class OFFManagerError(Exception):
    pass

class OFFSectionError(OFFManagerError):
    pass

class OFFManager(object):
    "Manage Amber OFF file types. Only for reading."
    SECTIONS = ['atoms','atomspertinfo','boundbox',
                'childsequence','connect','connectivity',
                'hierarchy','name','positions','residueconnect',
                'residues','residuesPdbSequenceNumber', 'solventcap',
                'velocities']
                
    def __init__(self, offFile=None, offString=None, *args, **kwargs):
        """
        Read Amber OFF files. Two options are accepted for initialization:
            - Give a path to a valid object file in :attr:`offFile`
            - Give a string with an object file content in :attr:`offString`

        :arg str offFile: Path to Amber OFF file to read.
        :arg str offString: String with the content of an already read offFile.
                
        """
        self.tmpfile = None
        
        if offFile:
            if osp.exists(offFile): self.off = open(offFile, 'r').read()
            else: raise OFFManagerError, "Object File %s not found"%offFile
        elif offString:
            self.off = offString
        else:
            raise OFFManagerError, "offFile or offString are needed for initializing instance."

    def __iterOff(self):
        "Returns the file as an iterable list"
        return iter(self.off.split('\n'))

    def getResidue(self, res, skipH=False):
        """
        Fetch residue in off and return a :py:class:`~pyMDMix.containers.Residue` instance
        aontaining also atomic information.

        :arg str unit: Residue name which should correspond to a valid unit in off file.
        :arg bool skipH: Skip hydrogen atom information. Default:False.

        :Returns: Residue instance info.
        :rtype: :py:class:`~pyMDMix.containers.Residue` or **False** if unit not found.
        """
        return Residue(name=res, atoms=self.getAtoms(res, skipH=skipH),
                        connectivity=self.getConnectivity(res), xyz=self.getCoords(res))

    def getCoords(self, unit):
        """
        Fetch positions information for unit selected.
        This is section !entry.UNIT.unit.positions in off file.

        :arg str unit:  Unit name

        :Returns: Coordinates of unit atoms.
        :rtype: :class:`numpy.ndarray` of floats with size *Nx3*
        """
        import numpy as npy
        positions = self.readOffSection(unit, 'positions')
        xyz = npy.array([l.split() for l in positions], dtype=npy.float32)
        return xyz
    
    def getConnectivity(self, unit):
        """
        Fetch connectivity table for unit selected.
        This is section !entry.UNIT.unit.connectivity in off file.

        :arg str unit: Unit name

        :Returns: List with bonded pair indexes verbosely. Example: ((1,2),(1,3),(3,1),(2,1)...)
        :rtype: list
        """
        try: connectinfo = self.readOffSection(unit, 'connectivity')
        except: return []
        connectivity = []
        for line in connectinfo:
            pair = map(int, line.split())[:2]
            connectivity.append(tuple(pair))
            pair.reverse()
            connectivity.append(tuple(pair))
        return connectivity

    def getAtoms(self, unit, skipH=False):
        """
        Fetch atomic information for the unit selected.
        Will return a dictionary with atom names and types.

        :arg str unit: Unit to search.
        :arg bool skipH: If atom is Hydrogen, skip it.

        :Returns: List of :py:class:`~pyMDMix.containers.Atom` instances containing id, name, atomtype, element and charge information.
        :rtype: list
        """
        atomlist = []
        atominfo = self.readOffSection(unit,'atoms')
        for line in atominfo:
                line = line.split()
                element = int(line[6])
                if skipH and element == 1: continue

                id = int(line[5])
                name = line[0].split('"')[1]
                type = line[1].split('"')[1]
                charge = float(line[7])

                d = {'name':name, 'id':id, 'type':type, 'charge':charge, 'element':element}
                atomlist.append(Atom(**d))
        return atomlist

    def getUnits(self):
        "Return the list of units int he object file."
        off = self.__iterOff()
        # Skip first line and store names until '!' is found again
        off.next()
        units = []
        while 1:
          line = off.next().strip()
          if "!" in line: break
          units.append(line.split("\"")[1])
        del off
        return units

    def hasUnit(self, unitname):
        "Return True if the OFF file has the unit with name 'unitname'."
        return unitname in self.getUnits()

    def isParameter(self, unit):
        """
        Check if unit *unit* is a parameter unit inside OFF.

        :arg str unit: Name of the unit

        :Returns: True if its a parameter unit. False otherwise.
        """
        line = self.__find('!entry.%s.'%unit)
        if line:
            if 'parm' in line[0]: return True
        return False

    def __find(self, expr, returnlines=True):
        """
        Find expresion in off file. Will match ANY occurrence.
        Wildcards are accepted. Will use fnmatch module from python.

        :Arguments:
            *expr*
                Any string that should be present in the off file.
                Example::
                    '!entry.*.parm' will match all sections for parameter units.

            *returnlines*
                Bool. If True, return the matching lines.
        :Returns:
            *match*
                If there was a match, will return::
                    - Matching lines if *returnline* is True.
                    - Number of matches (integer) if *returnlist* is False.
                If there was no match, will return False
        """
        import fnmatch
        off = list(self.__iterOff())
        match = fnmatch.filter(off, '*%s*'%expr)

        if not match: return False
        if returnlines: return match
        else: return len(match)

    def getResidueList(self, unit, unique=True):
        """
        Get a list of residue names for the *unit* chosen.

        :arg str unit: Unitname to search.
        :arg bool unique: If True, return a list with unique names. If False, the complete list of names will
                be returned.
        :Returns: List with residuenames inside unit *unit*.
        """
        section = self.readOffSection(unit, 'residues')
        resnamelist = []
        for line in section:
            resnamelist.append(line.split()[0].replace('"',''))
        if unique:
            return list(set(resnamelist))
        else:
            return resnamelist

    def readOffSection(self, unit, section, with_header=False):
        """
        Parse the object file and read a whole section for the unit selected.

        :arg str unit: Unit name to search.
        :arg str section: OFF file section name. Example: *residues* section will correspond to "!entry.UNITNAME.unit.residues .."  part of the file.
        :arg bool with_header: If True, return output with heading line of the section.

        :Returns: Content of the section unitl next '!entry' is found. Returned with our without the heading line depending
                on the value of *with_header*
        :rtype: list of strings
        """
        off = self.__iterOff()
        search = '!entry.'+unit+'.unit.'+section
        line = off.next()
        while line and not search in line:
            line = off.next()
        if not line: raise OFFSectionError, "Section %s for unit %s not in file."%(section, unit)

        out = []
        if with_header: out.append(line)
        line = off.next()
        while not '!entry' in line:
            out.append(line.strip())
            line = off.next()
        del off
        return out

    def getNumRes(self, unit, residue):
        "Count the number of residues with name *residue* inside unit *unit*"
        res = self.getResidueList(unit, unique=False)
        n = res.count(residue)
        return n

    def getNumAtoms(self, unit, residue, atomname):
        # First the number of residues
        nRes = self.getNumRes(unit, residue)
        section = self.readOffSection(residue, 'atoms')
        n = 0
        for line in section:
            if atomname in line: n+= 1
        return n*nRes

    def getBoxDimensions(self, unit):
        "Get box dimension information from the object file for :attr:`self.boxunit`"
        boxsection = self.readOffSection(unit, 'boundbox')
        return map(float, boxsection[2:])

    def getVolume(self, unit):
        "Get volume information from the object file for :attr:`self.boxunit`"
        boxdims = self.getBoxDimensions(unit)
        return reduce(lambda x,y: x*y, boxdims)

    def write(self, outname):
        """
        Write OFF content to disk file.

        :arg str outname: File name to store the Object File.

        :return: Absolute path to the file.
        """
        open(outname, 'w').write(self.off)
        return osp.abspath(outname)

    def writeTmp(self):
        """
        Write OFF content to temprorary disk file.
        Filepath is also stored in :attr:`tempfile`
        
        :return: Absolute path to the file.
        """
        import tempfile
        outname = tempfile.mktemp()
        self.tmpfile = outname
        open(outname, 'w').write(self.off)
        return outname

    def cleanTmp(self):
        """
        Remove temporary file if created with :meth:`writeTmp`
        """
        import tools as T
        if self.tmpfile:
            return T.tryRemove(self.tmpfile)

import Biskit.test as BT
import tools as T

class Test(BT.BiskitTest):
    """Test"""

    def test_OFFManager(self):
        """OFFManager test"""

        f_in = T.testRoot('solvents','ETAWAT20.off')
        
        self.m = OFFManager(offFile=f_in)
        self.assertEqual(self.m.getUnits(), ['ETA','ETAWAT20','WAT'])
        self.assertEqual(self.m.getResidueList('ETAWAT20',unique=True), ['WAT','ETA'])
        self.assertAlmostEqual(self.m.getVolume('ETAWAT20'), 7988.43038, 4)
        self.assertFalse(self.m.getVolume('ETA'))
        self.assertFalse(self.m.isParameter('ETA'))

        res = self.m.getResidue('ETA')
        self.assertEqual(res.charge, 0)
        self.assertEqual(res.connectivity[0], (1,2))
#        self.assertEqual( r, 42) ## from 'int-testparam = 42' in settings.cfg

#
#    def cleanUp(self):
#        T.tryRemove( self.f_out )

if __name__ == '__main__':
    BT.localTest()
