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

import logging
from operator import itemgetter
from itertools import groupby
import numpy as npy

import Biskit as bi


# Amber residue names sometimes not well identified by PDBModel
AMBER_RESNAMES = ['CYM','CYX','ASH','LYN','GLH']
ION_RESNAMES = ['CIO','Cl-','Cs+','IB','K+','Li+','MG2','Na+','Rb+'] # from ions94.lib

class SolvatedPDBError(Exception):
    pass

class SolvatedPDB(bi.PDBModel):
    """
    Subclass of Biskit.PDBModel to specifically deal with proteins or nucleic acids in solvent mixtures.
    """
    def __init__(self, pdb=None, solvent=None, extraResidues=[], *args, **kwargs):
        bi.PDBModel.__init__(self, source=pdb)
        self.log = logging.getLogger("SolvatedPDB")
        self.extraResidues = extraResidues
        self.soluteMask = None
        self.solventMask= None
        self.probeMasks = {}
        self.resMasks = {}
        if pdb: self.setSolvent(solvent)
        
#        self.fixNumbering()
#        self.maskSoluteSolvent()
        
    def __getstate__(self):
        d = self.__dict__.copy()
        del d['log']
        return d

    def __setstate__(self, d):
        d['log'] = logging.getLogger("SolvatedPDB")
        self.__dict__.update(d)

    def getNumResidues(self):
        "Return a counter dictionary with the number of residues byname"
        from collections import Counter
        c = Counter()
        res = self.atom2resProfile("residue_name")
        for r in res: c[r] += 1
        return c

    def setExtraResidues(self, extrares):
        """
        Add residue names which should be considered as part of the solute in the system.
        :arg list extrares: Name of residues to consider as solute
        """
        self.extraResidues += extrares
        self.setSoluteSolventMask()

    def setSolvent(self, solvname=None):
        """
        Set solvent for current system. Used for identification of solvent residues.
        If *solvname* is **None**, will try to identify the corresponding solvent mixture from the database
        by comparing residue names in the solvent mask.

        :arg str solvname: Solvent name to add. Will fetch all information from the solvent database.
        """
        from Solvents import SolventManager
        man = SolventManager()
        if not solvname:
            # Try to identify solvent from pdb composition
            self.setSoluteSolventMask()
            resnames = npy.unique(npy.array(self['residue_name'])[self.solventMask])
            resnames = set(resnames) - set(ION_RESNAMES)  # Remove known ions
            dbsolvents = [man.getSolvent(s) for s in man.listSolvents()]
            solvresnames = dict([(s.name, map(repr,s.residues)) for s in dbsolvents])
            possible = []
            for sol, res in solvresnames.iteritems():
                if set(res) == resnames: possible.append(sol)

            if not possible:
                self.log.debug("Solvent not identified: %s"%self.source)
                self.solvent = None
            else:
                if len(possible) > 1:
                    self.log.debug("More than one possible solvent identified for SolvatedPDB: %s. Using the first one: %s"%possible, possible[0])
                else:
                    self.log.debug("Identified pdb solvation as %s solvent box."%possible[0])
                self.solvent = man.getSolvent(possible[0])

        else:
            self.solvent = man.getSolvent(solvname)

    def setSoluteSolventMask(self):
        "When missing residues, only accept as part of the system those present in acceptList or in EXTRARESIDUES file.\
         Waters are always removed as this is intended for generating a refernce pdb and for masking issues."
        self.soluteMask = self.maskProtein() + self.maskNA() + \
                        self.maskFrom('residue_name', self.extraResidues) + self.maskFrom('residue_name', AMBER_RESNAMES)
        self.solventMask = ~self.soluteMask

    def fixNumbering(self, force=False):
        """
        This method fixes numbering of residues and atoms if they are greater
        than the maximum allowed by the column width (9990 for residues and 99999 for atoms).

        :arg bool force: If *True* force renumbering even maximum number is not reached.
        """
        if self.lenResidues() > 9999 or force:
            self['resnum'] = npy.arange(self.lenResidues()) + 1
            self['residue_number'] = self.res2atomProfile('resnum')
        if len(self) > 99999 or force:
            self['serial_number'] = npy.arange(len(self)) + 1
    
    def __prepareSolventResMasks(self):
        "Prepare mask for each solvent residue name"
        for res in self.solvent.residues:
            self.resMasks[res.name] = []
            resids = npy.unique(self['residue_number'][self.maskFrom('residue_name', res.name)])
            [self.resMasks[res.name].append(self.maskFrom('residue_number', rn)) for rn in resids]

    def __prepareSolventProbeMasks(self):
        "Get mask of atoms of interest per probe"
        for probe in self.solvent.probes:
            self.probeMasks[probe.name] = self.maskFrom('residue_name', probe.residue.name) * self.maskFrom('name', probe.atoms)
        return True
    
    def iterResidues(self, residuename):
        "Iterate over requested residues coordinates"
        if not self.resMasks: self.__prepareSolventResMasks()
        masks = self.resMasks.get(residuename)
        if not npy.any(masks): raise SolvatedPDBError, "Invalid residue name %s. Residue not in solvent %s."%(residuename,self.solvent.name)
        for m in masks:
                yield self.xyz[m]

    def getProbeCoords(self, probename):
        "Obtain all coordinates for probe *probename*"
        if not self.probeMasks: self.__prepareSolventProbeMasks()
        masks = self.probeMasks.get(probename)
        if npy.any(masks): 
            return self.xyz[masks]
        else:
            # Not normal probe, may be a COM probe?
            if not probename in self.solvent.comprobes.keys(): raise SolvatedPDBError, "Invalid probe name %s. Probe not in solvent %s."%(probename, self.solvent.name)
            else:
                # Its a com probe, determine residue and 
                # Fetch COM coordinates
                res = self.solvent.comprobes[probename].name
                return npy.array([coord.mean(axis=0) for coord in self.iterResidues(res)])
    
#    def checkAll(self):
##        if self.checkMissingResidues():
##            self.log.warn("Identified some residues that will be missing during reference pdb generation or automask identification: %s. Ignore this message if the system is solvated with organic molecules."%self.missingResidueList)
#        if self.checkMissingResidues():
#            self.log.warn("Identified some residues that will be missing during reference pdb generation or automask identification: %s. \
#                        If you want to keep them, add 'extrares' field under config file 'PROJECT' section. If you want to let them out ignore this message."%self.missingResidueList)
#        if self.checkWaters():
#            self.log.info("PDB has structural waters. They will be removed for generating the reference PDB")
#        if self.checkChains():
#            self.log.info("Protein has multiple chains. Beaware to check the automatic trajectory alignment is correct when runing the analysis.")
#        #if self.checkCloseAtomContacts():
#        #    self.log.warn("Close contacts (<1.5A) found in the input protein. Consider doing a small minimization before proceding.")
#
##    def extraResiduesGiven(self):
##        "If 'EXTRARESIDUES' text file exists, grab the list of residue names that should be considered as part of the system and saved in the reference pdb."
##        if exists('EXTRARESIDUES'):
##            keepResidueList = open('EXTRARESIDUES').readlines()[0].split(',')
##            self.extraResGiven = True
##            self.extraResList = keepResidueList
##            return True
##        return False
#
#    def checkChains(self):
#        "Take only the protein and check wether more than 1 chain is present."
#        pdb = self.getSolute()
#        nChains = pdb.lenChains()
#        if nChains > 1:
#            self.hasMultipleChains = True
#            self.nChains = nChains
#            return nChains
#        else:
#            return False

    def chainResIDs(self):
        "Return the residue numbers for each chain"
        protein = self.getSolute()
        startIDs = npy.take(protein['residue_number'],protein.chainIndex())
        lastID = protein['residue_number'][-1]
        ids = map(lambda x,y: range(x,y), startIDs.tolist(), startIDs.tolist()[1:]+[lastID])
        # Due to range function, last residue number is missing
        ids[-1].append(lastID)
        return ids

    def __maskToResIds(self, residueMask):
        "Given a residue number mask in form: '1-150,160-200' obtain all the residue ids"
        # Split first the commas
        resIds = []
        groups = residueMask.split(',')
        if len(groups) > 1:
            for group in groups:
                if '-' in group:    # range in the group
                    first, last = map(int, group.split('-'))
                    resIds += range(first, last+1)
                else:               # only one residue in the group
                    resIds += [int(group)]
        else:
            if '-' in groups:
                first, last = map(int, groups[0].split('-'))
                resIds += range(first, last+1)
            else:               # only one residue in the group
                    resIds += [int(groups)]
        return resIds

    def getBBMaskSelectedRes(self, residueMask):
        selectedResIds = self.__maskToResIds(residueMask)
        mask = self.maskFrom('residue_number', selectedResIds)
        mask *= self.maskBB()
        return npy.array(self['serial_number'])[mask]

    def getHAMaskSelectedRes(self, residueMask):
        selectedResIds = self.__maskToResIds(residueMask)
        mask = self.maskFrom('residue_number', selectedResIds)
        mask *= self.maskHeavy()
        return npy.array(self['serial_number'])[mask]

    def getBBAutoMaskAtomIds(self, extraResidueList = None):
        if extraResidueList:
            self.addAcceptedResidues(extraResidueList)
        if not npy.any(self.soluteMask): self.setSoluteSolventMask()
        mask = self.soluteMask*self.maskBB()
        return npy.array(self['serial_number'])[mask]

    def getHAAutoMaskAtomIds(self, extraResidueList = None):
        if extraResidueList:
            self.addAcceptedResidues(extraResidueList)
        if not npy.any(self.soluteMask): self.setSoluteSolventMask()
        mask = self.soluteMask*self.maskHeavy()
        return npy.array(self['serial_number'])[mask]

    def getAutoMaskResIds(self, extraResidueList = None):
        "Obtain AMBER mask format, e.g. '1-125,130-140', from automask"
        if extraResidueList:
            self.addAcceptedResidues(extraResidueList)
        if not npy.any(self.soluteMask): self.setSoluteSolventMask()
        resIds = npy.unique(npy.array(self['residue_number'])[self.soluteMask])
        resGroups = [map(itemgetter(1),g) for k,g in groupby(enumerate(resIds), lambda (i,x):i-x)]
        outMask = ''
        nGroups = len(resGroups)
        ng = 1
        for group in resGroups:
            if len(group) > 1:
                outMask += '%i-%i'%(group[0], group[-1])
            else:   # only one number
                outMask += '%i'%group[0]
            if ng < nGroups: outMask +=','
            ng += 1

        return outMask

    def checkMissingResidues(self):
        """ Given a pdb file, check if besides normal protein residues there are
        some extra ones that will be missing if applying normal maskProtein.
        Will compare a normal Protein mask with the total number
        of residues in the pdbfile that should be prepared for runing in the void."""
            
        pdb = self.removeWaters()
        total = set(pdb['residue_name'])
        pmask = pdb.maskProtein() + pdb.maskNA() + pdb.maskFrom('residue_name', self.extraResidues) + \
                                    pdb.maskFrom('residue_name', AMBER_RESNAMES)
        prot = pdb.compress(pmask)
        pres = set(prot['residue_name'])

        extra = total - pres - set(ION_RESNAMES) # remove ion names
        if extra:
            self.hasMissingResidues = True
            self.missingResidueList = list(extra)
            return True
        else:
            self.hasMissingResidues = False
            return False

    def getSolute(self):
        "Returns Prot or RNA or DNA (if present) in the system"
        if not npy.any(self.soluteMask): self.setSoluteSolventMask()
        return self.compress(self.soluteMask)

    def removeWaters(self):
        "Return a PDB instance without the waters"
        return self.compress(~self.maskH2O())

    def clone( self ):
        """
        Clone SolvatedPDB into new independent object
        """
        s = self.take( self.atomRange() )
        s.__dict__.update(self.__dict__)
        return s
    
#    def checkWaters(self):
#        "Detect if waters are present"
#        maskWats = self.pdb.maskH2O()
#        if npy.any(maskWats):
#            self.hasWats = True
#            return True
#        else:
#            self.hasWats = False
#            return False
#
#    def checkCloseAtomContacts(self, d=1.5):
#        "Check for close contacts (<d A)"
#        dm = npy.zeros((len(self.pdb), len(self.pdb)))
#        for i in range(len(self.pdb)):
#            dm[i,:] = npy.sqrt(npy.sum((self.pdb.xyz - self.pdb.xyz[i,:])**2,axis=1))
#        dm[dm == 0] = 999
#        closeContacts = npy.where(dm < d)
#        if closeContacts:
#            self.hasCloseContacts = True
#            self.closeContacts = closeContacts
#            self.dm = dm
#            return True
#        else: 
#            self.hasCloseContacts = False
#            return False
#
#    def writeReferencePDB(self, outname):
#        "Translate only protein to origin stripping waters/ions and residues not specified previously"
#        pdb = self.getSolute()
#        pdb.xyz -= pdb.center()
#        pdb.writePdb(outname, amber=1)
#        self.refPdb = pdb

if __name__ == '__main__':
    from Config import TEST_DIR
    import os


    