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
__date__ ="$06-feb-2014 15:36:29$"

"""
This module contains several classes for information storage and cross linking information.
"""

class Atom(object):
    "Simple container for atomic information gathered in the OFF file: name, type, element, charge"
    def __init__(self, id, name, type, element, charge, *args, **kwargs):
        self.id = id    #: Integer. Atom ID.
        self.name = name    #: Atom name
        self.type = type    #: Atom AMBER TYPE
        self.element = element  #: Element. Integer.
        self.charge = charge    #: Partial charge

    def __repr__(self):
        return self.name

    def __str__(self):
        s = "ID: {id}\nATOM NAME: {name}\nATOM TYPE: {type}\nELEMENT: {element}\NCHARGE: {chargestr}"
        d = self.__dict__
        d.update({'chargestr':'%.4f'%self.charge})
        return s.format(d)

    def __eq__(self, other):
        "Compare this to other by name"
        if isinstance(other, Atom): other=other.name
        return self.name == other

class Residue(object):
    """
    Simple container for whole residue unit information gathered in the OFF file.
    Basically to sotre atomic information and possibly masks for later quick identify atom
    positions.
    """
    def __init__(self, name, atoms, connectivity, xyz, *args, **kwargs):
        self.name = name    #: Name of the residue
        self.atoms = atoms  #: List of :class:`Atom` instances that belong to residue
        self.connectivity = connectivity    #: Tuple with connectivity information
        self.xyz = xyz                      #: XYZ coordinates in a numpy array Nx3

        # Calculate total residue charge
        #: Total charge of the residue
        self.charge = round(reduce(lambda x,y: x+y, [a.charge for a in self.atoms]), 4) # round to 4 decimals

        # Build some maps
        #: Dictionary mapping atom ids to :class:`Atom` instances
        self.atids = dict([(a.id,a) for a in self.atoms])

        #: Dictionary mapping atom names to :class:`Atom` instances
        self.atnames = dict([(a.name,a) for a in self.atoms])

    def __repr__(self):
        return self.name

    def __str__(self):
        s = "RESIDUE NAME: {name}\nATOMS: {atoms}"
        return s.format(**self.__dict__)
    
    def __eq__(self, other):
        "Compare this to other by name"
        if isinstance(other, Residue): other=other.name
        return self.name == other
    
    def center(self):
        return self.xyz.mean(axis=0)
    
    def setxyz(self, xyz):
        "Set new coordinates"
        if xyz.shape == self.xyz.shape:
            self.xyz = xyz
        else:
            self.xyz = False

class Probe(object):
    """
    Container for probe information. This object will store information about a particular probe linking
    atom names, residues and chemical types with probabilities. Will also contain a mask for the residue.
    """
    def __init__(self, name, residue, atoms, type, probability):
        self.name = name        #: Name of the probe as given in :attr:`Solvent.probesmap`
        self.residue = residue  #: :class:`Residue` instance with corresponding residue information
        self.atoms = atoms      #: Atom name list
        self.type = type        #: Chemical type
        self.p = probability    #: Probability of probe to be found in grid voxel volume
        self.mask = self.__createMask()

    def __eq__(self, other):
        "Compare this to other by name"
        if isinstance(other, Probe): other=other.name
        return self.name == other
    
    def __repr__(self):
        return self.name

    def __createMask(self):
        "Create a mask for selected atoms in the residue. This will match the atom order."
        # Set a mask for the atoms in residue that we are interested in
        import numpy as npy
        mask = npy.zeros(len(self.residue.atoms))
        for i, at in enumerate(self.residue.atoms):
            if at.name in self.atoms: mask[i] = 1
        return mask.astype(bool)
    
    def istype(self, type):
        "Check if current replica has assigned type *type*"
        return type in self.type
    


if __name__ == "__main__":
    print "Hello World"
