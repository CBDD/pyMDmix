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

__author__="daniel"
__date__ ="$Nov 15, 2012 6:20:24 PM$"

import os
import time
import logging
import numpy as npy
from scipy.cluster import hierarchy
from scipy.spatial import distance

from GridsManager import Grid
import tools as T

class HotSpot(object):
    "Contain a set of points defining a hotspot and its properties: volume, chemical type, etc..."
    def __init__(self, coords, energies, energymethod='volume', centroid='min',
                probe=False,index=0, spacing=[0.5,0.5,0.5], info='',**kwargs):
        """
        HotSpot (HS) object. It may contain multiple coordinates and energies for the same hotspot.From *coords*, the minimum energy coordinate will be selected as HS centroid.
        Moreover, the extension in x,y,z coordinates and the volume will be calcualted. From *energies*, 
        a mean or averaged energy will be considered the total energy of the hotspot (depending on *energymethod*).
        
        :arg coords: 3D coordinates of the hotspot (may be 1 or more).
        :type coords: :class:`npy.array`
        :arg energies: 1D array of energies for each coordinate.
        :type coords: :class:`npy.array`
        :arg str energymethod: Algorithm to compute the hotspot energy from the energy list::
                                            - min:   minimmum energy as HS total energy
                                            - avg:   boltzmann average as the HS energy
                                            - volume:volume average as the HS energy
        :arg str centroid: Method for centroid calculation. Options: *min* for centroid in minimum energy coordinate,
                            *mean* for centroid in average coordinate (weighted by energy).
        :arg str probe: Probe name identifying hotspot nature
        :arg int index: identification number if needed (e.g. cluster id)
        :arg spacing: Spacing of the grid. Used to calculate the volume. 1x3 array.
        :type spacing: :class:`npy.array` or list
        :arg str info: String with some description or information about the HP if desired.
        """
        self.coord = None
        self.energy = None
        self.extension = None
        self.volume = None
        self.npoints = 0
        self.sphereindex = 0
        self.info = info
        self.centroid = centroid

        if not isinstance(coords, npy.ndarray): coords = npy.array(coords)
        if not isinstance(energies, npy.ndarray): energies = npy.array(energies)

        self.coordList = coords
        self.energyList = energies
        self.index = index
        self.probe = probe

        # Check energymethod is valid
        valmethods = ['min','avg','volume']
        if energymethod in valmethods:
            self.energymethod = energymethod
        else:
            raise AttributeError, "Wrong energy averaging method. Valid methods are: %s"%valmethods
        
        # Calculate volume element from spacing
        self.spacing = spacing
        if len(spacing) == 3:
            self.volelement = npy.prod(spacing)
        else:
            self.log.warn("spacing should be a 3 number array (one distance for each dimension)")

        self.setup()
        self.calcSphereSimilarity()

    def setIndex(self, i):
        self.index = i
    
    def setProbe(self, probe):
        self.probe = probe
        
    def calcSphereSimilarity(self):
        """
        Considering the volume and the extension, return an value of sphere ressemblance
        sphereindex value ranges from 1 (complete sphere) to any value > 1 giving idea of sphere deformation.
        e.g. 1.3 = 30% sphere deformation

        A minimum volume is required for this index to take some meaning as very small hotspots
        are prone to yield higher deformations. Take in consideration the volume before discarding hotspots!
        """
        # Get the theoretical radius a sphere should have for the hotspot volume
        # Remember the sphere volume formula:    V = (4/3)*pi*r^3
        expectedRadius = ((3*self.volume)/(4*npy.pi))**(1/3.)
        expectedDiameter = expectedRadius*2 + self.spacing

        # Calculate relation extension/expectedRadius
        self.sphereRatios = self.extension/expectedDiameter # times 2 because extension would be a diameter
        self.sphereindex = npy.linalg.norm(self.sphereRatios-1) + 1

    def writeHotSpotPDB(self, outpdbfile, onlycenter=False, resname='HOT', atom='C'):
        """
        Write a PDB file containing current hotspot with info in a header and
        atoms in the body.

        :arg str outpdbfile:       Filename to write pdb
        :arg bool onlycenter:          If True, write only centroid and total energy.all coordinates and energies in hotspot.
                                    If False, write all coordinates and energies in hotspot.
        :arg str atom:       Atom for defining hotspot coordinates.
            
        :Returns: **True** if PDB file was correctly saved. **False** otherwise.
        """
        head, body = self.getHotSpotPDBstr(onlycenter=onlycenter, resname=resname, atom=atom)
        with open(outpdbfile,'w') as k:
            k.write(head+'\n')
            k.write(body)
        return True

    def getHotSpotPDBstr(self, onlycenter=False, resname=False, atom=False, occupancy=False, bfactor=False):
        """
        Prepare a PDB string containing current hotspot with info in a header string and
        atoms in a body string.

        :Args bool onlycenter: If True, write only center coordinate in hotspot and total energy.
                                    If False, write all hotspot points and partial energies.
        :args str resname: Residue name for defining hotspot coordinates.
        :args str atom: Atom name for defining hotspot coordinates.
        :args float occupancy: Value to assign in occupancy column. By default, sphereindex will be written.
        :args float bfactor: Value to assign in bfactor column. By default, energy will be written.
        
        :Returns: Tuple (header, body) with hotspot info in header and atoms in body.
        """
        # Write a PDB for each point found
        head=[]
        pstr='ATOM  %5d %4s %3s %s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f'
        out = []
        
        # define resname and atom
        if not resname: resname = self.probe.split('_')[0][:3]
        if not atom: 
            if '_' in self.probe: atom = self.probe.split('_')[1][:4]
            else: atom = '{:>3s}'.format(self.probe.strip())
        
        if self.index: id = self.index
        else: id = 1

        if not onlycenter and self.coordList.ndim >1:
            for i,coord in enumerate(self.coordList):
                out.append(pstr%(i+1,atom,resname,'A',id, coord[0], coord[1], coord[2], self.meanradius, self.energyList[i]))
        else:
            # Only one point or centroid to be use
            coord = self.coord
            out.append(pstr%(1,atom,resname,'A',id, coord[0], coord[1], coord[2], self.sphereindex, self.energy))

        head.append("CLUSTER ID:%i RADIUS:%.2f VOLUME:%.2f SPHEREINDEX:%.2f ENERGY:%.2f"%(id,self.meanradius,self.volume,self.sphereindex,self.energy))
        return ('\n'.join(head), '\n'.join(out))

    def setup(self):
        "Compute centroid and total energy"
        if self.coordList.ndim == 1:
            # Only one coord and energy
            self.coord = self.coordList
            self.energy = self.energyList
            self.volume = self.volelement
            self.extension = self.spacing
            self.meanradius = npy.mean(self.extension)/2.
            self.npoints = 1
        else:
            RT = 0.001986*300 # kcal/mol
            # More dimensions, meaning multiple coords and energies
            if self.centroid.lower() == 'min':
                # The centroid coord will be the coordinate with minimum energy
#                minenergyIndex = npy.where(self.energyList == self.energyList.min())[0]
#                self.coord = self.coordList[minenergyIndex,:][0]
                self.coord = self.coordList[npy.argmin(self.energyList),:]
            else:
                # Centroid in average coordinates
                probabilities = npy.exp(self.energyList/-RT)
                probabilities = probabilities/probabilities.sum()
                self.coord = npy.average(self.coordList, axis=0, weights=probabilities)
                #print probabilities
            
            # Other properties
            self.npoints = len(self.coordList)
            self.volume = self.npoints*self.volelement
            # introduce spacing in extension calculation to include uncertainty and avoid bad behaviour with small number of points
            self.extension = self.coordList.max(axis=0) - self.coordList.min(axis=0) + self.spacing 
            self.meanradius = npy.mean(self.extension)/2.

            # Calculate energy according to energymethod
            if self.energymethod == 'min':
                self.energy = self.energyList.min()
            elif self.energymethod == 'avg':
                # boltzmann average
                exps = npy.exp(-self.energyList/RT)
                avg = (self.energyList*exps)/npy.sum(exps)
                self.energy = avg
            elif self.energymethod == 'volume':
                # consider volume elements
                # same as averaging probabilities then converting back to energy
                self.energy = -RT*npy.log(npy.exp(self.energyList/-RT).mean())
            else:
                self.energy = None
                
    def __cmp__(self, other):
        if self.energy < other.energy:  # compare name value (should be unique)
            return -1
        elif self.energy > other.energy:
            return 1
        else: return 0              # should mean it's the same instance

    def __sub__(self, other):
        "Overloads substraction method. Will return distance between the two hotspots compared OR to coordinates in tuple, list or npy.array format"
        if isinstance(other, npy.ndarray) or isinstance(other, list) or isinstance(other, tuple):
            if len(other) == 3: coord = npy.array(other)
            else: return False
        else:
            coord = other.coord
        return npy.sqrt(((self.coord-coord)**2).sum())

    def __repr__(self):
        return "HotSpot ID:%d, E:%.2f C:%s Vol:%.2f Extension:%s SphereIndex: %.2f"%(self.index, self.energy,self.coord,self.volume, self.extension, self.sphereindex)
    

class HotSpotSet(object):
    "Base for a Set of HotSpots of the same chemical type"
    def __init__(self, probe='', name='', info='', *args, **kwargs):
        """
        Set of hotspots and operations to manipulate them.
        It is a good idea to identify the set by means of a 'name' and give some description in 'info'.
        'probe' references the probe grid over which the hotspotset was created.
        """
        self.log = logging.getLogger("ProjectManager.HotSpotSet (%s)"%name)
        self.probe = probe
        self.name = name
        self.hotspots = []
        self.nhotspots = 0
        self.dm = None
        self.linkage = False
        self.clusterIndexes = False
        self.nclusters = 0
        self.info = info
        self.kwargs=kwargs
        

    def __repr__(self):
        return "HotSpotSet Name: %s, Probe: %s. %i Hotspots."%(self.name, self.probe,len(self.hotspots))

    def __iter__(self):
        for h in self.hotspots:
            yield h

    def __getitem__(self, i):
        return self.hotspots[i]

    def __getstate__(self):
        d = self.__dict__.copy()
        del d['log']
        return d

    def __setstate__(self, d):
        d['log'] = logging.getLogger("ProjectManager.HotSpotSet (%s)"%d['name'])
        self.__dict__.update(d)


    def setInfo(self, info):
        self.info = info

    def updateDistanceMatrix(self):
        "Updates condensed distance matrix between hotspots coords"
        self.dm = distance.pdist([h.coord for h in self.hotspots])
#        self.clusterHotSpots()

    def getHSbyID(self, id):
        """
        Return hotspots with corresponding ids.
        
        :arg list id:   Id(s) of hotspots to retrieve from the set.
        :returns: List of hotspots
        :rtype: list of :class:`HotSpot` instances
        """
        return [h for h in self.hotspots if h.index in id]
    
    def getCondensedDMatrix(self):
        return self.dm
    
    def getSquareDMatrix(self):
        return distance.squareform(self.dm)

    def addHotSpots(self, HotSpotList):
        "Add a hotspot or a list of hotspots and update distance matrix between them"
        if isinstance(HotSpotList, list):  self.hotspots += HotSpotList
        else: self.hotspots.append(HotSpotList)
        self.nhotspots = len(self.hotspots)
#        if self.nhotspots > 1: self.updateDistanceMatrix()
        self.hotspots.sort()    # Keep ascending energy order with 1st hotspots the lower energy one

    def writeHotSpotsPDB(self, outpdbfile, onlycenter=False, atom=False):
        # Write a PDB for each point found
        head=[]
        out = []
        for i, hp in enumerate(self.hotspots):
            if not hp.index: hp.index = i+1
            h, b = hp.getHotSpotPDBstr(onlycenter=onlycenter, atom=atom)
            head += [h]
            out += [b]
        with open(outpdbfile,'w') as k:
            k.write('\n'.join(head))
            k.write('\n')
            k.write('\n'.join(out))
        self.log.info("DONE writting PDB %s with hotspots"%outpdbfile)

    def writeOneHotSpotPDB(self, outpdbfile, hotspot, onlycenter=False, atom=False):
        """
        Write a PDB file containing only the chosen hotspot
        """
        if hotspot.writeHotSpotPDB(outpdbfile, onlycenter=onlycenter, atom=atom):
            self.log.info("DONE writting PDB %s with hotspot information"%outpdbfile)
        else:
            return False

    def clusterHotSpots(self, cutDistance=2.5):
        # obtain linkage matrix
        self.updateDistanceMatrix()
        self.linkage = hierarchy.linkage(self.dm, method='average')
        self.clusterIndexes = hierarchy.fcluster(self.linkage, t=cutDistance, criterion='distance')
        self.nclusters = npy.unique(self.clusterIndexes).size
        return self.clusterIndexes

    def plotClusters(self, colorThreshold=2.5):
        "Plot clusters obtained by clusterHotSpots. Clusters"
        import pylab
        if not npy.any(self.linkage): self.clusterHotSpots()
        hierarchy.dendrogram(self.linkage, color_threshold=colorThreshold)
        pylab.show()

    def pruneByVolume(self, minvolume):
        """Prunning by minimmum volume. Remove hotspots which do not have a minimum volume.
        
        Args:
            minvolume   (float)     Minimum volume a hotspot should have to be considered.
        
        Returns:
            HotSpotSet containing the new hotspots after prunning
        """
        selectedHS = []
        for hp in self.hotspots:
            if hp.volume >= minvolume: selectedHS.append(hp)
        newSet = HotSpotSet(self.probe, self.name)
        newSet.addHotSpots(selectedHS)
        return newSet
    
    def pruneByEnergy(self, maxenergy):
        """Prunning by maximum mean energy. Remove hotspots which exceed this energy value.
        
        :arg float maxenergy: Maximum energy a hotspot should have to be kept
        :returns: HotSpotSet containing the new hotspots after prunning
        """
        selectedHS = []
        for hp in self.hotspots:
            if hp.energy <= maxenergy: selectedHS.append(hp)
        newSet = HotSpotSet(self.probe, self.name)
        newSet.addHotSpots(selectedHS)
        return newSet

    def pruneByShape(self, maximumSphereIndex):
        """Prunning by shape. All hotspots should have a SphereIndex <= maximumSphereIndex.
        
        Args:
            maximumSphereIndex  (float)     SphereIndex value all hotspot should not trapass to be considered.
                                            Values range from 1 to any. Being 1 a perfect sphere.
        
        Returns:
            HotSpotSet containing the new hotspots after prunning
        """
        selectedHS = []
        for hp in self.hotspots:
            if hp.sphereindex <= maximumSphereIndex: selectedHS.append(hp)
        newSet = HotSpotSet(self.probe, self.name)
        newSet.addHotSpots(selectedHS)
        return newSet

    def pruneByNpoints(self, minimumNumberPoints):
        """Prune by number of points. Remove hotspots that contain a number of poitns < minimumNumberPoints.

        Args:
            minimumNumberPoints     (int)   Number of points a hotspot should contain to be considered.

        Returns:
            HotSpotSet containing the new hotspots after prunning
        """
        selectedHS = []
        for hp in self.hotspots:
            if hp.npoints >= minimumNumberPoints: selectedHS.append(hp)
        newSet = HotSpotSet(self.probe, self.name)
        newSet.addHotSpots(selectedHS)
        return newSet

    def getHotSpotNearCoord(self, coordinate, maxdist=2):
        """
        Return the hotspot closest to the coordinate 'coordinate'.

        Args:
            coordinate      (npy.array)     3D coordinate [x,y,z]
            maxdist         (float)         If no HS within this distance, False is returned.

        Returns:
            HotSpots closest to this coordinate.
            FALSE if distance exceeds maxdist
        """
        distances = npy.array([hspot - coordinate for hspot in self.hotspots])
        mindist = distances.min()
        minindex = npy.where(distances == mindist)[0]
        if len(minindex) > 1:
            e = [self.hotspots[i].energy for i in minindex]
            best = npy.where(e == e.min())[0]
            minindex=minindex[best]
            self.log.warn("Two hotspots found at same minimum distance to coordinate %s. Returning the lower energy one."%(coordinate))
        if mindist <= maxdist:
            return self.hotspots[minindex]
        else:
            return False

    def getHotSpotNearProteinAtom(self, pdb, atomnumber, maxdist=2):
        """
        Return the closest hotspot to pdb atom id.

        Args:
            pdb     (PDBModel)      PDB with the system
            atomnumber   (int)      Serialnumber of atom
            maxdist      (float)    Maximum allowed distance to atom.
                                    If no HS within this distance, False is returned.

        Returns:
            HotSpot which is closer to the atomnumber
            FALSE if distance exceeds maxdist
        """
        import Biskit as bi
        if not isinstance(pdb, bi.PDBModel):
            raise AttributeError, "pdb argument should be a Biskit.PDBModel instance"

        mask = pdb['serial_number'] == atomnumber
        atomcoord = pdb.xyz[mask][0]
        return self.getHotSpotNearCoord(atomcoord, maxdist=maxdist)

    def getSphericalFraction(self, maxsphereindex=1.6):
        """
        Calculate the fraction of hotspots in the current hotspotset that have a sphereindex <= maxsphereindex.
        This function is helpful for deciding if the clustering and hotspot creation was correctly done. If this value is too low,
        it means many hotspots do not have a sphere-like shape (not even close) meaning that probably they extend too much over the
        protein surface. It is possible to have disseminated hotspots but not common. 
        
        So a big fraction means we have well localized hotposts with sphere-like shapes. A low fraction the contrary
        and the user should re-consider rebuilding the hotspotset with a lower cutvalue (for instance) or decrease 
        the cutdistance in the clustering algorithm to separate better.
        
        Args:
            maxsphereindex      (float)     sphereindex to use for fraction calculation. 
        
        Returns:
            fraction    (float)     Fraction of hotspots in the hotspotset which sphereindex is <= maxsphereindex.
        """
        si = npy.array([h.sphereindex for h in self.hotspots])
        return npy.sum(si <= maxsphereindex)/float(si.size)

    def reduceClusters(self, update=True):
        "Calculate a centroid for the clusters identified and average the energy\
        Substitute clusters for centroids if update=True. Return the list of new hotspots otherwise."
        if not npy.any(self.clusterIndexes):
#            self.log.warn("Clustering not calculated.")
            self.clusterHotSpots()
#            return self.hotspots
        if self.nhotspots != self.nclusters:
            newhotspotslist = []
            for i in range(self.nclusters):
                i += 1
                cluster = self.clusterIndexes == i
                if npy.sum(cluster == i) > 1:
                    clustHP = npy.array(self.hotspots)[cluster]
                    newhotspotslist.append(HPCluster(clustHP.tolist(), clusterID=i))
                else:
                    newhotspotslist.append(npy.array(self.hotspots)[cluster][0])

            # Return only or update values?
            if update:
                self.hotspots = newhotspotslist
                self.nhotspots = len(self.hotspots)
#                self.updateDistanceMatrix()
                self.clusterHotSpots()
            else:
                return newhotspotslist

    def dump(self, file):
        "Write pickle to file *file*"
        T.dump(self, file)

class HPCluster(HotSpot):
    def __init__(self, hotspotlist, clusterID=None):
        self.hotspots = hotspotlist
        if not isinstance(hotspotlist, list): self.hotspots = [self.hotspots]
        self.index = clusterID or 1
        self.setuphp()
    
    def setuphp(self):
        "Calculate centroid and average energy and STD"
        self.coord = npy.array([h.coord for h in self.hotspots]).mean(axis=0)
        elist = [h.energy for h in self.hotspots]
        self.energy = npy.mean(elist)
        self.var = npy.var(elist)
        self.std = npy.std(elist)
        self.probes = [h.probe for h in self.hotspots]
        self.probe = self.probes[npy.argsort(elist)[0]]
        self.coordList = npy.array([h.coord for h in self.hotspots])
        self.energyList = npy.array([h.energy for h in self.hotspots])
        self.sphereindex = 1
        self.volume = npy.array([h.volume for h in self.hotspots]).mean()

    def addHotSpots(self, HotSpots):
        if not isinstance(HotSpots, list): self.hotspots = [HotSpots]
        self.hotspots += HotSpots
        self.setuphp()

class HotSpotMultipleSetError(Exception):
    pass

class HotSpotMultipleSet(HotSpotSet):
    """
    Set containing multiple HotSpot Sets (from different chemical types) and 
    methods to obtain a consensus Hotspot Set (with best hotspots from any chemical type).
    """
    def __init__(self, hsetlist=False, probelist=False, probetypesmap={}, name='', info='', clustdistance=1.5, *args, **kwargs):
        """
        :args list hsetlist: List of Hot spot sets (HSets) to add.
        :args list typelist: List of names to identify types corresponding to the hsetlist. Must have same order. If not given, will try to determine type from the HotSpotSet.probe attribute.
        
        """
        HotSpotSet.__init__(self, name=name, info=info, *args, **kwargs)
        self.combinedhset = False
        self.clustdistance = clustdistance
        self.setProbeTypesMap(probetypesmap)
        if hsetlist: self.addHSets(hsetlist, probelist)
    
    def setProbeTypesMap(self, probetypesdict):
        """
        Dictionary to identify names of probes with chemical types
        """
        assert isinstance(probetypesdict, dict)
        self.probestotypes = probetypesdict
    
    def addHSets(self, hsetlist, probelist=False):
        "Add Hot Spot Sets (HSets) to the current Multiple Set"
        if hsetlist:
            if not isinstance(hsetlist, list): hsetlist = [hsetlist]
            for hi, hset in enumerate(hsetlist):
#                if not isinstance(hset, HotSpotSet): 
#                    raise HotSpotMultipleSetError, "Error adding %s. This should be a HotSpotSet object."%hset
                
                # Assign type from typelist if present or keep from HotSpotSet.probe attribute
                if probelist: hset.probe = probelist[hi]
                for i in range(hset.nhotspots): hset[i].probe = hset.probe
                
                # Add hotspots to current set
                self.addHotSpots(hset.hotspots)

    def writeCombinedHSetPDB(self, outpdb, combinedHset=False, proberesmap={}, probeatommap={}, onlycenter=False, **kwargs):
        """
        Write a PDB with the hotspots from a combined hotspot set obtained after combine method is called.
        Will write relative probability of the selected HS when combining HSets in occupancy column. Energy will be writen in bfactos column.
        
        :arg str outpdb: File name to save
        :arg combinedHSet: Hot Spot set with combined hotspots from different probes/types.
        :type combinedHset: :class:`HotSpotSet`
        :arg dict proberesmap: Map probe names in combinedHset with residue names you wish to appear in the PDB. Default: HOT.
        :arg dict probeatommap: Map probe names to atom names to print in the pdb. Default: C.
        :arg bool onlycenter: Write only the hotspot minimum energy point. Write all hotspot points if False (default).
        """
        combinedHset = combinedHset or self.combinedhset
        if not combinedHset: raise HotSpotMultipleSetError, "No combined HSet saved or given as argument"
        self.sortCombinedHSet(combinedHset)
        
         # Write a PDB for each point found
        head=[]
        out = []
        for i, hp in enumerate(combinedHset.hotspots):
            res = False
            atom = False
            if proberesmap: res = proberesmap.get(hp.probe)
            if probeatommap: atom = probeatommap.get(hp.probe)
            h, b = hp.getHotSpotPDBstr(onlycenter=onlycenter, resname=res, atom=atom, occupancy=hp.probability)
            head += [h]
            out += [b]
        with open(outpdb,'w') as k:
            k.write('\n'.join(head))
            k.write('\n')
            k.write('\n'.join(out))
        self.log.info("DONE writting PDB %s with combined hotspot set"%outpdb)

    def sortCombinedHSet(self, combinedHSet=False):
        """Sort hotspot set by increasing energy values and renumber indexes accordingly"""
        combinedHSet = combinedHSet or self.combinedhset
        if not combinedHSet: return
        combinedHSet.hotspots.sort()
        [h.setIndex(i+1) for i,h in enumerate(combinedHSet.hotspots)]
        
    def combine(self, clustdistance=None):
        """
        Combine Hotspots from different HSets in self.hsets list. Cluster them and give a score
        to the output hotspot depending on the number of partners, and chemical types.
        """
        # Cluster the hotspots
        clustdistance = clustdistance or self.clustdistance
        self.clusterHotSpots(cutDistance=clustdistance)
        RT = 0.001986*300 # kcal/mol
        resultHSet = HotSpotSet(**self.kwargs)
        
        # Work inside each cluster to select who will stay and give all information
        for i in range(1, self.nclusters+1):
            clustids = npy.where(self.clusterIndexes == i)[0]
            clusths = [self.hotspots[ci] for ci in clustids]
            energies = npy.array([c.energy for c in clusths])
            probabilities = npy.exp(energies/-RT)  # this returns number of times observed is greater than expected
#            volumes = npy.array([c.volume for c in clusths])
#            probes = npy.array([c.probe for c in clusths])
#            nhs = npy.unique(probes).size
            
            # Get the hotspot(s) with higher relative probability
            relp = probabilities / probabilities.sum()
            argsort = npy.argsort(relp)[::-1]
            relp = relp[argsort]
#            probes = probes[argsort]
            energies = energies[argsort]
#            volumes = volumes[argsort]
            clusths = npy.array(clusths)[argsort].tolist()
            for i,c in enumerate(clusths):
                c.probability = relp[i]
            keephp = clusths[0]
            hp = HotSpot(keephp.coordList, keephp.energyList, probe=keephp.probe,**self.kwargs)
            hp.others = clusths[1:]
            hp.probability = keephp.probability
            resultHSet.addHotSpots(hp)
        
        self.combinedhset = resultHSet
        self.sortCombinedHSet()
        return resultHSet
        
class CreateHotSpotSet(object):
    """
    Basic class for different hotspots construction algorithms.
    """
    def __init__(self, grid, info='', name='', *args, **kwargs):
        """Constructor will always need a grid to extract hotspots.
        Args:
            grid    (str or Grid)   Filename to load (str) or Grid Instance.

        Kwargs:
            Will depend on the implementation. Kwargs will be used in setup process.
            And usually will contain parameters to setup analysis.

        """
        self.log = logging.getLogger("CreateHotSpotSet")
        self.kwargs= kwargs
        
        if isinstance(grid, str) and os.path.exists(grid):
            self.grid = Grid(grid)
        elif isinstance(grid, Grid):
            self.grid = grid
        else:
            raise AttributeError, "grid argument should be an existing filename or a Grid instance"

        self.hotspotset = None
        self.info = info
        self.name = name
        self.setup(**kwargs)

    def setup(self, **kwargs):
        "Set parameters needed for running the algorithm"
        pass
    
    def calculate(self):
        "Method to extract hotspot set from the grid"
        pass

    def getSet(self):
        return self.hotspotset
    
class CreateBySpheres(CreateHotSpotSet):
    """
    Visit minimum energy points in the grid and get a sphere hotspot with center in the min and a given radius. Only hotspot with an energy below
    certain cutoff will be considered. 
    
    The algorithm iteratively searches the minimum energy point MIN in the grid. Then removes all surounding points at CANCELRADIUS angstroms around.
    It establishes a hotspot in MIN with radius VALSRADIUS angstroms. This hotspot is only finally considered if the total energy is below ENERGYCUT.
    """
    def setup(self, cutvalue=-0.4, valsradius=2.0, cancelRadius=2.0, protValue=1, rejectvalues=0, **kwargs):
        """
        protValue are the values for points occupied by protein (initial zeros in counts).
        If the grid is already corrected by StandardState DG, this 1 will be modified. That's why we should modify this argument in those cases.
        If not average, return the minimum value as hotspot value
        
        :args float cutvalue: Energy value to tell appart hotspots
        :args float valsradius: Distance in angstroms. Considrer surrounding values as hotspot.
        :args float cancelRadius: Radius in angstroms to tell apart different hotspots centers.
        :args float protValue: Any value used to mask the protein occupied points. All points corresponding to this value will be ignored.
        :args float rejectvalues: All grid points between MIN and VALSRADIUS with energy higher to this value will not be considered as part of the hotspot.
    
        """
        self.protValue = protValue
        self.cutvalue = cutvalue
#        self.convolvPointEnergy = convolvPointEnergy
        self.cancelRadius = cancelRadius
        self.valsradius= valsradius
        self.rejectvalues = rejectvalues
        
#        self.average = average
        self.log.debug("CreateBySpheres SETUP. energycut:%.3f valsradius:%.2f cancelRadius: %.2f protValue:%.2f"%(cutvalue,
                        valsradius,cancelRadius,protValue))

    def calculate(self):
        # Mask ones in energy grid corresponding to protein volume points
        # Extend grid data with zeros 7 points in each direction
        searchgrid = self.grid.copy()
        results = []
        index=0
        enestop = self.cutvalue +0.2 # Increase margin for loop stop. Do not miss hotspots with energy below energycut because the loop stops early
        henergy = -99999
        while henergy < enestop:
            
            # fetch minimum energy point in searchgrid
            mins = searchgrid.getMinIndex()
            # If more than one match, just take one and leave the rest for other loops
            if len(mins.shape) > 1: min = tuple(mins[0,:])
            else: min = tuple(mins)
            # Finally set all values around min to 999 to not consider the place again
            searchgrid.cancelPoints(min, self.cancelRadius, 999)
            
            # Fetch coordinates and energy values for the minimum and surrounding
            sphereindexes = self.grid.getRadialIndices(self.valsradius, min)
            coords = self.grid.origin+(sphereindexes*self.grid.delta)
            vals = npy.array([self.grid.data[tuple(i)] for i in sphereindexes])

            # Keep only negative valued coordinates
            maskneg = vals < self.rejectvalues
            if not npy.any(maskneg): continue
            vals = vals[maskneg]
            coords=coords[maskneg,:]

            # finally create a hotspot for these values and check if we keep it or not
            hp = HotSpot(coords, vals, probe=self.grid.probe, index=index, energymethod='volume',**self.kwargs)
            henergy = hp.energy
            if hp.energy <= self.cutvalue:
                index+=1
                hp.index=index
                results.append(hp)
                self.log.debug("Found hotspot: %s"%hp)#, surroundconvolv

        self.hotspotset = HotSpotSet(self.grid.probe, name=self.name, info=self.info)
        self.hotspotset.addHotSpots(results)

class createByCutoff(CreateHotSpotSet):
    """
    Creates hotspot set by taking all points below a energy cutoff value.
    A hotspot will be defined as a set of points within a cutdistance of each other
    """
    def setup(self, percentile=0.02, cutvalue=None, cutdistance=1.5, energymethod='volume', maskcutvalue=False, **kwargs):
        """
        Setup for calculation

        Args:
        :arg float percentile: Percentile of points to consider for automatic cutoff setup. The lower *percentile*% number of points will stablish the cutoff with which hotspots will be stablished.
        This is ignored if cutvalue is given.
        :arg float cutvalue: Values under or equal to this energy will be taken for hotspot building. Optional.
        :arg float cutdistance: Distance in Angstroms for deciding if a hotspot belongs to same cluster or not.
        :arg str energymethod: HS energy computing method. Same as :attr:`HotSpot.energymethod` .
        :arg float maskcutvalue: Remove all points with value *maskvalue* when calculating cutoff by percentile.

        """
        self.percentile = percentile
        self.cutvalue = cutvalue
        self.cutdistance = cutdistance
        self.energymethod = energymethod
        self.maskcutvalue = maskcutvalue
        self.log.debug("CreateByCutoff: Percentile %.2f. Cutdistance %.2f. Energymethod %s. Cutvalue %s. Maskcutvalue %s"%(percentile, cutdistance, energymethod, str(cutvalue),maskcutvalue))

    def calculate(self):
        """
        The process involves taking all points below certain *cutoff* and cluster them toghether to obtain single hotspots.
        *cutdistance* is the distance under which two neighbour points will be considered as part of the same hotspot.
        """
        data = self.grid.data
        probe = self.grid.probe
        # Stablish cutoff value if needed
        if not self.cutvalue:
            self.cutvalue = self.grid.getPercentileCutValue(self.percentile, self.maskcutvalue)
        self.log.info("Using energy cut value: %.2f"%self.cutvalue)
        
        # Identify indices for points < cutoffenergy
        maskpoints = data <= self.cutvalue
        indices = npy.vstack(npy.where(maskpoints)).T
        energies = npy.apply_along_axis(lambda x: self.grid.data[tuple(x)], 1, indices)

        # Transform indices to cartesian coordinates
        spacing = self.grid.delta
        origin = self.grid.origin
        coords = npy.apply_along_axis(lambda x: (x*spacing)+origin, 1, indices)

        # Cluster by cutdistance
        t0=time.time()
        self.log.debug("Clustering %d points..."%len(coords))
#        print "Building DM"
        self.dm = distance.pdist(coords)
#        print "Linkage and fcluster"
        self.linkage = hierarchy.linkage(self.dm)
        self.clusterIndexes = hierarchy.fcluster(self.linkage, t=self.cutdistance, criterion='distance')
        self.uClusterIds = npy.unique(self.clusterIndexes)
        self.nclusters = self.uClusterIds.size

#        print self.clusterIds, npy.min(self.clusterIds). npy.max(self.clusterIds)

        # Create a hotspot for each cluster
        results = []
        self.log.debug( "Found %d clusters"%self.nclusters)
        self.log.debug( "Time:",time.time()-t0)
        for clustid in self.uClusterIds:
            clustMask = self.clusterIndexes == clustid
            clustCoords = coords[clustMask,:]
            clustEnergies = energies[clustMask]
            results.append(HotSpot(clustCoords, clustEnergies, energymethod=self.energymethod,
                                        probe=self.grid.probe, index=clustid, spacing=spacing, **self.kwargs))

        # Renumber hotspots. 1st HS should be the lowest energy one
        results.sort()
        [h.setIndex(i+1) for i,h in enumerate(results)]
        self.hotspotset = HotSpotSet(self.grid.probe, name=self.name, info=self.info,**self.kwargs)
        self.hotspotset.addHotSpots(results)

#: Current hot spot set calculation methods
HS_CREATE_METHODS = {'spheres':CreateBySpheres, 'cutoff':createByCutoff}

class HotSpotsManager(object):
    "Operate over HotSpotSets and HotSpots: Creation from  energy grids and so on."
    def __init__(self):
        self.log = logging.getLogger("ProjectManager.HotSpotsManager")
        self.info = ''
#        self.setup()

    def __iter__(self):
        for h in self.hotspotSet.hotspots:
            yield h

    def __getitem__(self, i):
        return self.hotspotSet.hotspots[i]
    
    def createHotSpotSetFromGrid(self, grid, method, name='', info='', **kwargs):
        """
        From a Grid, extract hotspots.

        Args:
            method  (str)       Hotspot creation algorithm.
                                Available: %s
                                Each method requires different kwargs, fetch info in the setup method of each class.
            grid    (Grid instance) Grid over which to calcualte hotspotset
            name    (str)       Name to identify the hotspotset to be created.
            info    (str)       Optional. Extra info to attach to the hotspotset
            

        """%HS_CREATE_METHODS.keys()
        createhotspots = HS_CREATE_METHODS.get(method)(grid, info=info, name=name, **kwargs)
        self.log.info("Creating hotspot from grid %s. Method: %s."%(grid, method))
        createhotspots.calculate()
        self.hotspotSet = createhotspots.getSet()
        return self.hotspotSet

    def saveHotSpot(self, picklefileout, hotspot):
        "Save into pickle the hotspot"
        import cPickle
        cPickle.dump(hotspot, open(picklefileout, 'wb'))

    def saveHotSpotSet(self, picklefileout, hset=False):
        "Save into pickle the hotspotset found"
        import cPickle
        if hset:
            cPickle.dump(hset, open(picklefileout, 'wb'))
        elif not hset and self.hotspotSet:
            cPickle.dump(self.hotspotSet, open(picklefileout, 'wb'))
        else:
            self.log.error("No hotspotSet stored in HotSpotsManager and no HSet given. Cannot save anything.")
        
    def loadHotSpotSet(self, picklefilein):
        "Load previously saved hotspots from pickle file"
        import cPickle
        self.hotspotSet = cPickle.load(open(picklefilein, 'rb'))

    def writeHotSpotSetPDB(self, outpdbfile, hset=False, **kwargs):
        # Write a PDB for each point found
        if hset and isinstance(hset, HotSpotSet):
            self.hotspotSet = hset
        self.hotspotSet.writeHotSpotsPDB(outpdbfile, **kwargs)

    def getProteinAtomsAroundHotspot(self, hotspot, pdb, fromCentroid=True, radius=6, includeH=True, returnmask=False):
        """
        Get pdb indices for the atoms neighbouring the hotspot around a given distance. Distance will be calculated from the hotspot centroid or considering all points insed depending on fromCentroid value.
        Also it is optional to include the hydrogens (default) or remove them from the returned list.
        
        Args:
            hotspot     (HotSpot or HotSpotSet)     Fetch atoms around this hotspots
            pdb         (PDBModel or string)        What PDB to return atoms from. Should be a Biskit.PDBModel or a string pointing to a correct pdb file to read.
            fromCentroid    (bool)                  Calculate from centroid point or considering all points in the hotspot?
            radius          (float)                 Radius around we consider an atom as neighbour
            includeH        (bool)                  Include hydrogens in returen indices?
            returnMask      (bool)                  Return as mask instead of list of indices?
        
        Returns:
            indicesList     (list of ints)          Indices in the pdb corresponding to neighbouring atoms
            or
            indicesMask     (npy.ndarray of bools)  Mask for the PDBModel
        """
        import Biskit as bi
        import scipy
        if scipy.__version__ > "0.11.0":
            from scipy.spatial import KDTree as KDTree
        else:
            from scipy.spatial import KDTree
        
        # Load PDB if string
        # and get Hydrogens mask
        if isinstance(pdb, str): pdb = bi.PDBModel(pdb)
        hids = npy.where(pdb.maskH())[0]
        
        # First construct a KDTree with protein coordinates.
        tree = KDTree(pdb.xyz) # 3dims, 1 atom per bucket
        
        # Find neighbours
        # check first if hotspot has multiple coords or just one
        if hotspot.coordList.ndim == 1:
            dims = 1
        else:
            dims = 2

        if fromCentroid or dims == 1:
            # Calculate n'bours to one point
            rawids = tree.query_ball_point(hotspot.coord, radius)
            if not includeH:
                ids = [ri for ri in rawids if ri not in hids]   # check index is not a hydrogen atom 
            else:
                ids = rawids
        else:
            # Do a search for each point and get a unique set of n'bours
            rawids = tree.query_ball_point(hotspot.coordList, radius)
            if not includeH:
                ids = [ri for ri in rawids if ri not in hids]   # check index is not a hydrogen atom
            else:
                ids = rawids
        
        if returnmask:
            m = [i in ids for i in range(len(pdb.xyz))]
            return  m
        else:
            return ids

def createHotSpotsByCutoff(gridlist, percentile=0.02, cutoff=None, outprefix=None, onlycenter=False, **kwargs):
    """
    Create a combined hot spot set from all input grids using cutoff hot spot creation method.

    :arg list gridlist: List of :class:`GridsManager.Grid` instances.
    :arg float percentile: Lower percentile of points to take to automatically stablish a cutoff value for each grid.
    :arg float cutoff: Use this hard cutoff value and ignore percentile.
    :arg str outprefix: If given, will save a PDB file and a pickle file containing the hot spot set found.
    :arg bool onlycenter: Write only hotspots centers in the output pdb.
    
    :return: :class:`HotSpotMultipleSet` with results. If outprefix is given, will save two files with this prefix.
    """
    if not isinstance(gridlist, list): gridlist = [gridlist]
    hm = HotSpotsManager()
    multiset = HotSpotMultipleSet(**kwargs)
    hsets = []
    for grid in gridlist:
        gset = hm.createHotSpotSetFromGrid(grid, method='cutoff', percentile=percentile, clustdistance=1.5, cutvalue=cutoff, **kwargs)
#        gset.clusterHotSpots(1.5)
        #gset.reduceClusters(update=True)
        hsets.append(gset)
    multiset.addHSets(hsets)
    multiset.combine(1.5)
    if outprefix:
        multiset.writeCombinedHSetPDB(outprefix+'.pdb', onlycenter=onlycenter)
        multiset.combinedhset.dump(outprefix+'.hset')
    return multiset.combinedhset