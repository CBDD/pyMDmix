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
__date__ ="$Mar 15, 2014 6:24:25 PM$"

import os
import os.path as osp
import copy
import numpy as npy
import logging

import settings as S
import GridData 

class GridError(Exception):
    pass

class TypeError(GridError):
    pass

class ProbeError(GridError):
    pass

class BadFile(GridError):
    pass

class ValueError(GridError):
    pass

class Grid(GridData.GridData):
    "Grid information for a single atomtype/probe. It wraps GridData instance to add functionalities."
    def __init__(self, fname=False, probe="", type="", name="", description="", info="", *args, **kwargs):
        self.probe = probe
        self.headerinfo = info
        self.name = name
        self.type = type
        self.description = description
        self.typeheader = ''
        GridData.GridData.__init__(self, fname=fname, *args, **kwargs)
        if fname: 
            self.getProbeFromHeader()
            self.getTypeFromHeader()
            if self.probe == 'UNK':
                self.getProbeFromName(fname)
            if self.type == 'MDMIX_UNK':
                self.getTypeFromPath(fname)
        if type: self.setType(type)
        if probe: self.setProbe(probe)
    
    def setType(self, type, info=False):
        if type in S.GRIDTYPES:
            self.type = type
            self.setHeader(info=info)
        else:
            raise TypeError, "%s type not valid. Valid types are: %s"%(type, S.GRIDTYPES)

    def getPercentileCutValue(self, percentile, maskout=False):
        """
        Helper function. Will return the approximate cutvalue that will contain the percentile X number of points.
        Args:
            percentile  (float) Percentile value
            maskout     (float) If given, remove all points with this value from percentile calculation
        Returns:
            float. Cutvalue.
            
        Example:
            >> grid.getPercentileCutValue(0.1)
            >> -1.5
            >> 100*(npy.sum(grid.data < -1.5))/float(grid.data.size)
            >> 0.1
        """
        if maskout: d = npy.compress(self.data.flat!=maskout, self.data.flat)
        else: d=self.data.flat
        return npy.percentile(d, percentile)

    def getTypeFromHeader(self):
        "Look in grid file header to detect the type: MDMIX_DENS MDMIX_RAW MDMIX_CORR MDMIX_OTHER"       
        for gtype in S.GRIDTYPES:
            if self.header is not None and gtype in self.header:
                self.type = gtype
                return gtype
        # No type info found. Use unknown type MDMIX_UNK
        # And rebuild header
        self.type = 'MDMIX_UNK'
        self.setHeader()
        return self.type
    def getTypeFromPath(self, fname):
        """Fetch the grid type from grid path. It should be determined by the last folder where the grid is.
        Only called if getTypeFromHeader does not work.
        """
        partitioned_path = fname.split('/')
        grid_directory = partitioned_path[-2]
        path2type_map = {'dgrids': 'MDMIX_DENS', 'PROBE_AVG':'MDMIX_RAW_AVG', 'egrids':'MDMIX_RAW'} # modify if there is more
        if grid_directory in path2type_map:
            self.type = path2type_map[grid_directory]   
    def setProbe(self, probe):
        "Change probe name and adapt header info"
        self.probe = probe
        self.setHeader()
    
    def getProbeFromHeader(self):
        "Fetch probe info in header. It should be the second word in the line."
        if self.header is None:
            self.probe = 'UNK'
            return "UNK"

        self.probe = self.header.split()[1]
        return self.probe
        
    def getProbeFromName(self, fname):
        """Fetch the probe from grid name. It should be the last two words before the file extension.
        Only called if getProbeFromHeader does not work."""
        elements_withoutextension = os.path.splitext(fname)[0].split('_')
        self.probe = '_'.join(elements_withoutextension[-2:])
        return self.probe

    def setHeader(self, info=False):
        if not info: info=self.headerinfo
        else: self.headerinfo = info
        self.typeheader = self.type+' '+self.probe
        self.typeheader+= ' '+info
    
    def getMinIndex(self):
        """
        Get the grid indices for the minimum energy value(s).
        :returns: Array Nx3 with index of minimum energy points. N > 1 if minimum energy value in more than 1 point.
        :rtype: :class:`numpy.ndarray`
        """
        mins = npy.where(self.data == self.data.min())
        return npy.vstack(mins).T
        
    def writeDX(self, outname):
        "Overwrite gr.writeDX function to include file existance checking."
        GridData.GridData.writeDX(self, outname, self.typeheader)
        if not os.path.exists(outname):
            return False
        return True
    
    def writeXPLOR(self, outname):
        "Overwrite gr.writeXPLOR function to include file existance checking."
        GridData.GridData.writeXPLOR(self, outname, self.typeheader)
        if not os.path.exists(outname):
            return False
        return True

    def getIndexFunction(self):
        """Return a simplified function for converting coordinates to this grid indices
        ::
            >>> grid = NewGrid([10,10,10]) # Create grid with shape 10x10x10 centered at origin
            >>> fx = grid.getIndexFunction()
            >>> fx([1,2,1])  # find index at coordinate [1,2,1]
            (6, 8, 6)
            >>>

        """
        spacing = npy.array(self.delta)
        origin = npy.array(self.origin)
        shape = npy.array(self.shape)
        def getIndex(coords, spacing, origin, shape):
            ind = ((coords - origin) / spacing).astype(int)
            # Check that index obtained is inside the grid.
            if npy.any(ind < 0.) or npy.any(ind >= shape):
                return False
            else:
                return tuple(ind.astype(int).tolist())
        return lambda x: getIndex(x,spacing,origin,shape)

class NewGrid(Grid):
    def __init__(self, dimensions=None, origin=None, spacing=S.GRID_SPACING, probe='probe',type='MDMIX_UNK', *args, **kwargs):
        "Create a new grid instance specifying the dimensions in x,y,z direction in angstroms.\
        If origin=None, the origin will be 0,0,0 centered in the middle of the grid."
        if isinstance(dimensions, int) or isinstance(dimensions,float): dimensions = [dimensions]*3
        if not isinstance(spacing, list): spacing = [spacing]*3
        dimensions = npy.array(dimensions, dtype=float)
        if not npy.any(origin):
            origin = -npy.array(dimensions) / 2.
        shape = (dimensions/spacing).astype(int)
        Grid.__init__(self, fname=None, probe=probe, type=type, shape=shape,
                                    origin=origin, spacing=spacing, *args, **kwargs)
        

class GridSpace(object):
    "Set of Grid objects belonging to the same simulated structure. \
    Containts one or different chemical types and functionalities to work traversing all of them."
    kB = 1.38065E-23*6.022142E23/4.184/1000.0 # kcal/mol = R

    def __init__(self, GridList=None, spaceName=None, filterPositive=False, T=300., R=1.5, *args, **kwargs):
        "GridList: list of Grid instances or filepaths (also a list with 1 object is valid).\
        If GridList == None, then no action is done, waiting for a load command to restore from pickle.\
        If filterPositive is True, all values > 0 will be set to 0.\
        T is temperature. Used when doing boltzmann averages. Can be set later with setT().\
        Same for R, radii in angstroms to consider when fetching points. Set later with setR()."
        self.log = logging.getLogger("GridSpace(%s)"%spaceName)
        self.spaceName = spaceName
        self._defaultmode = 'point' # Calculation mode
        self._defT = T
        self._defR = R
        self.probelist = []
        self.probeMapping = {}
        if GridList:
            self.loadGrids(GridList, filterPositive)
            self.nGrids = len(GridList)
        else:
            self.log.debug("Empty GridSpace created, waiting for a pickle loading: GridSpace.load(picklefile)")

    def __getstate__(self):
        d = self.__dict__.copy()
        del d['log']
        return d

    def __setstate__(self, d):
        d['log'] = logging.getLogger("ProjectManager.GridSpace(%s)"%d['spaceName'])
        self.__dict__.update(d)

    def dump(self, filename):
        self.log.debug("Dumping GridSpace to pickle file: %s"%filename)
        import cPickle
        cPickle.dump(self, open(filename, 'wb'))
        self.log.debug("DONE")

    def load(self, filename):
        self.log.debug("Loading GridSpace from pickle file: %s"%filename)
        import cPickle
        self = cPickle.load(open(filename, 'rb'))
        self.log.debug("DONE")

    def addDegeneracy(self, probeDegeneracyDict):
        """Add more probes to the same dimension(s). Useful to redirect probes
        we are missing and assign to a similar chemical type: eg. redirect aromatics to hydrophobics.
        
        probeDegeneracyDict:    dictionary with the form {existingProbe: [newProbeName, newProbename,...],...}
        
        No return, just modification of self.probeMapping
        """
        self.__buildMapping()   # Rebuild clean mapping from loaded grids
        for probe, newprobes in probeDegeneracyDict.iteritems():
            for np in newprobes:
                if self.probeMapping.has_key(np): self.log.warn("Overwriting existing Mapping Key: %s"%np)
                self.probeMapping[np] = self.probeMapping[probe]
        self.__buildReverseMapping()
        return True

    def reduce(self, probes=None):
        """This method builds a new Grid instance with a reduced data. It is to get minimum value over all points from all grids in the 4th dimension (or the selected probes) and obtain a 3D grid.
        'probes' is a list of probe names we are interested in for the reduce operation, rejecting those not selected.
        This Grid returned bears two new attributes:
            - chemotype: a npy.array with same shape as data indicating the origind dimension where the point was extracted from
            - probeMapping: same as GridSpace.probeMapping to identify the chemotypes
        """
        self.log.info("Reducing GridSpace...")

        # Selection of dimensions corresponding to probe names in probe list
        # or all dims by default
        if probes:
            dims = []
            for p in probes:
                k = self.probeMapping.get(p)
                if k: [dims.append(i) for i in k]
                else:
                    self.log.error("(reduce) Invalid probe name: %s"%p)
                    raise ProbeError, "Invalid probe name: %s"%p
            dims = tuple(npy.unique(dims).tolist())
            gspace = self.gspace[:,:,:,dims]
            ndims = gspace.shape[-1]+1
        else:
            gspace = self.gspace    #all dims
            ndims = self.ndim

        # Proceed with reduction operation
        mins = npy.min(gspace, axis=3)
        minMask = []
        [minMask.append(gspace[:,:,:,d] == mins) for d in range(ndims)]
        index = npy.vstack(npy.where(minMask)).T
        reducedGridData = npy.zeros_like(gspace[:,:,:,0])
        reducedChemoType = npy.zeros_like(gspace[:,:,:,0])

        # This is the operation where all min points are fetched and
        # corresponding dims stored. This implementation will be very slow.
        # Should find faster ways :S TODO
        for i in index:
            reducedGridData[tuple(i[1:])] = gspace[tuple(i)]
            reducedChemoType[tuple(i[1:])] = i[0]

        # Finally build and return the new reducedGrid instance
        reducedGrid = self.container.copy()
        reducedGrid.update(reducedGridData)
        reducedGrid.chemotype = reducedChemoType
        reducedGrid.probeMapping = self.probeMapping
        reducedGrid.solvent = 'Reduced Grid'
        reducedGrid.probe = 'Reduced Grid'
        self.log.info("DONE")
        return reducedGrid

    def sum(self):
        """
        Obtain a single grid resulting from the sum of all the individuals.
        """
        reducedGridData = npy.sum(self.gspace, axis=3)
        reducedGrid = self.container.copy()
        reducedGrid.update(reducedGridData)
#        self.log.info("DONE")
        return reducedGrid

    def __buildMapping(self):
        "Build a dictionary mapping probe names with the 4D-array dimension\
        (1:1 or 1:N when multiple grids for same probe, eg: different replicas, same probe)"
        probes = [grid.probe for grid in self.grids]
        self.probelist = probes
        for i in range(self.ndim):
            if self.probeMapping.has_key(probes[i]):
                self.probeMapping[probes[i]] += [i]
            else:
                self.probeMapping[probes[i]] = [i]
        return True

    def __buildReverseMapping(self):
        # REVERSE DICTIONARY TO FIND TYPES AND NAMES FROM DIM POSITION
        # If degeneracy, same dimension index may correspond to different probes
        map_dict = copy.deepcopy(self.probeMapping)
        for n, d in self.probeMapping.iteritems():
            if isinstance(d, list):
                for i in d:
                    if not map_dict.has_key(i):
                        if isinstance(n, list):
                        	map_dict[i] = n
			else: map_dict[i] = [n]
                    else:
                        if isinstance(n, list): [map_dict[i].append(e) for e in n]
                        else: map_dict[i].append(n)
            else:
                i = d
                if not map_dict.has_key(i):
                    if isinstance(n, list):
                    	map_dict[i] = n
		    else: 
			map_dict[i] = [n]
                else:
                    if isinstance(n, list): [map_dict[i].append(e) for e in n]
                    else: map_dict[i].append(n)
        self.probeMapping = map_dict
        return True

    def loadGrids(self, GList, filterPositive=False):
        "Build a set of grids (a GridSpace) with all the grids in GList"
        self.grids = []
        self.log.debug("Loading grids and creating GridSpace...")
        # Check all members of GList are Grid instances or convert otherwise
        for i,g in enumerate(GList):
            if isinstance(g, GridData.GridData) or isinstance(g, Grid):
                self.grids.append(g)
            else:
                self.grids.append(Grid(g))

        # Set all positive values to zero if filterPositive
        if filterPositive:
            self.log.debug("Ignoring positive values in all grids (all set to zero) (filterPositive=True. If that's not the desired behaviour, use filterPositive=False.")
            for g in self.grids:
                g.data[g.data > 0] = 0

        # Set arguments about the grids and
        # Obtain same shape, origin for all of them
        self.ndim = len(self.grids)
        if self.ndim > 1:
            #Check if they differ in shape and origin to trim them
            shapes = npy.array([g.data.shape for g in self.grids])
            origins = npy.array([g.origin for g in self.grids])
            if npy.any(origins != origins[0]) or npy.any(shapes != shapes[0]):
                self.log.info("Different shapes and/or origins between all grids in Grid Space. Will try to trim them to normalize and match the space.")
                trimmed = trim(self.grids)
            else:
                trimmed = self.grids
        else: trimmed = self.grids
        self.sources = [g.source for g in self.grids]
        self.origin = trimmed[0].origin
        self.spacing = trimmed[0].delta
        self.shape = npy.array(trimmed[0].data.shape)
        self.size = trimmed[0].data.size

        # Create unique npy.array with all the grids (4D array)
        a = npy.zeros(list(self.shape)+[self.ndim], dtype='float')
        for i,g in enumerate(trimmed):
            a[:,:,:,i] = g.data
        self.gspace = a

        # Finally store one of the grid instances
        # to allow possible saving of DX data or easy Index-Cartesian conversions
        self.container = trimmed[0].copy()
#        self.container.data = []
        self.container.source = ""
        # Remove data from original grids to save memory
        for g in self.grids: del g.data
        self.log.debug("DONE")
        self.__buildMapping()
        self.__buildReverseMapping()
        return True
    
    def getValues(self, coord, name=None, ndim=None, cross=False, ownF=False, ignoreValue=False,
                    choose='min', r=None, mode=None):
        """Return grid value at point indicated by cartesian coordinates 'coord'.
         'name' or number of the dimension ('ndim') must be given to identify what grid to use.
         
         'name' can be either a List or a string. If a list is given, the value returned will be
         the 'choose' function of the two independent values of each element in the list. Same to 'ndim' with ints.
         'choose' options:
            - min   str       min of the values
            - max   str       max of the values
            - False bool      return all values
            
         Optionally, if 'cross' is True, values will be returned for all the dimensions (grids).
         'ownF' can be a function to operate over the data points extracted, as if it was mean, boltz or whatever.
         'ownF' overrides 'mode' and is always used over the values closer to 'r'.
         Different modes can be used:
            - point     Return the point value.
            - volmean   Return volume averaged energy (Hotspot energy) on radii 'r'.
            - avg       Return the mean of all values around 'r' angstroms from 'coord'.
            - boltz     Return the boltzman average of all values around 'r' angs from 'coord'.
            - min       Return MIN value of all values 'r' around 'coord'.
            - max       Return MAX value of all values 'r' around 'coord'.
        
        'ignoreValue' allows the specification of any value that should be ignored when computing the values.
        For instance, if we have protein occupied voxels with 999, those values will be ignored for computing 
        averages, min, etc... functions.
        """
        if not r: r = self._defR
        if not mode: mode = self._defaultmode

        if not ndim and not name and not cross:
#            raise ValueError, "name or atype or ndim argument must be provided."
            self.log.warn("Some atoms do not have a recognised type!")
            return 0    # Allow some atoms not having any type identified to not contribute and don't make this crash

        idx = self.toIndex(coord)
        if not idx: return False

        # Choose way of mapping the grids
        # Either give a number/s of the dimension
        # Or a list or str with name in the mapping
        # Or cross for returning all values
        if ndim:
            if isinstance(ndim, int):
                if ndim < self.dim:
                    ndim = [ndim]
                else:
                    ndim = False
            elif isinstance(ndim, list):
                if npy.any([n >= self.dim for n in ndim]):
                    ndim = False
                else:
                    ndim = ndim
            if not ndim:
                raise ValueError, "ndim out of maximum dimension"
        elif name:
            if isinstance(name, str):
                ndim = [self.probeMapping.get(name)]
            elif isinstance(name, list):
                ndim = [self.probeMapping.get(n) for n in name]
#                print name, ndim
            if npy.any([el == None for el in ndim]):
                raise ValueError, "%s name not in GridSpace mapping dict"%name
        elif cross:
            ndim = False
        else:
            raise ValueError, "valid ndim, name or atype argument must be given. Else choose 'cross'"

        # Choose options
        if choose == 'min':
            fchoose = npy.min
        elif choose == 'max':
            fchoose = npy.max
        elif choose == 'mean':
            fchoose = npy.mean
        else:
            # If other values, retur All
            fchoose = lambda x: x

        # WORK WITH PASSED FUNCTION
        if ownF:
            process = ownF
        # JUST THE POINT
        elif mode == 'point':
            if cross:
                p = self._vals(idx, ndim, radii=0)                
            else:
                p = npy.array([self._vals(idx, d, radii=0) for d in ndim])
            
            p = p.reshape((len(ndim),))
            if ignoreValue and p==ignoreValue: p = 0
            return p
            
        # CHOOSE PROCESSING MODE
        elif mode == 'avg':
            process = npy.mean
        elif mode == 'volmean':
            kBT=GridSpace.kB*self._defT
            def process(x):
                if x.sum() == 0 or not npy.any(x): return 0
                return -kBT*npy.log(npy.exp(x/-kBT).mean())
        elif mode == 'boltz':
            process = self._boltz
        elif mode == 'min':
            process = npy.min
        elif mode == 'max':
            process = npy.max
        else:
            raise ValueError, "Invalid mode, should be point, avg, boltz, min or max"

        # Get values
        values = [self._vals(idx, d, r) for d in ndim]

        # Process in cross or normal mode
        if cross:
            result = []
            for d in range(self.ndim):
                vals = values[0][:,:,:,d]
                if ignoreValue or type(ignoreValue) == int or type(ignoreValue) == float: # allow zeros as ignoredvalues
                    v = npy.ma.masked_equal(vals, ignoreValue)
                    v = v.flatten()
                    if npy.any(v.mask): v = v.compress(~v.mask)
                else:
                    v = vals
                result.append(process(v))
            return result
        else:
            if ignoreValue or type(ignoreValue) == int or type(ignoreValue) == float:
                result = []
                for val in values:
                    v = val.flatten()
                    v = npy.ma.masked_equal(v,ignoreValue)
                    if npy.any(v.mask): v = v.compress(~v.mask)
                    result.append(process(v))
            else:
                result = [process(val) for val in values]
#            print values, result, fchoose(result)
            return fchoose(result)

    def setMode(self, mode):
        "Set de default processing mode of the values obtained with XYZ. Default is 'point'.\
         See getXYZvalues for more information on the modes."

        if mode in ('point','avg','boltz','min','max','volmean'):
            self._defaultmode = mode
        else:
            print "Warning mode ",mode," is invalid."

    def setT(self, T):
        self._defT = T

    def setR(self, R):
        self._defR = R

    def getDimInfo(self, i=None):
        info = []
        if i != None:
            return self.probeMapping.get(i)

        for i in range(self.ndim):
            info.append(self.probeMapping[i])[0]

        return info

    def averageGrids(self, dims=False, names=False, mode='boltz'):
        "Return a new grid instance with the average of the grids specified by 'dims' or 'names',\
        if 'dims'/'names' is False, all dims are used.\
        It is possible to change averaging mode between:\
            - boltz - boltzman average\
            - min   - just minimum value (not really averaging)\
            - mean  - arithmetic mean."
        # Set averaging mode
        if mode == 'boltz':
            func = self._boltz
        elif mode == 'mean':
            func = npy.mean
        elif mode == 'min':
            func = npy.min
        else:
            self.log.error("Invalid averaging mode chosen (should be 'boltz','mean' or 'min').")
            self.log.error("Using default 'boltz'")

        if not dims and not names:
            select = False
        elif dims:
            dims = tuple(dims)
            if len(dims) <=1:
                self.log.warn("Cannot average a single grid!")
                return False
            select = dims
        elif names:
            # Identify dims for names given
            if not isinstance(names, list):
                self.log.warn("Cannot average a single grid!")
                return False
            dims = tuple([self.getDimInfo(n) for n in names])
            select = dims

        # Extract dimensions if needed
        if select: grids = self.extractDim(select)
        else: grids = self.gspace

        # Average along 4th-Dimension
        avgGridData = func(grids,axis=3)
        avgGrid = self.container.copy()
        avgGrid.update(avgGridData)
        return avgGrid

    def _boltz(self, x, axis=None):
        if not isinstance(x, npy.ndarray): x = npy.array(x)
        kBT=GridSpace.kB*self._defT
        expVals = npy.exp(x/-kBT)
        return npy.sum(x*expVals,axis=axis)/npy.sum(expVals,axis=axis)

    def __getitem__(self, xyz):
        return self.getXYZvalue(xyz, cross=True)

#    def _vals(self, index, dim, radii=0):
#        "Return value of the grid specified in dimension 'dim'. Return values of 'r' angstroms around\
#        if radii != 0. Finally, if dim is False, return the values for all the grids (CROSS)."
#        npoints = npy.round(radii/npy.array(self.spacing))
#
#        # Avoid negative indexes
#        maskout = (index - npoints) < 0
#        npoints[maskout] = npy.take(index, npy.where(maskout)[0])
#
#        # Avoid indexes that exceed grid dimensions
#        maskout = npoints >= self.shape
#        npoints[maskout] = npy.take((self.shape - 1), npy.where(maskout)[0])
#
#        i,j,k = index
#        a,b,c = npoints
#
#        if dim is not False: return self.gspace[i-a:i+a+1,j-b:j+b+1,k-c:k+c+1,dim]
#        else:   return self.gspace[i-a:i+a+1,j-b:j+b+1,k-c:k+c+1,:]

    def _vals(self, index, dim, radii=0):
        "Return value of the grid specified in dimension 'dim'. Return values of 'r' angstroms around\
        if radii != 0. Finally, if dim is False, return the values for all the grids (CROSS)."
        if radii: points = self.container.getRadialIndices(radii, point=index)
        else: points = [index]
        
        vals = []
        if dim is not False:
            space = self.gspace[:,:,:,dim]
            [vals.append(space[tuple(ind)]) for ind in points]
        else:
            for dim in range(self.gspace.shape[-1]):
                space = self.gspace[:,:,:,dim]
                [vals.append(space[tuple(ind)]) for ind in points]
        return npy.array(vals)
    
    def toIndex(self, xyz):
        return self.container.getIndex(xyz)

    def toCartesian(self, ijk):
        return self.container.getCartesian(ijk)

    def extractDim(self, dim):
        "Return npy.ndarray with the data in 'dim' dimension/s (if a tuple is given, multiple dimensions are extracted)"
        return self.gspace[:,:,:,dim]

# Auxiliary independent functions
def getEnergyValues(grid, coords, radius=0, temp=300.):
    """
    From input grid *grid* read value at coordinates in *coord* array. If *radius* parameter is given,
    return the average value at *coord* and all points *radius* angstroms around. *radius* can be a list of same length
    as coordinates for using different radius at each coordinate. By default, if *radius* is zero, only the value at *coord* will be returned.
    If no value is found, npy.nan will be returned.

    :arg grid: Grid to fetch values from
    :type grid: :class:`Grid`
    :arg coords: Array Nx3 with xyz coordinates
    :type coords: :class:`npy.ndarray`
    :arg radius: Radius around each coordinate to consider.
    :type radius: float or list of floats. If float, use same radius for all points. If list of floats, should be of shape Nx3.
    :arg float temp: Temperature for boltzmann relationship when averaging energy values.

    :return: numpy array with lenght N
    """
    g = Grid(grid)
    results = []
    if isinstance(radius, float): radius = [radius]*len(coords)

    RT = 0.0019857*temp
    for ci,coord in enumerate(coords):
        i = g.getIndex(coord)
        if npy.any(i):
            r = radius[ci]
            if r:
                vals = g.getSphereValues(coord, r)
                vals = vals[vals != 999.]
                if vals.size:
                    e = -RT*npy.log(npy.exp(vals/-RT).mean())
                else:
                    e = npy.nan
            else:
                e = g.data[i]
                if e == 999.: e = npy.nan
        else:
            e = npy.nan

        results.append(e)

    return npy.array(results)

def getEnergyFromTxtCoords(grid, txtfile, radius=0, temp=300.):
    """
    Using a txt file for defining coordinates and radius to extract energies from grid *grid*.
    """
    if osp.exists(txtfile): txt = open(txtfile, 'r').readlines()
    else: raise BadFile, "File name %s not found."%txtfile

    # Collect coordinates and radius
    coords = []
    rads = []
    for line in txt:
        line = map(float, line.strip().split())
        if len(line) == 3:
            x,y,z = line
            coords.append([x,y,z])
            rads.append(radius)
        elif len(line) == 4: # contains radius
            x,y,z,r = line
            coords.append([x,y,z])
            rads.append(r)
        else:
            raise BadFile, "File %s contains lines with wrong format. Expected 3 or 4 elements per line."%txtfile

    # Fetch results
    return getEnergyValues(grid, coords, rads, temp=temp)

def getEnergyFromPDBCoords(grid, pdbfile, forceradius=0, temp=300.):
    """
    Using a PDB file for defining coordinates and radius (in angstroms) to extract energies from grid *grid*.
    Radius will be defined at b-factor column in the PDB. A radius of **zero** will be considered as
    the point energy. If *forceradius* is defined diferent than zero, this value will override b-factor column.

    """
    import Biskit as bi
    if osp.exists(pdbfile): pdb = bi.PDBModel(pdbfile)
    else: raise BadFile, "File name %s not found."%pdbfile

    # Collect coordinates and radius
    coords = pdb.xyz
    if forceradius: rads = forceradius
    else: rads = pdb['temperature_factor']
 
    # Fetch results
    return getEnergyValues(grid, coords, rads, temp=temp)

def gridDifference(grid1, grid2, outname):
    "Save a grid with the difference in data from grid1-grid2 with name outname."
    g1=Grid(grid1)
    g2=Grid(grid2)
    if not npy.all(npy.array(g1.shape) == g2.shape) or not npy.all(npy.array(g1.origin) == g2.origin):
#        print g1.shape, g2.shape
#        print g1.origin, g2.origin
        g1 = g1.trim(g2)
    
    g1.data[g1.data>0] = 0
    g2.data[g2.data>0] = 0
    
    diffdata = g1.data - g2.data
    outg = g2.copy()
    outg.update(diffdata)
    outg.source="Difference grid %s - %s"%(grid1, grid2)
    outg.type="MDMIX_RAW"
    outg.probe="UNK"
    if 'xplor' in outname: outg.writeXPLOR(outname)
    else: outg.writeDX(outname)
    print "Saved %s with difference in grids (%s - %s)"%(outname, grid1, grid2)

def gridSum(grid1, grid2, outname):
    "Save a grid with the sum in data from grid1+grid2 with name outname."
    g1=Grid(grid1)
    g2=Grid(grid2)
    if not npy.all(npy.array(g1.shape) == g2.shape) or not npy.all(npy.array(g1.origin) == g2.origin):
#        print g1.shape, g2.shape
#        print g1.origin, g2.origin
        g1 = g1.trim(g2)
    
    sumdata = g1.data + g2.data
    outg = g2.copy()
    outg.update(sumdata)
    outg.source="Sum grid %s + %s"%(grid1, grid2)
    outg.type="MDMIX_RAW"
    outg.probe="UNK"
    if 'xplor' in outname: outg.writeXPLOR(outname)
    else: outg.writeDX(outname)
    print "Saved %s with sum of grids (%s + %s)"%(outname, grid1, grid2)

def trim(*Glist):
    """
    Trim all grids in Glist.
    
    :arg list Glist: List containing more than 1 Grid instances to crop.
    :returns: a list of trimmed grid instances
    """     
    if len(Glist)<2:
        if type(Glist[0]) is list:
          Glist = Glist[0]
        else:
          return "ERROR. Must provide at least 2 grids."
        
    deltaList = npy.array( [grid.delta for grid in Glist] )
    
    if npy.any(deltaList != deltaList[0]): 
      return "ERROR. All grids should have same spacing."
    
    delta = deltaList[0]
    originList = npy.array( [grid.origin for grid in Glist] )
    maxCoordList = npy.array( [grid.getCartesian(grid.data.shape) for grid in Glist] )
    
    #Get maximum Origin and Cartesian Coordinates for minimum shape
    newOrigin = originList.max(axis=0)
    newMaximum = maxCoordList.min(axis=0)
    newShape = ( newMaximum - newOrigin ) / delta
    newShape.astype('int')
    
    trimmedList = []
    
    for grid in Glist:
        newarray = Grid(shape=newShape, origin=newOrigin, spacing=delta)
        newarray.source = 'Trimmed array of %s'%(grid.source)
        #Copy data to the new array
        init = grid.getIndex(newarray.origin)
        newarray.update( grid.data[ init[0]:init[0]+newarray.data.shape[0], init[1]:init[1]+newarray.data.shape[1], init[2]:init[2]+newarray.data.shape[2] ].copy() )
        if hasattr(grid, 'type'): newarray.setType(grid.type)
        else: newarray.setType('MDMIX_UNK')
        if hasattr(grid, 'probe'): newarray.setProbe(grid.probe)
        else: newarray.setProbe('UNK_UNK')
        trimmedList.append(newarray)
        
    return trimmedList

def similarity(gridList, percentile=0.02, hardcutoff=False, comparepositive=False, ignoreValue=False):
    """
    Compute similarity index between grids in gridList. First grids are discretized to 1 and 0. 1 for points
    below certain cutoff (or above if comparepositive is True). The cutoff can be dynamic with a percentile
    calculation or stablished by the user as a hard cutoff.
    
    Will use a modified Tanimoto Index to evaluate similarity. This version takes 
    into acount neighbouring positions (expanded grid) to allow some flexibility 
    and reduce noise introduced by the grid spacing.
    
    SimIndex = (Nab + Na'b + Nab' + Na'b') / (Na + Nb + Na'b' + (Na'b/2) + (Nab'/2) - Nab)
    
    a' = expanded a without original points. Thus only expansion.
    b' = expanded b without 1s in b
    Na = Ones in original A
    Na'= Ones in expansion of A
    Nab = Ones common in A and B        
    Na'b = Ones in common in expansion of A and original B
    etc...
    
    :arg list gridList: List of grid instances or path to files to be compared.
    :arg float percentile: Percentile of points to convert to 1.
    :arg float hardcutoff: Ignore percentile calculation and use this hard cutoff to assign 1 and 0.
    :arg bool comparepositives: Compare negative or positive points? By default use negative tail.
    :arg float ignoreValue: Value to ignore during cutoff and comparison calculation. E.g. 999 masking excluded volume.
    
    :return: Redundant similarity matrix.
    """
    import itertools
    
    if not isinstance(gridList, list): raise AttributeError, "Expected gridList of type list. Got %s instead."%(type(gridList))
    
    # Parse each list element to check if it is a Grid instance or a file to be loaded
    grids = []
    for el in gridList:
        if isinstance(el, Grid): grids.append(el)
        elif isinstance(el, str):
            # Try to load as a file
            if os.path.exists(el): grids.append(Grid(el))
            else: raise AttributeError, "Element in list is not a valid file path or a Grid instance: %s"%el
        else:
            raise AttributeError, "Wrong list element type: %s %s"%(el, type(el))
    
    # Trim Grids if necessary
    grids = trim(grids)
    ngrids = len(grids)
    
    # Mask out ignoreValues if needed
    if ignoreValue:
        ignoreMasks = [g.data != ignoreValue for g in grids]
        
    # Invert data if comparepositive
    if comparepositive:
        if hardcutoff: hardcutoff *= -1
        for g in grids: g.data *= -1
    
    # Compute cutoffs if needed
    # Convert to zeros and ones
    if not hardcutoff:         
        cutoffs = [g.getPercentileCutValue(percentile, ignoreValue) for g in grids]
    else:
        cutoffs = [hardcutoff]*len(grids)
        
    for i,g in enumerate(grids):
        m = g.data <= cutoffs[i]
        g.data[m] = 1
        g.data[~m] = 0
        if ignoreValue: g.data[ignoreMasks[i]] = 0

    # Save an expanded version of each grid
    expanded_grids = []
    for g in grids:
        tmpg= g.copy()
        tmpg.expand(2)
        ones_idx = npy.vstack(npy.where(tmpg.data)).T
        for i in ones_idx:
            tmpg.cancelPoints(i, 0.5, 1.0)
        tmpg.contract(2)
        # Substract original ones to take only expansion
        # and save
        tmpg.data -= g.data
        expanded_grids.append(tmpg)
    
    data = npy.zeros((ngrids,ngrids))
    # Compare 
    for combi in itertools.combinations(npy.arange(ngrids), 2):
        i, j = combi
        gi, gj = grids[i].data, grids[j].data
        gi_, gj_ = expanded_grids[i].data, expanded_grids[j].data
        
        A = gi.sum()
        B = gj.sum()
        A_ = gi_.sum()
        B_ = gj_.sum()
        AB = (gi*gj).sum() # Common ones
        A_B = (gi_*gj).sum()
        AB_ = (gi*gj_).sum()
        A_B_ = (gi_*gj_).sum()
        
        data[i,j] = (AB + A_B + AB_ + A_B_)/float(A + B + A_B_ + (A_B/2.) + (AB_/2.) - AB)
    
    # Recompose full matrix
    # and return
    data += data.T
    data += npy.eye(ngrids) # Reconstruct diagonal with ones
    return data

if __name__ == "__main__":
    print "Testing GridSpace and Grid"
    
