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

import os
import sys
import time
import re
import copy

import numpy as npy

class GridDataError(Exception):
    pass

class BadFile(GridDataError):
    pass

class BadAttribute(GridDataError):
    pass

class GridData(object):
    """Class to manage grid files"""
    def __init__(self, fname=False, origin=None, shape=None, spacing=None, geomExtent=None,
                value_filling=0, margin=0, *args, **kwargs):
        """
        Constructor method.

        :arg str fname: path to existing grid file to be loaded.
        :arg origin: when creating an empty grid, specify origin coordinate
        :type origin: list of floats or numpy.array
        :arg shape: shape of new grid to be created in number of cell points.
        :type shape: list of ints or numpy.array
        :arg spacing: Grid spacing in angstroms. E.g. [0.5,0.5,0.5]
        :type spacing: list of floats or numpy.array.
        :arg geomExtent: distances to span in x,y,z directions. Considering the spacing, shape will be automatically determiend.
        :type geomExtent: numpy.array or list of floats
        """    
        self.data = False
        self.origin = origin
        self.delta = spacing
        self.atype = ''
        self.source = ''
        self.header = None

        # Load from file constructor
        if fname:
            if os.path.exists(fname):
              if '.xplor' in fname or '.cns' in fname: self.readXPLOR(fname)
              elif '.dx' in fname: self.readDX(fname)
              else: 
                raise BadFile, "Can not identify %s file format.It should contain 'xplor', 'cns' or 'dx' extension"%fname
            else:
                raise BadFile, "File %s not found."%fname   

        # Construct empty grid
        elif origin is not None:
            self.source = "User created Grid"
            # Check and convert to correct types
            if len(origin)==3:
                self.origin = npy.array(origin) - margin
            else:
                raise BadAttribute,"Invalid origin. Should be a tuple or list of 3 floats or numpy array"

            if spacing != None:
                if type(spacing) is int or type(spacing) is float:
                    self.delta = npy.array([spacing, spacing, spacing])
                elif type(spacing) is list or isinstance(spacing, npy.ndarray):
                    self.delta = npy.array(spacing)

                # Build grid from shape information
                if shape is not None:
                    self.shape = npy.array(shape).astype(int)

                elif geomExtent is not None:
                    if isinstance(geomExtent,float) or isinstance(geomExtent, int):
                        # Equidistant in all directions
                        geomExtent = npy.array([geomExtent, geomExtent, geomExtent]).astype(float)
                    elif isinstance(geomExtent, list) or isinstance(geomExtent, tuple) or isisntance(geomExtent, npy.ndarray):
                        geomExtent = npy.array(geomExtent).astype(float)
                    self.shape=(geomExtent+margin*2)/self.delta
                    self.shape= map(int, self.shape)
                else:
                    raise BadAttribute, "Constructing empty Grid: When spacing is given, origin+shape or origin+geomExtent are mandatory."

            else:
                # No spacing, to create grid, origin+shape+geometricExtent are needed, spacing will be adapted
                if shape is not None and geomExtent is not None:
                    self.shape = npy.array(shape).astype(int)
                    if isinstance(geomExtent,float) or isinstance(geomExtent, int):
                        # Equidistant in all directions
                        geomExtent = npy.array([geomExtent, geomExtent, geomExtent]).astype(float)
                    elif isinstance(geomExtent, list) or isinstance(geomExtent, tuple) or isisntance(geomExtent, npy.ndarray):
                        geomExtent = npy.array(geomExtent).astype(float)
                    self.delta = (geomExtent+margin*2)/self.shape
                else:
                    raise BadAttribute, "Constructing empty Grid: When no spacing is given, origin+shape+geomExtent are mandatory. Grid not created."

            # Build data array. dtype can be assigned through kwargs
            self.data = npy.zeros(self.shape, **kwargs) + value_filling

        else:
            raise BadAttribute, "To create a grid, a filename of an existing grid or an origin+shape or origin+extent must be given."

    def __getitem__(self, point):
        """Iterate Data giving cartesian coordinates"""
        index = self.getIndex(point)
        if index: return self.data[index]
        else: return False
    
    def __setitem__(self, key, value):
        """When used as self[a,b,c]=x for a,b,c as cartesian coordinates
        Thus we can assign grid values giving cartesian coordinates instead of indexes."""
        index = self.getIndex(key)
        self.data[index] = value
  
    def __repr__(self):
        if self.source =='': return "Empty Grid instance"
        return "Grid instance for %s"%(self.source)
  
  
    def __str__(self):
        if self.source =='': return "Empty Grid instance"
        return "Grid instance for %s"%(self.source)
    
    def __call__(self):
        """Returns information about the grid if already defined"""
        if (self.data !='' and self.origin !='' and self.delta !=''):
            print "Grid instance source:",self.source
            print "Data shape: ",self.data.shape," number of entries: ",self.data.size
            print "Grid Origin: ",self.origin
            print "Grid Spacing: ",self.delta
        else: print "Empty grid instance or data missing!"
  
    def readXPLOR(self, CNS):
        """
        Read XPLOR/CNS Formatted grids
        """
        #Open XPLOR/CNS Format Grid File
        # Unzip if zipped
        if 'gz' in CNS:
            import gzip
            op = lambda x: gzip.open(x, 'rb')
        else:
            op = lambda x: open(x, 'r')

        f = op(CNS)
        #Skip 1st line
        f.readline()
        #Header lines:
        head = f.readline().split()
        head = int(head[0])
        self.header = ''
        while head:
            self.header += f.readline()
            head-=1

        #Read grid format
        nx, minx, maxx, ny, miny, maxy, nz, minz, maxz = map(int, f.readline().split())
        dx, dy, dz, deltax, deltay, deltaz = map(float, f.readline().split())

        #Skip ZYX line
        f.readline()

        #Calculate size to read
        #XPLOR Format has ZYX format, X fast, Y medium, Z slow.
        #X data is sorted in 6 columns
        #Each data is 12 bits and each line is 1 bit more (newline charater)
        #This format is repeated ny times. Then a new Z block starts (one line to be skipped then).

        # Some programs write a new line every time an nx block finishes (ptraj)
        # and others do not write this nx block new line considering nx*ny as a single block (moe)
        # we need a reader that considers these two situations

        def readconsecutivenxny(f):
            nxy = ny*nx # bits per number
            # Calulate bits due to linebreaks
            lbreaks = nxy/6  # 6 columns
            if nxy%6: lbreaks+=1   # remainder
            size = nxy*12+lbreaks
            chunkstring = lambda b: [b[i:i+12] for i in range(0,len(b),12)]
            return npy.array(chunkstring(f.read(size).replace('\n','')), dtype=float)

        def readsplittednxny(f):
            nxy = ny*nx # bits per number
            # Calulate bits due to linebreaks
            lbreaks = (nx/6)*ny  # 6 columns
            if nx%6: lbreaks+=nx   # remainder
            size = nxy*12+lbreaks
            chunkstring = lambda b: [b[i:i+12] for i in range(0,len(b),12)]
            return npy.array(chunkstring(f.read(size).replace('\n','')), dtype=float)

        grid = npy.zeros([nx,ny,nz])
        for z in xrange(nz):
            f.readline() #skip block identifier line
            pos = f.tell()
            try:
                grid[:,:,z]=readconsecutivenxny(f).reshape(ny,nx).T
            except:
                f.seek(pos)
                grid[:,:,z]=readsplittednxny(f).reshape(ny,nx).T

        #Prepare standard descriptors:
        delta = [(dx/nx), (dy/ny), (dz/nz)]
        origin = [minx*delta[0],miny*delta[1],minz*delta[2]]

        #Assign attributes
        self.data = grid
        self.shape = npy.array([nx,ny,nz])
        self.origin = origin
        self.delta = delta
        self.source = CNS

    def readBinXPLOR(self, binxplor):
        """
        Read XPLOR grid in binary format.
        
        :arg str binxplor: Path to file.
        """
        import struct
        of = open(binxplor, 'rb')

        # Read header
        struct.unpack('s',of.read(1))
        hlen, =struct.unpack('i',of.read(4))
        header=struct.unpack('s'*hlen, of.read(hlen))
        print ''.join(header)

        # Write grid specifications
        maxd = [0,0,0]
        nx, minx, maxd[0], ny, miny, maxd[1],nz, minz, maxd[2] = struct.unpack('i'*9, of.read(4*9))
        dx, dy, dz,deltax, deltay, deltaz = struct.unpack('f'*6, of.read(4*6))

        grid = npy.zeros((nx,ny,nz), dtype=float)
        # Write data without floritures
        for z in xrange(nz):
            grid[:,:,z] = npy.array(struct.unpack('f'*nx*ny, of.read(4*nx*ny)), dtype=float).reshape(nx, ny)

        sentinel, = struct.unpack('i',of.read(4))    #end of file testimony
        of.close()
        if sentinel != -9999: print 'read error!'

        self.data = grid
        self.delta = npy.array([(dx/nx), (dy/ny), (dz/nz)])
        self.origin = npy.array([minx*self.delta[0],miny*self.delta[1],minz*self.delta[2]])
        self.shape = npy.array([nx,ny,nz])
        self.source = binxplor

    def readDX(self, DX):
        """
        Read grid data from a DX formated file.
        
        :parm str DX: Path to file to read.
        """
        # Unzip if zipped
        if 'gz' in DX:
            import gzip
            op = lambda x: gzip.open(x, 'rb')
        else:
            op = lambda x: open(x, 'r')

        f=op(DX)
        #read the header if any
        header=""
        l = f.readline()
        while l.startswith('#'):
            header= header + l
            l = f.readline()

        #read the grid size
        r=re.compile('\w+')
        gsize=r.findall(l) # Read first non hashed line
    #        print gsize
        gsize=[int(gsize[-3]),int(gsize[-2]),int(gsize[-1])]

        #read the origin of the system
        line=f.readline().split()
        origin=[float(line[-3]),float(line[-2]),float(line[-1])]

        #read grid space
        line=f.readline().split()
        deltax=[float(line[-3]),float(line[-2]),float(line[-1])]
        line=f.readline().split()
        deltay=[float(line[-3]),float(line[-2]),float(line[-1])]
        line=f.readline().split()
        deltaz=[float(line[-3]),float(line[-2]),float(line[-1])]

        #pay attention here, this assumes always orthogonal normalized space, but normally it should be ok
        delta=[deltax[0],deltay[1],deltaz[2]]

        #read the number of data
        f.readline()
        r=re.compile('\d+')
        n_entries=int(r.findall(f.readline())[2])

        #check correpondence with expected data points
        if(n_entries!=gsize[0]*gsize[1]*gsize[2]) : sys.exit("Error reading the file. The number of expected data points does not correspond to the number of labeled data points in the header.")

        #load data into numpy array
        grid = npy.fromstring(f.read(), sep=' ', dtype=float).reshape(gsize) #reshaping to fit grid format (it keeps Z fast, Y medium, X slow data organization)

        if grid.size != n_entries: sys.exit("Error reading the file. The number of expected data points does not correspond to the number of labeled data points in the header.")
        f.close()

        self.header = header
        self.data= grid
        self.shape = npy.array(self.data.shape)
        self.origin = origin
        self.delta = delta
        self.source = DX

    def loadData(self, pick):
        """
        DEPREACTED.Loads grid data from existing pickle. DEPREACTED.
        
        :parm str pick: Path to pickle file.
        """
        if os.path.exists(pick):
            grid = npy.load(pick)
            self.data = grid
            self.source = pick
        else: raise IOError, "File not found."

    def loadVar(self, pick):
        """DEPREACTED. Loads origin and delta information from a pickled file"""
        if os.path.exists(pick):
            V = npy.load(pick)
            origin, delta= V[:]
            self.origin = origin
            self.delta = delta
        else: raise IOError,"File not found."

    def writeDX(self, dxname, header=False, gzip=False):
        """
        Writes data into a DX Formated File. 
        
        :parm str dxname: File name to write.
        :parm str header: Custom string (one line without breaks) to be introduced at the top of the file for identification. If false a predefined string will be writen.
        :parm bool gzip: Write file in gzipped form.
        """
        if (self.data !='' and self.origin !='' and self.delta !=''):
            grid = self.data
            origin = self.origin
            delta = self.delta

            if gzip:
                import gzip
                op = lambda x: gzip.open(x, 'wb')
            else:
                op = lambda x: open(x, 'w')

            if not 'gz' in dxname and gzip:
                outname = dxname+'.gz'
            else:
                outname = dxname

            dxf = op(outname)
            if not header:
                header = ["#Grid generated with MDMix Grids.py module","#Time: %s"%time.ctime(),"#Source: %s"%self.source]
            elif header:
                header = '#'+header
                header = [header,"#Time: %s"%time.ctime(),"#Source: %s"%self.source]
            dxf.write('\n'.join(header))
            dxf.write("""\nobject 1 class gridpositions counts %i %i %i\norigin %.2f %.2f %.2f\ndelta %.5f   0   0\ndelta 0   %.5f   0\ndelta 0   0   %.5f\nobject 2 class gridconnections counts %i %i %i\nobject 3 class array type double rank 0 items %i data follows\n"""%(grid.shape[0],grid.shape[1],grid.shape[2],origin[0],origin[1],origin[2],delta[0],delta[1],delta[2],grid.shape[0],grid.shape[1],grid.shape[2],grid.size))

            count=0
            for data in grid.flat:
                if count==3: dxf.write('\n'); count=0
                dxf.write('%6.3f\t'%(float(data)))
                count+=1
            dxf.write('\n')
            dxf.close()

        else: raise GridError, "Data or parameters missing. Can't write DX File."
        
    def dump(self, dataname=None, varname=None):
        """
        DEPRECATED. This method uses by default the source name of self to generate the dumped file's names. For a specific name, it must be provided when calling the method.
        """
        if (self.data !='' and self.origin !='' and self.delta !=''):

          if not self.source:
              raise BadAttribute, "Must provide a name for the dumped files. Source attribute not defined."      

          if dataname:
            dname = dataname      #Grid data pickle name
          else:
            dname = str(os.path.splitext(self.source)[0])+"_grid.pick"   

          if varname:
            vdname = varname    #Grid variables pickle name
          else:
            vdname = str(os.path.splitext(self.source)[0])+"_grid_var.pick"

          npy.save(dname, self.data)
          V = [self.origin, self.delta]
          npy.save(vdname, V)
        else:
          raise BadAttribute, "Data or parameters missing. Check attributes."

    def writeXPLOR(self, xplorname, header="", gzip=False, binary=False):
        """
        Write data into XPLOR/CNS formated file.
        
        :parm str xplorname: File to write.
        :parm str header: Custom header line to add at the beggining of the file. If empty, a custom line will be written.
        :parm bool gzip: Write in gzipped form.
        :parm bool binary: Write in binary format to speed up reading.
        """
        grid = self.data
        origin = self.origin
        delta = self.delta
        if isinstance(delta, npy.ndarray):
            delta = delta.tolist()
        elif type(delta) is int or type(delta) is float:
            delta = map(float, [delta, delta, delta])
        nx, ny, nz = grid.shape

        #minx, miny, minz, = origin[0]/delta[0], origin[1]/delta[1], origin[2]/delta[2]
        #maxx, maxy, maxz = minx+nx-1, miny+ny-1, minz+nz-1
        mind = npy.round(npy.array(origin) / delta)
        maxd = npy.round(mind + grid.shape - 1)
        space = npy.array(delta)*grid.shape   #Calculate total system expansion (distance) in x, y, z.
        deltax, deltay, deltaz = 90., 90., 90.      #Assuming rectangular shape

        if not header: header = "Data Grid generated with MDMix.Grids"
        header = ["", "       1",header]

        if binary:
            import struct
            of = open(xplorname, 'wb')

            # Write header
            of.write(struct.pack('s',header[0]))
            hlen=len(header[2])
            of.write(struct.pack('i',hlen))
            of.write(struct.pack('s'*hlen, *tuple(header[2])))

            # Write grid specifications
            of.write(struct.pack('i'*9,*(nx, mind[0], maxd[0],
                                        ny, mind[1], maxd[1],
                                        nz, mind[2], maxd[2])))
            of.write(struct.pack('f'*6, space[0], space[1], space[2],deltax, deltay, deltaz))

            # Write data without floritures
            for z in xrange(nz):
                d = tuple(grid[:,:,z].flatten())
                dlen = len(d)
                of.write(struct.pack('f'*dlen, *d))

            of.write(struct.pack('i',-9999))    #end of file testimony
            of.close()

        else:

            if gzip:
                import gzip
                op = lambda x: gzip.open(x, 'wb')
            else:
                op = lambda x: open(x, 'w')

            # Give extension if not given by user
            if not 'gz' in xplorname and gzip:
                outname = xplorname+'.gz'
            else:
                outname = xplorname

            of = op(outname)

            #Header and grid format
            of.write("\n".join(header))
            of.write("""\n%8d%8d%8d%8d%8d%8d%8d%8d%8d\n%12.5E%12.5E%12.5E%12.5E%12.5E%12.5E\nZYX\n"""%(nx, 
                mind[0], maxd[0], ny, mind[1], maxd[1], nz, mind[2], maxd[2], space[0], space[1], space[2], 
                deltax, deltay, deltaz))

            chunkiter = lambda a: (a[i:i+6] for i in xrange(0, len(a), 6))
            fmt='%12.5E'
            #Write Grid data
            for z in xrange(nz):
                of.write("%8d\n"%(mind[2]))
                for ch in chunkiter(grid[:,:,z].flatten('F')):
                    s = fmt*len(ch)+'\n'
                    of.write(s%(tuple(ch)))
                mind[2]+=1
            of.write('%8i\n'%-9999)
            of.write('%12.5E %12.5E\n'%(grid.mean(),grid.std()))
            of.close()

    def averageData(self, contract=False):
        """ 
        Will average each grid point by the surrounding ones.
        
        :parm bool contract: Points in the edges will be removed, thus size of the grid will change
        by one point. Also origin coordinate will be modified accordingly.
        """
        calcgrid = npy.zeros((3,3,3))
        calcgrid += 1/27. 

        narr = copy.deepcopy(self.data)
        for x in xrange(1,narr.shape[0]-1):
            for y in xrange(1,narr.shape[1]-1):
                for z in xrange(1,narr.shape[2]-1):
                    sl = self.data[x-1:x+2,y-1:y+2,z-1:z+2]
                    narr[x,y,z] = (sl * calcgrid).sum()

        # Update data
        self.update(narr)

        if contract:
            # Remove borders. Origin is modified in contract method already.
            maxdelta = npy.array(self.delta).max()
            self.contract(maxdelta)

    def count2DG(self, expected, T=300., correction=1, maskvalue=999.):
        """
        Use invert Boltzmann equation for converting observed distribution to energies. Results in kcal/mol units.
        DGBind calulation:
        
        .. math::
            DG_{bind} = -R*T*ln(N_{i}/N_{0})
            
        where :math:`R` is gas constant, :math:`T` temperature, :math:`N_{i}` observed distribution and :math:`N_{0}` the
        expected distribution.
        
        :parm float expected: Expected number (N_0).
        :parm float T: Temperatue to use in Boltzmann expression.
        :parm float correction: Correction factor to apply to the expected number N0.
        :parm float maskvalue: Mask with this value point with zero density.
        
        :return: numpy.ndarray with energy converted values in kcal/mol.
        """
        Cgrid = self.data   # Grid count data
        Kb = 1.987/1000.    # Boltzmann Constant in kcal/mol units   
        N0 = expected * correction
        # Convert zeros to ones to avoid log(0)
        # later will restore values after applying formula
        maskZ = Cgrid == 0   
        Cgrid[maskZ] = 1
        DGgrid = - Kb * T * ( npy.log(Cgrid) - npy.log(N0) )
        DGgrid[maskZ] = maskvalue
        return DGgrid

    def update(self, data):
        """
        Updates data information.
        
        :parm numpy.ndarray data: Grid data to be updated.
        """
        self.data = data

    def getCartesian(self, indexes):
        """
        Returns cartesian coordinates of the index value.
        
        :parm list indexes: List, numpy array or tuple with size 3 identifying a point in the grid.
        
        :return: numpy array with length 3 with the corresponding artesian coordinates.
        """
        if len(indexes)==3:
          indexes = npy.array(indexes, dtype=float)
          spacing = self.delta
          origin = self.origin
          return (indexes * spacing) + origin
        else: raise BadAttribute, "indexes should be tuple or list of len 3."

    def getIndex(self, cartesian):
        """
        Returns grid indices for a given cartesian point (x,y,z)
        
        :parm list cartesian: List, numpy array or tuple with length 3 cartesian coordinates to be transformed to indices.
        
        :return: tuple with the corresponding indices (length 3).
        """
        spacing = self.delta
        origin = self.origin
        if type(cartesian) is not npy.array: cartesian = npy.array(cartesian)
        ind = npy.floor((cartesian - origin) / spacing)
        # Check that index obtained is inside the grid.
        if ( npy.any( ind < 0. ) or npy.any(ind >= self.data.shape ) ):    
            return False
        else:
            return tuple(ind.astype(int).tolist())

    def cancelPoints(self, point, cutoff, value=999 ):
        """
        DEPRECATED: Use setRadialValues directly.
        """
        self.setRadialValues(rmax=cutoff,rmin=0,value=value,point=point)

    def maskOut(self, value, over=False, include=True):
        "Return a masked array with values 'over' or under (if 'over' is false) 'value' masked.\
        If 'include' is True, 'value' will also be masked."
        if over:
            if include: func = npy.ma.masked_greater_equal
            else:       func = npy.ma.masked_greater
        else:
            if include: func = npy.ma.masked_less_equal
            else:       func = npy.ma.masked_less

        return func(self.data, value)

    def mergeDeleteProt(self, PDB, cutoff, value=999):
        """
        Method to "Delete" grid information on the points overlaping atoms of a PDB. It is actually a method to mask the overlapping points with a certain value.

        @param1 PDB:     protein
        @param1 type:    Biskit PDBModel instance or PDB File name

        @param2 cutoff: Cutoff allows to expand the assignement to N angstroms around that point.
        @param2 type;   int or float

        @param3 value:  Optionally one can provide the value to assign to those points (default = 999)
        @param3 type:   int
        """
        import Biskit
        if isinstance(PDB, Biskit.PDBModel):
            protein = PDB
        else:
            protein = Biskit.PDBModel(PDB)
        origin = self.origin
        delta = self.delta
        for atom in protein.xyz:
          idx = self.getIndex(atom)
          #a,b,c = atom
          #i, j, k = int(round((a-origin[0])/delta[0])), int(round((b-origin[1])/delta[1])), int(round((c-origin[2])/delta[2]))
          if idx != False: self.cancelPoints(idx, cutoff, value)

    def mergeConserveProt(self, PDB, cutoff, value=999):
        """
        Method to Conserve only grid information near the points overlaping atom's coordinates of a PDB.

        @param1 PDB:     protein
        @param1 type:    Biskit PDBModel instance or PDB File name

        @param2 cutoff: How far from the atom coordinate should keep grid information. In Angstroms.
        @param2 type;   int or float

        @param3 value:  Value to assign to the rejected points (default = 999)
        @param3 type:   int

        return data
        """
        import Biskit
        if isinstance(PDB, Biskit.PDBModel):
            protein = PDB
        else:
            protein = Biskit.PDBModel(PDB)

        origin = self.origin
        delta = self.delta
        _dataBackup = copy.deepcopy(self.data)
        self.data = npy.zeros(self.data.shape) == 0

        for atom in protein.xyz:
          a,b,c = atom
          i, j, k = int(round((a-origin[0])/delta[0])), int(round((b-origin[1])/delta[1])), int(round((c-origin[2])/delta[2]))
          self.cancelPoints([i,j,k], cutoff, False)

        _dataBackup[self.data] = value
        self.data = _dataBackup
        return _dataBackup

    def copy(self):
        """Returns a new instance as a deepcopy of self"""
        new_inst = copy.deepcopy(self)
        return new_inst

    def min_indx(self):
        """Returns index for the minimum value in data"""
        return npy.array( zip( *npy.where( self.data == self.data.min() ) ) )[0]

    def expand(self, buff,fill=0):
        """
        Expand the array size building a Zero filled Array adding buff Angstrom each side.
        Place data inside and modify origin information.
        """
        delta = self.delta
        size = self.data.shape
        buffering = [ int( buff / delta[0] ), int( buff / delta[1] ), int( buff / delta[2] ) ]

        G = npy.zeros( [ size[0] + (buffering[0]*2), size[1] + (buffering[1]*2), size[2] + (buffering[2]*2) ] )+fill
        G[ buffering[0]:size[0] + buffering[0], buffering[1]:size[1] + buffering[1], buffering[2]:size[2] + buffering[2] ] = self.data.copy()

        self.data = G.copy()
        self.origin = npy.array(self.origin) - buff

    def contract(self, buff):
        """ 
        Inverser operation to expand. 
        buff in angstroms.
        """
        delta = self.delta
        size = self.data.shape
        buffering = [ int( buff / delta[0] ), int( buff / delta[1] ), int( buff / delta[2] ) ]
        self.data = self.data[ buffering[0]:size[0]-buffering[0], buffering[1]:size[1]-buffering[1], buffering[2]:size[2]-buffering[2] ].copy()
        self.origin = npy.array(self.origin) + buff

    def trim(self, X):
        """
        Trim self array data to fit common data with grid X. 
        New shape and coordinates is common in the two arrays but data is copied from self
        """
        if npy.all(self.delta == X.delta):

            A = self

            Amaxind = npy.array( A.data.shape )     
            Xmaxind = npy.array( X.data.shape )

            #Get new Origin and new Shape
            newOrigin = npy.vstack((A.origin, X.origin)).max(axis=0)
            newMaximum = npy.vstack((A.getCartesian(Amaxind), X.getCartesian(Xmaxind))).min(axis=0)
            newShape = ( newMaximum - newOrigin ) / A.delta
            newShape.astype('int')

            newarray = GridData(shape=newShape, origin=newOrigin, spacing=A.delta)
            newarray.source = 'Trimmed array'

            #Copy self data to the new array
            init = A.getIndex(newarray.origin)
            newarray.update( A.data[ init[0]:init[0]+newarray.data.shape[0], init[1]:init[1]+newarray.data.shape[1], init[2]:init[2]+newarray.data.shape[2] ].copy() )

            return newarray

        else:
            raise BadAttribute, "Spacing is not equal. Can not trim."

    def takeSubGridBox(self, bot, top,forceCubic=0):
        """
        Return a grid instance containing a part of the main grid.

        arg1:   bot :   bottom Cartesian coordinates of the subgrid
        arg1    type:   tuple or list [x,y,z]

        arg2:   top:    top Cartesian coordinates of the subgrid
        arg2    type:   tuple or list [x,y,z]

        arg3:   forceCubic :    a boolean to know if we force  cubic grid to be returned
        arg3:   type:   int (0 or 1)
        """

        subOrigin = bot
        subOriginIndx = npy.array(self.getIndex(bot))
        subTopIndx = npy.array(self.getIndex(top))

        if forceCubic:
            diffIndx=subTopIndx-subOriginIndx
            maxDiff=npy.max(diffIndx)
            subTopIndx+=maxDiff-diffIndx

        tmp=self.data[ subOriginIndx[0]:subTopIndx[0]+1, subOriginIndx[1]:subTopIndx[1]+1, subOriginIndx[2]:subTopIndx[2]+1 ]
        subGrid = GridData(origin=subOrigin,spacing=self.delta,shape=tmp.shape)
        subGrid.data = npy.ones_like(tmp)
        subGrid.data[:] = tmp[:]
        subGrid.source = self.source + 'SUBGrid'

        return subGrid

    def takeSubGridPoint(self, point, distance):
        """
        Return a grid instance containing a part of the main grid.
        arg1:   point:    coordinates for the center of the box
        arg1    type:      tuple or list [x,y,z]

        arg2:   distance:   Distance of the box side size in Angstroms
        arg2    type:       tuple or list [dx, dy, dz] for rectangular box
        arg2    type:       int or float for cubic box
        """

        subOrigin = npy.array(point) - (distance/2.)
        subOriginIndx = self.getIndex(tuple(subOrigin))
        subTop = npy.array(point) + (distance/2.)
        subTopIndx = self.getIndex(tuple(subTop))
        data = self.data[ subOriginIndx[0]:subTopIndx[0]+1, subOriginIndx[1]:subTopIndx[1]+1, subOriginIndx[2]:subTopIndx[2]+1 ]
        subGrid = GridData(origin=subOrigin, spacing=self.delta, shape=data.shape)
        subGrid.data = data
        subGrid.source = self.source + 'SUBGrid'

        return subGrid

    def getRadialIndices(self, radius, point=None, minradius=False):
        "Return an array of indices belonging to a radius in angstroms considering the grid spacing.\
        If 'point' is given as indices for the center, indices outside the box will be removed.Also indices returned will be centered at that point.\
        If 'minradius' in angstroms is given, points not reaching that minimum will not be included, returning the radial indices"
#        spacing = npy.array(self.delta)
        spacing = npy.array(self.delta)
        span = npy.round(float(radius) / spacing).astype(int)
        indices = []
        for i in range(-span[0], span[0]+1):
            for j in range(-span[1], span[1]+1):
                for k in range(-span[2], span[2]+1):
                    r = npy.linalg.norm([i,j,k]*spacing)
                    if not minradius and r<=radius:indices.append((i,j,k))
                    if minradius:
                        if r<=radius and r>=minradius: indices.append((i,j,k))

        if not indices: return indices

        if point != None:
            indices = npy.array(indices) + point
            # Remove indices outside superior bound
            out = npy.where(indices - self.data.shape >= 0)[0]
            indices = npy.delete(indices, out, axis=0)

            # If there are still some indices left, try to remove indices outside inferior bound
            if npy.any(indices):
                out = npy.where(indices < 0)[0]
                indices = npy.delete(indices, out, axis=0)
        return indices

    def getSphereValues(self, center, r):
        """ This method returns an array containing the values around a certain
        cartesian point given a radius """

        idx = self.getIndex(center)
        indices = self.getRadialIndices(r, point=idx)
        return npy.array([self.data[tuple(i)] for i in indices])
        #return self.data[idx[0]-span:idx[0]+span+1, idx[1]-span:idx[1]+span+1, idx[2]-span:idx[2]+span+1]

    def setSphereValues(self, center, r, value):
        "Set all values around 'r' angstroms of coordinate 'center' to 'value'"
        idx = self.getIndex(center)
        indices = self.getRadialIndices(r, point=idx)
        for i in indices:
            self.data[tuple(i)] = value
        return True

    def getRadialValues(self, center, rmax, rmin):
        """ This method returns an array containing the values around a certain
        cartesian point from rmin to rmax radius """
        idx = self.getIndex(center)
        indices = self.getRadialIndices(rmax, point=idx, minradius=rmin)
        return npy.array([self.data[tuple(i)] for i in indices])

    def setRadialValues(self, rmax, rmin, value, center=None, point=None):
        "Set all values around 'center' or 'point' from rmin to rmax angstroms to 'value'"
        if center is None and point is None:
            raise AttributeError, "Calling setRadialValues requires a cooridnate ('center' arg) or a grid index ('point' arg)."

        if center != None: idx = self.getIndex(center)
        elif point != None: idx = point

        indices = self.getRadialIndices(rmax, point=idx, minradius=rmin)
        for i in indices:
            self.data[tuple(i)] = value
        return True


class GridFromPDB(GridData):
    def __init__(self, PDB, spacing, buff=0, value_filling=0, takeProtein=True):
	"""
	PDB	PDBModel or PDB file
	spacing grid spacing
	buff	buffer around atom coordinates in angstroms
	value_filling	initial value of the grid
	takeProtein	Boolean. If true, will try to compress only protein from PDB.
			If False, will take the whole PDB itself
	"""
        import Biskit
        if isinstance(PDB, Biskit.PDBModel):
            prot = PDB
        elif os.path.exists(PDB):
            pdb_model = Biskit.PDBModel(PDB)                    #Open pdb as Biskit PDBModel
            if takeProtein: prot = pdb_model.compress(pdb_model.maskProtein())  #Read only protein from PDBModel
	    else: prot = pdb_model
        else: 
            print "Couldn't locate PDB file: ",PDB
            return False        
        
        xyz = prot.xyz                                                  # get protein coordinates
        extension = ( xyz.max(axis=0) - xyz.min(axis=0) ) + buff * 2    # protein extension in x,y,z axis
        shape = (extension / spacing).astype(int)                       # New grid shape
        origin = xyz.min(axis=0) - buff               # grid min cartesian coordinates as origin 
        
        GridData.__init__(self, shape=shape, origin=origin.tolist(), spacing=spacing, value_filling=value_filling)
        self.source = 'Generated from PDB: '+str(PDB)


class GridFromMol2(GridData):
  def __init__(self, Mol2Model, shape, margin=0, value_filling=0):
        
        xyz = Mol2Model.xyz                                             #get molecule coordinates
        extension = ( xyz.max(axis=0)-xyz.min(axis=0) ) + margin*2      #protein extension in x,y,z axis
        spacing = extension/npy.array(shape)                             #New grid size with buff*2 angstroms added in each axis
        origin = xyz.min(axis=0) - margin            #grid min cartesian coordinates as origin build a zeros grid with 0.5 spacing
        #print "Applying ",buff," Angstroms of buffer each side."
        #print "Origin: ",origin
        #print "Spacing: ",spacing
        GridData.__init__(self,shape=shape,origin=origin,spacing=spacing,geomExtension=extension,value_filling=value_filling)
        self.source = 'Generated from Mol2 file: '+str(Mol2Model)

def get_grid(filename, buff=None):
    "DEPRECATED. BACKCOMPATIBILITY. Use Grid constructor directly."
    return GridData(filename,buff=buff)

def createFromPDB(PDB, spacing, buff=0, value_filling=0, takeProtein=True):
    "DEPRECATED. Mantained for back-compatibility. Use GridFromPDB class."
    return GridFromPDB(PDB,spacing,buff,value_filling,takeProtein)

def createFromMol2(Mol2Model, shape, margin=0, value_filling=0):
    "DEPRECATED. Mantained for back-compatibility. Use GridFromPDB class."
    return GridFromMol2(Mol2Model, shape, margin=0, value_filling=0)

def tanimotoIdentityIsoval(Grid1, Grid2, isovalue, up=True, expand = False, verbose = False):
    """
    This functions returns a Tanimoto coefficient to compare two grids ressemblance.
        Tanimoto = Nab / (Na + Nb + Nab)
        Nab = Ones in A and B
        Na = Ones in A
        Nb = Ones in B
    
    @arg1   Grid1:      Grid instance
    @arg2   Grid2:      Grid instance
    @arg3   isovalue:   cut value to transform to binary data
            type:       int or float
    @arg4   up:         if up is True, take grid values over isovalue, else take values under isovalue
            type:       bool
    @arg5   expand:     Expand the 1s this number of angstroms around. This way we avoid grid spacing problem.
            type:       float or int
    @arg6   verbose:    If True, it returns also Na, Nb, Nab, minimum value and cutoff value for each grid as dictionary
            type:       bool
    """
    
    # Grids are first trimmed to obtain the common part in space.    
    A, B = trim(Grid1, Grid2)
    A.source = Grid1.source
    B.source = Grid2.source
    
    # Mask TRUE (1) values under/over isovalue
    
    if up is True:
        AmaskTrue = A.data >= isovalue
        BmaskTrue = B.data >= isovalue
    else:
        AmaskTrue = A.data <= isovalue
        BmaskTrue = B.data <= isovalue   
    
    # EXPAND POSITIVE POINTS 0.5 ANGSTROMS AROUND
    if expand is True:
        ExpandA = A.copy()
        ExpandB = B.copy()
        ExpandA.data = AmaskTrue*1
        ExpandA.expand(2)   #Expand to avoid bad indexing
        ExpandB.data = BmaskTrue*1
        ExpandB.expand(2)
        for grid in (ExpandA, ExpandB):
            point_lst = npy.vstack(npy.nonzero(grid.data==1)).T
            for point in point_lst:
                grid.cancelPoints(point = point, cutoff = expand, value = 1)
        ExpandA.contract(2)
        ExpandB.contract(2)
        AmaskTrue = ExpandA.data == 1
        BmaskTrue = ExpandB.data == 1
    
    # Mix bitwise AND
    ABTrue = AmaskTrue * BmaskTrue
    
    # Tanimoto calculation
    Na, Nb, Nab = len(npy.nonzero(AmaskTrue)[0]), len(npy.nonzero(BmaskTrue)[0]), len(npy.nonzero(ABTrue)[0])
    Tc = float(Nab) / (Na + Nb - Nab)
    
    if verbose is False:
        return Tc
    else:
        return {'TanimotoCoeff':Tc, 'Na':Na, 'Nb':Nb, 'Nab':Nab, 'Amin': A.data.min(),'Bmin':B.data.min()}


def tanimotoIdentityCutoff(Grid1, Grid2, cutoff, up=False, expand = False, verbose = False):
    """
    This functions returns a Tanimoto coefficient to compare two grids ressemblance.
        Tanimoto = Nab / (Na + Nb + Nab)
        Nab = Ones in A and B
        Na = Ones in A
        Nb = Ones in B
    
    @arg1   Grid1:      Grid instance
    @arg2   Grid2:      Grid instance
    @arg3   isovalue:   cut value to transform to binary data
            type:       int or float
    @arg4   up:         if up is True, take grid values over isovalue, else take values under isovalue
            type:       bool
    @arg5   expand:     Expand the 1s 0.5 angstroms around. This way we avoid grid spacing problem.
            type:       bool
    @arg6   verbose:    If True, it returns also Na, Nb, Nab, minimum value and cutoff value for each grid as dictionary
            type:       bool
    """
    
    #Grids are first trimmed to obtain the common part in space.    
    A, B = trim(Grid1, Grid2)
    A.source = Grid1.source
    B.source = Grid2.source
    
    Amin = A.data.min()        
    Bmin = B.data.min()   
    
    #Flatten negative data, sort it and get the value of the CUT position of negative point as cutoff
    #Sort number of points for each grid and choose the lower one
    
    fsortA = npy.array( A.data[A.data<0].flat )
    fsortA.sort()
    fsortB = npy.array( B.data[B.data<0].flat )
    fsortB.sort()    
    
    cutIndex = []
    [ cutIndex.append( int( (sortedlist.shape[0]-1) * (cutoff/100.) )  ) for sortedlist in (fsortA, fsortB) ]
    
    cutIndex.sort()
    finalIndex = cutIndex[0]
    
    AcutValue = fsortA[finalIndex]
    BcutValue = fsortB[finalIndex]    
    
    #Mask TRUE (1) values under/over isovalue
    
    if up is True:
        AmaskTrue = A.data > AcutValue
        BmaskTrue = B.data > BcutValue
    
    else:
        AmaskTrue = A.data < AcutValue
        BmaskTrue = B.data < BcutValue  
    
    #EXPAND POSITIVE POINTS 0.5 ANGSTROMS AROUND
    if expand is True:
        ExpandA = A.copy()
        ExpandB = B.copy()
        ExpandA.data = AmaskTrue*1
        ExpandA.expand(2)   #Expand to avoid bad indexing
        ExpandB.data = BmaskTrue*1
        ExpandB.expand(2)
        for grid in (ExpandA, ExpandB):
            point_lst = npy.vstack(npy.nonzero(grid.data==1)).T
            for point in point_lst:
                grid.cancelPoints(point = point, cutoff = 0.5, value = 1)
        ExpandA.contract(2)
        ExpandB.contract(2)
        AmaskTrue = ExpandA.data == 1
        BmaskTrue = ExpandB.data == 1
    
    Na, Nb = AmaskTrue.sum(), BmaskTrue.sum()
    
    #Mix bitwise AND
    Nab = (AmaskTrue * BmaskTrue) *1
    Nab = Nab.sum()

    Tc = float(Nab) / (Na + Nb - Nab)
    
    if verbose is False:
        return Tc
    else:
        return {'TanimotoCoeff':Tc, 'Na':Na, 'Nb':Nb, 'Nab':Nab, 'Amin': Amin,'Bmin':Bmin}

def similarityIndexIsoval(Grid1, Grid2, isovalue, up = True, expansion = 1., verbose = False):
    """
    This functions returns a modified Tanimoto coefficient to compare two grids ressemblance.
    It takes in account the expansion in the calculation.
        a' = expanded a without original points. Thus only expansion.
        b' = expanded b without 1s in b
        SimIndex = (Nab + Na'b + Nab' + Na'b') / (Na + Nb + Na'b' + (Na'b/2) + (Nab'/2) - Nab)
        Na = Ones in original A
        Na'= Ones in expansion of A
        Nab = Ones common in A and B        
        Na'b = Ones in common in expansion of A and original B
        etc..
    
    @arg1   Grid1:      Grid instance
    @arg2   Grid2:      Grid instance
    @arg3   isovalue:   cut value to transform to binary data
            type:       int or float
    @arg4   up:         if up is True, take grid values over isovalue, else take values under isovalue
            type:       bool
    @arg5   expansion:  Amount of expansion of true points in Angstroms.This way we avoid grid spacing problem.
            type:       int or float (default = 1)
    @arg6   verbose:    If True, it returns a dictionary with index details.
            type:       bool
            
    """
    
    #Grids are first trimmed to obtain the common part in space.    
    A, B = trim(Grid1, Grid2)
    A.source = Grid1.source
    B.source = Grid2.source
    
    Amin = A.data.min()
    Bmin = B.data.min()
    
    #Mask TRUE (1) values under/over isovalue
    
    if up is True:
        Amask = A.data >= isovalue
        Bmask = B.data >= isovalue
    
    else:
        Amask = A.data <= isovalue
        Bmask = B.data <= isovalue
    
    #EXPAND POSITIVE POINTS 0.5 ANGSTROMS AROUND
    A_ = A.copy()
    B_ = B.copy()
    A_.data = Amask*1
    A_.expand(2)   #Expand to avoid bad indexing
    B_.data = Bmask*1
    B_.expand(2)
    for grid in (A_, B_):
        point_lst = npy.vstack(npy.nonzero(grid.data==1)).T
        for point in point_lst:
            grid.cancelPoints(point = point, cutoff = expansion, value = 1)
    A_.contract(2)
    B_.contract(2)
    
    A = Amask*1
    B = Bmask*1
    
    A_Mask = (A_.data - A)==1
    B_Mask = (B_.data - B)==1
    
    #Mix bitwise AND and sum to obtain number of 1s
    A  = float( A.sum() )
    B  = float( B.sum() )
    A_ = float( A_.data.sum() )
    B_ = float( B_.data.sum() )
    AB  = float( ( ( Amask*Bmask)*1  ).sum() )
    A_B = float( ( (A_Mask*Bmask)*1  ).sum() )
    AB_ = float( ( (Amask*B_Mask)*1  ).sum() )
    A_B_= float( ( (A_Mask*B_Mask)*1 ).sum() )
    
    #SimilarityIndex calculation
    SimIdx = ( AB + A_B + AB_ + A_B_) / (A + B + A_B_ + (A_B/2.) + (AB_/2.) - AB)
    
    if verbose is False:
        return SimIdx
   
    else:
        return {'SimilarityIndex':SimIdx, 'A':A, 'B':B, 'AB':AB, 'A_':A_, 'B_':B_, 'A_B':A_B, 'AB_':AB_, 'A_B_':A_B_, 'Amin': Amin, 'Bmin':Bmin}


def getCutValue(A, B, cutoff):
    """
    Complement function to trim grids.
    It returns the value for A and B grids that will be used as 
    cutoff % threshold to return approximately the same amount of points
    """
    #Flatten negative data, sort it and get the value of the CUT position of negative point as cutoff
    #Sort number of points for each grid and choose the lower one
    fsortA = npy.array( A.data[A.data<0].flat )
    fsortA.sort()
    fsortB = npy.array( B.data[B.data<0].flat )
    fsortB.sort()    
    
    cutIndex = []
    [ cutIndex.append( int( (sortedlist.shape[0]-1) * (cutoff/100.) )  ) for sortedlist in (fsortA, fsortB) ]
    
    cutIndex.sort()
    finalIndex = cutIndex[0]
    
    AcutValue = fsortA[finalIndex]
    BcutValue = fsortB[finalIndex]    

    return AcutValue, BcutValue


def similarityIndex2(Grid1, Grid2, cutoff= 5., verbose = False, trimfirst=True):
    """
    Obtain Similarity Comparison using SimIndex.
    Both grids will be trimmed to contain same number of points.
    
    cutoff (float)  -   Threshold of most negative points to be compares.
                        E.g. 5% of most negative points.
                        It will take into account the number of points independently
                        of the threshold value.
        
    trimfirst (bool)    If True, cutting value will be obtained before trimming.
                        If False, cutting value will be obtained after trimming (thus
                        we only consider values in the trimmed zone).
    
    """
    if not trimfirst: AcutValue, BcutValue = getCutValue(Grid1, Grid2, cutoff)
        
    A, B = trim(Grid1, Grid2)
    A.source = Grid1.source
    B.source = Grid2.source
    
    Amin = A.data.min()        
    Bmin = B.data.min()   
    
    if trimfirst: AcutValue, BcutValue = getCutValue(A, B, cutoff)    
    
    #Mask TRUE (1) values under cutoff
    Amask = A.data < AcutValue     
    Bmask = B.data < BcutValue
    
    #EXPAND POSITIVE POINTS 0.5 ANGSTROMS AROUND
    A_ = A.copy()
    B_ = B.copy()
    A_.data = Amask*1
    A_.expand(2)   #Expand to avoid bad indexing
    B_.data = Bmask*1
    B_.expand(2)
    for grid in (A_, B_):
        point_lst = npy.vstack(npy.nonzero(grid.data==1)).T
        for point in point_lst:
            grid.cancelPoints(point = point, cutoff = 0.5, value = 1)
    A_.contract(2)
    B_.contract(2)
    
    A = Amask*1
    B = Bmask*1
    
    A_Mask = (A_.data - A)==1
    B_Mask = (B_.data - B)==1
    
    
    #Mix bitwise AND and sum to obtain number of 1s
    A  = float( A.sum() )
    B  = float( B.sum() )
    A_ = float( A_.data.sum() )
    B_ = float( B_.data.sum() )
    
    AB  = float( ( ( Amask*Bmask)*1  ).sum() )
    A_B = float( ( (A_Mask*Bmask)*1  ).sum() )
    AB_ = float( ( (Amask*B_Mask)*1  ).sum() )
    A_B_= float( ( (A_Mask*B_Mask)*1 ).sum() )
    
    #SimilarityIndex calculation
    SimIdx = ( AB + A_B + AB_ + A_B_) / (A + B + A_B_ + (A_B/2.) + (AB_/2.) - AB)
    
    if verbose is False:
        return SimIdx
    
    else:
        return {'SimilarityIndex':SimIdx, 'A':A, 'B':B, 'AB':AB, 'A_':A_, 'B_':B_, 'A_B':A_B, 'AB_':AB_, 'A_B_':A_B_, 'Amin': Amin, 'Acut':AcutValue, 'Bmin':Bmin, 'Bcut':BcutValue}


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
        newarray = GridData(shape=newShape, origin=newOrigin, spacing=delta)
        newarray.source = 'Trimmed array of %s'%(grid.source)
        #Copy data to the new array
        init = grid.getIndex(newarray.origin)
        newarray.update( grid.data[ init[0]:init[0]+newarray.data.shape[0], init[1]:init[1]+newarray.data.shape[1], init[2]:init[2]+newarray.data.shape[2] ].copy() )
        trimmedList.append(newarray)
        
    return trimmedList


def cleanIndices(indices, point, shape):
    """
    Auxiliary function to remove indices outside boundaries.
    """
    indices = npy.array(indices) + point
    # Remove indices outside superior bound
    out = npy.where(indices - shape >= 0)[0]
    indices = npy.delete(indices, out, axis=0)
    # Remove indices outside inferior bound
    out = npy.where(indices < 0)[0]
    indices = npy.delete(indices, out, axis=0)
    return indices

def getRadialIndices(radius, spacing, minradius=False):
    "Return an array of indices belonging to a radius in angstroms considering the grid spacing.\
    If 'point' is given as indices for the center, indices outside the box will be removed.Also indices returned will be centered at that point.\
    If 'minradius' in angstroms is given, points not reaching that minimum will not be including, returning the radial indices"
    span = npy.round(float(radius) / spacing).astype(int)
    indices = []
    for i in range(-span[0], span[0]+1):
        for j in range(-span[1], span[1]+1):
            for k in range(-span[2], span[2]+1):
                r = npy.linalg.norm([i,j,k]*spacing)
                if not minradius and r<=radius:indices.append((i,j,k))
                if minradius:
                    if r<=radius and r>=minradius: indices.append((i,j,k))
    return indices


###TESTING

import Biskit.test as BT
import tools as T

class Test(BT.BiskitTest):
    """Test"""
    def prepare(self):
        self.f_out = T.tempfile.mkdtemp()
        self.dxf = T.testRoot('grids','ETA_CT.dx')
        self.xplorf = T.testRoot('grids','ETA_CT.xplor')

    def cleanUp(self):
        T.tryRemove(self.f_out, tree=1)

    def test_GridData_XPLOR(self):
        """Load XPLOR grid"""
        g = GridData(self.xplorf)
        self.assertAlmostEqual(g.data.mean(), 4.215793, 4)
        self.assertEqual(g.delta, [0.5,0.5,0.5])
        self.assertEqual(g.origin, [-33., -33., -33.])

    def test_GridData_DX(self):
        """Load DX grid"""
        g = GridData(self.dxf)
        self.assertAlmostEqual(g.data.mean(), 4.215793, 4)
        self.assertEqual(g.delta, [0.5,0.5,0.5])
        self.assertEqual(g.origin, [-33., -33., -33.])

    def test_GridData_writeDX(self):
        """Read XPLOR, write DX"""
        os.chdir(self.f_out)
        gxplor = GridData(self.xplorf)
        gdx = GridData(self.dxf)
        gxplor.writeDX('test.dx')
        gread = GridData('test.dx')
        self.assertAlmostEqual(gread.data[10,10,10], gdx.data[10,10,10])
        self.assertEqual(gdx.origin, gread.origin)

    def test_GridData_writeXPLOR(self):
        """Read DX, write XPLOR"""
        os.chdir(self.f_out)
        gxplor = GridData(self.xplorf)
        gdx = GridData(self.dxf)
        gdx.writeXPLOR('test.xplor')
        gread = GridData('test.xplor')
        self.assertAlmostEqual(gxplor.data[10,10,10], gread.data[10,10,10])
        self.assertEqual(gxplor.origin, gread.origin)


if __name__ == '__main__':
    BT.localTest()