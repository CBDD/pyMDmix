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

__author__="dalvarez"
__date__ ="$16-ene-2014 17:09:33$"

import os.path as osp
import settings as S
import Biskit as bi


class TrajFileError(Exception):
    pass

class TrajectoryError(Exception):
    pass

class TrajFile(object):
    def __init__(self, file, pdb, step=1, frameselection=[]):
        "pdb is a PDBModel, file should meet expected extensions (nc, netcdf, dcd, x)"
        self.fname = file
        self.pdb = pdb
        self.nframes = 0
        self.step=step
        self.frameselection = frameselection
        self.extension = osp.splitext(self.fname)[1].lstrip('.')
        if not self.extension in S.avail_trajext:
            raise TrajFileError, "Wrong extension %s. Expected extensions: %s"%(self.extension, ','.join(S.avail_trajext))
        self.loadFile()

    def __iter__(self):
        return self

    def loadFile(self):
        if self.extension in ('x','x.gz'):
            self.traj = bi.AmberCrdParser(self.fname, self.pdb, 1)
            self.traj.crd.readline()    #skip first line
            self.nextFunction = self.traj.nextFrame
        elif self.extension in ('nc','netcdf'):
            try:
                #import Scientific.IO.NetCDF as ncdf
                import scipy.io.netcdf as ncdf
                self.traj = ncdf.NetCDFFile(self.fname,'r').variables['coordinates']
                self.nextFunction = self.returnFrameFromNetcdf
            except ImportError:
                raise TrajFileError, "Can't read NetCDF trajectory"
        elif self.extension in ('dcd',):
            from NamdDCDParser import NamdDCDParser
            self.traj = NamdDCDParser(self.fname, self.pdb, box=1)
            self.nextFunction = self.traj.read_dcdstep

    def returnFrameFromNetcdf(self):
        return self.traj[self.nframes,:,:]

    def next(self):
        try:
            frame = self.nextFunction()
            self.nframes += 1
            if self.frameselection:
                while not self.nframes in self.frameselection:
                    try: 
                        frame = self.nextFunction()
                        self.nframes += 1
                    except:
                        raise StopIteration()
                return frame
            else:
                while self.nframes%self.step:   # Continue looping until matching spected step
                    frame = self.nextFunction()
                    self.nframes += 1
                return frame
        except:
            raise StopIteration()

    def close(self):
        if self.extension == 'dcd' or self.extemsion in ('nc','nectdf'): self.traj.close()
        else: pass

class Trajectory(object):
    """Allow reading trajectories from NAMD or AMBER"""
    def __init__(self, filelist, pdb, step=1, frameselection=[], **kwargs):
        """
        Trajectory parser for files in filelist. Files should have expected extension names:
        dcd for NAMD and nc, netcdf, x, x.gz for AMBER.
        
        :arg int step: Return frame every *step* frames.
        :arg list frameselection: List of integers selecting specific frames from each file.
        """
        self.files = filelist
        self.pdb = pdb
        self.step = step
        self.frameselection = frameselection
        self.nfiles = len(self.files)
        self.file_i = 0

    def __iter__(self):
        return self

    def next(self):
        "i is the index of the file in the files list"
        if self.file_i < self.nfiles:
            self.file = TrajFile(self.files[self.file_i], self.pdb, step=self.step, frameselection=self.frameselection)
            self.file_i += 1
            return self.file
        else:
            self.file_i = 0
            raise StopIteration()


import Biskit.test as BT
import tools as T

class Test(BT.BiskitTest):
    """Test"""

    def test_ReadNAMDDCD(self):
        """Read DCD trajectory"""
        dcd = T.testRoot('namd','min.dcd')
        pdb = bi.PDBModel(T.testRoot('namd','pep_WAT_WAT_1.pdb'))
        try:
            traj = TrajFile(dcd, pdb)
        except:
            traj = False
        self.assertTrue(traj)

    def test_readAmberAscii(self):
        """Read Amber Ascii trajectory"""
        x = T.testRoot('amber','traj.x')
        pdb = bi.PDBModel(T.testRoot('amber','pep_WAT_WAT_1.pdb'))
        try:
            traj = TrajFile(x, pdb)
        except:
            traj = False
        self.assertTrue(traj)

    def test_readAmberNetCDF(self):
        """Read Amber NetCDF trajectory"""
        nc = T.testRoot('amber','traj.nc')
        pdb = bi.PDBModel(T.testRoot('amber','pep_WAT_WAT_1.pdb'))
        try:
            traj = TrajFile(nc, pdb)
        except:
            traj = False
        self.assertTrue(traj)
#
#    def cleanUp(self):
#        if self.f_out: T.tryRemove( self.f_out, tree=1 )


if __name__ == '__main__':
    BT.localTest()

