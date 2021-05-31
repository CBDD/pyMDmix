import argparse
import os
import numpy
import pyMDMix
from pyMDMix import MDMixError
from pyMDMix.Commands.Command import Command
from pyMDMix.Projects import returnMDMixProject, returnMDMixProjectOrFail
class Tools(Command):
    def __init__(self):
        self.cmdstring = "tools"

    def create_parser(self, subparsers):
        util_parser = subparsers.add_parser("tools", help="Complementary tools")
        util_cmd = util_parser.add_subparsers(help='Tool commands', dest='util_command')

        #GET VALUES
        getvals_parser = util_cmd.add_parser('getvalues', help="Get energy values from grids at given coordinates.")
        getvals_parser.add_argument("-i", action="store", required=True, dest="ingrid", help="Input grid to query")
        getvals_parser.add_argument("-o", action="store", required=False, dest="outfile", help="Output file with results. Print to STDOUT if not given.")
        getvals_parser.add_argument("-xyz", action="store", dest="xyzfile", help="Coordinates text file. Each line is a coordinate (x y z) plus an optional radius (r) argument.")
        getvals_parser.add_argument("-pdb", action="store", dest="pdbfile", help="Define coordinates using a PDB file with optional radius in bfactor column")
        getvals_parser.add_argument("-r", action="store", dest="radius", type=float, help="Radius around each coordinate to average values")

        # SUBTRACT GRIDS
        diffgrids_parser = util_cmd.add_parser('diffgrids', help="Save new grid with the difference from other two (grid1 - grid2)")
        diffgrids_parser.add_argument("-g1", action="store", required=True, dest="grid1", help="Input grid file path (grid 1)")
        diffgrids_parser.add_argument("-g2", action="store", required=True, dest="grid2", help="Input grid file path (grid 2)")
        diffgrids_parser.add_argument("-o", action="store", required=True, dest="outfile", help="Output grid file name (DX format)")

        # SUM GRIDS
        sumgrids_parser = util_cmd.add_parser('sumgrids', help="Save new grid with the sum from other two (grid1 + grid2)")
        sumgrids_parser.add_argument("-g1", action="store", required=True, dest="grid1", help="Input grid file path (grid1)")
        sumgrids_parser.add_argument("-g2", action="store", required=True, dest="grid2", help="Input grid file path (grid 2)")
        sumgrids_parser.add_argument("-o", action="store", required=True, dest="outfile", help="Output grid file name (DX format)")

        # TRIM USING PDB
        trim_parser = util_cmd.add_parser('trim', help="Use a PDB file as reference to crop/trim a grid and save just a part of the space (useful to reduce sizes and save disk space). WARNING: There might be a slight shift of coordinates in the trimmed grid.")
        trim_parser.add_argument("-i", action="store", required=True, dest="ingrid", help="Input grid file to trim.")
        trim_parser.add_argument("-o", action="store", required=True, dest="outname", help="Out grid file name where to save trimmed grid.")
        trim_parser.add_argument("-pdb", action="store", required=True, dest="refpdb", help="PDB file used as reference to trim input grid. All coordinates will be used to select the space.")
        trim_parser.add_argument("--buff", action="store", dest="buff", default=8.0, help="Add this buffer distance (Angstroms) to minimum and maximum coordinates extracted from the PDB when trimming. DEFAULT: 8.0 Angstroms")

        # WRITE PROJECT TEMPLATE FILE
        write_template = util_cmd.add_parser('projecttemplate', help="Write template for project configuration")
        write_template.add_argument("-f", action="store", required=True, dest="filename", help="Filename to save")

        # EXTEND SIMULATIONS
        extend_parser = self.subparseronreplica(util_cmd, 'extendsim', help="Extend simulation length for selected replicas", extras=False)
        extend_parser.add_argument("-n",action="store",type=int, dest="extrananos",help="Extra number of nanoseconds to extend selected replicas", required=True)

    def action(self, parserargs):
        if parserargs.util_command == 'getvalues':
            grid = parserargs.ingrid
            out = parserargs.outfile
            rad = parserargs.radius or 0
            xyz = parserargs.xyzfile
            pdb = parserargs.pdbfile
            if not pdb and not xyz: raise MDMixError, "Either PDB or XYZ file required."
            if not os.path.exists(grid): raise MDMixError, "Grid file %s not found"%grid

            from pyMDMix.GridsManager import getEnergyFromTxtCoords, getEnergyFromPDBCoords

            if pdb:
                print "Fetching values at PDB coordinates"
                results = getEnergyFromPDBCoords(grid, pdb, forceradius=rad)
            else:
                print "Fetching values at XYZ file coordinates"
                results = getEnergyFromTxtCoords(grid, xyz, radius=rad)

            if out:
                print "Saving results to %s"%out
                numpy.savetxt(out, results, fmt="%.3f")
            else:
                print '\n'.join(['%.3f'%r for r in results])
            print "DONE"
        
        elif parserargs.util_command == 'diffgrids':
            grid1 = parserargs.grid1
            grid2 = parserargs.grid2
            out = parserargs.outfile
            
            for grid in (grid1, grid2):
                if not os.path.exists(grid): raise MDMixError, "Grid file %s not found"%grid

            from pyMDMix.GridsManager import gridDifference
            gridDifference(grid1, grid2, out)
        
        elif parserargs.util_command == 'sumgrids':
            grid1 = parserargs.grid1
            grid2 = parserargs.grid2
            out = parserargs.outfile
            
            for grid in (grid1, grid2):
                if not os.path.exists(grid): raise MDMixError, "Grid file %s not found"%grid

            from pyMDMix.GridsManager import gridSum
            gridSum(grid1, grid2, out)
        
        elif parserargs.util_command == 'trim':
            ingrid = parserargs.ingrid
            outname = parserargs.outname
            refpdb = parserargs.refpdb
            buffer = parserargs.buff
            
            if not os.path.exists(refpdb): raise MDMixError, "PDB file %s not found"%refpdb
            if not os.path.exists(ingrid): raise MDMixError, "Grid file %s not found"%ingrid
            
            from pyMDMix.GridsManager import trim, Grid, GridData
            smallgrid = GridData.createFromPDB(refpdb, spacing=0.5, buff=buffer, takeProtein=False)
            ingrid = Grid(ingrid)
            # trimmed = trim(smallgrid, ingrid)
            # trimmed[1].writeDX(outname)
            trim(smallgrid, ingrid)[1].writeDX(outname)
            print "DONE trimming"
        
        elif parserargs.util_command == 'projecttemplate':
            import shutil
            fname = parserargs.filename
            shutil.copy(pyMDMix.tools.templatesRoot('project.cfg'), fname)

        elif parserargs.util_command == 'extendsim':
            # EXTEND SIMULATIONS
            p = returnMDMixProjectOrFail(parserargs)
            replicas = self.fetchReplicaSelection(parserargs, p)
            if replicas:
                p.extendSimulations(replicas, parserargs.extrananos)
