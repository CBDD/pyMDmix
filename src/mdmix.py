#! /usr/bin/env python
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
__date__ ="$11-mar-2014 17:33:33$"

import argparse
import sys
import pyMDMix
import numpy as npy

class MDMixError(Exception):
    pass

def createParser():
    parser = argparse.ArgumentParser("mdmix")
    parser.add_argument("--log", action="store", dest="logfile", help="Logging file. Default: output to stdout")
    parser.add_argument("--debug", action="store_true", dest="debug", default=False, help="Print debugging info")
    subparsers = parser.add_subparsers(help='commands', dest='command')
    
    #CREATE
    create_parser = subparsers.add_parser('create', help='Create Project or Solvents')
    create_parser.add_argument("action", choices=("project", "replica", "solvents", "template"), help="Create new project, a new stand-alone replica or add solvent to the database. 'template' will save a project template for new project creation.")
    create_parser.add_argument("-n", action="store", dest='projname', help="Project name or Replica Name. Default: mdmix_project.", default='mdmix_project')
    create_parser.add_argument("-f", action="store", dest="file", help="Configuration file path or output file name in 'template' action. Project can be created empty. Replica will need a config file with MDSETTINGS section.")
    create_parser.add_argument("-top", action="store", dest="top", help="Amber Topology file for stand-alone solvated replica creation.")
    create_parser.add_argument("-crd", action="store", dest="crd", help="Amber CRD file for stand-alone solvated replica creation.")

    #INFO
    info_parser = subparsers.add_parser('info', help='Print information about the project or solventDB')
    info_parser.add_argument("infoselection", choices=('project', 'solvents',), help="Print info about the project or about available solvents")
   
    #ADD
    add_parser = subparsers.add_parser('add', help="Add new replicas, systems or create group of replicas in an existing project.")
    add_parser.add_argument("action", choices=('system', 'replicas', 'group'), help="Create new system, new replicas or create groups of replicas.")
    add_parser.add_argument("-f", action="store", dest="file", help="REPLICAS action & SYSTEM action: Configuration files. MANDATORY if action is 'replicas' or 'system'.")
    add_parser.add_argument("-gn", action="store", dest="groupname", help="GROUP action: Name of the group of replicas to be created.")
    add_parser.add_argument("-s", action="store", dest="selection", help="GROUP action: List replicanames to add to the group.", nargs='+')
    add_parser.add_argument("-sys", action="store", dest="sysname", help="REPLICAS action: When creating new replicas, specify the system name which should be prepared. System name should exist in the project. If not given, only one system should be present in the project.")
   
    #REMOVE
    remove_parser = subparsers.add_parser('remove', help="Remove groups from project. To remove systems or replicas, simply remove the forlders or system files from the project folder.")
    remove_parser.add_argument("action", choices=('group'), help="Remove group.")
    remove_parser.add_argument("-gn", action="store", dest="groupname", help="Name of the group to remove. Groupname must exists in project.", required=True)
   
    #QUEUE
    queue_parser = subparsers.add_parser('queue', help="Queue input files options.")
    queue_parser.add_argument("action", choices=('list','write'), help="LIST: Show installed queue system templates. WRITE: Write input files for all replicas in current project or for REPLICA in current folder.")
    queue_parser.add_argument("-n", action="store", dest="queuename", help="WRITE action: queue system to use. Mandatory.")
       
    def subparseronreplica(subprs, name, help, extras=True):
        sp = subprs.add_parser(name, help=help)
        sp.add_argument("mode", choices=('all','bysolvent','byname','group'), action="store", default="all", help="Perform selection of replicas based on solvent name, replica names or groups. If 'all', do action on all replicas.")
        sp.add_argument('-s', action="store", nargs='+', dest="selection", help="Selection list. If selecting 'bysolvent', list of solvent names is expected. If 'byname', list of replica names. If 'group', group name. Skip if 'all' is selected.")
        if extras:
#            sp.add_argument("--force", help="Force calculations when some parts where already done. Default: False",action="store_true", dest="force", default=False)
            sp.add_argument("-N", help="For action commands, list production steps to consider for the analysis using a colon separated range. Ex: 1:20 - first to 20th step.", default=False, nargs=1,dest="nanoselect")
            sp.add_argument("--step", help="Take snapshots every STEP number. Default:1.", action="store", default=1, type=int, dest="step")
            sp.add_argument("-C", type=int, help="Number of cpus to use for the action. If option not given, will use 1 CPU serial mode.", dest="ncpus", default=1, action="store")
        return sp
        
    # PLOT SIMULATION OUTPUT
    plot_parser = subparsers.add_parser('plot', help="Plotting command")
    plot_cmd = plot_parser.add_subparsers(help='Plotting options', dest='plot_command')
    amber_plot = subparseronreplica(plot_cmd,'ambermd', help="Plot Amber MD properties",extras=False)
    amber_plot.add_argument("-o","--out",action="store", dest="outname", help="Name of the output file. Extension should be .png, .pdf, .jpeg, .ps or .eps. Default: ambermdplot.pdf", default="ambermdplot.pdf")
    amber_plot.add_argument("-N", help="List production steps to be plotted using a colon separated range. Ex: 1:20 - first to 20th step.", default=False, nargs=1,dest="nanoselect")
    rmsd_plot = subparseronreplica(plot_cmd,'rmsd', help="Plot Backbone and Heavy atoms RMSD plot",extras=False)
    rmsd_plot.add_argument("-o","--out",action="store", dest="outname", help="Name of the output file. Extension should be .png, .pdf, .jpeg, .ps or .eps. Default: rmsdplot.pdf", default="rmsdplot.pdf")
    rmsd_plot.add_argument("-N", help="List production steps to be plotted using a colon separated range. Ex: 1:20 - first to 20th step.", default=False, nargs=1,dest="nanoselect")
    
    # ANALYSIS OPTIONS
#    anl_parser = parser.add_subparsers(help='Analysis commands', dest='anl_command')
    anl_parser = subparsers.add_parser("analyze", help="Several analysis tools to run on the replicas.")
    anl_cmd = anl_parser.add_subparsers(help='Analysis commands', dest='anl_command')

    # ANALYSIS: DENSITY
#    dens_commands = anl_parser.add_subparsers(help='Analysis commands', dest='anl_command')
    density_parser = subparseronreplica(anl_cmd, 'density', help="Calculate probe density grids from aligned trajectory")
    density_parser.add_argument("--probes", "-P",nargs='+', dest="probelist", help="Selection of probenames to calculate density for. If not given, all probes for the solvent will be selected.")
    density_parser.add_argument("-opref", action="store", dest="outprefix", help="Prefix for grids saved inside density folder for each replica. If this is given, automatic energy conversion will not work until you restablish expected names or explicitely give the prefix in energy command. Default: False")
    density_parser.add_argument("--onlycom", help="Density calculations ONLY for center of masses of each co-solvent. Probe list is ignored. Default: False",action="store_true", dest="onlycom", default=False)
    density_parser.add_argument("--com", help="Include center of mass to list of probes probes in density calculation. Default: False",action="store_true", dest="com", default=False)
    density_parser.add_argument("--ref", help="Reference PDB file to build grid shape over which counts will be added. By default, replica reference PDB will be used. This option is useful to build density grids on a subregion of the space or when replica trajectory was also aligned using --ref option.", action="store", dest="ref")

    # ANALYSIS: RESIDENCE
    residence_parser = subparseronreplica(anl_cmd, 'residence', help="Calculate hotspot residence dens from aligned trajectory")
    residence_parser.add_argument("--hpfile", "-hf", action="store", dest="hpfile", help="HotspotSet pickled file from 'analyze hotspots create' action.")
    residence_parser.add_argument("--hpid", "-id", action="store", type=int, dest="hpid", help="Hotspot ID in HPFILE to define region to study.")
    residence_parser.add_argument("--center", "-ce", nargs='+', dest="center", help="Sphere center to define region to study. Give a space separated 3 float list.")
    residence_parser.add_argument("--tol", "-t", default=0.5,type=float,dest="tolerance", help="Tolerance in angstroms around sphere center or hotspot coordinates to consider as occupied space. (default: 0.5)")
    

    # ANALYSYS: ALIGN
#    anl_commands = anl_parser.add_subparsers(help='Analysis commands', dest='anl_command')
    align_parser = subparseronreplica(anl_cmd, 'align', 'Align trajectory of selected replicas', extras=False)
    align_parser.add_argument("-N", help="List production steps to consider for alignment using a colon separated range. Ex: 1:20 - first to 20th step.", default=False, nargs=1,dest="nanoselect")
    align_parser.add_argument("-C", type=int, help="Number of cpus to use for the action. If option not given, will use 1 CPU serial mode.", dest="ncpus", default=1, action="store")
    align_parser.add_argument("--mask", action="store", dest="alignmask", help="Modify alignment mask defined when creating the replicas. By default the macromolecule will be automatically identified. Give a list with comma separated residue numbers or hyphen separated range. E.g. 10-100,120-240.")
    align_parser.add_argument("--ref", action="store", dest="ref", help="Path to reference PDB file. By default, pyMDmix generates one automatic reference pdb file which can be found inside each replica folder. This option will override it.")
    align_parser.add_argument("--only-write", action="store_true", default=False, dest="onlywrite", help="Only write ptraj input scripts BUT don't execute them. Useful when manual editing is needed. (default: False)")
    align_parser.add_argument("--only-exe", action="store_true", default=False, dest="onlyexe", help="Only execute existing ptra scripts, do not overwrite them. If scripts don't exist, this function will fail. (Default: False)")
    
    # In conversion from density to grids, allow also merging of densities between replicas before energy conversion
    energy_parser = subparseronreplica(anl_cmd, 'energy', help="Calculate free energy maps from density grids for selected probes and replicas", extras=False)
    energy_parser.add_argument("-nsnaps", action="store", dest="nsnaps", type=int, help="If given, use this number of snapshots for calculating the expected number instead of the total number. Useful when a subset of the trajectory is analyzed.")
    energy_parser.add_argument("-N", help="Specify production steps used for density grid calculation to correctly calculate the expected number when a subset of the trajectory was analyzed. Use a colon separated range. Ex: 1:20 - first to 20th step. If not given, all trajectory is considered.", default=False, nargs=1,dest="nanoselect")
    energy_parser.add_argument("--probes", "-P",nargs='+', dest="probelist", help="Selection of probenames to convert. If not given, all probes for the solvent will be converted.")
    energy_parser.add_argument("-noavg",action="store_false",dest="noavg",default="true", help="By default, densities for the diferent replicas will be  merged before energy conversion and a single replica-averaged energy map for each probe will be saved at project folder. \
                                                                                            To save separately each replica it's own energy grids, use -noavg. Energy grids will then be saved inside each replica folder independently.")
    energy_parser.add_argument("-ipref", action="store", dest="inprefix", help="If density grids were saved with a specific prefix, give it here so the program knows what density grids to take. Default: no prefix (predefined names).")
    energy_parser.add_argument("-opref",action="store",dest="outprefix",help="Prefix for output average grids. Default: no prefix (predefined names).")
    energy_parser.add_argument("-nodg0",action="store_false",dest="nodg0",default="true", help="Disable standard state correction. Ignore concentration issues")

    # ANALYSYS: Hotspots
    hs_parser = anl_cmd.add_parser("hotspots", help="Analysis tools on hot spots")
    hs_cmd = hs_parser.add_subparsers(help='Hot spots analysis commands', dest='hs_command')
    
    createhs = hs_cmd.add_parser("create", help="Create hot spots sets from selected grids using cutoff clustering method")
    createhs.add_argument("-i", action="store", dest="ingrids", required=True, nargs='+', help="List of grid files to use for hotspots creation. If multiple probes are given, a hot spot set with the minimum energy probes is obtained.")
    createhs.add_argument("-o", action="store", required=True, dest="outprefix", help="Output prefix. A PDB with the hotspots and a pickle file will be saved")
    createhs.add_argument("-centroid",action="store_true", dest="centroid", default="false", help="Calculate centroid as the average of all coordinates. Default: Minimum energy coordinate.")
    createhs.add_argument("--allpoints",action="store_false", dest="onlycenter", default="true", help="Write all points belonging to the hotspots in the output PDB instead of only the minimums with mean energy")
    createhs.add_argument("-p", action="store", dest="percentile", default=0.02, help="Percentile of points to use for establishing the cutoff value. Default: 0.02")
    createhs.add_argument("-x", action="store", dest="ignore", type=float, help="Ignore points with this value during percentile calculation. Useful when grids have values masking certain regions in space (like 999 to mask protein space).")
    createhs.add_argument("--hardcutoff", "-H", default=False, type=float, action="store", dest="hardcutoff", help="Stablish a hard energy cutoff instead of the adaptive one.")
    
    createhs_min = hs_cmd.add_parser("create_min", help="Create hot spots sets from selected grids using Minima Search method")
    createhs_min.add_argument("-i", action="store", dest="ingrids", required=True, nargs='+', help="List of grid files to use for hotspots creation. If multiple probes are given, a hot spot set with the minimum energy probes is obtained.")
    createhs_min.add_argument("-o", action="store", required=True, dest="outprefix", help="Output prefix. A PDB with the hotspots and a pickle file will be saved")
    createhs_min.add_argument("-cut", required=True, type=float, action="store", dest="cutoff", help="Stablish an energy cutoff to stop searching for more minima.")
    createhs_min.add_argument("--meanradius", type=float, action="store",default=0.5, dest="meanradius", help="Distance around the minimum which values will be used to calculate hotspot average value. DEFAULT: 0.5.")
    createhs_min.add_argument("--cancelradius", type=float, action="store",default=2.0, dest="cancelradius", help="Distance around the minimum that will be discarded in further iterations. DEFAULT: 2.0.")
    createhs_min.add_argument("--ignoreval", type=float, action="store",default=999, dest="ignoreval", help="Ignore values equal to this number (e.g. excluded volume with 999 value). DEFAULT: 999.")
    
    
    
    # UTILS
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
    extend_parser = subparseronreplica(util_cmd, 'extendsim', help="Extend simulation length for selected replicas", extras=False)
    extend_parser.add_argument("-n",action="store",type=int, dest="extrananos",help="Extra number of nanoseconds to extend selected replicas", required=True)


#    # ON THE FLY ALIGNMENT
#    onfly_parser = subparseronreplica('trackalign', 'Track on-runing simulation and run trajectory alignment for already completed nanoseconds')
#
#    # BFACTORS FOR ALL PROTEIN
#    bfact_parser = subparseronreplica('bfactors', help="Calculate protein bfactors for selected replicas and save results in PDB format and text file.")
#    bfact_parser.add_argument("--prefix", action="store", dest="prefix", default="bfactors", help="Prefix to add to output pdb and txt files (replica name will be automatically included). Default='bfactors'")
#
    return parser

def fetchReplicaSelection(parserargs, project):
    "Function that checks and returns replicas matching args in current project."
    if parserargs.mode == 'all':
        #return all replicas in project
        returnlist = project.replicas.values()
        print "Selected replica names: %s"%[r.name for r in returnlist]
        return returnlist
    elif parserargs.mode == 'bysolvent' or parserargs.mode == 'byname' or parserargs.mode == 'group':
        #Check a selection list was given
        if not parserargs.selection:
            raise AttributeError, "Selection (-s fag) is mandatory when 'byname', 'bysolvent' or 'group' mode is chosen."
        selection = parserargs.selection

        if parserargs.mode == 'bysolvent':
            #Check all solvents are present in project
            for s in selection:
                if s not in project.solventCounter.keys():
                    print >>sys.stderr, "Solvent %s not in project. Skipping."%s
            #Return only replicas matching solvents in selection list
            returnlist = []
            [returnlist.append(r) for r in project.replicas.values() if r.solvent in selection]
            
            if returnlist:
                print "Selected replica names: %s"%[r.name for r in returnlist]
                return returnlist
            else:
                raise MDMixError, "Replicas not found."

        if parserargs.mode == 'byname':
            #Return only replicas is name matches selection
            returnlist = []
            [returnlist.append(r) for r in project.replicas.values() if r.name in selection]

            if returnlist:
                print "Selected replica names: %s"%[r.name for r in returnlist]
                return returnlist
            else:
                raise MDMixError, "Replicas not found."
        
        if parserargs.mode == 'group':
            for s in selection:
                if s not in project.listGroups():
                    print >> sys.stderr, "Groupname %s not in current project. Skipping..."%s
                returnlist = []
                [returnlist.extend(project.getGroup(s)) for s in selection if s in project.listGroups()]
                
                if returnlist:
                    print "Selected replica names: %s"%[r.name for r in returnlist]
                    return returnlist
                else:
                    raise MDMixError, "Replicas not found"
    else:
        return False

def returnMDMixProject(parserargs):
    if parserargs.debug: level='DEBUG'
    else: level='INFO'
    pyMDMix.setLogger(level=level)
    try:
        p = pyMDMix.loadProject()
        return p
    except:
        return False

def returnMDMixProjectOrFail(parserargs):
        #When command is different to CREATE PROJECT or INFO, this program should be executed in project folder
        #Let's try to load a project or exit
        p = returnMDMixProject(parserargs)
        if not p: raise MDMixError, 'No project file found in current folder. Make sure you are in a pyMDMix project folder.'
        return p

def parseNumMask(mask):
    "Transform mask of form 1,2,5:10 to a list of ints"
    s = mask.split(',')
    out = []
    for el in s:
        if ':' in el:
            a, b = map(int, el.split(':'))
            out += npy.arange(a,b+1).tolist()
        else:
            out += [int(el)]
    out.sort()
    return out

def parsenanos(argparser):
    if not argparser.nanoselect: return False
    nanosel = parseNumMask(argparser.nanoselect[0])
    print "Selected steps: %s"%(', '.join(map(str, nanosel)))
    return nanosel

def createRootLogger(level, logFile=None):
    pyMDMix.setLogger(level, logFile)

def header():
    
    print """
    ==========================================================
    ||              pyMDMix User Interface                  ||
    ==========================================================
    ||  Author: Daniel Alvarez-Garcia                       ||
    ||  Version : %s                                     
    ==========================================================
    
    """%pyMDMix.__version__


def main():
    import time
    import os
    import sys
    t0 = time.time()
    parser = createParser()
    parserargs =  parser.parse_args()
    command = parserargs.command
    done = False

    # If logging file, redirect stdout and stderr to file
    if parserargs.logfile:
        sys.stderr = open(parserargs.logfile, 'w+')
        sys.stdout = sys.stderr
        print ' '.join(sys.argv)

    if parserargs.debug: level = 'DEBUG'
    else: level = 'INFO'
    createRootLogger(level,parserargs.logfile)

    #Study the different actions
    if command == 'create':
        if parserargs.action == 'project':
            #Create Project. Check config file exists and run creation.
            if parserargs.file:
                if os.path.exists(parserargs.file):
                    ######### ALL correct! Proceed with MDMix Project Creation
                    print "Creating PROJECT %s from config file: %s"%(parserargs.projname, parserargs.file)
                    p = pyMDMix.createProject(parserargs.file, parserargs.projname)
                    print "DONE"
                    #########
                else:
                    raise MDMixError, "File %s not found. Cannot create the project."%parserargs.file
            else:
                # CREATE empty project
                p = pyMDMix.Project(name=parserargs.projname)
                p.createProjectFolder()
                print "DONE"
                
        elif parserargs.action == 'replica':
            name = parserargs.projname
            file = parserargs.file
            top = parserargs.top
            crd = parserargs.crd
            if not file: raise MDMixError, "Input config file is needed with valid MDSETTINGS section for new replica creation from TOP and CRD files."
            if not top or not crd: raise MDMixError, "Amber Topology (-top) and Amber Crd files (-crd) are needed to create a solvated replica."
            sysname = os.path.splitext(os.path.basename(top))[0]
            solvatedsys = pyMDMix.SolvatedSystem(name=sysname,top=top, crd=crd)
            sets = pyMDMix.parseSettingsConfigFile(file, noSolvent=True) # Parse MDSETTINGS ignoring solvent info
            repl = solvatedsys+sets
            repl.setName(name=name)
            repl.createAll()
                        
        elif parserargs.action == 'solvents':
            # CREATE NEW SOLVENT IN THE DATABASE
            #Checking mandatory file option is given and exists
            if not parserargs.file:raise MDMixError, "Missing file in solvents create action. -f FILE is mandatory in this option."
            file = parserargs.file
            if not os.path.exists(file): raise MDMixError, "File not found: %s"%file
            man = pyMDMix.Solvents.SolventManager()
            solv = man.createSolvent(file)
            man.saveSolvent(solv)
            print "DONE"

        elif parserargs.action == 'template':
            import os
            import shutil
            # COPY PROJECT TEMPLATE
            f = pyMDMix.tools.templatesRoot('complete.cfg')
            file = parserargs.file or 'complete.cfg'
            shutil.copy(f,file)
            print "DONE. Project template saved: %s"%file

    elif command == 'info':
        if parserargs.infoselection == 'project':
            p = returnMDMixProjectOrFail(parserargs)
            print p.longdesc()
        elif parserargs.infoselection == 'solvents':#Actions on solvent database
            man = pyMDMix.Solvents.SolventManager()
            print man

    elif command == 'add':
        p = returnMDMixProjectOrFail(parserargs)
        action = parserargs.action
        if action == 'system' or action == 'replicas':
            file = parserargs.file
            if not file:
                raise MDMixError, "Configuration file needed to add new systems or replicas (use -f option)."
        if action == 'system':
            file = parserargs.file
            system = pyMDMix.parseSystemConfigFile(file)
            p.addNewSystems(system)
            print "DONE"
        elif action == 'replicas':
            file = parserargs.file
            sysname = parserargs.sysname
            if not sysname:
                avail = p.systems.keys()
                if len(avail) == 1: sysname = avail[0]
                else: raise MDMixError, "More than one system in current project. Choose which one you wish to prepare with -sys option. Available systems: %s"%avail
            else:
                if not sysname in p.systems.keys():
                    raise MDMixError, "Wrong system name. Project systems are: %s"%p.systems.keys()
            print "Creating replicas for system %s"%sysname
            settings = pyMDMix.parseSettingsConfigFile(file)
            p.createReplicas(sysname, settings)
            print "DONE"
        elif action == 'group':
            if not parserargs.groupname:
                raise MDMixError, "Groupname is required (-gn option)."
            if not parserargs.selection:
                raise MDMixError, "Selection list is mandatory for creating group %s (use -s option)."%parserargs.groupname
            p.createGroup(parserargs.groupname, parserargs.selection)
            print "DONE"

    elif command == 'remove':
        p = returnMDMixProjectOrFail(parserargs)
        if parserargs.action == 'group':
            name = parserargs.groupname
            r = p.removeGroup(name)
            if not r:
                raise MDMixError, "Group name %s does not exists. Project groups: %s"%p.listGroups()
            print "DONE"

    elif command == 'queue':
        if parserargs.action == 'list':
            q = pyMDMix.Queue.listQueueSystems()
            if q:
                print '\nInstalled queue system templates: %s\n'%(', '.join(q))
            else:
                print "\nNo queue templates found!\n"
        elif parserargs.action == 'write':
            qname = parserargs.queuename
            qlist = pyMDMix.Queue.listQueueSystems()
            if not qname:
                raise MDMixError, "Queue name is needed: give it with option (-n). Avilable queues: %s"%qlist
            if qname not in qlist:
                raise MDMixError, "Wrong queue name. Available queues: %s"%qlist
            p = returnMDMixProject(parserargs)
            if p: p.createQueueInputs(qname)
            else:
                # No project in current folder, try if there is a replica
                try:
                    r = pyMDMix.loadReplica()
                    r.createQueueInput(qname)
                except:
                    raise MDMixError, "Could not find any project file or replica file in current folder."

    #  PLOT COMMANDS
    elif command == 'plot':
        p = returnMDMixProjectOrFail(parserargs)
        import pyMDMix.Plotter as P
        plot = P.Plot()
        replicas = fetchReplicaSelection(parserargs, p)
        stepselection=parsenanos(parserargs)
        if parserargs.plot_command == 'ambermd':
            if replicas: plot.plotMDAmber(replicas, outfilename=parserargs.outname, selectedsteps=stepselection)
        elif parserargs.plot_command == 'rmsd':
            if replicas: plot.plotRMSDReplicas(replicas, outfilename=parserargs.outname, selectedsteps=stepselection)
                

    # ANALYSIS COMMANDS
    elif command == 'analyze':
        if parserargs.anl_command != 'hotspots': p = returnMDMixProjectOrFail(parserargs)
        if parserargs.anl_command == 'align':
            replicas = fetchReplicaSelection(parserargs, p)
            if replicas:
                nanosel = parsenanos(parserargs)
                ncpus = parserargs.ncpus
                mask = parserargs.alignmask
                onlywrite = parserargs.onlywrite
                onlyexe = parserargs.onlyexe
                reference = parserargs.ref
                
                if onlywrite or onlyexe:
                    # Options for partial execution selected
                    if onlywrite:
                        # Only write ptraj input scripts and do not execute them
                        p.alignReplicas(replicas, ncpus=ncpus, steps=nanosel, alignmask=mask, reference=reference, run=False)
                    else:
                        # Just execute existing alignment scripts
                        # Fail if scripts do not exist
                        p.alignReplicas(replicas, ncpus=ncpus, steps=nanosel, alignmask=mask, run=True, reference=reference, write=False)
                else:
                    p.alignReplicas(replicas, ncpus=ncpus, steps=nanosel, reference=reference, alignmask=mask)
                print "DONE"

        elif parserargs.anl_command == 'density':
            replicas = fetchReplicaSelection(parserargs, p)
            nanosel = parsenanos(parserargs)
            ncpus = parserargs.ncpus
            step = parserargs.step
            outprefix = parserargs.outprefix
            probelist = parserargs.probelist
            includeCOM=parserargs.com
            onlyCOM=parserargs.onlycom
            ref=parserargs.ref
            
            if replicas:
                anal = pyMDMix.Analysis.ActionsManager(ncpus=ncpus)
                anal.addReplicas(replicas)
                anal.addActions('DensityGrids')
                anal.prepareRun(probeselection=probelist, outprefix=outprefix, 
                                includeCOM=includeCOM, onlyCOM=onlyCOM, stepselection=nanosel,
                                reference=ref)
                anal.run(stepselection=nanosel, framestep=step)
            
            print "DONE"

        elif parserargs.anl_command == 'residence':
            replicas = fetchReplicaSelection(parserargs, p)
            nanosel = parsenanos(parserargs)
            ncpus = parserargs.ncpus
            step = parserargs.step
            hotspotfile = parserargs.hpfile
            hselection = parserargs.hpid
            center = parserargs.center
            tolerance = parserargs.tolerance
            hotspot = False

            if hotspotfile and hselection:
                # Read hotspots from pickle file and
                import cPickle
                inhset = cPickle.load(open(hotspotfile,'rb'))
                hotspot = inhset.getHSbyID(hselection)
                if not hotspot: raise MDMixError, "No hotspot selected."
                
            elif center:
                # Study residence at a hotspot defined by sphere and radius
                center = npy.array(center,dtype=float)
            else:
                raise MDMixError, "To do a residence analysis, the spot to study must be defined either with a sphere giving a center+tolerance \
                                    or by giving a hotspot pickled file and an ID identifying the hotspot to use."
                
            if replicas:
                anal = pyMDMix.Analysis.ActionsManager(ncpus=ncpus)
                anal.addReplicas(replicas)
                anal.addActions('Residence')
                anal.prepareRun(hotspot=hotspot, spherecenter=center, tolerance=tolerance, stepselection=nanosel)
                anal.run(stepselection=nanosel, framestep=step)
                anal.processResults()
            print "DONE"

        elif parserargs.anl_command == 'energy':
            replicas = fetchReplicaSelection(parserargs, p)
            if replicas:
                inprefix=parserargs.inprefix
                outprefix=parserargs.outprefix
                nsnaps = parserargs.nsnaps
                probelist = parserargs.probelist
                dg0 = parserargs.nodg0
                avg = parserargs.noavg
                nanosel = parsenanos(parserargs)

                from pyMDMix.Energy import EnergyConversion
                econv = EnergyConversion()
                econv.convert(replicas, probelist=probelist, average=avg, dg0correct=dg0,
                                    inprefix=inprefix, outprefix=outprefix, nsnaps=nsnaps, stepselection=nanosel)

        # HOTSPOTS SUBPARSER UNDER ANALYSIS
        elif parserargs.anl_command == 'hotspots':
            
            if parserargs.hs_command == 'create':
                ingrids = parserargs.ingrids
                outprefix = parserargs.outprefix
                percentile = parserargs.percentile
                hardcutoff = parserargs.hardcutoff
                centroid = parserargs.centroid
                maskcutvalue = parserargs.ignore
                
                for g in ingrids:
                    if not os.path.exists(g): raise MDMixError, "File %s not found."%g
                
                import pyMDMix.HotSpotsManager as HM
                if centroid: centroid = 'avg'
                else: centroid = 'min'
                HM.createHotSpotsByCutoff(ingrids, outprefix=outprefix, percentile=percentile, 
                                            cutoff=hardcutoff, centroid=centroid, 
                                            onlycenter=parserargs.onlycenter, maskcutvalue=maskcutvalue)
            elif parserargs.hs_command == 'create_min':
                ingrids = parserargs.ingrids
                outprefix = parserargs.outprefix
                meanradius = parserargs.meanradius
                cancelradius = parserargs.cancelradius
                ignorevalue = parserargs.ignoreval
                cutoff = parserargs.cutoff
                
                for g in ingrids:
                    if not os.path.exists(g): raise MDMixError, "File %s not found."%g
                
                import pyMDMix.HotSpotsManager as HM
                HM.createHotSpotsByMinSearch(ingrids, cutoff, outprefix=outprefix, meanradius=meanradius, 
                                            cancelRadius=cancelradius, protValue=ignorevalue)
                
    # UTILS COMMANDS
    elif command == 'tools':
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
                npy.savetxt(out, results, fmt="%.3f")
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
            trimmed = trim(smallgrid, ingrid)
            trimmed[1].writeDX(outname)
            print "DONE trimming"
        
        elif parserargs.util_command == 'projecttemplate':
            import shutil
            fname = parserargs.filename
            shutil.copy(pyMDMix.tools.templatesRoot('project.cfg'), fname)

        elif parserargs.util_command == 'extendsim':
            # EXTEND SIMULATIONS
            p = returnMDMixProjectOrFail(parserargs)
            replicas = fetchReplicaSelection(parserargs, p)
            if replicas:
                p.extendSimulations(replicas, parserargs.extrananos)

    sys.stdout.flush()
    print "Total execution time: %.3fs"%(time.time()-t0)

if __name__ == "__main__":
    try: 
        header()
        main()
    except KeyboardInterrupt: print "Forcing MDMix UI exit!"
    sys.exit(0)
