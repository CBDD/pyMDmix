import argparse
import os
import numpy
import pyMDMix
from pyMDMix import MDMixError
from pyMDMix.Commands.Command import Command
from pyMDMix.Projects import returnMDMixProject, returnMDMixProjectOrFail
class Analyze(Command):
    def __init__(self):
        self.cmdstring = "analyze"

    def create_parser(self, subparsers):
        anl_parser = subparsers.add_parser("analyze", help="Several analysis tools to run on the replicas.")
        anl_cmd = anl_parser.add_subparsers(help='Analysis commands', dest='anl_command')

        # ANALYSIS: DENSITY
    #    dens_commands = anl_parser.add_subparsers(help='Analysis commands', dest='anl_command')
        density_parser = self.subparseronreplica(anl_cmd, 'density', help="Calculate probe density grids from aligned trajectory")
        density_parser.add_argument("--probes", "-P",nargs='+', dest="probelist", help="Selection of probenames to calculate density for. If not given, all probes for the solvent will be selected.")
        density_parser.add_argument("-opref", action="store", dest="outprefix", help="Prefix for grids saved inside density folder for each replica. If this is given, automatic energy conversion will not work until you restablish expected names or explicitely give the prefix in energy command. Default: False")
        density_parser.add_argument("--onlycom", help="Density calculations ONLY for center of masses of each co-solvent. Probe list is ignored. Default: False",action="store_true", dest="onlycom", default=False)
        density_parser.add_argument("--com", help="Include center of mass to list of probes probes in density calculation. Default: False",action="store_true", dest="com", default=False)
        density_parser.add_argument("--ref", help="Reference PDB file to build grid shape over which counts will be added. By default, replica reference PDB will be used. This option is useful to build density grids on a subregion of the space or when replica trajectory was also aligned using --ref option.", action="store", dest="ref")

        # ANALYSIS: Cpptraj DENSITY
        cpp_density_parser = self.subparseronreplica(anl_cmd, 'cpp_density', help="Calculate probe density grids from aligned trajectory")
        cpp_density_parser.add_argument("--probes", "-P",nargs='+', dest="probelist", help="Selection of probenames to calculate density for. If not given, all probes for the solvent will be selected.")
        cpp_density_parser.add_argument("-opref", action="store", dest="outprefix", help="Prefix for grids saved inside density folder for each replica. If this is given, automatic energy conversion will not work until you restablish expected names or explicitely give the prefix in energy command. Default: False")
        cpp_density_parser.add_argument("--onlycom", help="Density calculations ONLY for center of masses of each co-solvent. Probe list is ignored. Default: False",action="store_true", dest="onlycom", default=False)
        cpp_density_parser.add_argument("--com", help="Include center of mass to list of probes probes in density calculation. Default: False",action="store_true", dest="com", default=False)
        cpp_density_parser.add_argument("--ref", help="Reference PDB file to build grid shape over which counts will be added. By default, replica reference PDB will be used. This option is useful to build density grids on a subregion of the space or when replica trajectory was also aligned using --ref option.", action="store", dest="ref")
        cpp_density_parser.add_argument("--only-write", action="store_true", default=False, dest="onlywrite", help="Only write ptraj input scripts BUT don't execute them. Useful when manual editing is needed. (default: False)")
        cpp_density_parser.add_argument("--only-exe", action="store_true", default=False, dest="onlyexe", help="Only execute existing ptra scripts, do not overwrite them. If scripts don't exist, this function will fail. (Default: False)")


        # ANALYSIS: RESIDENCE
        residence_parser = self.subparseronreplica(anl_cmd, 'residence', help="Calculate hotspot residence dens from aligned trajectory")
        residence_parser.add_argument("--hpfile", "-hf", action="store", dest="hpfile", help="HotspotSet pickled file from 'analyze hotspots create' action.")
        residence_parser.add_argument("--hpid", "-id", action="store", type=int, dest="hpid", help="Hotspot ID in HPFILE to define region to study.")
        residence_parser.add_argument("--center", "-ce", nargs='+', dest="center", help="Sphere center to define region to study. Give a space separated 3 float list.")
        residence_parser.add_argument("--tol", "-t", default=0.5,type=float,dest="tolerance", help="Tolerance in angstroms around sphere center or hotspot coordinates to consider as occupied space. (default: 0.5)")
        

        # ANALYSYS: ALIGN
    #    anl_commands = anl_parser.add_subparsers(help='Analysis commands', dest='anl_command')
        align_parser = self.subparseronreplica(anl_cmd, 'align', 'Align trajectory of selected replicas', extras=False)
        align_parser.add_argument("-N", help="List production steps to consider for alignment using a colon separated range. Ex: 1:20 - first to 20th step.", default=False, nargs=1,dest="nanoselect")
        align_parser.add_argument("-C", type=int, help="Number of cpus to use for the action. If option not given, will use 1 CPU serial mode.", dest="ncpus", default=1, action="store")
        align_parser.add_argument("--mask", action="store", dest="alignmask", help="Modify alignment mask defined when creating the replicas. By default the macromolecule will be automatically identified. Give a list with comma separated residue numbers or hyphen separated range. E.g. 10-100,120-240.")
        align_parser.add_argument("--ref", action="store", dest="ref", help="Path to reference PDB file. By default, pyMDmix generates one automatic reference pdb file which can be found inside each replica folder. This option will override it.")
        align_parser.add_argument("--only-write", action="store_true", default=False, dest="onlywrite", help="Only write ptraj input scripts BUT don't execute them. Useful when manual editing is needed. (default: False)")
        align_parser.add_argument("--only-exe", action="store_true", default=False, dest="onlyexe", help="Only execute existing ptra scripts, do not overwrite them. If scripts don't exist, this function will fail. (Default: False)")
        
        # In conversion from density to grids, allow also merging of densities between replicas before energy conversion
        energy_parser = self.subparseronreplica(anl_cmd, 'energy', help="Calculate free energy maps from density grids for selected probes and replicas", extras=False)
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

    def action(self, parserargs):
        if parserargs.anl_command != 'hotspots': p = returnMDMixProjectOrFail(parserargs)
        if parserargs.anl_command == 'align':
            replicas = self.fetchReplicaSelection(parserargs, p)
            if replicas:
                nanosel = self.parsenanos(parserargs)
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

        elif parserargs.anl_command == 'cpp_density':
            replicas = self.fetchReplicaSelection(parserargs, p)
            if replicas:
                args = {
                    "nanosel": self.parsenanos(parserargs),
                    "ncpus": parserargs.ncpus,
                    "steps": parserargs.step,
                    "outprefix": parserargs.outprefix,
                    "probelist": parserargs.probelist,
                    "includeCOM": parserargs.com,
                    "onlyCOM": parserargs.onlycom,
                    "ref": parserargs.ref,
                    "run": not parserargs.onlywrite,
                    "write": not parserargs.onlyexe,
                }
                p.calc_cppdensityReplicas(replicas, **args)
                print "DONE"

        elif parserargs.anl_command == 'density':
            replicas = self.fetchReplicaSelection(parserargs, p)
            nanosel = self.parsenanos(parserargs)
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
            replicas = self.fetchReplicaSelection(parserargs, p)
            nanosel = self.parsenanos(parserargs)
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
                center = numpy.array(center,dtype=float)
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
            replicas = self.fetchReplicaSelection(parserargs, p)
            if replicas:
                inprefix=parserargs.inprefix
                outprefix=parserargs.outprefix
                nsnaps = parserargs.nsnaps
                probelist = parserargs.probelist
                dg0 = parserargs.nodg0
                avg = parserargs.noavg
                nanosel = self.parsenanos(parserargs)

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
