import argparse
import pyMDMix
from pyMDMix import MDMixError
from pyMDMix.Commands.Command import Command
from pyMDMix.Projects import returnMDMixProject, returnMDMixProjectOrFail

class Create(Command):
    def __init__(self):
        self.cmdstring = "create"

    def create_parser(self, subparsers):
        create_parser = subparsers.add_parser('create', help='Create Project or Solvents')
        create_parser.add_argument("action", choices=("project", "replica", "solvents", "template"), help="Create new project, a new stand-alone replica or add solvent to the database. 'template' will save a project template for new project creation.")
        create_parser.add_argument("-n", action="store", dest='projname', help="Project name or Replica Name. Default: mdmix_project.", default='mdmix_project')
        create_parser.add_argument("-f", action="store", dest="file", help="Configuration file path or output file name in 'template' action. Project can be created empty. Replica will need a config file with MDSETTINGS section.")
        create_parser.add_argument("-top", action="store", dest="top", help="Amber Topology file for stand-alone solvated replica creation.")
        create_parser.add_argument("-crd", action="store", dest="crd", help="Amber CRD file for stand-alone solvated replica creation.")

    def action(self, parserargs):
        import os
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
