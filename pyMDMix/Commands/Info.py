import argparse
import pyMDMix
from pyMDMix.Commands.Command import Command
from pyMDMix.Projects import returnMDMixProject, returnMDMixProjectOrFail
class Info(Command):
    def __init__(self):
        self.cmdstring = "info"

    def create_parser(self, subparsers):
        info_parser = subparsers.add_parser('info', help='Print information about the project or solventDB')
        info_parser.add_argument("infoselection", choices=('project', 'solvents',), help="Print info about the project or about available solvents")

    def action(self, parserargs):
        if parserargs.infoselection == 'project':
            p = returnMDMixProjectOrFail(parserargs)
            print(p.longdesc())
        elif parserargs.infoselection == 'solvents':#Actions on solvent database
            man = pyMDMix.Solvents.SolventManager()
            print(man)