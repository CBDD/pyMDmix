import argparse
import pyMDMix
from pyMDMix import MDMixError
from pyMDMix.Commands.Command import Command
from pyMDMix.Projects import returnMDMixProject, returnMDMixProjectOrFail

class Add(Command):
    def __init__(self):
        self.cmdstring = "add"

    def create_parser(self, subparsers):
        add_parser = subparsers.add_parser('add', help="Add new replicas, systems or create group of replicas in an existing project.")
        add_parser.add_argument("action", choices=('system', 'replicas', 'group'), help="Create new system, new replicas or create groups of replicas.")
        add_parser.add_argument("-f", action="store", dest="file", help="REPLICAS action & SYSTEM action: Configuration files. MANDATORY if action is 'replicas' or 'system'.")
        add_parser.add_argument("-gn", action="store", dest="groupname", help="GROUP action: Name of the group of replicas to be created.")
        add_parser.add_argument("-s", action="store", dest="selection", help="GROUP action: List replicanames to add to the group.", nargs='+')
        add_parser.add_argument("-sys", action="store", dest="sysname", help="REPLICAS action: When creating new replicas, specify the system name which should be prepared. System name should exist in the project. If not given, only one system should be present in the project.")

    def action(self, parserargs):
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
