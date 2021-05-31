import argparse
import pyMDMix
from pyMDMix import MDMixError
from pyMDMix.Commands.Command import Command
from pyMDMix.Projects import returnMDMixProject, returnMDMixProjectOrFail
class Remove(Command):
    def __init__(self):
        self.cmdstring = "remove"

    def create_parser(self, subparsers):
        remove_parser = subparsers.add_parser('remove', help="Remove groups from project. To remove systems or replicas, simply remove the folders or system files from the project folder.")
        remove_parser.add_argument("action", choices=('group'), help="Remove group.")
        remove_parser.add_argument("-gn", action="store", dest="groupname", help="Name of the group to remove. Groupname must exists in project.", required=True)

    def action(self, parserargs):
        p = returnMDMixProjectOrFail(parserargs)
        if parserargs.action == 'group':
            name = parserargs.groupname
            r = p.removeGroup(name)
            if not r:
                raise MDMixError, "Group name %s does not exists. Project groups: %s"%p.listGroups()
            print "DONE"
