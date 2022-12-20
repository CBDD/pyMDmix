import argparse
import pyMDMix
from pyMDMix import MDMixError
from pyMDMix.Commands.Command import Command
from pyMDMix.Projects import returnMDMixProject, returnMDMixProjectOrFail
class Queue(Command):
    def __init__(self):
        self.cmdstring = "queue"

    def create_parser(self, subparsers):
        queue_parser = subparsers.add_parser('queue', help="Queue input files options.")
        queue_parser.add_argument("action", choices=('list','write'), help="LIST: Show installed queue system templates. WRITE: Write input files for all replicas in current project or for REPLICA in current folder.")
        queue_parser.add_argument("-n", action="store", dest="queuename", help="WRITE action: queue system to use. Mandatory.")
       

    def action(self, parserargs):
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

