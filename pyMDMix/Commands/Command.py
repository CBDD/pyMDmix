import argparse
import sys
from pyMDMix import MDMixError
from pyMDMix.tools import parseNumMask

class Command(object):
    def __init__(self):
        self.cmdstring = ""
        pass
    
    def create_parser(self, subparser):
        pass
    
    def action(self, parserargs):
        pass
    
    def subparseronreplica(self, subprs, name, help, extras=True):
        sp = subprs.add_parser(name, help=help)
        sp.add_argument("mode", choices=('all','bysolvent','byname','group'), action="store", default="all", help="Perform selection of replicas based on solvent name, replica names or groups. If 'all', do action on all replicas.")
        sp.add_argument('-s', action="store", nargs='+', dest="selection", help="Selection list. If selecting 'bysolvent', list of solvent names is expected. If 'byname', list of replica names. If 'group', group name. Skip if 'all' is selected.")
        if extras:
            sp.add_argument("-N", help="For action commands, list production steps to consider for the analysis using a colon separated range. Ex: 1:20 - first to 20th step.", default=False, nargs=1,dest="nanoselect")
            sp.add_argument("--step", help="Take snapshots every STEP number. Default:1.", action="store", default=1, type=int, dest="step")
            sp.add_argument("-C", type=int, help="Number of cpus to use for the action. If option not given, will use 1 CPU serial mode.", dest="ncpus", default=1, action="store")
        return sp

    def fetchReplicaSelection(self, parserargs, project):
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

    def parsenanos(self, argparser):
        if not argparser.nanoselect: return False
        nanosel = parseNumMask(argparser.nanoselect[0])
        print "Selected steps: %s"%(', '.join(map(str, nanosel)))
        return nanosel