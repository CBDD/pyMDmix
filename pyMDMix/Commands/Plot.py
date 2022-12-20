import argparse
import pyMDMix
from pyMDMix.Commands.Command import Command
from pyMDMix.Projects import returnMDMixProject, returnMDMixProjectOrFail
class Plot(Command):
    def __init__(self):
        self.cmdstring = "plot"

    def create_parser(self, subparsers):
        plot_parser = subparsers.add_parser('plot', help="Plotting command")
        plot_cmd = plot_parser.add_subparsers(help='Plotting options', dest='plot_command')
        amber_plot = self.subparseronreplica(plot_cmd,'ambermd', help="Plot Amber MD properties",extras=False)
        amber_plot.add_argument("-o","--out",action="store", dest="outname", help="Name of the output file. Extension should be .png, .pdf, .jpeg, .ps or .eps. Default: ambermdplot.pdf", default="ambermdplot.pdf")
        amber_plot.add_argument("-N", help="List production steps to be plotted using a colon separated range. Ex: 1:20 - first to 20th step.", default=False, nargs=1,dest="nanoselect")
        rmsd_plot = self.subparseronreplica(plot_cmd,'rmsd', help="Plot Backbone and Heavy atoms RMSD plot",extras=False)
        rmsd_plot.add_argument("-o","--out",action="store", dest="outname", help="Name of the output file. Extension should be .png, .pdf, .jpeg, .ps or .eps. Default: rmsdplot.pdf", default="rmsdplot.pdf")
        rmsd_plot.add_argument("-N", help="List production steps to be plotted using a colon separated range. Ex: 1:20 - first to 20th step.", default=False, nargs=1,dest="nanoselect")

    def action(self, parserargs):
        p = returnMDMixProjectOrFail(parserargs)
        import pyMDMix.Plotter as P
        plot = P.Plot()
        replicas = self.fetchReplicaSelection(parserargs, p)
        stepselection = self.parsenanos(parserargs)
        if parserargs.plot_command == 'ambermd':
            if replicas: plot.plotMDAmber(replicas, outfilename=parserargs.outname, selectedsteps=stepselection)
        elif parserargs.plot_command == 'rmsd':
            if replicas: plot.plotRMSDReplicas(replicas, outfilename=parserargs.outname, selectedsteps=stepselection)
                
