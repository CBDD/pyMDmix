#  ------------------------ pyMDMix ----------------------------------- 
#                  http://mdmix.sourceforge.net
#  -------------------------------------------------------------------- 
# 
#  Software for preparation, analysis and quality control
#  of solvent mixtures molecular dynamics.
# 
#  Copyright (C) 2014 daniel
# 
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
# 
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
# 
#  Please cite your use of pyMDMix in published work:
# 
#              TOBEPUBLISHED
# 
#  -------------------------------------------------------------------- 

__author__="daniel"
__date__ ="$Mar 19, 2014 6:34:36 PM$"

import os
import os.path as osp
import logging
import numpy as npy
import matplotlib as mp
import matplotlib.pyplot as plt
import matplotlib.cm as cm

ALL_AMBER_PROPS = ['Etot','EPtot', 'EKtot','TEMP(K)','RESTRAINT','VDWAALS','EELEC']

class PlotError(Exception):
    pass

class Plot(object):
    def __init__(self):
        self.log = logging.getLogger("Plot")
    
    def plotMDAmber(self, replicalist, properties=[], outfilename=None, selectedsteps=[], hideylabels=True, show=True, **kwargs):
        """
        Plot production properties like energy, volume, temperature, restrains... etc. Name in properties should match those of AMBER
        output file names (eg. VOLUME, ETOT, ...)
        """
        import re
        
        if not properties: properties=ALL_AMBER_PROPS
        if not isinstance(replicalist, list): replicalist = [replicalist]
        if not isinstance(properties, list): properties = [properties]
        
        self.log.info("Plotting replicas %s. Properties: %s"%(','.join(map(lambda x: x.name, replicalist)), ','.join(properties)))
        
        # Construct regexps
        regexps = {}
        for prop in properties:
            regexps[prop] = re.compile('%s\s+=\s+(-{0,1}\d+\.{0,1}\d+)'%re.escape(prop))
        
        # Fetch data
        rdata = {}
        for repl in replicalist:
            check = repl.getChecker()
            rdata[repl.name] = dict([(p,[]) for p in properties])
            steps = selectedsteps or range(1, repl.ntrajfiles+1)
            for f in steps:
                outfile = check.getProductionOutputFile(f)
                
                if not outfile:
                    raise PlotError, "Output file for step %i could not be found. Make sure all selected steps are finished: %s"%(f,steps)
                
                # fetch properties from file
                tmpdata = dict([(p,[]) for p in properties])
                for line in outfile.split('\n'):
                    for prop, exp in regexps.iteritems():
                        m = exp.search(line)
                        if m: tmpdata[prop].append(float(m.groups()[0]))
                
                # Remove two last digits from values (Meand and fluctuation values)
                # and extend global results
                for k, v in tmpdata.iteritems():
                    rdata[repl.name][k].extend(v[:-2])
                
        # Organize data by property instead of by replica
        plotdata = {}
        for r, data in rdata.iteritems():
            for prop, vals in data.iteritems():
                if not plotdata.has_key(prop): plotdata[prop] = {}
                plotdata[prop][r] = vals
                
        # Calc shape of the plot
        nplots = len(properties)
        fig, axes = plt.subplots(nplots, 1, sharex=True)
        for i,ax in enumerate(axes):
            prop = plotdata.keys()[i]
            data = plotdata[prop]
            if i == (len(plotdata.keys())-1): ax.set_xlabel('STEP')
            ax.set_ylabel(prop)
            if hideylabels: ax.set_yticklabels([])
            [ax.plot(d, label=r) for r,d in data.iteritems()]
        
        del rdata
        
        plt.legend(prop={'size':6})
        if outfilename: fig.savefig(outfilename)
        if show: plt.show()
                
    
    def plotRMSDReplicas(self, replicaList, outfilename=None, show=True, selectedsteps=[], *args, **kwargs):
        """Plot rmsd evolution for the replica and save the plot to outfilename.  Plot using pyplot. 
        X axis is time and Y axis RMSD. args and kwargs can be used to tune the plot attributes.
        TRAJECTORY MUST BE ALIGNED BEFORE PLOTTING.
        
        Args:
                replicaList     (list of ReplicaInfo)   Replicas to plot
                outfilename     (str)           Filename to store the plot. Use appropriate extension (pdf, png, jpg, etc) as pyplot will try 
                                                to recognize the output format from filename.
                show            (bool)          Plot directly to screen.
                stepselection   (list)          Number of steps to plot.
                        
        Returns:
                Function writes figure to file directly and it returns the Figure instance
                so the user can further modify it.
        """
        plotdata = self.fetchRMSDdata(replicaList, selectedsteps)
        fig = self.plotRMSDdata(plotdata, outfilename=outfilename, show=show, *args, **kwargs)
        return fig
    
    def fetchRMSDdata(self, replicaList, selectedsteps=[]):
        """
        Collect RMSD data for the replicas in *replicaList*. Replicas must be aligned.
        All files ending with *rmsd.out in alignment folder will be collected, joined and stored in a dictionary
        that can be plotted with plotRMSDdata.
        
        :arg list replicaList: List of Replica instances to plot
        
        :return: Dictionary with replica name and pairs x,y values for BB and HA data.
        """
        if not isinstance(replicaList, list): replicaList = [replicaList]
        import glob
        import re

        digit = re.compile('(\d+)')
        plotdata = {}
        for i, replica in enumerate(replicaList):
            # Check if replica is centered and fetch rmsd output files
            if replica.isAligned(selectedsteps):
                plotdata[replica.name] = {}
                folder = osp.join(replica.path,replica.alignfolder)
                
                # Work on BB RMSDs
                fnames = "*_bb_rmsd.out"
                filelist = glob.glob(folder+os.sep+fnames)
                if filelist: 
                    filelist.sort(key=lambda x: int(digit.search(os.path.basename(x)).groups()[0]))
                    self.log.debug("Sorted file list to join: %s"%filelist)
                    data = npy.hstack([npy.loadtxt(f)[:,1] for f in filelist])
                    data = npy.vstack((range(len(data)),data)).T

                    # Save data to disk inside centering folder
                    npy.savetxt(folder+os.sep+'all_bb_rmsd.dat', data, fmt='%.4f', delimiter='\t')
                    self.log.info("Joined all independent BB rmsd files to a unique all_bb_rmsd.dat file for replica %s"%replica.name)

                    # Plot
                    x = data[:,0]/(replica.prod_steps/float(replica.trajfrequency))
                    y = data[:,1]
                    plotdata[replica.name]['BB'] = [x,y]
                else:
                    self.log.error("No backbone rmsd file found with *_bb_rmsd.out termination in alignment folder.")
                # Do same for HA plot
                fnames = "*_ha_rmsd.out"
                filelist = glob.glob(folder+os.sep+fnames)
                if filelist:
                    filelist.sort(key=lambda x: int(digit.search(os.path.basename(x)).groups()[0]))
                    self.log.debug("Sorted file list to join: %s"%filelist)
                    data = npy.hstack([npy.loadtxt(f)[:,1] for f in filelist])
                    data = npy.vstack((range(len(data)),data)).T

                    # Save data to disk inside centering folder
                    npy.savetxt(folder+os.sep+'all_ha_rmsd.dat', data, fmt='%.4f', delimiter='\t')
                    self.log.info("Joined all independent HA rmsd files to a unique all_ha_rmsd.dat file for replica %s"%replica.name)

                    # Plot                    x = data[:,0]/(replica.prod_steps/float(replica.trajfrequency))
                    y = data[:,1]
                    plotdata[replica.name]['HA'] = [x,y]
                else:
                    self.log.error("No heavy atoms rmsd file found with *_ha_rmsd.out termination in alignment folder.")
            else:
                self.log.warn("Replica %s does not have trajectory aligned. Skipping plotting."%replica.name)

        return plotdata
    
    def plotRMSDdata(self, replicarmsd, outfilename=None, show=False, *args, **kwargs):
        """
        Plot data obtained with fetchRMSDdata.
        
        :arg dict replicarmsd: Dictionary containing for each replica name a pair BB and HA entries with the RMSD data. Format: {'replica':{'BB':[x, y], 'HA':[x,y]},}
        :arg str outfilename: If given, save figure with this file name. Make sure you use a matplotlib compatible suffix (e.g. png pdf jpg)
        :arg bool show: Show figure interactively?
        
        :return: Matplotlib figure
        """
        mp.rcParams['lines.linewidth'] = 0.0
        mp.rcParams['axes.linewidth'] = 0.5
        colorspace = cm.rainbow(npy.linspace(0,1,len(replicarmsd.keys())))
        fig, axes = plt.subplots(2,1,sharex=True)
        
        # Plot replicas BB and HA
        i = 0
        for replica, rmsdata in replicarmsd.iteritems():
            x,y = rmsdata['BB']
            axes[0].plot(x, y, color=colorspace[i], label=replica, linestyle='-', linewidth='1.5', *args, **kwargs)
            x,y = rmsdata['HA']
            axes[1].plot(x, y, color=colorspace[i], label=replica, linestyle='-', linewidth='1.5', *args, **kwargs)
            i += 1
        
        # Decorate
        axes[0].set_title("BB RMSD plot")
        axes[0].set_ylabel('RMSD ($\AA^2$)')
        axes[1].set_title("HA RMSD plot")
        axes[1].set_xlabel('Time (ns)')
        axes[1].set_ylabel('RMSD ($\AA^2$)')
                
        plt.legend(prop={'size':'small'}, loc="upper center", ncol=min(len(replicarmsd.keys()), 3))
        if outfilename: 
            self.log.info("Saving RMSD plot to file %s"%outfilename)
            fig.savefig(outfilename, *args, **kwargs)
        if show: plt.show()
        return fig
    
    def plotResidenceResults(self, results, outfilename=None, show=False, colormap={}, *args, **kwargs):
        """
        Plot occupancy action results and save the plot to outfilename.  Plot using pyplot. 
        One color per residuename. X axis is time and Y axis residue ID. args and kwargs can be used to set up the plot attributes.
        
        :arg dict|str results:  If dict is given, its assumed that it's the resulting dict from Occupancy Action. If its a string. It's assumed it is a textfile with the result printed from Occupancy Action.
        :arg str outfilename:   Filename to store the plot. Use appropriate extension (pdf, png, jpg, etc) as matplotlib will try to recognize the output format from filename.
        :arg bool show:         Show plot in screen? pyplot.show command.

        :Returns: Figure instance so the user can further modify it.
        """
        self.log.info("Plotting occupancy results...")
        from collections import Counter
        
        mp.rcParams['lines.linewidth'] = 0.0
        mp.rcParams['axes.linewidth'] = 0.5
                
        # If results is string, process filename to obtain same results dict format as Occupancy Action output.
        if isinstance(results, str):
            results = self.getDictFromOccupancyFile(results)
        elif not isinstance(results, dict):
            raise AttributeError, "results should be a string with a valid filename containing occupancy results or a dictionary with the results from occupancy action."
        
        # One color per residue name
        names = results['map'].keys()
        if 'NO_OCCUPANCY' in names: names.remove('NO_OCCUPANCY')
        if not colormap:
            colorspace = cm.rainbow(npy.linspace(0,1,len(names)))
            colormap = dict(zip(names, colorspace))
        
        # Build a reverse map ID to Name
        idToName = {}
        for n, ids in results['map'].iteritems():
            for i in ids: idToName[i] = n
        
        # Finally re-order data by name to create different series
        # Each name will contain a list (resid, frame)
        data = {}
        frames = results.keys()
        frames.remove('map')
        maxframe = npy.array(frames, dtype=int).max()
        maxres = 0
        for frame in sorted(results.keys()):
            if frame == 'map': continue
#            if frame > maxframe: maxframe = int(frame)
            for idval in results[frame]:
                if idval > maxres: maxres = int(idval)
                name = idToName.get(idval)
                if not data.has_key(name): data[name] = []
                data[name].append((idval,frame))
        
#        labels = map(idToName.get, data.keys())
#        colors = map(colormap.get, labels)
#        coloridmap = dict(zip(data.keys(),colors))
                
        # Start a figure
        fig = plt.figure(*args, **kwargs)
        plot = fig.add_subplot(111)
        fig.suptitle("Occupancy plot")
        plot.set_ylim(-10, maxres+10)
        plot.set_xlim(0, maxframe+1)
        plot.set_xlabel('Frame number')
        plot.set_ylabel('Residue ID')
        
        # Work on each series (residues)
        for name, resFrame in data.iteritems():
            if name == 'NO_RESIDENCE': continue
            y,x = npy.array(resFrame).T
            plot.scatter(x, y, s=50, lw=0.0, alpha=0.5, facecolor=colormap[name], 
                                        label=name, *args, **kwargs)
                                        
            # Add text for 2 most famous residues (those occupying a lot the site)
            c = Counter(y.tolist()).most_common(2)
            [plot.text(-200, famous[0], str(famous[0]), verticalalignment='center', horizontalalignment='left',
                    bbox=dict(facecolor=colormap[name], alpha=0.5), fontsize=8) for famous in c]

        # Draw legend outside plots
#        box = plot.get_position()
#        plot.set_position([box.x0, box.y0 + box.height * 0.1,
#                 box.width, box.height * 0.9])
#        
        plot.legend(scatterpoints=1,
                    loc='center', bbox_to_anchor=(0.5, 1.03),
                    ncol=4, fancybox=False, shadow=True, 
                    prop={'size':6})
        if outfilename: 
            fig.savefig(outfilename, *args, **kwargs)
            self.log.info("Saved residence plot: %s"%os.path.abspath(outfilename))
        if show: plt.show()
                
        return fig


import Biskit.test as BT
import tools as T

class Test(BT.BiskitTest):
    """Test"""
    def test_plotRMSD(self):
        """RMSD data plotting with matplotlib"""
        dummydata = {'repl1':{'BB':[[0,1,2,3],[1,2,3,4]], 'HA':[[0,1,2,3],[1.4,5.6,7.8,9]]}}
        plotter = Plot()
        fig = plotter.plotRMSDdata(dummydata)
        self.assertTrue(fig)
        

if __name__ == "__main__":
    print "Hello World"
