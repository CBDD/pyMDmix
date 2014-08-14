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
__date__ ="$12-mar-2014 16:45:27$"

import logging
import os
import os.path as osp

import settings as S
import tools as T

class Align(object):
    def __init__(self, replica, steps=[], nthreads=False, warn=True, waitend=True, **kwargs):
        """
        Constructor
        
        :arg replica: Replica to work with
        :type replica: :class:`~Replicas.Replica`
        :arg list steps: Selected steps to align. If empty, align all expected steps.
        :arg int nthreads: Number of processors to use.
        :arg bool warn: Print warning messages.
        :arg bool waitend: When runing commands, wait until all steps are complete before exiting
        
        :kwargs: All parameters for :meth:`writePtrajInput` can be modified with keyword arguments.
        """
        self.replica = replica
        self.log = logging.getLogger('Align')
        self.steps = steps or range(1, replica.ntrajfiles+1)
        self.nthreads = nthreads
        self.warn = warn
        self.waitend=waitend
        
        # Register to executor
        # change number of threads if number given (use EXECUTOR nthreads if False)
        if nthreads: 
            self.log.debug("Changing nthreads in Executor to %d"%nthreads)
            T.EXECUTOR.changeNthreads(nthreads)
        T.EXECUTOR.start()

        self.writePtrajInput(**kwargs)
        self.run()
        
        # Exit executor
        if waitend: T.EXECUTOR.waitJobCompletion()
        T.EXECUTOR.terminate()
        
    def __aligncmd(self, step):
        "Return ptraj execution command for step *step*"
        if not self.replica.checkProductionExtension([step])[step]:
            return False
        inf = self.replica.mdoutfiletemplate.format(step=step, extension='ptraj')
        outf= inf.replace('.ptraj','_ptraj.log')
        top = os.pardir+os.sep+self.replica.top
        cmd = S.AMBEREXE+os.sep+'ptraj {top} < {inf} >& {outf}'.format(top=top, inf=inf, outf=outf)
        path = osp.join(self.replica.path, self.replica.alignfolder)
        return cmd, path
        
    def run(self):
        "Run centering process for steps selected at instantiation."
        self.log.info("Running alignment of trajectory in replica %s"%self.replica.name)
        for i in self.steps:
            cmdpath = self.__aligncmd(i)
            if cmdpath: 
                cmd, path = cmdpath
                T.EXECUTOR.submitCmd(cmd=cmd, path=path)
            else:
                if self.warn: self.log.warn("Production step %i not found")

    def writePtrajInput(self, writermsd=True, writeavgpdb=True, alignmask=False, **kwargs):
        """
        Write ptraj input scripts for aligning the selected md production steps

        :arg list steps: list of steps to write input for. Should be a list of integers.
        """
        from Amber import AmberWriter
        outrmsd = False
        outpdbavg = False
        
        self.replica.go()
        if not osp.exists(self.replica.alignfolder): os.mkdir(self.replica.alignfolder)

        awriter = AmberWriter(self.replica)

        # Expected extension names in production folder
        exts = self.replica.checkProductionExtension(self.steps)

        # Use user alignmask?
        if alignmask: self.log.info("Using user defined alignment mask: %s"%alignmask)

        # Will consider mdpath is in previous dir
        prev = os.pardir+os.sep
        mdpath = prev+self.replica.mdfolder
        scripts=[]

        self.log.debug("Writing input for replica %s steps %s writermsd %s writeavgpdb %s"%(self.replica.name, self.steps, writermsd, writeavgpdb))

        for i in self.steps:
            # For each step, generate alignment script and write file
            ext = exts[i]
            if not ext:
                if self.warn: self.log.warn("Production trajectory file for step %i not found. Skipping alignment scripts writting for this step."%i)
                continue
            n = self.replica.mdoutfiletemplate.format(step=i, extension=ext)
            trajin = mdpath+os.sep+n
            trajout=n
            if ext in ('nc','netcdf','ncdf'): trajout+=' netcdf'
            elif ext in ('dcd'): trajout+=' dcd'

            # Protavg and rmsdout names
            if writermsd: outrmsd = 'md%i'%i
            if writeavgpdb: avgpdb = 'prot_avg_%i.pdb'%i

            # Create script
            script = awriter.getPtrajAlignScript(trajin=trajin, trajout=trajout, writermsd=writermsd,
                                    writeavgpdb=writeavgpdb, rmsdoutprefix=outrmsd, avgpdbout=avgpdb, alignmask=alignmask)

            # Write ptraj script
            oname = osp.join(self.replica.alignfolder, n.replace(ext, 'ptraj'))
            self.log.debug("Writing ptraj alignment script %s"%oname)
            scripts.append(oname)
            out = open(oname, 'w')
            out.write(script)
            out.close()

        T.BROWSER.goback()
        return scripts



if __name__ == "__main__":
    print "Hello World"
