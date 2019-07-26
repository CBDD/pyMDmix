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

"""
Queue input writer.
"""
__author__="dalvarez"
__date__ ="$16-mar-2014 22:54:03$"

import re
import os
import os.path as osp
import glob
import settings as S
import tools as T

class QueueInputError(Exception):
    pass

class QueueInputWriter(object):
    """
    Create Queue Input folloring template file stored in package directory or user mdmix directory.
    A template file should exist in either the package templates directory or the user mdmix directory
    with name QUEUENAME_queue_temp.txt and following expected standards.
    """
    def __init__(self, queue, replica=None, templatefile=None, **kwargs):
        """
        Use template with name *queuename*_queue_temp.txt from user's mdmix directory or from pyMDMix package
        template folder. Create input files for Replica *replica*.
        """
        self.replica = replica
        self.queue = queue
        self.templfile = templatefile

        if not self.templfile:
            # Find template file in users directory or package templates
            # User directory has priority
            usertemplates = glob.glob(S.USER_MDMIX_HOME+os.sep+'*_queue_temp.txt')
            packagetemplates = glob.glob(T.templatesRoot()+os.sep+'*_queue_temp.txt')
            qname = re.compile('(\w+)_queue_temp.txt')
            for f in usertemplates:
                n = qname.match(osp.basename(f))
                if n:
                    n = n.groups()[0]
                    if n == self.queue: self.templfile = f

            # If still no match, continue with package folder
            if not self.templfile:
                for f in packagetemplates:
                    n = qname.match(osp.basename(f))
                    if n:
                        n = n.groups()[0]
                        if n == self.queue: self.templfile = f

            # If still no template found, raise QueueInputError
            if not self.templfile: raise QueueInputError, "No template file for queue name %s"%self.queue

        self.template = open(self.templfile,'r').read()

    def __writeTemplate(self, fname, jobname, precommands, commands, postcommands,
                        next, queuename, last=False):
        "Write formated tempalte to file fname"
        if last:
            # Remove line containing {next}
            int = self.template.split('\n')
            templ = []
            for l in int:
                if not '{next}' in l: templ.append(l)
            templ = '\n'.join(templ)
        else:
            # Just use normal template
            templ = self.template
        open(fname,'w').write(templ.format(jobname=jobname,
                                    precommands=precommands, commands=commands,
                                    postcommands=postcommands, next=next, queuename=queuename))

    def write(self, replica=None, queuename='', **kwargs):
        "Write queue input files for replica *replica*"
        if not replica: replica=self.replica
        if not replica: raise QueueInputError, "Replica not set"

        # CWD
        cwd = T.BROWSER.getcwd()

        # Get writter
        if replica.mdProgram == 'AMBER':
            from Amber import AmberWriter
            writer = AmberWriter(replica)
        elif replica.mdProgram == 'NAMD':
            from NAMD import NAMDWriter
            writer = NAMDWriter(replica)
        elif replica.mdProgram == 'OPENMM':
            from OpenMM import OpenMMWriter
            writer = OpenMMWriter(replica)
        else:
            raise QueueInputError, "Replica has unknown mdProgram attribute: %s"%replica.mdProgram

        # start writing files
        # MINIMIZATION
        jobname = replica.name+'_min'
        precommands = ''
        commands = writer.getCommand('min')
        postcommands = 'cd ..'+os.sep+replica.eqfolder
        next = 'eq.q'
        replica.go()
        T.BROWSER.chdir(replica.minfolder)
        self.__writeTemplate('min.q', jobname, precommands, commands, postcommands, 
                                next, queuename=queuename)

        # EQUILIBRATION
        jobname = replica.name+'_eq'
        precommands = ''
        if replica.mdProgram in {'AMBER', 'OPENMM'}: neqsteps = 5
        else: neqsteps = 2
        commands = '\n'.join([writer.getCommand('eq',i) for i in range(1,neqsteps+1)])
        postcommands = 'cd ..'+os.sep+replica.mdfolder
        next = 'md1.q'
        replica.go()
        T.BROWSER.chdir(replica.eqfolder)
        self.__writeTemplate('eq.q', jobname, precommands, commands, postcommands,
                                next, queuename=queuename)

        # PRODUCTION STEPS
        replica.go()
        T.BROWSER.chdir(replica.mdfolder)
        maxstep = replica.ntrajfiles
        for step in range(1, maxstep+1):
            jobname = replica.name+'_md%i'%step
            precommands = ''
            commands = writer.getCommand('md',step)
            postcommands = ''
            next = 'md%i.q'%(step+1)
            last = step == maxstep
            self.__writeTemplate('md%i.q'%step, jobname, precommands, commands, 
                                postcommands, next, queuename=queuename,last=last)

        T.BROWSER.chdir(cwd)
        return True

### FUNCTIONS
def listQueueSystems():
    usertemplates = glob.glob(S.USER_MDMIX_HOME+os.sep+'*_queue_temp.txt')
    packagetemplates = glob.glob(T.templatesRoot()+os.sep+'*_queue_temp.txt')
    qname = re.compile('(\w+)_queue_temp.txt')
    qsys = []
    usermatches = [qname.match(osp.basename(f)) for f in usertemplates]
    [qsys.append(m.groups()[0]) for m in usermatches if m]
    packagemathes = [qname.match(osp.basename(f)) for f in packagetemplates]
    [qsys.append(m.groups()[0]) for m in packagemathes if m]
    return list(set(qsys))

if __name__ == "__main__":
    print "Hello World"
