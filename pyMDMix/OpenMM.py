##
## pyMDMix --- http://mdmix.sourceforge.net
## Software for preparation, analysis and quality control
## of solvent mixtures molecular dynamics
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 3 of the
## License, or any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
## General Public License for more details.
##
## You find a copy of the GNU General Public License in the file
## license.txt along with this program; if not, write to the Free
## Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
##
## Please cite your use of pyMDMix in published work:
##
##    TOBEPUBLISHED
##

__author__="daniel"
__date__ ="$Feb 10, 2014 11:36:41 AM$"

#import sys
import os
import os.path as osp
import logging
import string
import numpy as npy
#import MDMix.Config as config
#from MDMix.PDB import PDBManager

import tools as T
import settings as S

class OpenMMWriterError(Exception):
    pass

class OpenMMCheckError(Exception):
    pass

class OpenMMWriter(object):
    def __init__(self, replica=False):
        self.log = logging.getLogger("OpenMMWriter")
        self.replica = replica

        # Load template inputs
        self.loadConfig()

    def loadConfig(self, min=None, eq1=None, eq2=None, mdNPT=None, mdNVT=None, restr=None):
        """
        Load standard template files or user given template files for OpenMM MD configuration.
        All template files should contain expected fieldnames. Files are read and loaded for later substitution.

        :arg str min: Filepath for minimization template file.
        :arg str eq1: Filepath for equilibration first step template file.
        :arg str eq2: Filepath for equilibration second step template file.
        :arg str mdNVT: Filepath for production with NVT ensemble template file.
        :arg str mdNPT: Filepath for production with NPT ensemble template file.
        """
        if min: self.minT = string.Template(open(min,'r').read())
        else:   self.minT = string.Template(open(T.templatesRoot('openmm_min_temp.txt'), 'r').read())

        if eq1: self.eq1T = string.Template(open(eq1, 'r').read())
        else:   self.eq1T = string.Template(open(T.templatesRoot('openmm_eq_temp.txt'), 'r').read())

        if eq2: self.eq2T = string.Template(open(eq2, 'r').read())
        else:   self.eq2T = string.Template(open(T.templatesRoot('openmm_eq_temp.txt'), 'r').read())

        if mdNPT:   self.cpmd = string.Template(open(mdNPT, 'r').read())
        else:       self.cpmd = string.Template(open(T.templatesRoot('openmm_eq_npt_temp.txt'), 'r').read())


        if mdNVT:   self.cvmd = string.Template(open(mdNVT, 'r').read())
        else:       self.cvmd = string.Template(open(T.templatesRoot('openmm_md_temp.txt'), 'r').read())

    def getBoxFromCRD(self, crd):
        "Read box size from CRD file bottom line"
        boxline = open(crd, 'r').readlines()[-1]
        return npy.array( boxline.split()[0:3] , dtype='float32')


    def getCommand(self, process, step=False, replica=None):
        """
        Return a command string to execute for the given process and number of step.
        It takes into account if it needs restraints and the output file formats.

        ::
            getCommand('min')   # Return minimization execution command for current replica
            getCommand('eq',1)  # Get first step equilibration execution command
            getCommand('md',4)  # Get 4th step production execution command.

        :arg str process: Process for which to return the exe command
        :type process: string that should be: **min** for minimization or **eq** for equilibration or **md** for production

        :returns: string with execution command.
        """
        replica = replica or self.replica
        if not replica: raise OpenMMWriterError, "Replica not assigned."

        prevsep = os.pardir+os.sep
        top = osp.basename(replica.top)
        crd = osp.basename(replica.crd)
        extension = 'nc'

        command = False

        # Set reference file if restrained simulation
        if replica.hasRestraints:
            if replica.minimizationAsRef: ref = prevsep+replica.minfolder+os.sep+'min.rst'
            else: ref = prevsep+crd
        else:
            ref = False

        if process == 'min':
            command = S.OPENMM_EXE+' min.py %s %s min.rst '%(prevsep+top, prevsep+crd)
            return command

        elif process == 'eq':
            if not step: return False

            if step == 1:
                #First step
                eqfname = replica.eqoutfiletemplate.format(step=step,extension='')
                command = S.OPENMM_EXE+' eq1.py %s %smin.rst %srst '%(prevsep+top, prevsep+replica.minfolder+os.sep, eqfname)

            elif step > 1:
                eqfname = replica.eqoutfiletemplate.format(step=step,extension='')
                preveqfname = replica.eqoutfiletemplate.format(step=step-1,extension='')

                command = S.OPENMM_EXE+' eq%i.py %s %srst %srst '%(step, prevsep+top, preveqfname, eqfname)

            return command

        elif process == 'md':
            if not step: return False

            mdouttemplate=replica.mdoutfiletemplate.replace('.{extension}','')
            if step == 1:
                fname = mdouttemplate.format(step=1)
                command = S.OPENMM_EXE+' md.py %s %seq5.rst %s.rst %s.nc %s.log'%(prevsep+top, prevsep+replica.eqfolder+os.sep, fname, fname, fname)

                command = command.format(fname=fname)

            elif step > 1:
                prevfname=mdouttemplate.format(step=step-1)
                nextfname=mdouttemplate.format(step=step)
                command = S.OPENMM_EXE+' md.py %s %s.rst %s.rst %s.nc %s.log'%(prevsep+top, prevfname, nextfname, nextfname, nextfname)
                fname = nextfname

                command = command.format(nextfname=nextfname, prevfname=prevfname)

            return command

        else: pass

    def getReplicaCommands(self, replica=None):
        """
        get a string of commands to execute for running the MD for *replica*.
        It will contain the expected file names for input/output files and directory chamnge commands.

        :args replica: Replica to write execution commands for.
        :type replica: :class:`~Replicas.Replica`
        """
        replica = replica or self.replica
        if not replica: raise OpenMMWriterError, "Replica not assigned."

        # Set variables
        outcommands = []

        # MINIMIZATION
        outcommands.append('cd %s'%replica.minfolder)
        outcommands.append(self.getCommand('min'))

        # EQUILIBRATION
        outcommands.append('cd %s'%osp.join(os.pardir,replica.eqfolder))
        [outcommands.append(self.getCommand('eq',i)) for i in range(1,6)]

        # PRODUCTION
        outcommands.append('cd %s'%osp.join(os.pardir,replica.mdfolder))
        [outcommands.append(self.getCommand('md',i)) for i in range(1, replica.ntrajfiles+1)]

        return outcommands

    def writeCommands(self, replica=False, outfile='COMMANDS.sh'):
        "Write list of commands to run the MD into an output file."
        replica = replica or self.replica
        if not replica: raise OpenMMWriterError, "Replica not assigned."
        commands = self.getReplicaCommands(replica)
        open(outfile,'w').write('\n'.join(commands))
        ok = osp.exists(outfile)
        if ok: self.log.debug("Wrote commands file %s"%outfile)
        else: self.log.warn("COMMANDS.sh not writen!")
        return ok

    def writeMinInput(self, replica=False):
        """Write minimization input file for OpenMM in replica.minfolder
         :args replica: Replica instance. If False, use replica assigned in instantiation.
         :type replica: :class:`~Replicas.Replica`
        """
        replica = replica or self.replica
        if not replica: raise OpenMMWriterError, "Replica not assigned."
        
        T.BROWSER.gotoReplica(replica)
        
        restr = ''
        if replica.hasRestraints:
            if not replica.minimizationAsRef: restr = self.restr
            else: 
                self.log.warn('Use of Minimized structure as restraint reference is still not possible with OpenMM. Will use starting PRMCRD.')
                restr = self.restr
        
        formatdict = {'top':replica.top, 'crd':replica.crd, 'restraints':restr, 
                        'box':self.getBoxFromCRD(replica.crd).max(), 'timestep':int(replica.md_timestep),
                        'freq':replica.trajfrequency}
        formatdict['minsteps'] = replica.minsteps
        out = replica.minfolder+os.sep+'min.py'
        open(out,'w').write(self.minT.substitute(formatdict))
        exists = osp.exists(out)
        T.BROWSER.goback()
        
        return exists

    def writeEqInput(self, replica=False):
        """Write equilibration input file for OpenMM in replica.eqfolder
         Args:
            replica     (ReplicaInfo)   Replica instance
        """
        replica = replica or self.replica
        if not replica: raise OpenMMWriterError, "Replica not assigned."
        
        restr = ''
        if replica.hasRestraints: restr = self.restr
            
        T.BROWSER.gotoReplica(replica)
        formatdict = {'top':replica.top, 'crd':replica.crd, 'restraints':restr, 
                        'timestep':replica.md_timestep, 'freq':replica.trajfrequency}
        # EQUILIBRATION
        # It comprises 2 steps:
        # - first step (NV) = 500000 steps (1ns) for heating from 100K up to 300K
        # - second step (NPT) = 100000 steps (2ns) of constant pressure constant temp at max temp

        # FIRST STEP
        # Heating up the system from 100 to Final Temperature during 1ns
        formatdict['temp'] = replica.temp
        formatdict['eqinput'] = os.path.join(os.pardir, replica.minfolder, 'min')
        formatdict['eqoutput'] = 'eq1'
        formatdict['first_step'] = 0
        formatdict['final_step'] = replica.namd_heating_steps # 1ns (timestep=0.002ps)
        eq1out = replica.eqfolder+os.sep+'eq1.py'
        open(eq1out,'w').write(self.eq1T.substitute(formatdict))

        formatdict['eqinput'] = formatdict['eqoutput']
        formatdict['first_step'] = formatdict['final_step']

        # SECOND STEP
        # NPT equilibration for 1ns
        formatdict['eqoutput'] = 'eq2'
        formatdict['final_step'] = formatdict['first_step'] + replica.npt_eq_steps
        eq2out = replica.eqfolder+os.sep+'eq2.py'
        open(eq2out,'w').write(self.mdNPT.substitute(formatdict))
        exists = osp.exists(eq1out) and osp.exists(eq2out)
        
        T.BROWSER.goback()
        return exists

    def writeReplicaInput(self, replica=False):
        replica = replica or self.replica
        if not replica: raise AmberWriterError, "Replica not assigned."

        self.log.info("Writing AMBER simulation input files for replica %s ..."%(replica.name))
        cwd = T.BROWSER.cwd
        T.BROWSER.gotoReplica(replica)

        if not (osp.exists(replica.top) and osp.exists(replica.crd)): # and osp.exists(replica.pdb)):
            raise AmberWriterError, "Replica top or crd files not found in current folder: %s, %s"%(replica.top, replica.crd)

        substDict = {}

        # Substitute frequency, ioutfm flags
        substDict['ioutfm'] = replica.mdnetcdf
        substDict['iwrap'] = replica.iwrap
        substDict['freq'] = '%i'%replica.trajfrequency
        substDict['timestep'] = '%.3f'%(replica.md_timestep/1000.)

        # check iwrap 
        # if 1 and replica has restraints write check_com.sh
        iwrap = replica.iwrap == '1'

        mfield = ''

        # Write minimization input
        substDict['minsteps'] = replica.minsteps
        # only add restraining field if we want the starting structure to be restrained
        # otherwise, the minimization will have no restraints and will use this output for 
        # future restrains
        if replica.minimizationAsRef: substDict['maskfield'] = ''
        else: substDict['maskfield'] = mfield
        outf = replica.minfolder+os.sep+'min.py'
        self.log.debug("Writing: %s"%outf)
        open(outf,'w').write(self.minT.substitute(substDict))

        # Write equilibration input
        # System warming from 100 to Project Temperature
        # Increase of 50K every 200ps
        # Finally run 1ns at NPT to equilibrate density
        substDict['maskfield'] = mfield
        substDict['nsteps'] = replica.heating_steps
        substDict['tempi'] = replica.parm_heating_tempi
        substDict['tempf'] = substDict['tempi'] + 50
        outf = replica.eqfolder+os.sep+'eq1.py'
        self.log.debug("Writing: %s"%outf)
        open(outf,'w').write(self.eq1T.substitute(substDict))
        substDict['tempi'] = substDict['tempf']

        for i in range(2,5):
            substDict['tempf'] = substDict['tempi'] + 50
            outf=replica.eqfolder+os.sep+'eq%i.py'%i
            self.log.debug("Writing: %s"%outf)
            open(outf,'w').write(self.eq2T.substitute(substDict))
            substDict['tempi'] = substDict['tempf']

        substDict['tempf'] = replica.temp
        substDict['nsteps'] = replica.npt_eq_steps# 1ns of equilibration under NPT
        substDict['pres0'] = replica.npt_pressure   # Pressure for NPT simulation, default = 1 bar
        outf = replica.eqfolder+os.sep+'eq5.py'
        self.log.debug("Writing: %s"%outf)
        open(outf,'w').write(self.cpmd.substitute(substDict))

        # Write md input, 1ns each file and run under NVT conditions
        substDict['nsteps'] = replica.prod_steps # 1ns each file
        outf=replica.mdfolder+os.sep+'md.py'
        self.log.debug("Writing: %s"%outf)
        if replica.production_ensemble == 'NPT': prodfile = self.cpmd
        else: prodfile = self.cvmd
        open(outf,'w').write(prodfile.substitute(substDict))

        # If replica has restraints, write ptraj imaging scripts for restart files
        if replica.hasRestraints and iwrap:
            self.writeCheckCOM(replica.top, replica.crd, m.lstrip(':'), self.getPtrajImagingCommands())

        # Once done leave replica directory back to previous dir
        T.BROWSER.chdir(cwd)
        return True


class OpenMMCheck(object):
    """
    Class to control execution status of an OpenMM simulation process.
    Will check if MD output files are complete or expected trajectory files exist.

    """
    def __init__(self, replica=False, warn=True, **kwargs):
        """
        :args replica: Replica to study. Assign it now or later in each method call.
        :type replica: :class:`~Replicas.Replica`
        :args bool warn: Print warnings when a file is not found or is incomplete.
        """
        self.log = logging.getLogger("OpenMMCheck")
        self.replica = replica
        self.warn = warn

    def checkMinimization(self, replica=False):
        "Check if minimization run correctly"
        replica = replica or self.replica
        if not replica: raise OpenMMCheckError, "Replica not assigned."

        # Move to replica path if not yet there
        T.BROWSER.gotoReplica(replica)
        # Check mimimisation. Loof for 'Maximum number of minimization cycles reached.' string in output
        if not osp.exists(replica.minfolder+os.sep+'min.log'):
            if self.warn: self.log.warn("Minimization output file not found: %smin.out"%(replica.minfolder+os.sep))
            T.BROWSER.goback()
            return False
        minout = open(replica.minfolder+os.sep+'min.log','r').read()
        T.BROWSER.goback()
        return True

    def checkEquilibration(self, replica=False, stepselection=[], outextension='log'):
        """
        Check if equilibration run correctly.

        :arg replica: Replica under study. Default: Loaded replica at instantiation.
        :type replica: :class:`~Replicas.Replica`

        :args list stepselection: Equilibration file step number to check. If False, check all.

        Returns: True or False
        """
        replica = replica or self.replica
        if not replica: raise OpenMMCheckError, "Replica not assigned."
        if not isinstance(stepselection, list): stepselection=[stepselection]

        # Check all equilibration steps (look for 'Total CPU time:')
        selection = stepselection or range(1, 6)
        for i in selection:
            out = self.getEquilibrationOutputFile(i, replica, outextension)
            if not out:
                if self.warn: self.log.warn("Equilibration output file not found for step %i"%i)
                return False
        return True

    def getEquilibrationOutputFile(self, step, replica=False, outextension='log'):
        """
        Return the file content of the equilibration output file for step *step*. Return **False** if file not found.

        :arg int step: Step number to identify equilibration output file. Will use :attr:`Replica.eqoutfiletemplate` to match names.
        :arg replica: Replica under study. Default: Loaded replica at initialization.
        :type replica: :class:`~Replicas.Replica`
        :arg str outextension: Expected file extension for the MD output.

        :return: file content or **False**
        :rtype: str or bool
        """
        replica = replica or self.replica
        if not replica: raise OpenMMCheckError, "Replica not assigned."

        # Move to replica path if not yet there
        T.BROWSER.gotoReplica(replica)
        mdoutfile = replica.eqfolder+os.sep+replica.eqoutfiletemplate.format(step=step,extension=outextension)
        if not osp.exists(mdoutfile): out = False
        else: out = open(mdoutfile, 'r').read()
        T.BROWSER.goback()
        return out


    def getProductionOutputFile(self, step, replica=False, outextension='log'):
        """
        Return the file content of the production output file for step *step*. Return **False** if file not found.

        :arg int step: Step number to identify production output file. Will use :attr:`Replica.mdoutfiletemplate` to match names.
        :arg replica: Replica under study. Default: Loaded replica at initialization.
        :type replica: :class:`~Replicas.Replica`
        :arg str outextension: Expected file extension for the MD output.

        :return: file content or **False**
        :rtype: str or bool
        """
        replica = replica or self.replica
        if not replica: raise OpenMMCheckError, "Replica not assigned."

        # Move to replica path if not yet there
        T.BROWSER.gotoReplica(replica)
        mdoutfile = replica.mdfolder+os.sep+replica.mdoutfiletemplate.format(step=step,extension=outextension)
        if not osp.exists(mdoutfile): out = False
        else: out = open(mdoutfile, 'r').read()
        T.BROWSER.goback()
        return out

    def checkProduction(self, replica=False, stepselection=[], outextension='log'):
        """
        Check if production run correctly.

        :arg replica: Replica under study. Default: Loaded replica at instantiation.
        :type replica: :class:`~Replicas.Replica`

        :args list stepselection: Production file step number to check. If False, check all.

        Returns: True or False
        """
        replica = replica or self.replica
        if not replica: raise OpenMMCheckError, "Replica not assigned."
        if not isinstance(stepselection, list): stepselection=[stepselection]

        selection = stepselection or range(1, replica.ntrajfiles+1)
        for i in selection:
            out = self.getProductionOutputFile(i, replica, outextension)
            if not out:
                if self.warn: self.log.warn("Production output file not found for step %i"%i)
                return False
        return True

    def checkMD(self, replica=False, returnsteps=False):
        "Check in current folder min/ eq/ and md/ for N nanoseconds taken from project info\
         if 'returnsteps', don't evaluate and just return True or False for min, eq and md steps in a dictironary."
        replica = replica or self.replica
        if not replica: raise OpenMMCheckError, "Replica not assigned."

        stepsdone = {}
        stepsdone['min'] = self.checkMinimization(self.replica)
        stepsdone['eq'] = self.checkEquilibration(self.replica)
        stepsdone['md'] = self.checkProduction(self.replica)

        # Return steps?
        if returnsteps: return stepsdone

        # Evaluate to True or False if all is done or not respectively.
        if npy.sum([stepsdone.values()]) == 3:
            self.log.info("Simulation completed for replica %s"%replica.name)
            return True
        else:
            if self.warn: self.log.warn("Checking replica MD failed. Some steps could not pass the check: %s"%stepsdone)
            return False

    
import Biskit.test as BT

class Test(BT.BiskitTest):
    """Test"""

    def test_OpenMMWriter(self):
        """Create new replica and write MDinput"""
        from MDSettings import MDSettings
        from Systems import SolvatedSystem
        
        top = osp.join(T.testRoot('pep', 'pep.prmtop'))
        crd = osp.join(T.testRoot('pep', 'pep.prmcrd'))
        sys = SolvatedSystem(name='pep',top=top, crd=crd)
        settings = MDSettings(solvent='WAT',mdProgram='OpenMM',restrMode='HA', restrForce=0.1)
        
        self.testdir =  T.tempDir()
        self.r1 = sys+settings
        self.r1.setName('testOpenMM')
        
        T.BROWSER.chdir(self.testdir)
        
        # write replica folder and check methods of AmberWriter
        self.r1.createFolder()
        self.r1.createMDInput()
        writer = OpenMMWriter(self.r1)
        
        self.testdir += os.sep+'testOpenMM'
    
    def cleanUp(self):
        T.tryRemove( self.testdir, tree=1 )

if __name__ == "__main__":
    BT.localTest()
