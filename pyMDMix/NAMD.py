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

class NAMDWriterError(Exception):
    pass

class NAMDCheckError(Exception):
    pass

class NAMDWriter(object):
    def __init__(self, replica=False):
        self.log = logging.getLogger("NAMDWriter")
        self.replica = replica

        # Load template inputs
        self.loadConfig()

    def loadConfig(self, min=None, heating=None, eqNPT=None, mdNVT=None, restr=None):
        """
        Load standard template files or user given template files for NAMD MD configuration.
        All template files should contain expected fieldnames. Files are read and loaded for later substitution.

        :arg str min: Filepath for minimization template file.
        :arg str heating: Filepath for heating step template file.
        :arg str eqNPT: Filepath for equilibration with NPT ensemble template file.
        :arg str mdNVT: Filepath for production with NVT ensemble template file.
        :arg str restr: Retrains template. Section that can be added to any file to include restraints.
        """
        if min: self.minT = string.Template(open(min,'r').read())
        else:   self.minT = string.Template(open(T.templatesRoot('namd_min_temp.txt'), 'r').read())
        if heating: self.eqIncrT = string.Template(open(heating, 'r').read())
        else:       self.eqIncrT = string.Template(open(T.templatesRoot('namd_eq_nvt_incrT_temp.txt'), 'r').read())
        if eqNPT:   self.eqNPT = string.Template(open(eqNPT, 'r').read())
        else:       self.eqNPT = string.Template(open(T.templatesRoot('namd_eq_npt_temp.txt'), 'r').read())
        if mdNVT:   self.mdNVT = string.Template(open(mdNVT, 'r').read())
        else:       self.mdNVT = string.Template(open(T.templatesRoot('namd_md_temp.txt'), 'r').read())
        self.restr = ''
        if restr:   self.restr = open(restr, 'r').read()
        else:       self.restr = open(T.templatesRoot('namd_restr_temp.txt'), 'r').read()

    def getBoxFromCRD(self, crd):
        "Read box size from CRD file bottom line"
        boxline = open(crd, 'r').readlines()[-1]
        return npy.array( boxline.split()[0:3] , dtype='float32')

    def createReplicaRestraintPDB(self, replica=False, inputpdb=False):
        """
        Create in the current folder a pdb with restraint force in B-factor column for NAMD input.
        The reference is the starting CRD file if **inputpdb** is not given. Output name will be *restrains.pdb*.
        
        :args replica: Replica to work with. If **False**, work with loaded replica at instantiation.
        :type replica: :class:`~Replicas.Replica`
        
        :args str inputpdb: PDB file path to use as reference coordinates for the restraints. If **False**, use replica.pdb as reference.
        
        :return: True if file correctly saved or replica.restrMode == FREE. False otherwise.
        """        
        replica = replica or self.replica
        if not replica: raise NAMDWriterError, "Replica not assigned."
        
        replica.go()
        self.log.info("Creating restraints.pdb with reference positions and restraining forces at b-factor column")
        if inputpdb:
            from PDB import SolvatedPDB
            pdb = SolvatedPDB(inputpdb, extraResidues=replica.extraResidues)
        else:
            pdb = replica.system.getSolvatedPDB()

        restrainLevel = replica.restrMode
        mask = replica.restrMask
        force = replica.restrForce
        if restrainLevel == 'FREE':  return True

        if mask == 'AUTO' or not mask: # automask
            if restrainLevel == 'BB':  # BB
                atIds = pdb.getBBAutoMaskAtomIds()
            else:   # Heavy atoms
                atIds = pdb.getHAAutoMaskAtomIds()
        else:
            if restrainLevel == 'BB': # BB
                atIds = pdb.getBBMaskSelectedRes(mask)
            else:
                atIds = pdb.getHAMaskSelectedRes(mask)

        self.log.debug("Restrain mode: %s"%restrainLevel)
        self.log.debug("Restrain mask: %s"%mask)
        self.log.debug("Identified atoms for restrain application: %s"%atIds)
        self.log.debug("Force to apply: %.3f"%force)
        maskAtoms = pdb.maskFrom('serial_number', atIds.tolist())
        pdb['temperature_factor'] = npy.zeros(len(pdb))
        pdb['temperature_factor'][maskAtoms] = force
        pdb.writePdb('restrains.pdb')
        
        exists = osp.exists('restrains.pdb')
        if exists: self.log.debug("restraints.pdb correctly saved")
        replica.goback()
        return exists

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
        if not replica: raise NAMDWriterError, "Replica not assigned."

        if process == 'min':
            command = S.NAMD_EXE+' min.in &> min.out'
            return command

        elif process == 'eq':
            if not step: return False
            outname = replica.eqoutfiletemplate.format(step=step, extension='out')
            command = S.NAMD_EXE+' eq%i.in &> %s'%(step, outname)
            return command

        elif process == 'md':
            if not step: return False

            mdouttemplate=replica.mdoutfiletemplate.replace('.{extension}','')
            outname = mdouttemplate.format(step=step)
            command = S.NAMD_EXE+' %s.in &> %s.out'%(outname, outname)
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
        if not replica: raise NAMDWriterError, "Replica not assigned."

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
        if not replica: raise NAMDWriterError, "Replica not assigned."
        commands = self.getReplicaCommands(replica)
        open(outfile,'w').write('\n'.join(commands))
        ok = osp.exists(outfile)
        if ok: self.log.debug("Wrote commands file %s"%outfile)
        else: self.log.warn("COMMANDS.sh not writen!")
        return ok

    def writeMinInput(self, replica=False):
        """Write minimization input file for NAMD in replica.minfolder
         :args replica: Replica instance. If False, use replica assigned in instantiation.
         :type replica: :class:`~Replicas.Replica`
        """
        replica = replica or self.replica
        if not replica: raise NAMDWriterError, "Replica not assigned."
        
        T.BROWSER.gotoReplica(replica)
        
        restr = ''
        if replica.hasRestraints:
            if not replica.minimizationAsRef: restr = self.restr
            else: 
                self.log.warn('Use of Minimized structure as restraint reference is still not possible with NAMD. Will use starting PRMCRD.')
                restr = self.restr
        
        formatdict = {'top':replica.top, 'crd':replica.crd, 'restraints':restr, 
                        'box':self.getBoxFromCRD(replica.crd).max(), 'timestep':int(replica.md_timestep),
                        'freq':replica.trajfrequency}
        formatdict['minsteps'] = replica.minsteps
        out = replica.minfolder+os.sep+'min.in'
        open(out,'w').write(self.minT.substitute(formatdict))
        exists = osp.exists(out)
        T.BROWSER.goback()
        
        return exists

    def writeEqInput(self, replica=False):
        """Write equilibration input file for NAMD in replica.eqfolder
         Args:
            replica     (ReplicaInfo)   Replica instance
        """
        replica = replica or self.replica
        if not replica: raise NAMDWriterError, "Replica not assigned."
        
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
        eq1out = replica.eqfolder+os.sep+'eq1.in'
        open(eq1out,'w').write(self.eqIncrT.substitute(formatdict))

        formatdict['eqinput'] = formatdict['eqoutput']
        formatdict['first_step'] = formatdict['final_step']

        # SECOND STEP
        # NPT equilibration for 1ns
        formatdict['eqoutput'] = 'eq2'
        formatdict['final_step'] = formatdict['first_step'] + replica.npt_eq_steps
        eq2out = replica.eqfolder+os.sep+'eq2.in'
        open(eq2out,'w').write(self.eqNPT.substitute(formatdict))
        exists = osp.exists(eq1out) and osp.exists(eq2out)
        
        T.BROWSER.goback()
        return exists

    def writeMDInput(self, replica=False):
        """Write production stage input file for NAMD in :attr:`~pyMDMixReplicas.Replica.mdfolder`
         Args:
            replica     (ReplicaInfo)   Replica instance
        """
        replica = replica or self.replica
        if not replica: raise NAMDWriterError, "Replica not assigned."
        
        T.BROWSER.gotoReplica(replica)
        
        restr = ''
        if replica.hasRestraints: restr = self.restr
        
        # PRODUCTION
        # Prepare md configuration files for each trajectory file
        substDict = {'top':replica.top, 'crd':replica.crd, 'restraints':restr,
                        'timestep':replica.md_timestep, 'freq':replica.trajfrequency,
                        'temp':replica.temp}
        # Set first step which is the sum of the minimization and the equilibration steps
        # Set also first input name (here pre_out) that refers to the files of the last run
        # (coordinates, velocity and box size). That on the first step of the production
        # refer to the equilibration files outputed in the 'eq/' folder.
        substDict['first_step'] = replica.namd_heating_steps + replica.npt_eq_steps
        substDict['prev_out'] = os.path.join(os.pardir,replica.eqfolder,'eq2')

        # Proceed with the file generation
        outfiles = []
        for nfile in range(1, replica.ntrajfiles+1):
            # Use template for output formats. Remove extension and trailing dot.
            step_outname = replica.mdoutfiletemplate.format(step=nfile, extension='')[:-1]
            substDict['post_out'] = step_outname
            substDict['final_step'] = substDict['first_step'] + replica.prod_steps

            outname = replica.mdfolder+os.sep+'%s.in'%step_outname
            outfiles.append(outname)
            open(outname, 'w').write(self.mdNVT.substitute(substDict))

            substDict['first_step'] = substDict['final_step']
            substDict['prev_out'] = substDict['post_out']
        
        exists = npy.all([osp.exists(f) for f in outfiles])
        T.BROWSER.goback()
        
        return exists
        
    def writeReplicaInput(self, replica=False):
        replica = replica or self.replica
        if not replica: raise NAMDWriterError, "Replica not assigned."
        
        T.BROWSER.gotoReplica(replica)
        self.log.info("Writing NAMD simulation input files for replica %s"%replica.name)

        if not (osp.exists(replica.top) and osp.exists(replica.crd) and osp.exists(replica.pdb)):
            m = "Replica top, crd or pdb files does not exist in the folder: %s, %s, %s"%(replica.top, replica.crd, replica.pdb)
            raise NAMDWriterError, m

        # Create restrained pdb if needed
        if replica.hasRestraints:
            if not self.createReplicaRestraintPDB(replica):
                raise NAMDWriterError, "Could not save restrain.pdb for replica %s"%replica.name

        # Write inputs
        minok = self.writeMinInput(replica)
        eqok = self.writeEqInput(replica)
        mdok = self.writeMDInput(replica)

        if not (minok and eqok and mdok): 
            raise NAMDWriterError, "MD input not generated for replica %s"%replica.name
        
        self.log.info("MD Input OK")
        T.BROWSER.goHome()
        return True
        
class NAMDCheck(object):
    """
    Class to control execution status of an NAMD simulation process.
    Will check if MD output files are complete or expected trajectory files exist.

    """
    def __init__(self, replica=False, warn=True, **kwargs):
        """
        :args replica: Replica to study. Assign it now or later in each method call.
        :type replica: :class:`~Replicas.Replica`
        :args bool warn: Print warnings when a file is not found or is incomplete.
        """
        self.log = logging.getLogger("NAMDCheck")
        self.replica = replica
        self.warn = warn

    def checkMinimization(self, replica=False):
        "Check if minimization run correctly"
        replica = replica or self.replica
        if not replica: raise NAMDCheckError, "Replica not assigned."

        # Move to replica path if not yet there
        T.BROWSER.gotoReplica(replica)
        # Check mimimisation. Loof for 'Maximum number of minimization cycles reached.' string in output
        if not osp.exists(replica.minfolder+os.sep+'min.out'):
            if self.warn: self.log.warn("Minimization output file not found: %smin.out"%(replica.minfolder+os.sep))
            T.BROWSER.goback()
            return False
        minout = open(replica.minfolder+os.sep+'min.out','r').read()
        if not S.NAMD_FILE_COMPLETE in minout:
            if self.warn: self.log.warn("Minimization step not completed! CWD:%s"%os.getcwd())
            T.BROWSER.goback()
            return False
        T.BROWSER.goback()
        return True

    def checkEquilibration(self, replica=False, stepselection=[], outextension='out'):
        """
        Check if equilibration run correctly.

        :arg replica: Replica under study. Default: Loaded replica at instantiation.
        :type replica: :class:`~Replicas.Replica`

        :args list stepselection: Equilibration file step number to check. If False, check all.

        Returns: True or False
        """
        replica = replica or self.replica
        if not replica: raise NAMDCheckError, "Replica not assigned."
        if not isinstance(stepselection, list): stepselection=[stepselection]

        # Check all equilibration steps (look for 'Total CPU time:')
        selection = stepselection or range(1, 6)
        for i in selection:
            out = self.getEquilibrationOutputFile(i, replica, outextension)
            if not out:
                if self.warn: self.log.warn("Equilibration output file not found for step %i"%i)
                return False
            else:
                if not S.NAMD_FILE_COMPLETE in out:
                    if self.warn: self.log.warn("MD equilibration step %i not completed or errors arised"%i)
                    return False
        return True

    def getEquilibrationOutputFile(self, step, replica=False, outextension='out'):
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
        if not replica: raise NAMDCheckError, "Replica not assigned."

        # Move to replica path if not yet there
        T.BROWSER.gotoReplica(replica)
        mdoutfile = replica.eqfolder+os.sep+replica.eqoutfiletemplate.format(step=step,extension=outextension)
        if not osp.exists(mdoutfile): out = False
        else: out = open(mdoutfile, 'r').read()
        T.BROWSER.goback()
        return out


    def getProductionOutputFile(self, step, replica=False, outextension='out'):
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
        if not replica: raise NAMDCheckError, "Replica not assigned."

        # Move to replica path if not yet there
        T.BROWSER.gotoReplica(replica)
        mdoutfile = replica.mdfolder+os.sep+replica.mdoutfiletemplate.format(step=step,extension=outextension)
        if not osp.exists(mdoutfile): out = False
        else: out = open(mdoutfile, 'r').read()
        T.BROWSER.goback()
        return out

    def checkProduction(self, replica=False, stepselection=[], outextension='out'):
        """
        Check if production run correctly.

        :arg replica: Replica under study. Default: Loaded replica at instantiation.
        :type replica: :class:`~Replicas.Replica`

        :args list stepselection: Production file step number to check. If False, check all.

        Returns: True or False
        """
        replica = replica or self.replica
        if not replica: raise NAMDCheckError, "Replica not assigned."
        if not isinstance(stepselection, list): stepselection=[stepselection]

        selection = stepselection or range(1, replica.ntrajfiles+1)
        for i in selection:
            out = self.getProductionOutputFile(i, replica, outextension)
            if not out:
                if self.warn: self.log.warn("Production output file not found for step %i"%i)
                return False
            else:
                if not S.NAMD_FILE_COMPLETE in out:
                    if self.warn: self.log.warn("MD production step %i not completed or errors arised"%i)
                    return False
        return True

    def checkMD(self, replica=False, returnsteps=False):
        "Check in current folder min/ eq/ and md/ for N nanoseconds taken from project info\
         if 'returnsteps', don't evaluate and just return True or False for min, eq and md steps in a dictironary."
        replica = replica or self.replica
        if not replica: raise NAMDCheckError, "Replica not assigned."

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

    def getSimVolume(self, replica=False, step=False, boxextension=False):
        """
        Fetch simulation volume information from restart files. 
        
        :arg Replica replica: Replica to study. If false, will take replica loaded in initalization.
        :arg int step: Step to fetch volume for. If False, will identify last completed production step and use that one.
        :arg str boxextension: Extension for the output file containing the restart information. DEFAULT: xsc.
        
        :return float Volume: Simulation volume.
        """
        replica = replica or self.replica
        if not replica: raise NAMDCheckError, "Replica not assigned."
        boxextension = boxextension or 'xsc'
        
        # Work on step. If not given, fetch last completed production step.
        step = step or replica.lastCompletedProductionStep()
        
        # Fetch file and read last line to get box side length
        fname = replica.mdoutfiletemplate.format(step=step, extension=boxextension)
        fname = osp.join(replica.path, replica.mdfolder, fname)
        if not os.path.exists(fname):
            self.log.error("No file found with name %s to fetch box volume in DG0 penalty calculation. Returning no penalty..."%fname)
            return False
        
        # Read file and fetch 3 vectors to calculate volume
        box = map(float, open(fname,'r').readlines()[2].strip().split())
        vec_a = npy.array(box[1:4]).astype(float)
        vec_b = npy.array(box[4:7]).astype(float)
        vec_c = npy.array(box[7:10]).astype(float)
        vol = npy.linalg.norm(a)*npy.linalg.norm(b)*npy.linalg.norm(c)
        if b[0] != 0: # Assume we have an orthorombic box
            vol *= 0.77
        return vol

import Biskit.test as BT

class Test(BT.BiskitTest):
    """Test"""

    def test_NAMDWriter(self):
        """Create new replica and write MDinput"""
        from MDSettings import MDSettings
        from Systems import SolvatedSystem
        
        top = osp.join(T.testRoot('pep', 'pep.prmtop'))
        crd = osp.join(T.testRoot('pep', 'pep.prmcrd'))
        sys = SolvatedSystem(name='pep',top=top, crd=crd)
        settings = MDSettings(solvent='WAT',mdProgram='NAMD',restrMode='HA', restrForce=0.1)
        
        self.testdir =  T.tempDir()
        self.r1 = sys+settings
        self.r1.setName('testNAMD')
        
        T.BROWSER.chdir(self.testdir)
        
        # write replica folder and check methods of AmberWriter
        self.r1.createFolder()
        self.r1.createMDInput()
        writer = NAMDWriter(self.r1)
        
#        self.assertEqual(writer.getAmberRestrMask(), ':1-8 & !@H=')
#        self.assertTrue(writer.writeCommands())
#        self.assertTrue(writer.writeReplicaInput())
#        self.assertEqual(writer.getReplicaCommands(), checkCommands)
        
        self.testdir += os.sep+'testNAMD'
    
    def cleanUp(self):
        T.tryRemove( self.testdir, tree=1 )

if __name__ == "__main__":
    BT.localTest()
