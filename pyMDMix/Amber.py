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

import os
import os.path as osp
import re
import sys
import time
import logging
import fnmatch
import string
import tempfile
import subprocess as sub

import numpy as npy

import tools as T
import settings as S
import Solvents

# ERRORS DEFINITION
class AmberWriterError(Exception):
    pass

class AmberCheckError(Exception):
    pass

class AmberCreateSystemError(Exception):
    pass

class Leaper(object):
    """Class to handle tLeap I/O"""
    def __init__(self,  objectFileName=None, extraff=[], cwd=None):
        self.tleap = sub.Popen(S.AMBEREXE+os.sep+"tleap -f - ",  shell=True,   stdin = sub.PIPE,  stdout=sub.PIPE,  stderr=sub.PIPE, cwd=cwd)
        self._in = self.tleap.stdin
        self._out = self.tleap.stdout
        self._err = self.tleap.stderr
        self._setterminator = 'mdmix_terminator_command = "mdmix_terminator_command"'
        self._terminator = 'desc mdmix_terminator_command'
        self.command(self._setterminator)
        if objectFileName: self.command('loadOff %s'%(objectFileName))

        # Load default and required forcefields
        ff = extraff
        for f in ff:
            if 'frcmod' in f: self.command('loadAmberParams %s'%f)
            else: self.command('source %s'%f)
        self.initialized = True


    def command(self,  command):
        """
        Command must be a string corresponding to a correct Leap command.
        Returns the output of leap for the given command.
        tLeap is a bit tricky to be controlled by stdin, it terminates when it reaches EOF.
        Thus we use a fake command to detect the last line before EOF.
        """
        self._in.write(command+'\n')
        logging.getLogger("AmberTleap").debug(command)
        self._in.write(self._terminator+'\n')
        return self.flush()
        #self._in.write('lastCommandOut\n')
        #out = []
        #while self.tleap.returncode is None:
            #line = self._out.readline().strip()
            #if 'Fatal Error!' in line: break
            #if 'ERROR: syntax error' in line: break
            #if 'Error from the parser: syntax error' in line: break
            #if 'Exiting LEaP' in line:
            #    logging.getLogger("AmberTleap").debug(line)
            #    break    

    def flush(self):
        out = []
        done = False
        while not done:
        #for line in self._out.readline().splitlines():
            #line = line.strip()
            #if not line: break
            line = self._out.readline().strip()
            if 'mdmix_terminator_command' in line:
                logging.getLogger("AmberTleap").debug(line)    
                done = True
                break
            if self.tleap.returncode is not None:
                break
            if 'Fatal Error!' in line:
                break
            if 'Exitting LEaP' in line:
                break
            if line:
                out.append(line)
                logging.getLogger("AmberTleap").debug(line)
        if not done:
            logging.getLogger("AmberTleap").warning("LEaP may have failed")
        return out
    
    def status(self):
        return self.tleap.returncode

    def close(self):
        #if self.tleap.returncode is None:
        #    self.command("quit")
        #return self.tleap.wait()
        self.tleap.terminate()


class AmberCreateSystem(object):
    def __init__(self, replica=False, FFlist=[], informative=True):
        """
        :arg bool informative: Log in INFO level and information on the process. Log only at DEBUG level.
        """
        # Set logger and controls
        self.log = logging.getLogger("AmberCreateSystem")
        self.replica = False
        self.leap = None

        # Load default and user FF lists
        self.FFlist = []
        [self.addFF(f) for f in S.DEF_AMBER_FF]
        [self.addFF(ff) for ff in FFlist]

        # Set replica
        self.workOnReplica(replica)

        # Pipe all information to debug mode if informative=False
        if not informative:
            self.log.info = self.log.debug
        
    def workOnReplica(self, replica):
        if replica:
            self.replica = replica
            [self.addFF(f) for f in replica.FF]

    def listFF(self):
        "List available forcefield files in Amber directory"
        dir = osp.join(S.AMBERHOME,'dat','leap','cmd')
        return fnmatch.filter(os.listdir(dir), 'leaprc.*')
    
    def listFrcmod(self):
        "List available forcefield files in Amber directory"
        dir = os.path.join(S.AMBERHOME,'dat','leap','parm')
        return fnmatch.filter(os.listdir(dir), '*frcmod*')

    def initLeap(self):
        "Initialize leap interface loading chosen forcefields and solvents library"
        self.leap = Leaper()
        self.log.debug("Initializing Leap...")
        for ff in self.FFlist:
            self.log.debug(ff)
            if 'leaprc' in ff: self.leap.command("source %s"%ff)
            elif 'frcmod' in ff: self.leap.command("loadamberparams %s"%ff)
            self.log.debug("Added FF: %s"%ff)
        self.log.debug("Leap initialized...")

    def checkFF(self, ffname):
        """
        Will check if filename exists as given by the user or try to locate it inside AMBERHOME directory.
        :return: absolute path to file if file exists and is found. **False** otherwise (file not found).
        """
        if osp.exists(ffname):
            ffname = osp.abspath(ffname)
            return ffname
        else:
            # If it does not exists, try to search it in amber directory
            path1 = osp.join(S.AMBERHOME,'dat','leap','parm')
            path2 = osp.join(S.AMBERHOME,'dat','leap','cmd')
            okffname = False
            for p in (path1, path2):
                for root, dirname, files in os.walk(p):
                    match = fnmatch.filter(files, '*%s*'%ffname)
                    if match and len(match) == 1:
                        okffname = osp.join(root, match[0])
                    elif match and ffname in match:
                        okffname = osp.join(root, ffname)
                    elif match:
                        self.log.warn("Several files found matching extra FF name '%s' found: %s. Using '%s'."%(ffname, match, match[0]))
                        okffname = osp.join(root, match[0])

            if okffname: return okffname
            else: return False
            
    def addFF(self, ffname):
        """
        Set forcefield and forcefield modifications to be loaded when executing Leap.
        Will previously check if filename exists as given by the user or try to locate it inside AMBERHOME directory

        :return: **True** if correctly loaded. **False** otherwise (file not found).
        """
        okffname = self.checkFF(ffname)
        if okffname:
            self.FFlist.append(okffname)
            self.log.info("Using Forcefield or FRCMOD file: %s"%okffname)
            self.FFlist = list(set(self.FFlist))
            return True
        else:
            self.log.warn("%s forcefield or frcmod not found. Will not be loaded!"%ffname)
            self.log.warn("FORCEFIELD LEAPRC FOUND: %s"%self.listFF())
            self.log.warn("FORCEFIELD FRCMOD FOUND: %s"%self.listFrcmod())
            return False

    def loadOff(self, objectFile):
        if not self.leap: self.initLeap()
        self.leap.command("loadOff %s"%objectFile)
        self.amberOFF = objectFile
        self.log.debug("Loaded off file %s into Leap"%objectFile)
        return True
    
    def ambpdb(self, top, crd, outpdb):
        # Save pdb from top and crd
        self.log.info("Creating PDB file %s from top and crd files '%s' and '%s'"%(outpdb, top, crd))
        null = os.devnull
#        null = osp.splitext(top)[0]+'_ambpdb.log'
        #ampdb = osp.join(S.AMBEREXE,S.AMBER_AMBPDB)
        # original line:
        #command = '(%s -p %s < %s > %s) 1> %s 2> %s'%(S.AMBER_AMBPDB, top, crd, outpdb,null,null)
        # debug line:
        command = '(%s -p %s < %s > %s) '%(S.AMBER_AMBPDB, top, crd, outpdb)
        self.log.debug(command)
        proc = sub.Popen(command, shell=True)
        exit_code = proc.wait()
        if exit_code: # Exit different to zero means error
            self.log.error("Could not save PDB from TOP and CRD. Check %s file in project folder for the reason."%null)
            raise AmberCreateSystemError, "Could not save PDB from TOP and CRD. Check leap.log file in project folder."
        return True

    def saveAmberParm(self, unit, top, crd, leapHandler=None):
        if leapHandler: self.leap = leapHandler   # Useful when working on same unit to not init leap each time (like solvate and then save)
        elif not self.leap: self.initLeap()
        self.log.debug("Saving unit '%s' TOP and CRD files: '%s', '%s'..."%(unit,top,crd))
        self.log.debug(self.leap.status())
        err = self.leap.command("saveAmberParm %s %s %s"%(unit, top, crd))
        time.sleep(1)
        if osp.exists(top) and osp.exists(crd) and os.stat(top).st_size and os.stat(crd).st_size:
            # Check files exists and do not have zero size!
            self.log.info("TOP and CRD files correctly created for unit '%s'."%unit)
            return True
        else:
            self.log.error("Error saving TOP '%s' and CRD '%s' file for unit '%s' in Object File '%s'."%(top, crd, unit, self.amberOFF))
            self.log.error('Most of the times this error is due to missing parameters.')
            self.log.error('CHECK that you correctly provided all needed forcefields in EXTRAFF flag when creating the project and that you correctly prepared the system Object File.')
            self.log.error("This is Leap Error Message:\n%s"%('\n'.join(err)))
            raise IOError("Error saving TOP '%s' and CRD '%s' file for unit '%s' in Object File '%s'."%(top, crd, unit, self.amberOFF))
        
    def leapCharge(self, unit):
        "return charge of object unit"
        if not self.leap: self.initLeap()
        out = self.leap.command("charge %s"%unit)
        self.log.debug(out)
        return float(out[0].split()[-1])
        
    def neutralizeWNaCl(self, unit):
        "Neutralize in leap using NaCl. unit should be a solvated system."
        if not self.leap: self.initLeap()
        self.log.debug("Neutralizing with NaCl unit %s"%unit)
        if self.replica: watmodel = Solvents.getSolvent(self.replica.solvent).watermodel
        else: watmodel = S.DEF_AMBER_WATBOX
        waterff = {
            'TIP3P': 'leaprc.water.tip3p',
            'TIP4P': 'leaprc.water.tip3pew',
            'SPC': 'leaprc.water.spce'
        }
        if waterff.has_key(watmodel):
            self.log.debug("Loading ions params: %s"%waterff[watmodel])
            self.leap.command("source %s"%waterff[watmodel])
        else:
            self.log.debug("No extra ions params loaded")
            
        self.leap.command("addIons %s Na+ 0"%unit)
        self.leap.command("addIons %s Cl- 0"%unit)
        charge = self.leapCharge(unit)
        if charge < -0.001 or charge > 0.001:
            self.log.warn("Could not neutralize box in Leap. Charge: %.4f"%charge)
            return False
        
        self.log.info("Neutralize NaCl OK")    
        return True

    def neutralizeIonicBox(self, unit, negativeres, positiveres):
        "Special neutralization of solvent mixtures which have positive and negative residues"
        if not self.leap: self.initLeap()
        charge = self.leapCharge(unit)
        desc = self.leap.command("desc %s"%unit)[5:]
        self.log.debug("Neutralize ionic box unit %s, negativeres: %s, positiveres:%s"%(unit, negativeres, positiveres))

        # Get first residue number
        digit = re.compile('\d+')
        first = int(digit.search(desc[0]).group())

        NEGids = npy.where([negativeres in e for e in desc])[0]
        POSids = npy.where([positiveres in e for e in desc])[0]
        NEGids = NEGids.tolist()
        POSids = POSids.tolist()
        [l.sort() for l in (POSids, NEGids)]

        # remove needed solvents to reach charge 0
        if charge > 0:
            for ion in range(int(charge)):
                idx = POSids.pop()
                self.leap.command("remove %s %s.%d"%(unit,unit,idx+first))
        elif charge < 0:
            for ion in range(abs(int(charge))):
                idx = NEGids.pop()
                self.leap.command("remove %s %s.%d"%(unit,unit,idx+first))
        finalcharge = self.leapCharge(unit)
        if finalcharge != 0:
            self.log.warn("Could not neutralize ionic box in Leap. Charge: %.2f. Will try adding NaCl..."%finalcharge)
            if self.neutralizeWNaCl(unit):
                self.log.warn("Neutralization succeded with NaCl. Continuing...")
                return True
            else:
                self.log.warn("No sucess neutralizing with NaCl.. might give problems during simulation.")
                return False
        
        self.log.debug("Neutralize ionic box OK")
        return True
    
    def __saveTempOff(self):
        tmp = tempfile.mktemp(prefix='solv')
        self.solvent.writeOff(tmp)
        self.log.debug("Write temporary off file:%s"%tmp)
        return tmp

    def loadSolventLeap(self, solvent, leapHandle=False):
        # Get solvent object file from DB. Store it temporary to load it in leap.
        if isinstance(solvent, str): self.solvent = Solvents.getSolvent(solvent)
        else: self.solvent = solvent
        if not self.solvent:
            self.log.error("Solvent Name %s not found in the database"%solvent)
            raise NameError("Solvent Name %s not found in the database"%solvent)

        tmpOff = self.__saveTempOff()
        if not leapHandle: 
            if not self.leap: self.initLeap()
            leapHandle = self.leap
            
        leapHandle.command('loadOff %s'%tmpOff)
        return self.solvent.boxunit

    def solvateOrganic(self, unit, solvent=False, buffer=False, cubicBox=False):
        """
        Solvate unit in solvent box. Use saveAmberParm to save top and crd files.
        
        :arg str unit: tLeap unit name to be solvated.
        :arg str solvent: Solvent mixture name to use for solvating the system. If not given, will be taken from self.replica.
        :arg float buffer: Buffer distance to add in solvate command. By default will take settings.DEF_AMBER_BUFFER.
        :arg bool cubicBox: Use a cubic box (TRUE) instead of an orthorombic box (FALSE). Default: Take from settings.DEF_AMBER_CUBICBOX.
        """
        if not solvent and self.replica: solvent = self.replica.solvent
        if not buffer: buffer = S.DEF_AMBER_BUFFER
        if not cubicBox: cubicBox = S.DEF_AMBER_CUBICBOX == 1 # Try to take info from settings
        if isinstance(solvent, str): solvent = Solvents.getSolvent(solvent)
        self.log.debug("SolvateOrganic: unit %s solvent %s buffer %.3f cubicBox %s"%(unit, solvent, buffer, cubicBox))
        
        # Load solvent into leap
        solventBox = self.loadSolventLeap(solvent)

        # Solvate command depends on cubic or orthorombic box choice
        if cubicBox: solvate_cmd = 'solvateBox'
        else: solvate_cmd = 'solvateOct'

        # add box of solvents
        # and neutralize charge
        self.leap.command("%s %s %s %.2f iso 1"%(solvate_cmd, unit,solventBox,buffer))

        if solvent.isIonic():        # Only in Ionic Box
            negative_res = [r.name for r in solvent.residues if r.charge < 0]
            positive_res = [r.name for r in solvent.residues if r.charge > 0]
            if negative_res and positive_res:
                corrected = self.neutralizeIonicBox(unit, negative_res[0], positive_res[0])
            else:
                # Just one of the partners is present, neutralize with ions
                corrected = self.neutralizeWNaCl(unit)
        else:
            corrected = self.neutralizeWNaCl(unit)

        if not corrected:
            self.log.error("Could not neutralize system. Fix it manually. Exiting...")
            raise RuntimeError("Could not neutralize system. Fix it manually.")

        self.log.debug("solvateOrganic OK")
        return True
    
    def solvateWater(self, unit, outPrefix, watbox=False, buffer=False, cubicBox=False):
        """
        Take unit and solvate in tip3pbox water box.
        
        :arg str unit: tLeap unit name to be solvated.
        :arg str outPrefix: Prefix to output prmtop and prmcrd files.
        :arg str watbox: Water model to use. Default will be taken from settings.DEF_AMBER_WATBOX.
        :arg float buffer: Buffer distance to add in solvate command. By default will take settings.DEF_AMBER_BUFFER.
        :arg bool cubicBox: Use a cubic box (TRUE) instead of an orthorombic box (FALSE). Default: Take from settings.DEF_AMBER_CUBICBOX.
        """
        if not watbox: watbox = S.DEF_AMBER_WATBOX+'BOX'
        if not buffer: buffer = S.DEF_AMBER_BUFFER
        if not cubicBox: cubicBox = S.DEF_AMBER_CUBICBOX == 1
        self.log.debug("SolvateWater unit %s watbox %s buffer %s cubicBox %s"%(unit, watbox, buffer, cubicBox))

        # output names
        top = outPrefix+'.top'
        crd = outPrefix+'.crd'

        # Solvate command depends on cubic or orthorombic box choice
        if cubicBox: solvate_cmd = 'solvateBox'
        else: solvate_cmd = 'solvateOct'

        # add box of water
        # and neutralize charge with NaCl
        if not self.leap: self.initLeap()
        outmessage = self.leap.command("%s %s %s %.2f iso 1"%(solvate_cmd, unit,watbox, buffer))
        if self.neutralizeWNaCl(unit):
            # System is neutralized
            # write topology and crd files
            if self.leapSaveAmberParm(unit, top, crd): 
                return True
            
        self.log.error("Error solvating in water. Leap error:\n%s"%(top,crd,'\n'.join(outmessage)))
        raise RuntimeError("Error solvating in water. Leap error:\n%s"%(top,crd,'\n'.join(outmessage)))

    def createOFF(self, outprefix, inpdb, unitname='sys', extraff=[], clean=True, **kwargs):
        """
        Create an amber object file from a pdb file.
        If *clean* is **True** will automatically try to clean the PDB: build SS bonds and add cappings.
        Check AmberPDBCleaner for more options that can be passed throught **kwargs** parameters to cleanPDB method.

        :arg str outprefix: Prefix name for cleaned PDB and Object File created.
        :arg str unitname: Unitname to assign to the unit created with leap. Default: 'sys'
        :arg list extraff: Extra forcefield parameters or modifications to be loaded in tLeap for correct recognition of non-standard residues.
        """
        import Biskit as bi
        from AutoPrepare import AmberPDBCleaner
        
        if isinstance(inpdb, str) and osp.exists(inpdb):
            inpdb = bi.PDBModel(inpdb)
        elif not isinstance(inpdb, bi.PDBModel):
            raise AmberCreateSystemError, "createOff requires a pdb file path or Biskit.PDBModel parameter as inpdb parameter"

        # Check/add FF
        for ff in extraff:
            self.addFF(ff)

        # Remove hydrogens before proceeding as many times fails occurr with hydrogen naming
        inpdb = inpdb.compress(~inpdb.maskH())

        cleaner = AmberPDBCleaner(inpdb, verbose=True, **kwargs)
        cleanpdb = cleaner.cleanPDB(**kwargs)

        if clean:
            towrite = cleanpdb
        else:
            towrite = inpdb

        towrite.writePdb(outprefix+'.pdb')
            
        # Init leap and save off
        if not self.leap: self.initLeap()
        self.leap.command('%s = loadpdb %s'%(unitname,outprefix+'.pdb'))
        err = self.leap.command('check %s'%unitname)
        if not 'Unit is OK' in err[-1]:
            raise AmberCreateSystemError, "Errors encountered when loading pdb file into Leap. Please check leap.log file for more info. Cannot automatically create the OFF file"

        else:
            # Unit is OK, let's check if there are warnings apart from charge.
            allerr = '\n'.join(err)
            if 'check:  Warnings:' in allerr:
                self.log.warn("Warnings arised during unit check. Read leap.log file for details. Charges and close contact warnings should not considered.")

            # Build ss bonds
            [self.leap.command(c) for c in cleaner.leap_ss]

            # Finally save the OFF file
            self.leap.command('saveOff %s %s'%(unitname, outprefix+'.lib'))
            self.leap.close()

            return True


class AmberCheck(object):
    """
    Class to control execution status of an AMBER simulation process.
    Will check if MD output files are complete or expected trajectory files exist.
    
    """
    def __init__(self, replica=False, warn=True, **kwargs):
        """
        :args replica: Replica to study. Assign it now or later in each method call.
        :type replica: :class:`~Replicas.SolvatedReplica`
        :args bool warn: Print warnings when a file is not found or is incomplete.
        """
        self.log = logging.getLogger("AmberCheck")
        self.replica = replica
        self.warn = warn

    def checkMinimization(self, replica=False):
        "Check if minimization run correctly"
        replica = replica or self.replica
        if not replica: raise AmberCheckError, "Replica not assigned."

        # Move to replica path if not yet there
        T.BROWSER.gotoReplica(replica)
        # Check mimimisation. Loof for 'Maximum number of minimization cycles reached.' string in output
        if not osp.exists(replica.minfolder+os.sep+'min.out'):
            if self.warn: self.log.warn("Minimization output file not found: %smin.out"%(replica.minfolder+os.sep))
            T.BROWSER.goback()
            return False
        minout = open(replica.minfolder+os.sep+'min.out','r').read()
        if not S.AMBER_MIN_COMPLETE in minout:
            if self.warn: self.log.warn("Minimization step not completed! CWD:%s"%os.getcwd())
            T.BROWSER.goback()
            return False
        T.BROWSER.goback()
        self.log.debug("Minimization complete (replica %s)"%replica.name)
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
        if not replica: raise AmberCheckError, "Replica not assigned."
        if not isinstance(stepselection, list): stepselection=[stepselection]

        # Check all equilibration steps (look for 'Total CPU time:')
        selection = stepselection or range(1, 6)
        for i in selection:
            out = self.getEquilibrationOutputFile(i, replica, outextension)
            if not out:
                if self.warn: self.log.warn("Equilibration output file not found for step %i"%i)
                return False
            else:
                if not S.AMBER_MD_COMPLETE in out:
                    if self.warn: self.log.warn("MD equilibration step %i not completed or errors arised"%i)
                    return False
        self.log.debug("Equilibration complete (replica %s, step_selection %s)"%(replica.name,stepselection))
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
        if not replica: raise AmberCheckError, "Replica not assigned."

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
        if not replica: raise AmberCheckError, "Replica not assigned."

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
        if not replica: raise AmberCheckError, "Replica not assigned."
        if not isinstance(stepselection, list): stepselection=[stepselection]

        selection = stepselection or range(1, replica.ntrajfiles+1)
        for i in selection:
            out = self.getProductionOutputFile(i, replica, outextension)
            if not out:
                if self.warn: self.log.warn("Production output file not found for step %i"%i)
                return False
            else:
                if not S.AMBER_MD_COMPLETE in out:
                    if self.warn: self.log.warn("MD production step %i not completed or errors arised"%i)
                    return False
        
        self.log.debug("Production complete (replica %s step_selection %s)"%(replica.name, stepselection))
        return True    

    def checkMD(self, replica=False, returnsteps=False):
        "Check in current folder min/ eq/ and md/ for N nanoseconds taken from project info\
         if 'returnsteps', don't evaluate and just return True or False for min, eq and md steps in a dictironary."
        replica = replica or self.replica
        if not replica: raise AmberCheckError, "Replica not assigned."

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
        :arg str boxextension: Extension for the output file containing the restart information. DEFAULT: rst.
        
        :return float Volume: Simulation volume.
        """
        replica = replica or self.replica
        if not replica: raise AmberCheckError, "Replica not assigned."
        
        boxextension = boxextension or 'rst'
        
        # Work on step. If not given, fetch last completed production step.
        step = step or replica.lastCompletedProductionStep()
        
        # Fetch rst file and read last line to get box side length and angle
        fname = replica.mdoutfiletemplate.format(step=step, extension=boxextension)
        fname = osp.join(replica.path, replica.mdfolder, fname)
        if not os.path.exists(fname):
            self.log.error("No file found with name %s to fetch box volume in DG0 penalty calculation. Returning no penalty..."%fname)
            return False
        box = map(float, open(fname,'r').readlines()[-1].strip().split())
        vol = box[0]*box[1]*box[2]
        
        if box[3] != 90.0: vol *= 0.77 # orthorombic volume correction
        return vol
        

class AmberWriter(object):
    "Write input files to run the simulations with AMBER for each replica"
    def __init__(self, replica=None):
        self.log = logging.getLogger("AmberWriter")
        self.replica = replica
        self.replnatoms = 0
        self.imaging = None
        self.nlines = 0
        self.loadConfig()

    def loadConfig(self, min=None, eq1=None, eq2=None, mdNVT=None, mdNPT=None, restr=None):
        """
        Load standard template files or user given template files for AMBER MD configuration.
        All template files should contain expected fieldnames. Files are read and loaded for later substitution.

        :arg str min: Filepath for minimization template file.
        :arg str eq1: Filepath for equilibration first step template file.
        :arg str eq2: Filepath for equilibration second step template file.
        :arg str mdNVT: Filepath for production with NVT ensemble template file.
        :arg str mdNPT: Filepath for production with NPT ensemble template file.
        """
        if min: self.minT = string.Template(open(min,'r').read())
        else:   self.minT = string.Template(open(T.templatesRoot('amber_min_templ.txt'), 'r').read())
        if eq1: self.eq1T = string.Template(open(eq1, 'r').read())
        else:       self.eq1T = string.Template(open(T.templatesRoot('amber_eq1_templ.txt'), 'r').read())
        if eq2:   self.eq2T = string.Template(open(eq2, 'r').read())
        else:       self.eq2T = string.Template(open(T.templatesRoot('amber_eq2_templ.txt'), 'r').read())
        if mdNVT:   self.cvmd = string.Template(open(mdNVT, 'r').read())
        else:       self.cvmd = string.Template(open(T.templatesRoot('amber_cv_md_templ.txt'), 'r').read())
        if mdNPT:   self.cpmd = string.Template(open(mdNVT, 'r').read())
        else:       self.cpmd = string.Template(open(T.templatesRoot('amber_cp_md_templ.txt'), 'r').read())
        if restr: self.restrT = string.Template(open(restr, 'r').read().rstrip())
        else: self.restrT = string.Template(open(T.templatesRoot('amber_restr_templ.txt'), 'r').read().rstrip())


#    def __replicaatoms(self, replica):
#        import Biskit as bi
#        return len(bi.PDBModel(osp.join(replica.path, replica.pdb)))

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
        if not replica: raise AmberWriterError, "Replica not assigned."

        prevsep = os.pardir+os.sep
        top = osp.basename(replica.top)
        crd = osp.basename(replica.crd)
        if replica.mdnetcdf == '1': extension = 'nc'
        else: extension = 'x'

        command = False

        # Set reference file if restrained simulation
        if replica.hasRestraints:
            if replica.minimizationAsRef: ref = prevsep+replica.minfolder+os.sep+'min.rst'
            else: ref = prevsep+crd
        else:
            ref = False
            
        # Check iwrap option used
        # if 1 and replica has restraints, use check_com.sh
        # do not use it otherwise
        iwrap = replica.iwrap == '1'
            
        if process == 'min':
            command = S.AMBER_MIN_EXE+' -O -i min.in -o min.out -p %s -c %s -r min.rst'%(prevsep+top, prevsep+crd)
            if ref and not replica.minimizationAsRef: command += ' -ref %s'%ref
            return command

        elif process == 'eq':
            if not step: return False
            if step == 1:
                #First step
                eqfname = replica.eqoutfiletemplate.format(step=step,extension='')
                command = S.AMBER_PROD_EXE+' -O -i eq1.in -o %sout -p %s -c %smin.rst -r %srst -x %s%s'%(eqfname, prevsep+top,
                                                prevsep+replica.minfolder+os.sep, eqfname, eqfname, extension)

            elif step > 1:
                eqfname = replica.eqoutfiletemplate.format(step=step,extension='')
                preveqfname = replica.eqoutfiletemplate.format(step=step-1,extension='')
                command =  S.AMBER_PROD_EXE+' -O -i eq%i.in -o %sout -p %s -c %srst -r %srst -x %s%s'% \
                                                            (step, eqfname, prevsep+top, preveqfname, eqfname, eqfname, extension)
            if ref:
                # Add reference flag and trajectory imaging of restart file
                command += ' -ref %s\n'%(ref)
                if iwrap: command += 'sh %scheck_com.sh %srst &> %simage.log'%(prevsep, eqfname, eqfname)

            return command
            
        elif process == 'md':
            if not step: return False

            mdouttemplate=replica.mdoutfiletemplate.replace('.{extension}','')
            if step == 1:
                fname = mdouttemplate.format(step=1)
                command = S.AMBER_PROD_EXE+' -O -i md.in -o {fname}.out -p %s -c %seq5.rst -r {fname}.rst -x {fname}.%s'%(prevsep+top, prevsep+replica.eqfolder+os.sep, extension)
                if ref:
                    command += ' -ref %s\n'%(ref)
                    if iwrap: command += 'sh %scheck_com.sh %s.rst &> %s.image.log'%(prevsep, fname, fname)

                command = command.format(fname=fname)

            elif step > 1:
                prevfname=mdouttemplate.format(step=step-1)
                nextfname=mdouttemplate.format(step=step)
                command = S.AMBER_PROD_EXE+' -O -i md.in -o {nextfname}.out -p %s -c {prevfname}.rst -r {nextfname}.rst -x {nextfname}.%s'%(prevsep+top, extension)
                fname = nextfname
                if ref:
                    command += ' -ref %s\n'%(ref)
                    if iwrap: command += 'sh %scheck_com.sh %s.rst &> %s.image.log'%(prevsep, fname, fname)

                command = command.format(nextfname=nextfname, prevfname=prevfname)

            return command
        
        else: pass

    def getReplicaCommands(self, replica=None):
        """
        get a list of commands to execute for running the MD for *replica*.
        It will contain the expected file names for input/output files and directory chamnge commands.

        :args replica: Replica to write execution commands for.
        :type replica: :class:`~Replicas.Replica`
        
        :returns: List of strings with commands to be executed.
        """
        replica = replica or self.replica
        if not replica: raise AmberWriterError, "Replica not assigned."

        if not (replica.top and replica.crd):
            raise AmberWriterError, "Replica %s does not have PRMTOP and PRMCRD files set."%replica.name

#        S.gotoReplica(replica)
#        self.log.info("Writing shell commands to run the MD for replica %s. Adapt it to your queue needs.\
#                    Filenames match the expected output names for later analysis."%replica.name)
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
        if not replica: raise AmberWriterError, "Replica not assigned."
        T.BROWSER.gotoReplica(replica)
        commands = self.getReplicaCommands(replica)
        open(outfile,'w').write('\n'.join(commands))
        if osp.exists(outfile):
            self.log.debug("Wrote commands file %s"%outfile)
            T.BROWSER.goback()
            return True
        T.BROWSER.goback()
        return False
    
    def getAmberRestrMask(self, replica=False):
        """
        Get a string with the amber formated restrain mask.
        If replica.restrMask is 'AUTO', calculate mask from solute residue ids.
        
        :args replica: Replica to obtain mask for. If False, use replica loaded at instantiation.
        :type replica: :class:`~Replicas.Replica`
        
        :returns: string with mask (E.g. ':1-124@C,N,O' for HA restrains) or **False** if
        replica has FREE restrain mode.
        """
        replica = replica or self.replica
        if not replica: raise AmberWriterError, "Replica not assigned."

        if replica.restrMode == 'FREE': return False

        if not replica.restrMask or replica.restrMask.upper() == 'AUTO':
            # Obtain mask from residue ids in solute
            if not replica.system:
                raise AmberWriterError, "Replica System not set. Can not generate mask."
            syspdb = replica.system.getSolvatedPDB()
            if not syspdb:
                raise AmberWriteError, "Error creating SolvatedPDB from System in replica %s"%replica.name
            
            out = ':'+syspdb.getAutoMaskResIds()
        else:
            out = ':'+replica.restrMask
        
        # Add restraining mode ending
        if replica.restrMode == 'BB':
            out += ' & @C,N,O'
        elif replica.restrMode == 'HA':
            out += ' & !@H='
        else:
            return
        self.log.debug("Mask: %s"%out)
        return out
    
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

        # Before starting writting, we need to obtain masks according to the option chosen
        if replica.hasRestraints:
            m = self.getAmberRestrMask()
            self.log.info("Restrain mask: %s"%m)
            mfield = self.restrT.substitute({'mask':m,'force':replica.restrForce})
        else:
            mfield = ''
            
        # Write minimization input
        substDict['minsteps'] = replica.minsteps
        # only add restraining field if we want the starting structure to be restrained
        # otherwise, the minimization will have no restraints and will use this output for 
        # future restrains
        if replica.minimizationAsRef: substDict['maskfield'] = ''
        else: substDict['maskfield'] = mfield
        outf = replica.minfolder+os.sep+'min.in'
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
        outf = replica.eqfolder+os.sep+'eq1.in'
        self.log.debug("Writing: %s"%outf)
        open(outf,'w').write(self.eq1T.substitute(substDict))
        substDict['tempi'] = substDict['tempf']

        for i in range(2,5):
            substDict['tempf'] = substDict['tempi'] + 50
            outf=replica.eqfolder+os.sep+'eq%i.in'%i
            self.log.debug("Writing: %s"%outf)
            open(outf,'w').write(self.eq2T.substitute(substDict))
            substDict['tempi'] = substDict['tempf']

        substDict['tempf'] = replica.temp
        substDict['nsteps'] = replica.npt_eq_steps# 1ns of equilibration under NPT
        substDict['pres0'] = replica.npt_pressure   # Pressure for NPT simulation, default = 1 bar
        outf = replica.eqfolder+os.sep+'eq5.in'
        self.log.debug("Writing: %s"%outf)
        open(outf,'w').write(self.cpmd.substitute(substDict))
        
        # Write md input, 1ns each file and run under NVT conditions
        substDict['nsteps'] = replica.prod_steps # 1ns each file
        outf=replica.mdfolder+os.sep+'md.in'
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

    def writeCheckCOM(self, top, crd, mask, imagingCommands):
        template = string.Template(open(T.templatesRoot('check_com_temp.txt'),'r').read())
        d = {'top':top, 'mask':mask, 'crd':crd, 'imaging':imagingCommands, 'ptraj':S.AMBER_PTRAJ}
        open('check_com.sh','w').write(template.substitute(d))

#    def writeRestartImaging(self, replica=False):
#        replica = replica or self.replica
#        if not replica: raise AmberWriterError, "Replica not assigned."
#
#        self.log.info("Writing AMBER restart file alignment to reference for each step (replica %s)..."%(replica.name))
#        cwd = T.BROWSER.cwd
#        T.BROWSER.gotoReplica(replica)
#
#        # Equilibration from 2nd to 5th step
#        outf = replica.eqfolder
#        for i in range(1,6):
#            outf=replica.eqfolder+os.sep+'alignptraj_%i.ptraj'%i
#            self.log.debug("Writing: %s"%outf)
#            trajin = replica.eqoutfiletemplate.format(step=i,extension='rst_back')
#            trajout = trajin.replace('_back','')+' restart'
#            script = self.getPtrajAlignScript(trajin=trajin, trajout=trajout, writermsd=False, writeavgpdb=False)
#            open(outf,'w').write(script)
#
#        for i in range(1,replica.ntrajfiles+1):
#            outf=replica.mdfolder+os.sep+'alignptraj_%i.ptraj'%i
#            self.log.debug("Writing: %s"%outf)
#            trajin = replica.mdoutfiletemplate.format(step=i,extension='rst_back')
#            trajout = trajin.replace('_back','')+' restart'
#            script = self.getPtrajAlignScript(trajin=trajin, trajout=trajout, writermsd=False, writeavgpdb=False)
#            open(outf,'w').write(script)
#
#        T.BROWSER.chdir(cwd)
#        return True

    def getAutoMask(self, pdb=None, extraResidues=[]):
        "From a PDB file or replica ref pdb if pdb is None, get an automask of form: 1-100,105-200"
        from PDB import SolvatedPDB
        pdb = pdb or osp.join(self.replica.path,self.replica.ref)
        extraResidues = extraResidues or self.replica.system.extraResList
        pdb = SolvatedPDB(pdb, extraResidues=extraResidues).getSolute()
        chainres = pdb.chainResIDs()
        chainres = T.simplifyNestedList(chainres,[])
        return T.numListToMask(chainres)

    def getPtrajAlignScript(self, trajin, trajout, reference=False, replica=False,
                            writermsd=True, writeavgpdb=True, rmsdoutprefix='rmsd', 
                            avgpdbout='prot_avg.pdb', alignmask=False):
        """
        Obtain a PTRAJ script for aligning the selected trajectory file (trajin) 
        and save in selected output name (trajout). All details if not given as argument will be taken
        from the replica information.
        
        :arg str trajin: Input trajectory file path to be aligned
        :arg str trajout: Output file trajectory path.
        :arg str reference: Path to reference structure for the alignment.
        :arg replica: Replica to fetch information from. If not given and AmberWriter was instantiated with
        a replica, information will be taken from that one.
        :arg bool writermsd: Output rmsd files with prefix *rmsdoutprefix*.
        :arg str rmsdoutprefix: Prefix for output files with RMSD information.
        :arg bool writeavgpdb: Write an average PDB of the protein.
        :arg str avgpdbout: Output pdb file name for the average structure.
        :arg str alignmask: If given, use this residue mask to identify over which residues the 
        trajectory should be aligned to. If not given, this information will be taken from the replic       
        
        """
        template = "trajin {trajin}\nreference {ref}\n{imaging}\n"

        if replica:
            self.replica=replica
            self.imaging = self.getPtrajImagingCommands()
        if not self.imaging:
            self.imaging = self.getPtrajImagingCommands()

        # Setting up reference pdb
        # if replica does not have referencePDB attribute set, use replica.pdb as reference
        # and issue a warning
        if not self.replica.ref:
            if self.warn: self.log.warn("Replica %s does not have a reference PDB assigned, will use replica system pdb %s"%self.replica.pdb)
            self.replica.ref = self.replica.pdb

        ref = reference or os.pardir+os.sep+self.replica.ref

        # Extract number of residues in the solute
        # and build a global mask for the solute
        
        if alignmask:
            refres = alignmask
        elif not self.replica.alignMask or self.replica.alignMask.lower() == 'auto':
            refres = self.getAutoMask()
        else:
            refres = self.replica.alignMask

        self.log.debug("Alignment residue mask: %s"%refres)
        allrefres = self.getAutoMask()

        # Add rmsd and averaging if requested
        # Do rmsd over backbone to align to reference

        if writermsd:
            # Take rmsd of all protein or region
            s = 'rms reference ":{refres}@CA,C,N,O" out {rmsdout}_bb_rmsd.out\n'
            s+= 'rms reference nofit ":{refres} & !@H=" out {rmsdout}_ha_rmsd.out\n'
            self.log.debug("Including rmsd output in trajectory centering script: %s"%s)
            template += s
        else:
            # Just rms but dont write
            s = 'rms reference ":{refres} & @CA,C,N,O"\n'
            template += s

        if writeavgpdb:
            s="average {protavg} :{allrefres} pdb\n"
            self.log.debug("Including protein average output in trajectory centering script: %s"%s)
            template+=s

        # Add trajout command
        template += "trajout {trajout}\n"

        return template.format(trajin=trajin, trajout=trajout, imaging=self.imaging,ref=ref, refres=refres,
                            allrefres=allrefres, protavg=avgpdbout, rmsdout=rmsdoutprefix)

    def getPtrajImagingCommands(self, system=None, pdbFile=None, extraResidues=[]):
        "Usually, when more than one chain is present in the system, amber writes the trajectory splitting the system. \
        This command is useful to obtain the proper ptraj commands to fix this."
        if system:
            pdb = system.getSolvatedPDB()
        elif pdbFile:
            from PDB import SolvatedPDB
            if extraResidues:
                if not isinstance(extraResidues, list): extraResidues = [extraResidues]
            pdb = SolvatedPDB(pdbFile, extraResidues=extraResidues)
        elif self.replica:
            pdb = self.replica.getPDB()
        else:
            raise AmberWriterError, "A PDB is needed to produce ptraj imaging strings. Either pass a pdb path or a SolvatedSystem instance or initialize AmberWriter with a Replica as parameter."

        template = "center :%i-%i mass origin\nimage :* origin center byres familiar\n"
        command = ""

        # Get solute
        pdb = pdb.getSolute()

        # Check PDB chains
        resIdsPerChain = pdb.chainResIDs()
        first = resIdsPerChain[0][0]
        for chain in resIdsPerChain:
            command += template%(first,chain[-1])

        return command

import Biskit.test as BT

class Test(BT.BiskitTest):
    """Test"""
    def test_tleap(self):
        """Test correct operation of tLeap"""
        leap = Leaper(extraff='leaprc.ff99SB')
        leap.command("savepdb LEU tmpleu.pdb")
        leap.close()
        self.testdir = 'tmpleu.pdb'
        self.assertTrue(osp.exists('tmpleu.pdb'))
    
    def test_AmberWriter(self):
        """Create new replica and write MDinput"""
        from Replicas import Replica
        
        top = osp.join(T.testRoot(), 'pep','pep.prmtop')
        crd = osp.join(T.testRoot(), 'pep','pep.prmcrd')
#        checkCommands = [l.strip() for l in open(osp.join(T.testRoot(),
#                                'pep','COMMANDS_test.sh'),'r').readlines()]
        
        self.testdir =  T.tempDir()
        self.r1 = Replica(name='testReplica', top=top, crd=crd, 
                            restrMode='HA', restrForce=0.1,
                            minimizationAsRef=1)
        T.BROWSER.chdir(self.testdir)
        
        # write replica folder and check methods of AmberWriter
        self.r1.createFolder()
        writer = AmberWriter(self.r1)
        
        self.assertEqual(writer.getAmberRestrMask(), ':1-8 & !@H=')
        self.assertTrue(writer.writeCommands())
        self.assertTrue(writer.writeReplicaInput())
        self.assertTrue(writer.getReplicaCommands())
        self.assertTrue(writer.getPtrajAlignScript('md1.nc','md1.nc'))
        
        self.testdir += os.sep+'testReplica'
    
    def test_AmberWriter_ImagingCommand(self):
        w = AmberWriter()
        pdb = osp.join(T.testRoot(), 'pep','pep.pdb')
        aligncommand=w.getPtrajImagingCommands(pdbFile=pdb)
        self.assertEqual(aligncommand, 'center :1-8 mass origin\nimage :* origin center byres familiar\n')
        self.testdir=''
    
    def cleanUp(self):
        T.tryRemove( self.testdir, tree=1 )

if __name__ == "__main__":
    BT.localTest()

