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
__date__ ="$08-abr-2014 0:22:58$"

import os
import os.path as osp
import logging

import numpy as npy

import settings as S
import tools as T
from Replicas import Replica
from GridsManager import Grid, GridSpace

class EnergyConversionError(Exception):
    pass

class EnergyConversion(object):
    def __init__(self):
        self.log = logging.getLogger('EnergyConversion')

    def averageBoltzmann(self, gridlist, temp=300):
        """
        Given a list of grids, return the boltzmann average grid

        :arg list gridlist: List of paths to grid files or :class:`Grid` instances.
        :arg float temp: Temperature for boltzmann expression

        :return: :class:`Grid` instance with average grid
        """
        if not isinstance(gridlist, list): gridlist = [gridlist]

        if len(gridlist) == 1:
            avgGrid = Grid(gridlist[0])   # Not doing average!
            self.log.warn("(averageBoltzmann) Only one replica for averaging! Return grid unmodified: %s"%gridlist[0])
        else:
            gspace = GridSpace()
            gspace.setT(temp)
            gspace.loadGrids(gridlist)
            avgGrid = gspace.averageGrids()
        return avgGrid

    def calcPDBRatio(self, pdb, unit, refunit='WAT'):
        """
        Calculate the correction ratio to apply to expected number before energy conversion.
        The expected number is calculated from the pre-equilibrated solvent box. In simulation setup, the real
        number of each residue finally added to the system may vary. For this reason, the expected number must be adjusted
        to match the real concentration.

        One way of doing this is to calculate the relative proportion of *unit* residue with respect to water in the
        solvated PDB and compare it with the pre-equilibrated ratio.

        :return: float with ratio refunit/unit
        """
        if isinstance(pdb, str):
            from PDB import SolvatedPDB
            pdb = SolvatedPDB(pdb)

        numres = pdb.getNumResidues()
        unitnum = numres.get(unit)
        refnum = numres.get(refunit)
        if not unitnum or not refnum: raise EnergyConversionError, "Unit %s or refunit %s not found in pdb %s"%(unit, refunit, pdb.source)
        return refnum/float(unitnum)

    def calcReplicaExpectedValue(self, replica, probe, numsnaps=False, gridspacing=S.GRID_SPACING):
        "Return the expected density value per grid voxel for the probe given. Will consider the number of snapshots in replica instance"
        # Voxel volume
        if (isinstance(gridspacing, list) or isinstance(gridspacing, npy.ndarray)) and len(gridspacing) == 3: voxel = npy.prod(gridspacing)
        else: voxel=gridspacing**3

        # Calculate number of snapshots if not given
        if not numsnaps:
            numsnaps = replica.ntrajfiles*(replica.prod_steps/float(replica.trajfrequency))

        # Fetch unit corresponding to probe and count residues in pdb
        solv = replica.getSolvent()
        if probe in solv.comprobes.keys(): unit = solv.comprobes[probe].name
        elif probe in solv.probelist: unit = solv.getProbeByName(probe).residue.name
        else:
            raise EnergyConvertError, "Probe %s not found for solvent in replica %s"%(probe, replica.name)

        # Calculate ratio-correction factor
        if unit == 'WAT' and len(solv.residues) > 1:
            solv.residues.remove(unit)
            refunit = solv.residues[0].name
            self.log.debug("Setting residue %s as reference for ratio-correcting water probe"%refunit)
        else:
            refunit = 'WAT'

        self.log.debug("Unit: %s RefUnit: %s"%(unit, refunit))
        expectedratio = solv.getNumRes(refunit) / float(solv.getNumRes(unit))
        observedratio = self.calcPDBRatio(replica.getPDB(), unit, refunit)
        ratio_correction = expectedratio/observedratio
        self.log.debug("Expected ratio: %.4f Observed ratio: %.4f Correctionfactor: %.4f"%(expectedratio, observedratio, ratio_correction))

        # Finally calculate expected number
        expectedVal = solv.getProbeProbability(probe)*numsnaps*ratio_correction
        self.log.debug("Replica %s probe %s. ExpectedValue: %.3f, snapshots: %i, Voxel: %.3f"%(replica.name, probe, expectedVal, numsnaps, voxel))
        return expectedVal

    def calcDG0correction(self, replicaList, unit, **kwargs):
        """
        Calculate the standard dg penalty by volume method.
        If several replicas are given, a mean penalty is calculated.

        Args:
            ReplicaList     (list|ReplicaInfo)      List of ReplicaInfo instances or a single instance.
                                                    If a list is given, the mean penalty is calculated.
            unit            (str)                   Residue name to consider (usually the residue name for the probe in study).

        :kwargs: Parameters of calcReplicaBoxVolume accepted

        Returns:
            dg0 penalty to be added up to the dg grid.
        """
        if not isinstance(replicaList, list): replicaList = [replicaList]

        volumes = []
        numres = []
        temp = []
        for replica in replicaList:
            simvolume = self.calcReplicaBoxVolume(replica, **kwargs)
            if not simvolume:
                self.log.warning("Replica %s simulation volume not found, applying zero penalty!"%replica.name)
                return 0
            numUnitRes = replica.getPDB().getNumResidues()[unit]
            volumes.append(simvolume)
            numres.append(numUnitRes)
            temp.append(replica.temp)

        meanvol = npy.mean(volumes)
        meannum = npy.mean(numres)
        RT = npy.mean(temp) * 0.001987
        return -RT*npy.log((float(meanvol)/meannum)/1660.5)

    def calcReplicaBoxVolume(self, replica, step=False, boxextension=False, **kwargs):
        """Calculate Volume of simulation. Will be done with last production nanosecond volume information.

        :args int step: Nano file to use for fetching box information. If False, Take last nano.
        :args str boxextension: File extension for restart files to search for box information.
            If **False**, will take default values (NAMD: xsc, AMBER: rst)

        """
        step = step or replica.ntrajfiles
        # Fetch last nanosecond of production volume information
        if replica.mdProgram == 'AMBER':
            # Fetch rst file and read last line to get box side length
            boxextension = boxextension or 'rst'
            fname = replica.mdoutfiletemplate.format(step=step, extension=boxextension)
            fname = osp.join(replica.path, replica.mdfolder, fname)
            if not os.path.exists(fname):
                self.log.error("No file found with name %s to fetch box volume in DG0 penalty calculation. Returning no penalty..."%fname)
                return False
            box = map(float, open(fname,'r').readlines()[-1].strip().split())
            if box[3] == 90.0:
                cubic = True
            else:
                cubic = False
            vol = box[0]*box[1]*box[2]
        elif replica.mdProgram == 'NAMD':
            # Fetch xsc file and read first length as cube side length
            boxextension = boxextension or 'xsc'
            fname = replica.mdoutfiletemplate.format(step=step, extension=boxextension)
            fname = osp.join(replica.path, replica.mdfolder, fname)
            if not os.path.exists(fname):
                self.log.error("No file found with name %s to fetch box volume in DG0 penalty calculation. Returning no penalty..."%fname)
                return False
            box = map(float, open(fname,'r').readlines()[2].strip().split())
            a = box[1]
            b = box[4]
            if a!=b: cubic = True # all sides same size, cubic box
            else: cubic = False   # if different, orthorombic box
            vol = a**3

        if not cubic: vol *= 0.77 # orthorombic volume correction
        return vol

    def count2DGProbe(self, grid, numsnaps, probename=None, temp=300., correction=1.):
        """
        Given a name of a probe, convert the grid to energies using probability of the probe.
        Will return a Grid instance.
        """
        from Solvents import SolventManager as SM

        # Load grid if not already done
        if isinstance(grid, str): grid = Grid(grid)
        if not probename: probename = grid.probe

        # Fetch probability
        sm = SM()
        solv = sm.fetchSolventByProbe(probename)
        if not solv: raise EnergyConversionError, "Probename %s not found."%probename
        prob = solv.getProbeProbability(probename)

        # Convert
        grid.update(grid.count2DG(prob*numsnaps, T=temp, correction=correction))
        grid.setType('MDMIX_RAW')
        grid.setProbe(probename)
        return grid

    def convert(self, replicalist, probelist=[], average=False, dg0correct=True, inprefix=None,
                    outprefix=None, nsnaps=None, outpath="PROBE_AVG", protvalue=999.,**kwargs):
        """
        Convert density grids to energies.
        If *average* is **True**, density grids will be added up before conversion. Only one grid will be saved inside *outpath* folder in project home directory.
        If average is **False**, one energy grid will be saved for each replica independently inside each replica energy grid folder.

        :arg list replicalist: List of replicas to fetch density grids to convert.
        :arg list probelist: List of probes to convert. If empty, convert all available.
        :arg bool average: Return an average energy grid from all the input densities.
        :arg bool dg0correct: Apply standard state correction.
        :arg str outprefix: Prepend prefix to output grid.
        :arg int nsnaps: Number of snapshots with wich density grid was constructed. If None, use expected number of snapshots according to replica simulation settings.
        :arg float protvalue: Identify zeros in density grid and substitute energy at these regions with *protvalue*. Usually zeros occur at protein occupied regions.
        """
        # check types
        if not isinstance(replicalist, list): replicalist = [replicalist]
        solvents = []
        probes = []
        for r in replicalist:
            if not isinstance(r, Replica): raise EnergyConversionError, "Wrong type %s. Expected Replica type."%type(r)
            solvents.append(r.solvent)
            probes.extend(r.getProbes())

        # Replica list of names
        replnames = [r.name for r in replicalist]
        if len(replnames) == 1: average = False # Only one replica, convert to energy only for that replica

        # If averaging, all replicas must be of same solvent
        if len(set(solvents)) > 1 and average:
            raise EnergyConversionError, "Cannot average grids from replicas run with different solvents: %s"%replnames

        # Check probes
        if probelist and not set(probelist) < set(probes):
            raise EnergyConversionError, "Some selected probes are not found in replica list: %s"%list((set(probelist) - set(probes)))
        if not probelist: probelist = list(set(probes))

        # set empty outprefix if not given
        if not outprefix: outprefix = ''
        if not inprefix: inprefix = ''

        # Go Home and create outputpath if missing
        T.BROWSER.goHome()

        # Fetch all density grids for all replicas
        allgrids = dict((r.name,r.getGridsByType('MDMIX_DENS', prefix=inprefix)) for r in replicalist)
        self.log.debug("allgrids: %s"%allgrids)

        # Will proceed different if average is requested
        allSaved = True
        # AVERAGE MODE
        if average:
            self.log.info("Converting and averaging grids for replicas %s"%replnames)
            for probe in probelist:
                # Check if probe is present in all replicas
                numok = sum([probe in v['MDMIX_DENS'].keys() for v in allgrids.values()])
                if numok != len(replicalist):
                    self.log.warn("Skipping probe %s. Not found in all replicas."%probe)
                    continue

                # All Ok, get grids
                self.log.info("Probe %s..."%probe)
                grids = [v['MDMIX_DENS'][probe] for v in allgrids.values()]
                self.log.debug("Grids: %s"%[g.source for g in grids])

                # Calc expected values for each replica and add them up
                expectedvals = [self.calcReplicaExpectedValue(r, probe=probe, numsnaps=nsnaps) for r in replicalist]
                sumexpectval = npy.sum(expectedvals)

                # Averaging if more than 1 replica
                if len(replicalist) == 1:
                    sumDensityGrid = grids[0]   # Not doing average!
                    self.log.warn("(sumDensitiesAndConvert) Only one replica for averaging! Can not do average of the grids. Saving same grid with avg output name. Replica: %s"%replicalist[0].name)
                    t = replicalist[0].temp
                else:
                    gspace = GridSpace()
                    t = npy.mean([repl.temp for repl in replicalist]) # Average temperature for energy conversion
#                    gspace.setT(t)
                    gspace.loadGrids(grids)
                    sumDensityGrid = gspace.sum()

                self.log.info("Expected value: %.3f, Temperature: %.2f"%(sumexpectval, t))
                # Convert sumDensityGrid to free energies and apply penalty
                RT = t* 0.001987
                
                # Get a mask of zero values (probably correspoding to protein positions
                maskzeros = sumDensityGrid.data == 0
                sumDensityGrid.data[maskzeros] = 1 # avoid zero-divisions
                dgData = -RT*npy.log(sumDensityGrid.data / float(sumexpectval))
                dgData[maskzeros] = protvalue # set back value to identify zeros

                suffix = '_DG'
                if dg0correct:
                    solv = replicalist[0].getSolvent()
                    if probe in solv.comprobes.keys(): unit = solv.getProbeByName(probe).name
                    else: unit = solv.getProbeByName(probe).residue.name
                    if unit != 'WAT' and unit !='HOH':
                        dgcorrection = self.calcDG0correction(replicalist, unit=unit) # Mean DG0 correction for all replicas
                        dgData += dgcorrection
                        dgData[maskzeros] = protvalue # set back value to identify zeros
                        suffix = '_DG0'
                        self.log.info("Applying DG0 standard state correction %.3f"%dgcorrection)

                sumDensityGrid.update(dgData)

                # Saving to output directory
                if not osp.exists(outpath): os.mkdir(outpath)
                outname = outpath+os.sep+outprefix+probe+suffix+'.dx'
                sumDensityGrid.setType('MDMIX_RAW_AVG')
                sumDensityGrid.setProbe(probe)
                self.log.info("Saving grid %s..."%outname)
                sumDensityGrid.writeDX(outname)

                if osp.exists(outname):
                    self.log.info("Energy averaged grid %s saved. OK."%outname)
                else:
                    self.log.warn("Could not save %s grid. Check reason!"%outname)
                    allSaved = False
        
        #########################
        # INDEPENDENT SAVE MODE #
        #########################
        else:
            # Save grids independently inside each replica
            for replica in replicalist:
                self.log.info("Converting grids for replica %s"%replica.name)
                grids = allgrids[replica.name]['MDMIX_DENS']
                for probe, g in grids.iteritems():
                    if probe in probelist:
                        self.log.info("Probe %s..."%probe)
                        # Convert counts to DG
                        expectedval = self.calcReplicaExpectedValue(replica, probe, numsnaps=nsnaps)
                        dgData = g.count2DG(expectedval, T=replica.temp, maskvalue=protvalue)
                        suffix = '_DG'
                        if dg0correct:
                            solv = replica.getSolvent()
                            if probe in solv.comprobes.keys(): unit = solv.getProbeByName(probe).name
                            else: unit = solv.getProbeByName(probe).residue.name
                            if unit != 'WAT' and unit !='HOH':
                                dgcorrection = self.calcDG0correction(replica, unit=unit) # Mean DG0 correction for all replicas
                                dgData += dgcorrection
                                suffix='_DG0'
                            
                        g.update(dgData)
                        g.setType('MDMIX_RAW')
                        g.setProbe(probe)
                        outpath = replica.energypath
                        if not osp.exists(outpath): os.mkdir(outpath)
                        outname = outpath+os.sep+outprefix+probe+suffix+'.dx'
                        g.writeDX(outname)
                        if osp.exists(outname):
                            self.log.info("Energy grid %s saved. OK."%outname)
                        else:
                            self.log.warn("Could not save %s grid. Check reason!"%outname)
                            allSaved = False
        if allSaved:
            self.log.info("Done with all probes for replicas %s"%replnames)
            return True
        else:
            self.log.warn("Some averaged grids could not be saved for replicas %s"%replnames)
            return False


if __name__ == "__main__":
    print "Hello World"
