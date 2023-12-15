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

"""
This module provides all classes to automatically prepare PDB files for running a simulation.
If departing from a completely new PDB file, the process would involve the following tasks:
    1 - Remove double occupancy residues, cap N and C terminus for the selected chain. 
    2 - Using PDB2PQR server find optimal protonation states and add hydrogens accordingly.
    3 - Save the protonated and clean PDB as a Amber Object File (OFF).
    
"""

__author__="daniel"
__date__ ="$Feb 10, 2014 4:09:06 PM$"

import re
import string
import time
import logging
import os.path as osp

import numpy as npy

import Biskit as bi
import Biskit.mathUtils as MU
from Biskit.PDBParser import PDBParserError
from Biskit.PDBParseFile import PDBParseFile

import tools as T
import settings as S

try:
    import mechanize
    MECHANIZE_AVAIL = True
except:
    MECHANIZE_AVAIL = False

# Some PDB2PQR class related errors
class PDB2PQRError(Exception):
    pass
class ConnectionError(PDB2PQRError):
    pass
class FormChange(PDB2PQRError):
    pass

class PDB2PQRInterface(object):
    """
    Interface to PDB2PQR webserver for automatically protonate a PDB according to propKa predictions.
    As the server does not provides any non-website interface, it requires :mod:`mechanize` module to work.
    If website changes, this module should be adapted accordingly...
    """
    #: PDB2PQR Webserver
    PDB2PQRWEBSITE = "http://nbcr-222.ucsd.edu/pdb2pqr_1.8/"
    
    #: Form control names
    KNOWNCONTROLS = set(['PDBSOURCE', 'PDBID', 'PDB', 'FF', 'USERFF', 'USERNAMES', 
                        'FFOUT', 'DEBUMP', 'OPT', 'PROPKA', 'PH', 'LIGANDCHECK', 
                        'LIGAND', 'INPUT', 'CHAIN', 'WHITESPACE', 'TYPEMAP', 'NEUTRALN', 
                        'NEUTRALC', None])
    
    MANDATORYCONTROLS = set(['PDBSOURCE','PDBID','PDB','FF','FFOUT','PROPKA','PH'])
    
    def __init__(self):
        self.log = logging.getLogger("PDB2PQRInterface")
        
        # Prepare browser
        self.br = mechanize.Browser()
#        self.br.set_all_readonly(False)    # allow everything to be written to
        self.br.set_handle_robots(False)   # ignore robots
        self.br.set_handle_refresh(False)  # can sometimes hang without this
        self.br.addheaders = [('User-agent', 'Firefox')]
        
        # Open webserver url and check form
        self.setConnection()
        self.checkSite()

    def setConnection(self):
        """
        Connect to webserver.
        :raise PDB2PQRError: Fail in connection.
        """
        try:
            self._html = self.br.open(PDB2PQRInterface.PDB2PQRWEBSITE)
        except Exception, err:
            raise ConnectionError, err
        
        # If correct, select first form
        self.br.select_form(nr=0)
    
    def checkSite(self):
        "Control wether the site has changed the form. If so, exit as everything will fail."
        # Check that html website is similar to the content when this module was written
        import difflib
        checkhtml = open(T.testRoot('pdb2pqr_html.txt'),'r').read()
        webhtml = self._html.read()
        ratio = difflib.SequenceMatcher(lambda x: x == " ", checkhtml, webhtml).ratio()        
        if ratio < 0.9:
            self.log.warn('Website seems different since this module was coded. Similarity: %3f'%ratio)
        
        # Check form has mandatory controls
        controls = set([c.name for c in self.br.form.controls])
        if controls != PDB2PQRInterface.KNOWNCONTROLS:
            # Something has changed!
            # First check if all mandatory controls are there
            self.log.warning("Changes in webserver form.")
            if PDB2PQRInterface.MANDATORYCONTROLS < controls: return True
            else:
                raise FormChange, "Missing controls in form: ", PDB2PQRInterface.MANDATORYCONTROLS - controls            
        else:
            return True

    def __uploadFile(self, pdbfile):
        "Set pdbsource to upload and put file in the form"
        self.br.form.find_control('PDBSOURCE').set(True, 'UPLOAD')
        self.br.form.add_file(open(pdbfile), 'text/plain', pdbfile, name="PDB")
        self.br.form.find_control('PDB').disabled = False

    def __setPDBID(self, pdbid):
        "Sett pdbsource to PDBID and puts id in the textfield"
        self.br.form.find_control('PDBSOURCE').set(True, 'ID')
        self.br.form['PDBID'] = pdbid

    def protonatePDB(self, pdbfile=False, pdbid=False, twait=2, tries=999):
        """
        Main function to
        """
        import fnmatch
        
        if pdbfile:
            self.__uploadFile(pdbfile)
        elif pdbid:
            self.__setPDBID(pdbid)
        else:
            return False

        # Check amber ff and naming
        # CHECK PROPKA CALCULATION ALSO
        self.br.form.find_control('FF').set(True, 'amber')
        self.br.form.find_control('FFOUT').set(True, 'amber')
        self.br.form.find_control('PROPKA').items[0].selected=True
        self.br.form.find_control('PH').disabled=False
        self.br.form.find_control('INPUT').items[0].selected=False   # Don't prepare APBS inputs
        self.br.form.find_control('CHAIN').items[0].selected=True    # Add or keep chain names
        self.br.form.find_control('WHITESPACE').items[0].selected=True # Add whitespaces between all columns

        # Submit and fetch results link
        resp = self.br.submit()
        results_url = None
        for l in self.br.links():
            if l.text == 'click here':
                results_url = l.url
        return self.fetchResults(results_url, twait=twait, tries=tries)
        
    def fetchResults(self, url, twait=3, tries=999):
        "Download PQR file from results url"
        self.br.open(url)
        # Loop until 'complete' Message appears
        match = re.compile('Message: (\w+)')
        while tries:
            response = self.br.reload()
            
            # Fetch Message
            status = match.search(response.read())
            if status: status = status.groups()[0]
            else:
                raise PDB2PQRError, "Error scraping reuslts website. Check website: %s."%self.br.geturl()
            
            if status == 'complete':
                print "Done"
                done = True
                break
            elif status == 'running':
                tries -= 1
                done = False
                print "Running... %d\r"%tries,
            else:
                raise PDB2PQRError, "Error in job execution. Check website: %s."%self.br.geturl()
        
        print
        
        if not done:
            raise PDB2PQRError, "Job execution timed out. Check later website for results: %s."%self.br.geturl()
            
        # If link was found download pqr file
        pqrout = None
        for l in self.br.links():
            if '.pqr' in l.text:
                import tempfile
                pqrout = tempfile.mktemp()
                self.br.retrieve(l.url, pqrout)

        if not pqrout: raise PDB2PQRError, "Unable to download PQR result file."
        #Return pqr as a pdbmodel
        return PQRParseFile(pqrout).getModel()

class PQRParseFile( PDBParseFile ):
    """
    Wraps :mod:`Biskit.PDBParseFile` and adapts it to read PQR formatted files
    in which the columns are always space separated and it doesn't matter the column number.
    """
    def __init__(self, pqr_in, *args, **kwargs):
        """
        Init PDBParseFile, read *pqr_in* PQR file and return a PDBModel.
        """
        PDBParseFile.__init__(self, *args, **kwargs)
        self.inpqr = pqr_in
    
    def getModel(self):
        return self.parse2new(self.inpqr)

    def update( self, model, source, skipRes=None, updateMissing=0, force=0,
                headPatterns=[]):
        """
        Update empty or missing fields of model from the source. The
        model will be connected to the source via model.source.
        Profiles that are derived from the source are labeled 'changed'=0.
        The same holds for coordinates (xyzChanged=0).
        However, existing profiles or coordinates or fields remain untouched.

        @param model: existing model
        @type  model: PDBModel
        @param source: source PDB file
        @type  source: str
        @param skipRes: list residue names that should not be parsed
        @type  skipRes: [ str ]
        @param updateMissing: ignored
        @type  updateMissing: 1|0
        @param headPatterns: [(putIntoKey, regex)] extract given REMARKS
        @type  headPatterns: [(str, str)]

        @raise PDBParserError - if something is wrong with the source file
        """

        try:
            # atoms and/or coordinates need to be updated from PDB
            if force or self.needsUpdate( model ):

                atoms, xyz, info = self.__collectAll( source, skipRes,
                                                        headPatterns )
#                print atoms["serial_number"]
#                for atom in atoms:
#                    print atom
                keys = MU.union( atoms.keys(),  self.DEFAULTS.keys() )

                for k in keys:

                    if model.atoms.get( k, default=0, update=False ) in \
                            (0,None):

                        dflt = self.DEFAULTS.get( k, None )
                        model.atoms.set(k, atoms.get(k, dflt), changed=0 )

                if model.xyz is None:
                    model.xyz = xyz
                    model.xyzChanged = 0

                model._resIndex  =None
                model._chainIndex=None

                model.fileName = model.fileName or source

                model.pdbCode = model.pdbCode or info.get('pdb_code', None) or \
                                self.idFromName( model.fileName)

                model.info.update( info )

        except:
            msg = self._PDBParseFile__xplorAtomIndicesTest( source ) or ' '
            raise PDBParserError('Cannot read ' + str(source) + ' as PQR\n'\
                           '\ERROR: ' + T.lastError() + msg)

        model.setSource( source )
        
    def __collectAll( self, fname, skipRes=None, headPatterns=[] ):
        """
        Parse ATOM/HETATM lines from PDB. Collect coordinates plus
        dictionaries with the other pdb records of each atom.
        REMARK, HEADER, etc. lines are ignored.

        Some changes are made to the dictionary from PDBFile.readline()::
            - the 'position' entry (with the coordinates) is removed
            - leading and trailing spaces are removed from 'name' ..
            - .. but a 'name_original' entry keeps the old name with spaces
            - a 'type' entry is added. Its value is 'ATOM' or 'HETATM'
            - a 'after_ter' entry is added. Its value is 1, if atom is
              preceeded by a 'TER' line, otherwise 0
            - empty 'element' entries are filled with the first non-number
              letter from the atom 'name'

        @param fname: name of pdb file
        @type  fname: str
        @param skipRes: list with residue names that should be skipped
        @type  skipRes: list of str

        @return: tuple of (1) dictionary of profiles
                 and (2) xyz array N x 3
        @rtype: ( list, array )
        """
        xyz   = []

        aProfs = {}

        info = {}

        in_header = True
        headPatterns = headPatterns or self.RE_REMARKS
        patterns = [ (key, re.compile(ex)) for key,ex in headPatterns ]

        for k in bi.PDBModel.PDB_KEYS:
            aProfs[k] = list()

        f = open( fname ,'r' )

        try:
            line, i = ('',''), 0

            while line[0] <> 'END' and line[0] <> 'ENDMDL':

                i += 1
                try:
                    line = f.readline().split()
                except ValueError, what:
                    self.log.add('Warning: Error parsing line %i of %s' %
                                 (i, T.stripFilename( fname )) )
                    self.log.add('\tError: '+str(what) )
                    continue

                if not line: break
                
                ## header handling
                if in_header and line[0] == 'HEADER':
                    info.update( self._PDBParseFile__parseHeader( line ) )

                if in_header and line[0] == 'REMARK':
                    info.update( self._PDBParseFile__parseRemark( line, patterns ) )


                ## preserve position of TER records
#                print line
                newChain = line[0] == 'TER'
                if newChain:
                    line = f.readline().split()

                if (line[0] in ['ATOM','HETATM'] ):

                    if in_header: in_header = False  ## switch off HEADER parsing
                    if len(line) == 11: # contains chain name
                        a = {'serial_number': int(line[1]),
                                 'name': line[2],
                                 'alternate': '',
                                 'residue_name': string.strip(line[3]),
                                 'chain_id': string.strip(line[4]),
                                 'residue_number': int(line[5]),
                                 'insertion_code': '',
                                 'position': map(float,line[6:9]),
                                 'occupancy': 1.0,
                                 'temperature_factor': 0.0,
                                 'segment_id': '',
                                 'element': '',
                                 'charge': line[9]}
                    else:   # without chain name
                        a = {'serial_number': int(line[1]),
                                 'name': line[2],
                                 'alternate': '',
                                 'residue_name': string.strip(line[3]),
                                 'chain_id': '',
                                 'residue_number': int(line[4]),
                                 'insertion_code': '',
                                 'position': map(float,line[5:8]),
                                 'occupancy': 1.0,
                                 'temperature_factor': 0.0,
                                 'segment_id': '',
                                 'element': '',
                                 'charge': line[8]}
#                    print line
#                    print a
                    if skipRes and a['residue_name'] in skipRes:
                        continue

                    a['name_original'] = a['name']
                    a['name'] = a['name'].strip()

                    a['type'] = line[0]

                    if newChain:
                        a['after_ter'] = 1
                    else:
                        a['after_ter'] = 0

                    if a['element'] == '':
                        a['element'] = self._PDBParseFile__firstLetter( a['name'] )

#                    print a
                    xyz.append( a['position'] )
                    del( a['position'])
                    
                    for k, v in a.items():
                        aProfs[k].append( v )

        except:
            raise PDBParserError("Error parsing file "+fname+": " \
                                 + T.lastError())
        try:
            f.close()
        except:
            pass

        if len( xyz ) == 0:
            raise PDBParserError("Error parsing file "+fname+": "+
                            "Couldn't find any atoms.")

        return aProfs, npy.array( xyz, npy.float64 ), info

class AmberPDBCleaner(bi.AmberParmBuilder):
    "Modification of AmberParmBuilder to save an objectfile. Use AmberParmBuilder cleaning capabilities."
    def __init__(self, pdb, *args, **kwargs):
        bi.AmberParmBuilder.__init__(self, model=pdb, *args, **kwargs)
        #self.log = logging.getLogger("AmberPDBCleaner")

    def capACE( self, model, chain ):
        """
        Cap N-terminal of given chain.

        :arg model: model
        :type  model: PDBMode
        :arg chain: index of chain to be capped
        :type  chain: int
        """
        self.log.add('Capping N-terminal of chain %i.' % chain )
        m_ace = bi.PDBModel( self.F_ace_cap )

        chains_before = None
        chains_after = None
        if chain > 0 : chains_before = model.takeChains( range(chain), breaks=1 )
        if chain < model.lenChains(1) - 1: chains_after  = model.takeChains( range(chain+1, model.lenChains(1)),
                                          breaks=1 )
        m_chain       = model.takeChains( [chain], breaks=1 )

        m_term  = m_chain.resModels()[0]

        ## we need 3 atoms for superposition, CB might mess things up but
        ## could help if there is no HN
        if 'HN' in m_term.atomNames():
            m_ace.remove( ['CB'] )

        ## rename overhanging residue in cap PDB
        for a in m_ace:
            if a['residue_name'] != 'ACE':
                a['residue_name'] = m_term.atoms['residue_name'][0]
            else:
                a['residue_number'] = m_term.atoms['residue_number'][0]-1
                a['chain_id']       = m_term.atoms['chain_id'][0]
                a['segment_id']     = m_term.atoms['segment_id'][0]

        ## fit cap onto first residue of chain
        m_ace = m_ace.magicFit( m_term )

        ## concat cap on chain
        ## merge chains
        m_chain = m_ace.resModels()[0].concat( m_chain )
        m_chain.mergeChains(0, renumberAtoms=True)

        ## re-assemble whole model
        assemble_list = []
        if chains_before: assemble_list.append(chains_before)
        assemble_list.append(m_chain)
        if chains_after: assemble_list.append(chains_after)

        if len(assemble_list) > 1:
            return assemble_list[0].concat(*assemble_list[1:])
        else:
            return assemble_list[0]
#        return chains_before.concat( m_chain, chains_after )


    def capNME( self, model, chain ):
        """
        Cap C-terminal of given chain.

        :arg model: model
        :type  model: PDBMode
        :arg chain: index of chain to be capped
        :type  chain: int
        """
        self.log.add('Capping C-terminal of chain %i.' % chain )
        m_nme   = bi.PDBModel( self.F_nme_cap )

        chains_before = None
        chains_after = None
        if chain > 0 : chains_before = model.takeChains( range(chain), breaks=1 )
        if chain < model.lenChains(1) - 1: chains_after  = model.takeChains( range(chain+1, model.lenChains(1)),
                                          breaks=1 )
        m_chain       = model.takeChains( [chain], breaks=1 )

        m_term  = m_chain.resModels()[-1]

        ## rename overhanging residue in cap PDB, renumber cap residue
        for a in m_nme:
            if a['residue_name'] != 'NME':
                a['residue_name'] = m_term.atoms['residue_name'][0]
            else:
                a['residue_number'] = m_term.atoms['residue_number'][0]+1
                a['chain_id']       = m_term.atoms['chain_id'][0]
                a['segment_id']     = m_term.atoms['segment_id'][0]

        ## chain should not have any terminal O after capping
        m_chain.remove( ['OXT'] )

        ## fit cap onto last residue of chain
        m_nme = m_nme.magicFit( m_term )

        ## concat cap on chain
        m_chain = m_chain.concat( m_nme.resModels()[-1] )
        m_chain.mergeChains(0, renumberAtoms=True)
        
        ## should be obsolete now
        if getattr( m_chain, '_PDBModel__terAtoms', []) != []:
            m_chain._PDBModel__terAtoms = [ len( m_chain ) - 1 ]

        ## re-assemble whole model
        assemble_list = []
        if chains_before: assemble_list.append(chains_before)
        assemble_list.append(m_chain)
        if chains_after: assemble_list.append(chains_after)

        if len(assemble_list) > 1:
            return assemble_list[0].concat(*assemble_list[1:])
        else:
            return assemble_list[0]

    def cleanPDB(self, hetatm=1, keepwaters=1, cap=1, capN=[], capC=[], **kw):
        """
        Try to save a clean object file that can be used in simulation from the loaded pdb
        when AmberOFFBuilder was instantiated. Modified method from parmSolvated

        @param f_out: target file for Amber OFF file
        @type  f_out: str
        @param outparm: write topology, coordinates and pdb from saved objectfile.
                        This can be used to check the OFF file was correctly generated.
        @type   outparm: bool

        @param hetatm: keep hetero atoms (default: 1)
        @type  hetatm: 1|0
        @param cap: put ACE and NME capping residue on chain breaks
                    (default: 0)
        @type  cap: 1|0
        @param capN: indices of chains that should get ACE cap (default: [])
        @type  capN: [int]
        @param capC: indices of chains that should get NME cap (default: [])
        @type  capC: [int]
        @param fmod: list of files with amber parameter modifications
                    (to be loaded into leap with loadAmberParams) (default:[])
        @type  fmod: [str]
        @param fprep: list of files with amber residue definitions
                    (to be loaded into leap with loadAmberPrep) (default: [])
        @type  fprep: [str]

        @raise IOError:
        """
        try:
            if self.verbose: self.log.add( 'Cleaning PDB file for Amber:' )
            self.m.setXyz(self.m.xyz - self.m.center()) # center object
            
            if keepwaters:
                wats = self.m.compress(self.m.maskH2O())


            m = self.m.clone()
            m = m.compress(~m.maskH2O())
            m.xplor2amber()

            if cap:
                if m.chainBreaks():
                    end_broken = m.atom2chainIndices( m.chainBreaks() )
                    capC = MU.union( capC, end_broken ) or [0]
                    capN = MU.union( capN, npy.array( end_broken ) + 1 ) or [0]

            for i in capN:
                if self.verbose:
                    self.log.add( 'Adding ACE cap to chain %i' % i )
                if cap: m = self.capACE( m, i )
            
            for i in capC:
                if self.verbose:
                    self.log.add( 'Adding NME cap to chain %i' % i )
                if cap: m = self.capNME( m, i )

            m.renumberResidues( addChainId=1 )  ## again, to accomodate capping

            ss = self._AmberParmBuilder__ssBonds( m, cutoff=4. )
            self._AmberParmBuilder__cys2cyx( m, ss )
            self.leap_ss  = self._AmberParmBuilder__fLines( self.ss_bond, ss )
            if self.verbose:
                self.log.add('Found %i disulfide bonds: %s' % (len(ss),str(ss)))

#            if self.verbose:
#                self.log.add( 'writing cleaned PDB to %s'  % self.clean_pdb )
#            m.writePdb( self.clean_pdb, ter=3, amber=1 )

            if keepwaters: 
                import scipy.spatial as spat
                # Add only waters that do not clash
                dmat = spat.distance_matrix(wats.xyz, m.xyz)
                clashes = npy.where(dmat < 1.5)[0]
                if npy.any(clashes):
                    clashres = wats.atom2resIndices(clashes)
                    # Remove residues which have clashes
                    wats.removeRes(clashres)                
                m = m.concat(wats)
                m.renumberResidues( addChainId=1 )
                
            self.m = m
            
            return m

        except IOError, why:
            raise IOError, why

class AutoPrepareError(Exception):
    pass

class AutoPrepare(object):
    """
    Giving an input PDB, automatically create an amber object file (OFF) that can be simulated.
    Will protonate the system using PDB2PQR server and add capping NME and ACE at each terminus
    and save the resulting correctly protonated pdb for user inspection. 
    
    If everything is correct, the user can come back and finally create the object file.
    """
    def __init__(self, pdb=False, chains=[], protonate=False, pdbid=None, *args, **kwargs):
        "pdb is initial pdbmodel to work with (optional)"
        if pdb:
            self.pdb=bi.PDBModel(pdb)
        else:
            self.pdb=False
        self.pdbid=pdbid
        self.chains=chains

        # If PDBID requested, always execute protonate
        if pdbid:
            self.fetchPDBandProtonate(pdbid, **kwargs)
        elif self.pdb:
            if protonate:
                self.protonatePDB(pdb=self.pdb, **kwargs)
        else:
            pass
        
        # Finally do cleaning if self.pdb was set
        if self.pdb: self.cleanPDB(pdb=None, chains=self.chains, **kwargs)


    def cleanPDB(self, pdb=None, chains=[], hetatm=True, keepwaters=True, cap=True, capC=[], capN=[], **kwargs):
        """
        TODO Document better
        
        - Rename residues according to the protonation state (HID HIE or HIP, etc.)
        - Rename CYS involved in disulfide bonds to CYX.
        - Remove all Hydrogens to let tLeap add them.
        - Cap with NME and ACE the terminus of a PDBModel.
        - Remove or keep the original waters. If original waters are kept, those clashing with the cappings will 
         still be removed.
         
         :return: Cleaned pdb
         :rtype: PDBModel
        """
        if not pdb: pdb=self.pdb
        if chains:
            # Take only chains specified
            pdb = pdb.takeChains(chains)
        cleaner = AmberPDBCleaner(pdb, verbose=True)
        if chains and not capC: capC=range(len(chains))
        if chains and not capN: capN=capC
        self.pdb = cleaner.cleanPDB(hetatm=hetatm, keepwaters=keepwaters, cap=cap, capC=capC, capN=capN, **kwargs)
               
    def protonatePDB(self, pdb=None, tries=50, twait=3, **kwargs):
        """
        Use PDB2PQR webserver to protonate the PDB following Propka predictions.
        PDB can be either a local file or a PDBID.
        
        TODO document
        """
        if pdb:
            if isinstance(pdb, bi.PDBModel):
                file = pdb.source.original()
            elif isinstance(pdb, str):
                file = pdb
            else:
                raise AutoPrepareError, 'pdb argument must be a filepath or a PDBModel with valid source files.'
            
            if not osp.exists(file): raise AutoPrepareError, 'pdbfile %s not found.'%fi
        elif self.pdb:
            file = self.pdb.source.original()
            if not osp.exists(file): raise AutoPrepareError, 'pdbfile %s not found.'%fi
        else:
            raise AutoPrepareError, 'protonatePDB needs a pdb. Set a PDB with setPdb() or give as argument.'
        
        self.b = PDB2PQRInterface()
        self.pdb = self.b.protonatePDB(pdbfile=file,twait=twait,tries=tries, **kwargs)
    
    def fetchPDBandProtonate(self, pdbid, tries=50, twait=5, **kwargs):
        "Fetch a pdb by PDBID and protonate using PDB2PQR server"
        self.b = PDB2PQRInterface()
        self.pdb = self.b.protonatePDB(pdbid=pdbid,twait=twait,tries=tries, **kwargs)
    
    def setPdb(self, pdb):
        """
        Set PDB to prepare.

        :arg str|PDBModel pdb: PDB file path or PDBModel.
        :raises AutoPrepareError: If argument of wrong type.
        """
        if isinstance(pdb, str) and osp.exists(pdb): self.pdb = bi.PDBModel(pdb)
        elif isinstance(pdb, bi.PDBModel): self.pdb = pdb
        else: raise AutoPrepareError, "setPdb argument should be a pdb filepath or a PDBModel"

    def getPdb(self):
        return self.pdb
    
    def savePdb(self, outname):
        "Save current PDB"
        self.pdb.writePdb(outname)
    
    def saveOFF(self, outname, inpdb=None, unitname='sys', extraff=[]):
        "From a PDBModel or File, load into tLeap and save as ObjectFile."
        from Amber import AmberCreateSystem
               
        if not inpdb and self.pdb: inpdb = self.pdb
        else: raise AutoPrepareError, "Input needed."

        prepare = AmberCreateSystem()
        prepare.createOFF(outname, inpdb, extraff=extraff, **kwargs)
        outoff = outname+'.lib'

        return osp.abspath(outoff)
          
# AUXILIAR FUNCTIONS
def prepareOFF(pdbinput, outfname, ff=[]):
    """Prepare a Amber Object File Format (OFF) file from a PDB. 
    :arg str pdbinput: Path to existing PDB file
    :arg str outfname: Output filename
    :arg list ff: Forcefield / Frcmod files to load. If empty, use defaults.
    """
    pass
          
          
# TESTING
import test as BT

class Test(BT.BiskitTest):
    """Test"""
    TAGS = [BT.FUTURE]

    def comparePDB(self, a, b, numatoms=30):
        "Compare first *numatoms* of two pdb models for matching name, serial_number and coordinates"
        self.assertEqual(a['name'][:numatoms], b['name'][:numatoms])
        self.assertEqual(a['serial_number'][:numatoms].tolist(), b['serial_number'][:numatoms].tolist())
        self.assertEqual(a.xyz[:numatoms,].tolist(), b.xyz[:numatoms,].tolist())

    def test_PQRParser(self):
        """PQRParser test. Read a PQR file"""
        f_in = T.testRoot('autoprepare', '1yer.pqr')
        f_check = T.testRoot('autoprepare', '1yer_h.pdb')
        self.checkpdb = bi.PDBModel(f_check)
        self.pqr = PQRParseFile(f_in).getModel()
#        self.pqr.writePdb(f_check)
        self.comparePDB(self.pqr, self.checkpdb)
       
    def test_PDB2PQR_byFile(self):
        """Test PDB2PQR interface. Protonate Uploaded PDB file."""
        f_in = T.testRoot('autoprepare', '1yer.pdb')
        pqr_check = PQRParseFile(osp.join(T.testRoot(), '1yer.pqr')).getModel()
        self.b = PDB2PQRInterface()
        pqr_out = self.b.protonatePDB(pdbfile=f_in,twait=2)
        self.comparePDB(pqr_out, pqr_check)
        
    def test_PDB2PQR_byID(self):
        """Test PDB2PQR interface. Protonate PDB selected by ID. """
        pqr_check = PQRParseFile(T.testRoot('autoprepare', '1yer.pqr')).getModel()
        self.b = PDB2PQRInterface()
        pqr_out = self.b.protonatePDB(pdbid='1yer',twait=2)
        self.comparePDB(pqr_out, pqr_check)

    def test_autoprepare_1(self):
        """Test main AutoPrepare class. Prepare existing pdb file."""
        f_in = T.testRoot('autoprepare', '1yer.pdb')
        f_check = T.testRoot('autoprepare', '1yer_h.pdb')
        self.prepare = AutoPrepare(pdb=f_in, protonate=True)
        self.checkpdb = bi.PDBModel(f_check)
        self.comparePDB(self.prepare.getPdb(), self.checkpdb)

    def test_autoprepare_2(self):
        """Test main AutoPrepare class. Prepare PDBID 1a0o chains A, B == [0,1]"""
        f_check = T.testRoot('autoprepare', '1yer_h.pdb')
        self.prepare = AutoPrepare(pdbid='1yer', chains=[0,1])
        self.checkpdb = bi.PDBModel(f_check)
#        self.savePdb(osp.join(T.testRoot(),'1a0o_prepared.pdb'))
        self.comparePDB(self.prepare.getPdb(), self.checkpdb)

#    def cleanUp(self):
#        T.tryRemove( self.f_out )

if __name__ == "__main__":
    #BT.localTest()
    Test().autoprepare_2()
    
