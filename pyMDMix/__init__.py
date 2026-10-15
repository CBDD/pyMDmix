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
##              TOBEPUBLISHED^
##
##  --------------------------------------------------------------------

__author__="dalvarez"
__date__ ="$19-ene-2014 18:43:19$"
__version__="0.2.8"
from . import settings as S # A basic logger is initiated here
from . import tools

browser = tools.BROWSER

from .Projects import Project, createProject, loadProject
from .Replicas import Replica, loadReplica
from .MDSettings import MDSettings, parseSettingsConfigFile
from .Systems import System, SolvatedSystem, parseSystemConfigFile, loadSystem
from .Solvents import SolventManager, getSolvent
from .GridsManager import Grid, NewGrid
from .PDB import SolvatedPDB
from . import QueueWriting as Queue
from . import Analysis


import atexit
atexit.register(tools.EXECUTOR.terminate)

#from PDB import SolvatedSystem
from .settings import setLogger

class MDMixError(Exception):
    pass
