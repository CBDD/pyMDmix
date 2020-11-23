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

"""
Settings (Biskit module wrapped)
================================

All this work has been adapted from original Biskit.settings work by
Raik. All merits to him!

This module provides global settings as fields. Throughout MDMix a
(environment-dependent) parameter such as, e.g., ssh_bin can be addressed as:

  >>> import pyMDMix.settings as S
  >>> bin = S.ssh_bin

However, since a user should not be required to hack python modules,
ssh_bin is not actually defined in settings.py. Instead, the value is
taken from C{~/.mdmix/settings.cfg} -- which should have an entry
like C{ssh_bin=/bin/ssh # comment}. If this entry (or the config file)
is not found, settings.py uses the default value from
C{MDMix/data/defaults/settings.cfg}.

If missing, the user configuration file C{~/.mdmix/settings.cfg} is
created automatically during the startup of MDMix (i.e. for any
import). The auto-generated file only contains parameters for which
the default values don't seem to work (invalid paths or binaries).

See L{Biskit.SettingsManager}

Summary for pyMDMix users
-------------------------
  If you want to change a pyMDMix parameter, do so in C{~/.mdmix/settings.cfg}

Summary for MDMix developpers
------------------------------
  Check original Biskit.settings for more info.

Summary for all
---------------
  !Dont't touch C{settings.py}!
"""

import os
import sys
import logging
import user
import os.path as osp
import tools as T
import SettingsParser as P

VERSION="0.1"

# First of all check if logging is activated or activate it!
# Check if a root logger has been created or run one with basic configuration
LOGLEVEL = {'INFO':logging.INFO, 'INFO_TIME':logging.INFO, 'DEBUG':logging.DEBUG, 'WARN':logging.WARN, 'ERROR':logging.ERROR}
rootlog = logging.getLogger()
if not rootlog.handlers:
    # Build a root logger with INFO level and custom format
    logging.basicConfig(level=logging.INFO)
    rootlog = logging.getLogger()
    rootlog.handlers[0].setFormatter(T.LogFormatter())

# USER MDMIX HOME
# change this to mdmix package data folder
#USER_MDMIX_HOME = osp.join(user.home, '.mdmix')
USER_MDMIX_HOME = osp.join(osp.dirname(__file__), 'data')

# Parse settings from package defaults and user defined
__CFG_DEFAULT = osp.join(T.dataRoot('defaults'),'settings.cfg')
__CFG_USER    = osp.join(USER_MDMIX_HOME,'settings.cfg')
try:
    m = P.SettingsManager(__CFG_DEFAULT, __CFG_USER, createmissing=True  )
    m.updateNamespace( locals() )
except Exception, why:
    raise P.SettingsError, 'Error importing pyMDMix settings: %s'%why

## REPLICA DEFAULTS PATHS
CFG_MD_DEFAULT = osp.join(T.dataRoot('defaults'),'md-settings.cfg')
CFG_MD_USER    = osp.join(USER_MDMIX_HOME,'md-settings.cfg')
try:
    m = P.SettingsManager(CFG_MD_DEFAULT, CFG_MD_USER, createmissing=True  )
    m.updateNamespace( {} )
except Exception, why:
    raise P.SettingsError, 'Error importing pyMDMix MD default settings: %s'%why


##
## Create some settings on the fly
##

# Directories of package files
ROOT = T.projectRoot()
DATAROOT = T.dataRoot()
TEMPLATE_DIR = T.dataRoot('templates')
SOLVENTS_DIR = T.dataRoot('solventlib')
TEST_DIR =  T.dataRoot('test')

#STEPDICT = {"DIRECTORIES":False,"PLACING":False,"SOLVATING":False,"MDINPUT":False,
#            "QUEUESCRIPTS":False, "READYTOSIM":False, "SIMDONE":False, "CENTERING":False, "CENTERSCRIPTS":False,
#            "DENSITYSCRIPTS":False, "DENSITYGRIDS":False,"RAWSCRIPTS":False,"RAWGRIDS":False,"RAWDGAVERAGING":False,
#            "RAWDGNOAVG":False, "CORRECTEDGRIDS":False}

# CHECK AMBER INSTALLATION
AMBERHOME = os.environ.get('AMBERHOME') or sys.exit('AMBERHOME env variable not defined!')

# Fetch executables folder
if osp.exists(AMBERHOME+os.sep+'exe'):
    AMBEREXE = AMBERHOME+os.sep+'exe'
elif osp.exists(AMBERHOME+os.sep+'bin'):
    AMBEREXE = AMBERHOME+os.sep+'bin'
else:
    AMBEREXE = os.environ.get('AMBEREXE') or sys.exit('AMBEREXE env variable not defined and executables folder not found automatically. Define AMBEREXE pointing to the folder containing executable files.')

# Check ptraj, ambpdb and tleap are found
for prog in (AMBER_TLEAP, AMBER_PTRAJ, AMBER_AMBPDB):
    if not osp.exists(AMBEREXE+os.sep+prog): sys.exit("%s not found"%prog)

# Check default forcefields are found to be loaded with tleap
__ffok = False
__defff = set(DEF_AMBER_FF)
for root, dirs, files in os.walk(AMBERHOME):
    if set(files) & __defff: __ffok = True

if not __ffok: rootlog.warn("Default Forcefields to load with tLeap (%s) not found under AMBERHOME directory!"%(",".join(DEF_AMBER_FF)))


# WORK ON SOLVENT DATABASE
# Check if SOLVENTS.db file exists in users .mdmix folder and set that one to default
# Else use package database
USER_SOLVENT_DB = osp.join(USER_MDMIX_HOME,'SOLVENTS.db')
SOLVENTDB = SOLVENTS_DIR+os.sep+'SOLVENTS.db'

# Check for pycuda installation
# Minimum version 2013!
try:
    import pycuda
    PYCUDA_VERSION = pycuda.VERSION[0]
    if PYCUDA_VERSION >= 2013: PYCUDA_AVAIL = True
    else: PYCUDA_AVAIL = False
except ImportError:
    PYCUDA_AVAIL = False

# Check if we can use CUDA
# User should define an env variable USE_CUDA = True
USE_CUDA = os.environ.get('USE_CUDA')
if USE_CUDA: USE_CUDA = USE_CUDA.lower().strip() == 'true'

# NETCDF
# Check if netcdf IO is available
NETCDF_AVAIL = False
try:
    import Scientific.IO.NetCDF as nc
    NETCDF_AVAIL=True
except ImportError:
    pass


## CREATE A GLOBAL BROWSER ALL CLASSES CAN ACCESS
from Browser import Browser
global BROWSER
BROWSER = Browser()

if DEBUG==1:
    import atexit
    atexit.register(T.traceback_plus)


# Set logging
def setLogger(level=None, logFile=None):
    "Build root logger"
    import sys
    import logging
    from tools import LogFormatter

    if not level: level = 'INFO'
    if level.upper() not in LOGLEVEL.keys():
        raise KeyError, "%s level not valid. Valid names: %s"%(level.upper(),LOGLEVEL.keys())

    rootlog = logging.getLogger()
    rootlog.setLevel(level.upper())

    if logFile:
        handle = logging.FileHandler(logFile, mode='w')
    else:
        handle = logging.StreamHandler()

    rootlog.handlers = []
    rootlog.addHandler(handle)
    rootlog.handlers[0].setFormatter(LogFormatter())

    if level.upper() != 'DEBUG':
        sys.tracebacklimit = 0


######################
## clean up name space

del T, user
del __CFG_DEFAULT, __CFG_USER, m
del __defff, __ffok, root, dirs, files
del rootlog, os, osp, logging, testparam, Browser

if __name__ == '__main__':
    for k, v in locals().iteritems():
        print k, v